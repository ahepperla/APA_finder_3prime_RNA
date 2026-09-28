"""The one-pass alignment scan must reproduce the per-source reference passes."""

from __future__ import annotations

import json
from itertools import groupby
from pathlib import Path

import pysam
import pytest
from alignment_builders import (
    aligned_segment,
    sorted_copy,
    to_cram,
    write_alignment,
    write_reference,
)

from pacusage.cli import _iter_observation_file, main
from pacusage.errors import PacusageError
from pacusage.evidence import (
    extract_evidence,
    extract_splice_continuations,
    write_bedgraphs,
    write_evidence,
    write_splice_continuations,
)
from pacusage.scan import (
    _scan_paired_end,
    _SourceAccumulator,
    _SpliceAccumulator,
    read_scan_manifest,
    scan_alignment,
    write_scan,
)
from pacusage.tableio import write_tsv

CONTIGS = [("chr1", 5000), ("chr2", 5000), ("1", 5000), ("MT", 5000), ("chrM", 5000), ("0", 5000)]
CHR1, CHR2, ONE, MT, CHRM, ZERO = range(len(CONTIGS))
# "MT" becomes the excluded "chrM" only after the raw-name exclusion check,
# and "0" becomes "chr3", which sorts differently from its raw name.
ALIASES = {"1": "chr1", "MT": "chrM", "0": "chr3"}
FILTER_SETS = {
    "default": {
        "min_mapq": 20,
        "require_unique": True,
        "require_proper_pair": True,
        "exclude_duplicates": True,
        "excluded_contigs": ["chrM"],
        "contig_aliases": ALIASES,
    },
    "permissive": {
        "min_mapq": 0,
        "require_unique": False,
        "require_proper_pair": False,
        "exclude_duplicates": False,
        "excluded_contigs": [],
        "contig_aliases": {},
    },
    "mixed": {
        "min_mapq": 30,
        "require_unique": True,
        "require_proper_pair": False,
        "exclude_duplicates": False,
        "excluded_contigs": ["chr2"],
        "contig_aliases": {"1": "chr1"},
    },
}
SE_SOURCES = ("polyA_junction", "read_3p", "read_5p")
PE_SOURCES = ("fragment_3p", "polyA_junction", "read_3p", "read_5p")


def single_end_records() -> list[pysam.AlignedSegment]:
    r = aligned_segment
    return [
        r("plain", 0, CHR1, 100, "30M"),
        r("same_end", 0, CHR1, 100, "30M"),
        r("reverse", 16, CHR1, 200, "30M"),
        r("low_mapq", 0, CHR1, 250, "30M", mapq=5),
        r("poly_a", 0, CHR1, 300, "30M10S", right="A" * 10),
        r("poly_t_reverse", 16, CHR1, 400, "10S30M", left="T" * 10),
        r("six_with_five_a", 0, CHR1, 500, "30M6S", right="AAAAAC"),
        r("five_a_too_short", 0, CHR1, 600, "30M5S", right="AAAAA"),
        r("ten_with_seven_a", 0, CHR1, 700, "30M10S", right="AAAAAAACCC"),
        r("ten_with_eight_a", 0, CHR1, 750, "30M10S", right="CCAAAAAAAA"),
        r("both_clips", 0, CHR1, 800, "5S30M8S", left="TTTTT", right="A" * 8),
        r("both_clips_reverse", 16, CHR1, 850, "8S30M5S", left="T" * 8, right="AAAAA"),
        r("duplicate", 1024, CHR1, 100, "30M"),
        r("spliced", 0, CHR1, 900, "10M100N20M"),
        r("twice_spliced_reverse", 16, CHR1, 1200, "15M50N15M40N10M"),
        r("deletion", 0, CHR1, 1500, "10M5D20M"),
        r("match_mismatch", 0, CHR1, 1600, "10=2X18="),
        r("insertion", 16, CHR1, 1700, "10M3I20M"),
        r("no_sequence", 0, CHR1, 1800, "30M10S", sequence=False),
        r("unmapped", 4, CHR1, 1900, "30M"),
        r("secondary", 256, CHR1, 2000, "30M"),
        r("supplementary", 2048, CHR1, 2100, "30M"),
        r("qc_failed", 512, CHR1, 2200, "30M"),
        r("mapq_25", 0, CHR1, 2450, "30M", mapq=25),
        r("excluded", 0, CHRM, 100, "30M"),
        r("alias_of_excluded", 0, MT, 100, "30M"),
        r("multimapped", 0, CHR1, 2500, "30M", nh=2),
        r("unique_tag", 0, CHR1, 2600, "30M", nh=1),
        r("alias_merge", 0, ONE, 100, "30M"),
        r("alias_sorts_differently", 16, ZERO, 300, "30M"),
        r("other_contig", 0, CHR2, 100, "30M"),
        r("other_contig_poly_t", 16, CHR2, 200, "6S30M", left="TTTTTT"),
        r("spliced_other_contig", 16, CHR2, 400, "20M30N10M"),
        r("spliced_alias", 0, ONE, 900, "10M100N20M"),
    ]


def paired_end_records() -> list[pysam.AlignedSegment]:
    r = aligned_segment
    groups = [
        [r("fr", 99, CHR1, 100, "30M"), r("fr", 147, CHR1, 250, "30M")],
        # Mate 2 is listed first in some groups: in coordinate-sorted input it
        # often precedes read 1, which must still be the selected read.
        [r("fr_same", 147, CHR1, 250, "30M"), r("fr_same", 99, CHR1, 100, "30M")],
        [r("rf", 163, CHR1, 300, "30M"), r("rf", 83, CHR1, 400, "30M")],
        [r("mate_past_read1", 99, CHR1, 500, "30M"), r("mate_past_read1", 147, CHR1, 520, "40M")],
        [r("improper", 97, CHR1, 600, "30M"), r("improper", 145, CHR1, 700, "30M")],
        [r("orphan", 99, CHR1, 800, "30M")],
        [r("unmapped_mate", 73, CHR1, 900, "30M"), r("unmapped_mate", 133, CHR1, 900, "30M")],
        [
            r("interchromosomal", 97, CHR1, 1000, "30M"),
            r("interchromosomal", 145, CHR2, 1000, "30M"),
        ],
        [
            r("extras", 99, CHR1, 1100, "30M"),
            r("extras", 147, CHR1, 1200, "30M"),
            r("extras", 355, CHR2, 50, "30M"),
            r("extras", 2195, CHR1, 3000, "20M"),
        ],
        [
            r("three_primaries", 99, CHR1, 1250, "30M"),
            r("three_primaries", 147, CHR1, 1300, "30M"),
            r("three_primaries", 147, CHR1, 1310, "30M"),
        ],
        [r("both_read1", 99, CHR1, 1350, "30M"), r("both_read1", 83, CHR1, 1400, "30M")],
        [
            r("low_mapq_mate", 99, CHR1, 1450, "30M"),
            r("low_mapq_mate", 147, CHR1, 1500, "30M", mapq=5),
        ],
        [r("duplicate", 1123, CHR1, 100, "30M"), r("duplicate", 1171, CHR1, 250, "30M")],
        [r("qc_failed", 611, CHR1, 1550, "30M"), r("qc_failed", 147, CHR1, 1600, "30M")],
        [r("multimapped", 99, CHR1, 1650, "30M", nh=2), r("multimapped", 147, CHR1, 1700, "30M")],
        [
            r("same_edge", 99, CHR1, 1800, "10M100N20M"),
            r("same_edge", 147, CHR1, 1800, "10M100N20M"),
        ],
        [r("two_edges", 99, CHR1, 2000, "10M50N20M"), r("two_edges", 147, CHR1, 2100, "10M60N20M")],
        [
            r("poly_a", 99, CHR1, 2500, "30M10S", right="A" * 10),
            r("poly_a", 147, CHR1, 2600, "30M"),
        ],
        [
            r("poly_t_rf", 163, CHR1, 2700, "30M"),
            r("poly_t_rf", 83, CHR1, 2800, "10S30M", left="T" * 10),
        ],
        [
            r("no_sequence", 99, CHR1, 2900, "30M5S", sequence=False),
            r("no_sequence", 147, CHR1, 3000, "30M"),
        ],
        [r("alias_pair", 99, ONE, 100, "30M"), r("alias_pair", 147, ONE, 250, "30M")],
        [r("alias_of_excluded", 99, MT, 100, "30M"), r("alias_of_excluded", 147, MT, 200, "30M")],
        [r("excluded", 99, CHRM, 100, "30M"), r("excluded", 147, CHRM, 200, "30M")],
        [
            r("alias_sorts_differently", 163, ZERO, 300, "30M"),
            r("alias_sorts_differently", 83, ZERO, 400, "30M"),
        ],
    ]
    return [record for group in groups for record in group]


def assert_matches_reference(scan, alignment, reference, layout, strandedness, sources, filters):
    assert sorted(scan.evidence) == sorted(sources)
    for source in sources:
        observations, filtering = extract_evidence(
            "S", alignment, reference, layout, strandedness, source, **filters
        )
        assert list(scan.evidence[source].observations()) == observations, source
        # The counter order is published in fragment_filtering.tsv.
        assert list(scan.evidence[source].filtering.items()) == list(filtering.items()), source
    continuations, splice_filtering = extract_splice_continuations(
        "S", alignment, reference, layout, strandedness, **filters
    )
    assert scan.splice is not None
    assert scan.splice[0] == continuations
    assert list(scan.splice[1].items()) == list(splice_filtering.items())


@pytest.mark.parametrize("filters", FILTER_SETS.values(), ids=FILTER_SETS)
@pytest.mark.parametrize("strandedness", ["forward", "reverse"])
def test_single_end_scan_matches_reference_passes(tmp_path, strandedness, filters) -> None:
    reference = write_reference(tmp_path / "genome.fa", CONTIGS)
    bam = write_alignment(tmp_path / "reads.bam", CONTIGS, single_end_records())
    scan = scan_alignment("S", bam, reference, "SE", strandedness, SE_SOURCES, True, **filters)
    assert_matches_reference(scan, bam, reference, "SE", strandedness, SE_SOURCES, filters)


@pytest.mark.parametrize("filters", FILTER_SETS.values(), ids=FILTER_SETS)
@pytest.mark.parametrize("strandedness", ["forward", "reverse"])
def test_paired_end_scan_matches_reference_passes(tmp_path, strandedness, filters) -> None:
    reference = write_reference(tmp_path / "genome.fa", CONTIGS)
    # Prepared alignments are coordinate-sorted, so mates arrive in either order.
    unsorted = write_alignment(tmp_path / "pairs.bam", CONTIGS, paired_end_records())
    bam = sorted_copy(unsorted, tmp_path / "pairs.sorted.bam")
    scan = scan_alignment("S", bam, reference, "PE", strandedness, PE_SOURCES, True, **filters)
    assert_matches_reference(scan, bam, reference, "PE", strandedness, PE_SOURCES, filters)


def test_single_end_scan_values_follow_the_documented_rules(tmp_path) -> None:
    reference = write_reference(tmp_path / "genome.fa", CONTIGS)
    bam = write_alignment(tmp_path / "reads.bam", CONTIGS, single_end_records())
    scan = scan_alignment(
        "S", bam, reference, "SE", "forward", SE_SOURCES, True, **FILTER_SETS["default"]
    )
    read_3p = {
        (item.contig, item.strand, item.coordinate): (item.count, item.poly_a_clip_count)
        for item in scan.evidence["read_3p"].observations()
    }
    # Two chr1 reads and one aliased from "1" end at 130; the duplicate is filtered.
    assert read_3p[("chr1", "+", 130)] == (3, 0)
    # The raw "MT" contig is not excluded; its alias is "chrM". The raw chrM
    # read at the same coordinate is excluded, so the count stays 1.
    assert read_3p[("chrM", "+", 130)] == (1, 0)
    # Interbase ends: 30 aligned bases from 300 end at 330, with a poly(A) clip.
    assert read_3p[("chr1", "+", 330)] == (1, 1)
    contigs = [item.contig for item in scan.evidence["read_3p"].observations()]
    assert contigs[0] == "chr1" and contigs[-1] == "chrM" and "chr3" in contigs
    poly_a = [
        (item.contig, item.strand, item.coordinate)
        for item in scan.evidence["polyA_junction"].observations()
    ]
    # Soft clips consume no reference, and a reverse read's 3' end is its start.
    assert poly_a == [
        ("chr1", "+", 330),
        ("chr1", "+", 530),
        ("chr1", "+", 780),
        ("chr1", "+", 830),
        ("chr1", "-", 400),
        ("chr1", "-", 850),
        ("chr2", "-", 200),
    ]
    filtering = scan.evidence["polyA_junction"].filtering
    # 34 records, 8 filtered: 26 pass, 7 with a poly(A)-like clip.
    assert filtering["no_poly_a_clip"] == 19
    assert filtering["accepted_fragments"] == 7
    for reason in ("low_mapq", "duplicate", "unmapped", "secondary", "supplementary"):
        assert filtering[reason] == 1, reason
    assert filtering["qc_failed"] == 1
    assert filtering["excluded_contig"] == 1
    assert filtering["multimapped"] == 1
    continuations, splice_filtering = scan.splice
    spliced = {
        (item.contig, item.strand, item.upstream_start, item.upstream_end): item.count
        for item in continuations
    }
    # The chr1 read and its "1" alias share one edge, merged under chr1.
    assert spliced[("chr1", "+", 900, 910)] == 2
    # One edge each for three reads, and two for the twice-spliced read.
    assert splice_filtering["splice_direct_edges"] == 5


def test_paired_end_scan_values_follow_the_documented_rules(tmp_path) -> None:
    reference = write_reference(tmp_path / "genome.fa", CONTIGS)
    bam = write_alignment(tmp_path / "pairs.bam", CONTIGS, paired_end_records())
    scan = scan_alignment(
        "S", bam, reference, "PE", "forward", PE_SOURCES, True, **FILTER_SETS["default"]
    )
    fragment = {
        (item.contig, item.strand, item.coordinate): item.count
        for item in scan.evidence["fragment_3p"].observations()
    }
    # Two chr1 pairs and one aliased pair end at 280; mate 2 extends past read 1.
    assert fragment[("chr1", "+", 280)] == 3
    assert fragment[("chr1", "+", 560)] == 1
    # Read 1 is reverse in an RF pair: the transcript runs on "-" and the
    # fragment's 3' end is its leftmost start.
    assert fragment[("chr1", "-", 300)] == 1
    read_3p = {
        (item.contig, item.strand, item.coordinate)
        for item in scan.evidence["read_3p"].observations()
    }
    assert ("chr1", "+", 130) in read_3p and ("chr1", "-", 400) in read_3p
    filtering = scan.evidence["read_3p"].filtering
    assert filtering["orphan_or_multiple_primary"] == 2
    assert filtering["mate_designation"] == 1
    assert filtering["interchromosomal"] == 1
    assert filtering["improper_pair"] == 1
    assert filtering["unmapped"] == 1
    continuations, splice_filtering = scan.splice
    # Both mates cross the same intron; the fragment counts it once.
    same_edge = [item for item in continuations if item.upstream_start == 1800]
    assert [item.count for item in same_edge] == [1]
    assert splice_filtering["splice_direct_edges"] == 3


class _CollatedStream:
    """Stands in for a name-collated alignment, with records in a chosen order."""

    def __init__(self, records: list[pysam.AlignedSegment]) -> None:
        self.records = records

    def fetch(self, until_eof: bool = True):
        return iter(self.records)


def test_paired_end_scan_selects_read1_by_flag_not_order(tmp_path) -> None:
    # samtools collate emits read 1 first, but read 1 must be chosen by its
    # flag, so a collate that emits mate 2 first gives the same evidence.
    bam = write_alignment(tmp_path / "pairs.bam", CONTIGS, paired_end_records())
    with pysam.AlignmentFile(str(bam)) as handle:
        records = list(handle.fetch(until_eof=True))
    groups = [list(group) for _, group in groupby(records, lambda r: r.query_name)]
    read1_first = [r for group in groups for r in sorted(group, key=lambda r: not r.is_read1)]
    mate2_first = [r for group in groups for r in sorted(group, key=lambda r: r.is_read1)]
    results = []
    for records in (read1_first, mate2_first):
        accumulators = [_SourceAccumulator(source) for source in PE_SOURCES]
        splice = _SpliceAccumulator()
        filters = (20, True, True, {"chrM"})
        _scan_paired_end(_CollatedStream(records), "forward", filters, True, accumulators, splice)
        results.append(
            (
                [(item.counts, item.poly_a, list(item.filtering.items())) for item in accumulators],
                splice.counts,
                list(splice.filtering.items()),
            )
        )
    assert results[0] == results[1]
    assert results[0][0][0][0]  # fragment_3p counted fragments


def test_paired_end_scan_collates_once_with_the_reference(tmp_path, monkeypatch) -> None:
    reference = write_reference(tmp_path / "genome.fa", CONTIGS)
    bam = write_alignment(tmp_path / "pairs.bam", CONTIGS, paired_end_records())
    calls = []
    collate = pysam.collate

    def counting_collate(*args, **kwargs):
        calls.append(args)
        return collate(*args, **kwargs)

    monkeypatch.setattr(pysam, "collate", counting_collate)
    scan_alignment("S", bam, reference, "PE", "forward", PE_SOURCES, True)
    assert len(calls) == 1
    assert calls[0][calls[0].index("--reference") + 1] == str(reference)


def test_paired_end_scan_reads_cram_like_bam(tmp_path, monkeypatch) -> None:
    # Keep htslib from looking references up over the network.
    monkeypatch.setenv("REF_PATH", str(tmp_path / "no-reference-cache"))
    reference = write_reference(tmp_path / "genome.fa", CONTIGS)
    bam = write_alignment(tmp_path / "pairs.bam", CONTIGS, paired_end_records())
    cram = to_cram(bam, reference, tmp_path / "pairs.cram")
    from_bam = scan_alignment("S", bam, reference, "PE", "reverse", PE_SOURCES, True)
    from_cram = scan_alignment("S", cram, reference, "PE", "reverse", PE_SOURCES, True)
    for source in PE_SOURCES:
        assert list(from_cram.evidence[source].observations()) == list(
            from_bam.evidence[source].observations()
        )
        assert from_cram.evidence[source].filtering == from_bam.evidence[source].filtering
    assert from_cram.splice == from_bam.splice


def test_scan_rejects_fragment_ends_for_single_end_data(tmp_path) -> None:
    reference = write_reference(tmp_path / "genome.fa", CONTIGS)
    bam = write_alignment(tmp_path / "reads.bam", CONTIGS, single_end_records())
    with pytest.raises(Exception, match="fragment_3p requires paired-end data"):
        scan_alignment("S", bam, reference, "SE", "forward", ["fragment_3p"], False)


def test_written_scan_round_trips(tmp_path) -> None:
    reference = write_reference(tmp_path / "genome.fa", CONTIGS)
    bam = write_alignment(tmp_path / "reads.bam", CONTIGS, single_end_records())
    for splice in (True, False):
        scan = scan_alignment("S", bam, reference, "SE", "forward", SE_SOURCES, splice)
        directory = tmp_path / f"splice-{splice}"
        directory.mkdir()
        manifest_path = write_scan(scan, directory / "S.scan")
        manifest = read_scan_manifest(manifest_path)
        identity = (manifest.sample_id, manifest.layout, manifest.strandedness)
        assert identity == ("S", "SE", "forward")
        for source in SE_SOURCES:
            parquet, filtering = manifest.sources[source]
            assert parquet == directory / f"S.scan.{source}.parquet"
            assert list(filtering.items()) == list(scan.evidence[source].filtering.items())
            assert list(_iter_observation_file(parquet, require_sorted=True)) == list(
                scan.evidence[source].observations()
            )
        splice_file = directory / "S.scan.splice_continuations.tsv.gz"
        if splice:
            assert manifest.splice == (splice_file, scan.splice[1])
        else:
            assert manifest.splice is None
            assert not splice_file.exists()


def test_scan_manifest_reports_missing_files(tmp_path) -> None:
    reference = write_reference(tmp_path / "genome.fa", CONTIGS)
    bam = write_alignment(tmp_path / "reads.bam", CONTIGS, single_end_records())
    scan = scan_alignment("S", bam, reference, "SE", "forward", SE_SOURCES, True)
    manifest = write_scan(scan, tmp_path / "S.scan")
    (tmp_path / "S.scan.read_5p.parquet").unlink()
    with pytest.raises(PacusageError, match="names missing files: S.scan.read_5p.parquet"):
        read_scan_manifest(manifest)


def _cli_inputs(tmp_path: Path, readthrough: bool, endpoint_model: str) -> dict[str, Path]:
    reference = write_reference(tmp_path / "genome.fa", CONTIGS)
    annotation = tmp_path / "genes.gtf"
    annotation.write_text(
        'chr1\ttest\texon\t101\t130\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    params = tmp_path / f"params-{readthrough}.json"
    params.write_text(
        json.dumps(
            {
                "input": str(tmp_path / "samples.tsv"),
                "assembly": "test",
                "fasta": str(reference),
                "gtf": str(annotation),
                "excluded_contigs": ["chrM"],
                "constitutive_readthrough_filter": readthrough,
                "calibration_min_genes": 1,
                "pac_min_sample_count": 1,
            }
        )
    )
    resolution = tmp_path / f"S.{endpoint_model}.resolution.json"
    resolution.write_text(
        json.dumps(
            {
                "sample_id": "S",
                "library_profile": "generic_3prime",
                "layout": "SE",
                "strandedness": "forward",
                "evidence_source": "auto",
                "endpoint_model": endpoint_model,
            }
        )
    )
    ends = tmp_path / "ends.tsv"
    arguments = ["calibration-reference", "--annotation", str(annotation), "--output", str(ends)]
    assert main(arguments) == 0
    bam = write_alignment(tmp_path / "S.bam", CONTIGS, single_end_records())
    return {
        "reference": reference,
        "params": params,
        "resolution": resolution,
        "ends": ends,
        "bam": bam,
    }


def _scan_cli(inputs: dict[str, Path], directory: Path) -> Path:
    directory.mkdir()
    assert (
        main(
            [
                "scan-alignment",
                "--sample-id",
                "S",
                "--alignment",
                str(inputs["bam"]),
                "--resolution",
                str(inputs["resolution"]),
                "--reference",
                str(inputs["reference"]),
                "--transcript-ends",
                str(inputs["ends"]),
                "--params",
                str(inputs["params"]),
                "--calibration",
                str(directory / "S.calibration.json"),
                "--output-prefix",
                str(directory / "S.scan"),
            ]
        )
        == 0
    )
    return directory / "S.scan.filtering.json"


def _run_resolution(path: Path, source: str, endpoint_model: str, strandedness: str) -> Path:
    path.write_text(
        json.dumps(
            {
                "evidence_source": source,
                "endpoint_model": endpoint_model,
                "samples": {"S": {"layout": "SE", "strandedness": strandedness}},
            }
        )
    )
    return path


OUTPUTS = {
    "--tsv": "S.3prime_evidence.tsv.gz",
    "--parquet": "S.3prime_evidence.parquet",
    "--plus-track": "S.plus.3prime_evidence.bedGraph.gz",
    "--minus-track": "S.minus.3prime_evidence.bedGraph.gz",
    "--splice-continuations": "S.splice_continuations.tsv.gz",
    "--qc": "S.fragment_filtering.tsv",
}


def _extract_cli(manifest: Path, run_resolution: Path, params: Path, directory: Path) -> list[str]:
    directory.mkdir()
    arguments = [
        "extract-evidence",
        "--sample-id",
        "S",
        "--scan",
        str(manifest),
        "--resolution",
        str(run_resolution),
        "--params",
        str(params),
    ]
    for option, name in OUTPUTS.items():
        arguments += [option, str(directory / name)]
    return arguments


def _legacy_extract(inputs, source, endpoint_model, readthrough, directory: Path) -> None:
    """The extract-evidence command as it read the alignment before the scan."""
    directory.mkdir()
    filters = {
        "min_mapq": 20,
        "require_unique": True,
        "require_proper_pair": True,
        "exclude_duplicates": True,
        "excluded_contigs": ["chrM"],
        "contig_aliases": {},
    }
    observations, qc = extract_evidence(
        "S", inputs["bam"], inputs["reference"], "SE", "forward", source, **filters
    )
    write_evidence(observations, directory / OUTPUTS["--tsv"], directory / OUTPUTS["--parquet"])
    write_bedgraphs(
        observations, directory / OUTPUTS["--plus-track"], directory / OUTPUTS["--minus-track"]
    )
    if endpoint_model == "proximal_tag" and readthrough:
        continuations, splice_qc = extract_splice_continuations(
            "S", inputs["bam"], inputs["reference"], "SE", "forward", **filters
        )
        splice_qc["splice_filter_enabled"] = True
    else:
        continuations = []
        splice_qc = {
            "splice_filter_enabled": False,
            "splice_records_examined": 0,
            "splice_accepted_fragments": 0,
            "splice_direct_edges": 0,
            "splice_unique_continuations": 0,
        }
    write_splice_continuations(continuations, directory / OUTPUTS["--splice-continuations"])
    qc.update(splice_qc)
    write_tsv([qc], directory / OUTPUTS["--qc"])


@pytest.mark.parametrize("readthrough", [True, False])
@pytest.mark.parametrize("endpoint_model", ["proximal_tag", "exact_boundary"])
@pytest.mark.parametrize("source", SE_SOURCES)
def test_extract_from_scan_writes_the_legacy_outputs(
    tmp_path, source, endpoint_model, readthrough
) -> None:
    inputs = _cli_inputs(tmp_path, readthrough, "auto")
    manifest = _scan_cli(inputs, tmp_path / "scan")
    run_resolution = _run_resolution(tmp_path / "run.json", source, endpoint_model, "forward")
    assert main(_extract_cli(manifest, run_resolution, inputs["params"], tmp_path / "new")) == 0
    _legacy_extract(inputs, source, endpoint_model, readthrough, tmp_path / "legacy")
    for name in OUTPUTS.values():
        new, legacy = tmp_path / "new" / name, tmp_path / "legacy" / name
        assert new.read_bytes() == legacy.read_bytes(), name


def test_extract_rejects_a_scan_with_other_strandedness(tmp_path, capsys) -> None:
    inputs = _cli_inputs(tmp_path, True, "auto")
    manifest = _scan_cli(inputs, tmp_path / "scan")
    run_resolution = _run_resolution(tmp_path / "run.json", "read_3p", "exact_boundary", "reverse")
    with pytest.raises(SystemExit) as error:
        main(_extract_cli(manifest, run_resolution, inputs["params"], tmp_path / "new"))
    assert error.value.code == 2
    assert "was scanned as SE/forward, but the run resolution says SE/reverse" in (
        capsys.readouterr().err
    )


def test_extract_rejects_a_scan_without_needed_splice_continuations(tmp_path, capsys) -> None:
    # An exact endpoint model at scan time skips splice continuations.
    inputs = _cli_inputs(tmp_path, True, "exact_boundary")
    manifest = _scan_cli(inputs, tmp_path / "scan")
    assert read_scan_manifest(manifest).splice is None
    run_resolution = _run_resolution(tmp_path / "run.json", "read_3p", "proximal_tag", "forward")
    with pytest.raises(SystemExit) as error:
        main(_extract_cli(manifest, run_resolution, inputs["params"], tmp_path / "new"))
    assert error.value.code == 2
    assert "did not collect splice continuations" in capsys.readouterr().err
