import pandas as pd
import pysam
from alignment_builders import aligned_segment, write_alignment, write_reference

from pacusage.cli import _iter_observations
from pacusage.evidence import (
    direct_splice_continuations,
    fragment_boundary,
    infer_strandedness,
    is_poly_a_like,
    read_boundary,
    terminal_soft_clip,
    transcript_strand,
    write_evidence,
)
from pacusage.models import EvidenceObservation
from pacusage.reference import GenomicFeature
from pacusage.scan import scan_alignment


def record(flag: int, start: int, cigar: list[tuple[int, int]], sequence: str = "A" * 30):
    value = pysam.AlignedSegment()
    value.query_name = "read"
    value.flag = flag
    value.reference_id = 0
    value.reference_start = start
    value.mapping_quality = 60
    value.cigartuples = cigar
    value.query_sequence = sequence
    return value


def test_read_edges_use_interbase_coordinates() -> None:
    plus = record(0, 100, [(0, 20)])
    minus = record(16, 100, [(0, 20)])
    assert transcript_strand(plus, "forward") == "+"
    assert transcript_strand(minus, "forward") == "-"
    assert read_boundary(plus, "+", "read_3p") == 120
    assert read_boundary(minus, "-", "read_3p") == 100
    assert read_boundary(plus, "+", "read_5p") == 100
    assert read_boundary(minus, "-", "read_5p") == 120


def test_fragment_edges_ignore_introns_and_soft_clips() -> None:
    first = record(65, 100, [(4, 5), (0, 10), (3, 100), (0, 10)])
    second = record(129, 300, [(0, 20), (4, 5)])
    assert fragment_boundary(first, second, "+") == 320
    assert fragment_boundary(first, second, "-") == 100


def test_direct_splice_continuations_are_transcript_oriented() -> None:
    plus = record(0, 100, [(0, 10), (3, 90), (0, 15), (3, 20), (0, 5)])
    minus = record(16, 100, [(0, 10), (3, 90), (0, 15), (3, 20), (0, 5)])
    assert direct_splice_continuations(plus, "+") == [
        (100, 110, 200, 215),
        (200, 215, 235, 240),
    ]
    assert direct_splice_continuations(minus, "-") == [
        (200, 215, 100, 110),
        (235, 240, 200, 215),
    ]


def test_terminal_poly_a_clip_on_both_strands() -> None:
    plus = record(0, 100, [(0, 20), (4, 10)], "C" * 20 + "A" * 10)
    minus = record(16, 100, [(4, 10), (0, 20)], "T" * 10 + "G" * 20)
    assert terminal_soft_clip(plus, "+") == "A" * 10
    assert terminal_soft_clip(minus, "-") == "A" * 10
    assert is_poly_a_like(terminal_soft_clip(minus, "-"))


def test_paired_end_extraction_counts_one_fragment(tmp_path) -> None:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n")
    pysam.faidx(str(fasta))
    bam = tmp_path / "pairs.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    first = record(99, 100, [(0, 30)])
    first.query_name = "pair1"
    first.next_reference_id = 0
    first.next_reference_start = 200
    second = record(147, 200, [(0, 30)])
    second.query_name = "pair1"
    second.next_reference_id = 0
    second.next_reference_start = 100
    with pysam.AlignmentFile(bam, "wb", header=header) as handle:
        handle.write(first)
        handle.write(second)
    scan = scan_alignment(
        "sample",
        bam,
        fasta,
        "PE",
        "forward",
        ["fragment_3p"],
        False,
        min_mapq=0,
    )
    observations = list(scan.evidence["fragment_3p"].observations())
    assert [(item.coordinate, item.count) for item in observations] == [(230, 1)]
    assert scan.evidence["fragment_3p"].filtering["accepted_fragments"] == 1


def test_paired_splice_continuations_deduplicate_mates(tmp_path) -> None:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n")
    pysam.faidx(str(fasta))
    bam = tmp_path / "pairs.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    first = record(99, 100, [(0, 10), (3, 90), (0, 10)], "A" * 20)
    first.query_name = "pair1"
    first.next_reference_id = 0
    first.next_reference_start = 100
    second = record(147, 100, [(0, 10), (3, 90), (0, 10)], "A" * 20)
    second.query_name = "pair1"
    second.next_reference_id = 0
    second.next_reference_start = 100
    with pysam.AlignmentFile(bam, "wb", header=header) as handle:
        handle.write(first)
        handle.write(second)
    scan = scan_alignment(
        "sample",
        bam,
        fasta,
        "PE",
        "forward",
        ["read_3p"],  # Need at least one source
        splice_continuations=True,
        min_mapq=0,
    )
    continuations, qc = scan.splice
    assert [
        (
            item.upstream_start,
            item.upstream_end,
            item.downstream_start,
            item.downstream_end,
            item.count,
        )
        for item in continuations
    ] == [(100, 110, 200, 210, 1)]
    assert qc["splice_direct_edges"] == 1


def test_evidence_parquet_is_written_in_batches(tmp_path) -> None:
    observations = [
        EvidenceObservation("sample", "chr1", "+", 100, 3),
        EvidenceObservation("sample", "chr1", "-", 200, 4, 2),
    ]
    tsv = tmp_path / "evidence.tsv.gz"
    parquet = tmp_path / "evidence.parquet"
    write_evidence(observations, tsv, parquet)
    frame = pd.read_parquet(parquet)
    assert frame[["coordinate", "count"]].to_records(index=False).tolist() == [
        (100, 3),
        (200, 4),
    ]


def test_parquet_evidence_streams_in_global_coordinate_order(tmp_path) -> None:
    first = [
        EvidenceObservation("a", "chr1", "+", 100, 3),
        EvidenceObservation("a", "chr2", "+", 50, 2),
    ]
    second = [
        EvidenceObservation("b", "chr1", "+", 90, 4),
        EvidenceObservation("b", "chr1", "-", 200, 1),
    ]
    paths = []
    for name, observations in (("a", first), ("b", second)):
        tsv = tmp_path / f"{name}.tsv.gz"
        parquet = tmp_path / f"{name}.parquet"
        write_evidence(observations, tsv, parquet)
        paths.append(str(parquet))
    merged = list(_iter_observations(paths, merge_sorted=True))
    assert [
        (row.contig, row.strand, row.coordinate, row.sample_id) for row in merged
    ] == [
        ("chr1", "+", 90, "b"),
        ("chr1", "+", 100, "a"),
        ("chr1", "-", 200, "b"),
        ("chr2", "+", 50, "a"),
    ]


def test_strandedness_inference_applies_contig_aliases(tmp_path) -> None:
    # Reads on contig "1"; the annotation names it "chr1".
    contigs = [("1", 5000)]
    reference = write_reference(tmp_path / "genome.fa", contigs)
    reads = [aligned_segment(f"r{index}", 0, 0, 100 + index, "30M") for index in range(20)]
    bam = write_alignment(tmp_path / "reads.bam", contigs, reads)
    exons = [GenomicFeature("chr1", 50, 400, "+", "exon", "g1", "t1")]
    settings = {
        "minimum_informative": 10,
        "maximum_sampled": 100,
        "decision_fraction": 0.8,
        "random_seed": 1,
    }
    aliased = infer_strandedness(bam, reference, exons, contig_aliases={"1": "chr1"}, **settings)
    assert (aliased["informative_fragments"], aliased["inferred_strandedness"]) == (20, "forward")
    assert infer_strandedness(bam, reference, exons, **settings)["informative_fragments"] == 0


def test_strandedness_inference_skips_excluded_contigs_and_their_aliases(tmp_path) -> None:
    # Ten chr1 reads match their exon's strand. Twenty mitochondrial reads, on
    # "M" (aliased to chrM) and on "MT", oppose theirs, so counting them turns
    # a forward call into an ambiguous one.
    contigs = [("chr1", 5000), ("M", 5000), ("MT", 5000)]
    reference = write_reference(tmp_path / "genome.fa", contigs)
    reads = [
        *[aligned_segment(f"chr1_{index}", 0, 0, 100 + index, "30M") for index in range(10)],
        *[aligned_segment(f"m_{index}", 0, 1, 100 + index, "30M") for index in range(10)],
        *[aligned_segment(f"mt_{index}", 0, 2, 100 + index, "30M") for index in range(10)],
    ]
    bam = write_alignment(tmp_path / "reads.bam", contigs, reads)
    exons = [
        GenomicFeature("chr1", 50, 400, "+", "exon", "g1", "t1"),
        GenomicFeature("chrM", 50, 400, "-", "exon", "g2", "t2"),
        GenomicFeature("MT", 50, 400, "-", "exon", "g3", "t3"),
    ]
    settings = {
        "minimum_informative": 10,
        "maximum_sampled": 100,
        "decision_fraction": 0.8,
        "random_seed": 1,
        "contig_aliases": {"M": "chrM"},
    }
    excluded = infer_strandedness(
        bam, reference, exons, excluded_contigs=["chrM", "MT"], **settings
    )
    assert (excluded["informative_fragments"], excluded["inferred_strandedness"]) == (10, "forward")
    kept = infer_strandedness(bam, reference, exons, excluded_contigs=[], **settings)
    assert (kept["informative_fragments"], kept["reverse_count"]) == (30, 20)
    assert kept["inferred_strandedness"] == "ambiguous"
