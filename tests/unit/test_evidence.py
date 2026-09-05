import pandas as pd
import pysam

from pacusage.cli import _iter_observations
from pacusage.evidence import (
    extract_evidence,
    fragment_boundary,
    is_poly_a_like,
    read_boundary,
    terminal_soft_clip,
    transcript_strand,
    write_evidence,
)
from pacusage.models import EvidenceObservation


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
    observations, qc = extract_evidence(
        "sample",
        bam,
        fasta,
        "PE",
        "forward",
        "fragment_3p",
        min_mapq=0,
    )
    assert [(item.coordinate, item.count) for item in observations] == [(230, 1)]
    assert qc["accepted_fragments"] == 1


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
