"""Synthetic references and alignments shared by the scan and calibration tests."""

from __future__ import annotations

from pathlib import Path

import pysam


def write_reference(path: Path, contigs: list[tuple[str, int]]) -> Path:
    with path.open("w") as handle:
        for name, length in contigs:
            handle.write(f">{name}\n{('ACGT' * (length // 4 + 1))[:length]}\n")
    pysam.faidx(str(path))
    return path


def aligned_segment(
    name: str,
    flag: int,
    reference_id: int,
    start: int,
    cigar: str,
    *,
    left: str = "",
    right: str = "",
    sequence: bool = True,
    mapq: int = 60,
    nh: int | None = None,
    mate: tuple[int, int] | None = None,
) -> pysam.AlignedSegment:
    """One record; ``left`` and ``right`` fill its soft clips, the rest is C."""
    record = pysam.AlignedSegment()
    record.query_name = name
    record.flag = flag
    record.reference_id = reference_id
    record.reference_start = start
    record.mapping_quality = mapq
    record.cigarstring = cigar
    if sequence:
        length = record.infer_query_length()
        record.query_sequence = left + "C" * (length - len(left) - len(right)) + right
        record.query_qualities = pysam.qualitystring_to_array("I" * length)
    if nh is not None:
        record.set_tag("NH", nh)
    if mate is not None:
        record.next_reference_id, record.next_reference_start = mate
    return record


def write_alignment(
    path: Path,
    contigs: list[tuple[str, int]],
    records: list[pysam.AlignedSegment],
    sort_order: str = "unsorted",
) -> Path:
    header = {
        "HD": {"VN": "1.6", "SO": sort_order},
        "SQ": [{"SN": name, "LN": length} for name, length in contigs],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for record in records:
            handle.write(record)
    return path


def sorted_copy(source: Path, destination: Path) -> Path:
    pysam.sort("-o", str(destination), str(source), catch_stdout=False)
    return destination


def to_cram(bam: Path, reference: Path, cram: Path) -> Path:
    pysam.view("-C", "-T", str(reference), "-o", str(cram), str(bam), catch_stdout=False)
    return cram
