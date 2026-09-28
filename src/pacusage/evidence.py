"""Strand-aware extraction of aggregate 3-prime evidence from BAM/CRAM."""

from __future__ import annotations

import random
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import asdict
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from .errors import PacusageError
from .models import EvidenceObservation, SpliceContinuation
from .reference import GenomicFeature
from .tableio import open_text, write_tsv

REFERENCE_CONSUMING = {0, 2, 3, 7, 8}


def transcript_strand(record: pysam.AlignedSegment, strandedness: str) -> str:
    alignment_strand = "-" if record.is_reverse else "+"
    if strandedness == "forward":
        return alignment_strand
    if strandedness == "reverse":
        return "+" if alignment_strand == "-" else "-"
    raise PacusageError("Strandedness must be resolved to forward or reverse before extraction.")


def read_boundary(record: pysam.AlignedSegment, strand: str, evidence_source: str) -> int:
    if record.reference_start is None or record.reference_end is None:
        raise PacusageError(f"Read {record.query_name!r} has no aligned reference boundary.")
    if evidence_source in {"read_3p", "polyA_junction"}:
        return record.reference_end if strand == "+" else record.reference_start
    if evidence_source == "read_5p":
        return record.reference_start if strand == "+" else record.reference_end
    raise PacusageError(f"read_boundary does not support {evidence_source!r}.")


def fragment_boundary(
    first: pysam.AlignedSegment, second: pysam.AlignedSegment, strand: str
) -> int:
    starts = [first.reference_start, second.reference_start]
    ends = [first.reference_end, second.reference_end]
    if any(value is None for value in starts + ends):
        raise PacusageError(f"Pair {first.query_name!r} has an undefined aligned boundary.")
    return max(ends) if strand == "+" else min(starts)


def direct_splice_continuations(
    record: pysam.AlignedSegment,
    strand: str,
) -> list[tuple[int, int, int, int]]:
    """Return transcript-oriented exon-block pairs joined by direct CIGAR N operations."""
    if record.reference_start is None or not record.cigartuples:
        return []
    reference_position = record.reference_start
    block_start = reference_position
    blocks: list[tuple[int, int]] = []
    for operation, length in record.cigartuples:
        if operation == 3:
            if block_start < reference_position:
                blocks.append((block_start, reference_position))
            reference_position += length
            block_start = reference_position
        elif operation in REFERENCE_CONSUMING:
            reference_position += length
    if block_start < reference_position:
        blocks.append((block_start, reference_position))
    pairs = zip(blocks, blocks[1:], strict=False)
    if strand == "+":
        return [(*left, *right) for left, right in pairs]
    if strand == "-":
        return [(*right, *left) for left, right in pairs]
    raise PacusageError(f"Splice continuation strand must be '+' or '-', not {strand!r}.")


def terminal_soft_clip(record: pysam.AlignedSegment, strand: str) -> str:
    """Return the soft clip adjacent to the transcript-oriented 3-prime edge."""
    if not record.cigartuples or not record.query_sequence:
        return ""
    if strand == "+":
        operation, length = record.cigartuples[-1]
        return record.query_sequence[-length:].upper() if operation == 4 else ""
    operation, length = record.cigartuples[0]
    if operation != 4:
        return ""
    return reverse_complement(record.query_sequence[:length].upper())


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def is_poly_a_like(sequence: str, minimum_length: int = 6, minimum_fraction: float = 0.8) -> bool:
    return (
        len(sequence) >= minimum_length and sequence.count("A") / len(sequence) >= minimum_fraction
    )


def record_filter_reason(
    record: pysam.AlignedSegment,
    min_mapq: int,
    require_unique: bool,
    exclude_duplicates: bool,
    excluded_reference_ids: set[int] | frozenset[int],
) -> str | None:
    if record.is_unmapped:
        return "unmapped"
    if record.is_secondary:
        return "secondary"
    if record.is_supplementary:
        return "supplementary"
    if record.is_qcfail:
        return "qc_failed"
    if exclude_duplicates and record.is_duplicate:
        return "duplicate"
    if record.mapping_quality < min_mapq:
        return "low_mapq"
    if record.reference_id in excluded_reference_ids:
        return "excluded_contig"
    if require_unique and record.has_tag("NH") and int(record.get_tag("NH")) != 1:
        return "multimapped"
    return None


def pair_filter_reason(
    first: pysam.AlignedSegment,
    second: pysam.AlignedSegment,
    min_mapq: int,
    require_unique: bool,
    require_proper_pair: bool,
    exclude_duplicates: bool,
    excluded_reference_ids: set[int] | frozenset[int],
) -> str | None:
    if first.query_name != second.query_name:
        return "query_name_mismatch"
    if first.is_read1 == second.is_read1 or first.is_read2 == second.is_read2:
        return "mate_designation"
    for record in (first, second):
        reason = record_filter_reason(
            record, min_mapq, require_unique, exclude_duplicates, excluded_reference_ids
        )
        if reason:
            return reason
    if first.reference_id != second.reference_id:
        return "interchromosomal"
    if require_proper_pair and not (first.is_proper_pair and second.is_proper_pair):
        return "improper_pair"
    return None


def query_name_groups(
    records: Iterable[pysam.AlignedSegment],
) -> Iterator[list[pysam.AlignedSegment]]:
    current_name: str | None = None
    current: list[pysam.AlignedSegment] = []
    for record in records:
        if current_name is not None and record.query_name != current_name:
            yield current
            current = []
        current_name = record.query_name
        current.append(record)
    if current:
        yield current


EVIDENCE_FIELDS = [
    "sample_id",
    "contig",
    "strand",
    "coordinate",
    "count",
    "poly_a_clip_count",
    "evidence_source",
]
EVIDENCE_SCHEMA = pa.schema(
    [
        ("sample_id", pa.string()),
        ("contig", pa.string()),
        ("strand", pa.string()),
        ("coordinate", pa.int64()),
        ("count", pa.int64()),
        ("poly_a_clip_count", pa.int64()),
        ("evidence_source", pa.string()),
    ]
)
EVIDENCE_PARQUET_BATCH = 100_000


def write_evidence(
    observations: list[EvidenceObservation],
    tsv_path: str | Path,
    parquet_path: str | Path,
) -> None:
    write_tsv((asdict(observation) for observation in observations), tsv_path, EVIDENCE_FIELDS)
    write_evidence_parquet(observations, parquet_path)


def write_evidence_parquet(
    observations: Iterable[EvidenceObservation],
    path: str | Path,
) -> None:
    """Write evidence in fixed-size row groups, streaming the observations."""
    with pq.ParquetWriter(path, EVIDENCE_SCHEMA) as writer:
        rows: list[dict[str, Any]] = []
        for observation in observations:
            rows.append(asdict(observation))
            if len(rows) == EVIDENCE_PARQUET_BATCH:
                writer.write_table(pa.Table.from_pylist(rows, schema=EVIDENCE_SCHEMA))
                rows = []
        if rows:
            writer.write_table(pa.Table.from_pylist(rows, schema=EVIDENCE_SCHEMA))


def write_splice_continuations(
    continuations: list[SpliceContinuation],
    path: str | Path,
) -> None:
    fields = [
        "sample_id",
        "contig",
        "strand",
        "upstream_start",
        "upstream_end",
        "downstream_start",
        "downstream_end",
        "count",
    ]
    write_tsv((asdict(continuation) for continuation in continuations), path, fields)


def write_bedgraphs(
    observations: list[EvidenceObservation],
    plus_path: str | Path,
    minus_path: str | Path,
) -> None:
    for strand, path in (("+", plus_path), ("-", minus_path)):
        with open_text(path, "wt") as handle:
            for item in observations:
                if item.strand == strand:
                    handle.write(
                        f"{item.contig}\t{item.coordinate}\t{item.coordinate + 1}\t{item.count}\n"
                    )


class ExonBinIndex:
    """Compact interval bins for diagnostic strandedness overlap queries."""

    def __init__(self, features: Iterable[GenomicFeature], bin_size: int = 10000):
        self.bin_size = bin_size
        self.bins: dict[tuple[str, int], list[GenomicFeature]] = defaultdict(list)
        for feature in features:
            if feature.feature_type != "exon":
                continue
            for bin_number in range(feature.start // bin_size, (feature.end - 1) // bin_size + 1):
                self.bins[(feature.contig, bin_number)].append(feature)

    def query(self, contig: str, start: int, end: int) -> list[GenomicFeature]:
        found: dict[tuple[str, str], GenomicFeature] = {}
        for bin_number in range(start // self.bin_size, max(start, end - 1) // self.bin_size + 1):
            for feature in self.bins.get((contig, bin_number), []):
                if feature.start < end and start < feature.end:
                    found[(feature.gene_id, feature.strand)] = feature
        return list(found.values())


def infer_strandedness(
    alignment: str | Path,
    reference: str | Path,
    exons: Iterable[GenomicFeature],
    minimum_informative: int,
    maximum_sampled: int,
    decision_fraction: float,
    random_seed: int,
    contig_aliases: dict[str, str] | None = None,
) -> dict[str, Any]:
    index = ExonBinIndex(exons)
    alignment = Path(alignment)
    mode = "rc" if alignment.suffix.lower() == ".cram" else "rb"
    reservoir: list[tuple[bool, str]] = []
    eligible = 0
    rng = random.Random(random_seed)
    aliases = contig_aliases or {}
    with pysam.AlignmentFile(str(alignment), mode, reference_filename=str(reference)) as handle:
        for record in handle.fetch(until_eof=True):
            if (
                record.is_unmapped
                or record.is_secondary
                or record.is_supplementary
                or (record.is_paired and not record.is_read1)
            ):
                continue
            query_contig = aliases.get(record.reference_name, record.reference_name)
            overlaps = index.query(
                query_contig, record.reference_start, record.reference_end
            )
            genes = {(feature.gene_id, feature.strand) for feature in overlaps}
            if len(genes) != 1:
                continue
            gene_strand = next(iter(genes))[1]
            eligible += 1
            value = (record.is_reverse, gene_strand)
            if len(reservoir) < maximum_sampled:
                reservoir.append(value)
            else:
                replacement = rng.randrange(eligible)
                if replacement < maximum_sampled:
                    reservoir[replacement] = value
    forward = sum(
        1 for reverse, gene_strand in reservoir if ("-" if reverse else "+") == gene_strand
    )
    reverse = len(reservoir) - forward
    total = len(reservoir)
    forward_fraction = forward / total if total else 0.0
    reverse_fraction = reverse / total if total else 0.0
    if total < minimum_informative:
        inferred = "ambiguous"
        reason = f"only {total} informative fragments; need {minimum_informative}"
    elif forward_fraction >= decision_fraction:
        inferred = "forward"
        reason = "forward fraction passed threshold"
    elif reverse_fraction >= decision_fraction:
        inferred = "reverse"
        reason = "reverse fraction passed threshold"
    else:
        inferred = "ambiguous"
        reason = "neither orientation passed the decision threshold"
    return {
        "informative_fragments": total,
        "forward_count": forward,
        "reverse_count": reverse,
        "forward_fraction": forward_fraction,
        "reverse_fraction": reverse_fraction,
        "inferred_strandedness": inferred,
        "reason": reason,
    }
