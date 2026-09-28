"""One pass over an alignment for every candidate evidence source.

``scan_alignment`` reads an alignment once, with at most one name-collate for
paired-end data, collecting every candidate evidence source, direct splice
continuations, and per-source filtering counters. Counters keep first-seen
order because ``fragment_filtering.tsv`` publishes them in that order.
"""

from __future__ import annotations

import json
import tempfile
from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam

from .errors import PacusageError
from .evidence import (
    direct_splice_continuations,
    fragment_boundary,
    is_poly_a_like,
    pair_filter_reason,
    query_name_groups,
    read_boundary,
    record_filter_reason,
    terminal_soft_clip,
    transcript_strand,
    write_evidence_parquet,
    write_splice_continuations,
)
from .models import EvidenceObservation, SpliceContinuation
from .tableio import write_json

EVIDENCE_SOURCES = ("fragment_3p", "polyA_junction", "read_3p", "read_5p")


@dataclass
class SourceEvidence:
    """One source's evidence, keyed by (contig, strand) until it is written.

    Each group holds coordinate -> fragments, and coordinate -> fragments with
    a poly(A)-like terminal clip, where only nonzero clip counts are stored.
    """

    sample_id: str
    source: str
    groups: dict[tuple[str, str], tuple[dict[int, int], dict[int, int]]]
    filtering: dict[str, Any]

    def observations(self) -> Iterator[EvidenceObservation]:
        """Yield observations sorted by contig, strand, and coordinate."""
        for contig, strand in sorted(self.groups):
            counts, poly_a = self.groups[(contig, strand)]
            for coordinate in sorted(counts):
                yield EvidenceObservation(
                    sample_id=self.sample_id,
                    contig=contig,
                    strand=strand,
                    coordinate=coordinate,
                    count=counts[coordinate],
                    poly_a_clip_count=poly_a.get(coordinate, 0),
                    evidence_source=self.source,
                )


@dataclass
class AlignmentScan:
    sample_id: str
    layout: str
    strandedness: str
    evidence: dict[str, SourceEvidence]
    splice: tuple[list[SpliceContinuation], dict[str, Any]] | None


@dataclass(frozen=True)
class ScanManifest:
    """A written scan: per-source parquet files and ordered counters."""

    sample_id: str
    layout: str
    strandedness: str
    sources: dict[str, tuple[Path, dict[str, Any]]]
    splice: tuple[Path, dict[str, Any]] | None


class _SourceAccumulator:
    __slots__ = ("source", "counts", "poly_a", "filtering")

    def __init__(self, source: str) -> None:
        self.source = source
        # Keyed by (reference_id, strand); aliases are applied when finished.
        self.counts: dict[tuple[int, str], dict[int, int]] = {}
        self.poly_a: dict[tuple[int, str], dict[int, int]] = {}
        self.filtering: Counter[str] = Counter()

    def add(self, key: tuple[int, str], coordinate: int, poly_a_like: bool) -> None:
        bucket = self.counts.get(key)
        if bucket is None:
            bucket = self.counts[key] = {}
        bucket[coordinate] = bucket.get(coordinate, 0) + 1
        if poly_a_like:
            clipped = self.poly_a.get(key)
            if clipped is None:
                clipped = self.poly_a[key] = {}
            clipped[coordinate] = clipped.get(coordinate, 0) + 1
        self.filtering["accepted_fragments"] += 1


class _SpliceAccumulator:
    __slots__ = ("counts", "filtering")

    def __init__(self) -> None:
        self.counts: dict[tuple[int, str, int, int, int, int], int] = {}
        self.filtering: Counter[str] = Counter()

    def add(self, key: tuple[int, str], edges: Iterable[tuple[int, int, int, int]]) -> None:
        edges = list(edges)
        for edge in edges:
            full = (*key, *edge)
            self.counts[full] = self.counts.get(full, 0) + 1
        self.filtering["splice_accepted_fragments"] += 1
        self.filtering["splice_direct_edges"] += len(edges)


def scan_alignment(
    sample_id: str,
    alignment: str | Path,
    reference: str | Path,
    layout: str,
    strandedness: str,
    sources: Iterable[str],
    splice_continuations: bool,
    min_mapq: int = 20,
    require_unique: bool = True,
    require_proper_pair: bool = True,
    exclude_duplicates: bool = True,
    excluded_contigs: Iterable[str] = (),
    contig_aliases: dict[str, str] | None = None,
    threads: int = 1,
) -> AlignmentScan:
    if layout not in {"SE", "PE"}:
        raise PacusageError(f"Sample {sample_id}: layout must resolve to SE or PE, not {layout!r}.")
    if strandedness not in {"forward", "reverse"}:
        raise PacusageError(f"Sample {sample_id}: strandedness must resolve to forward or reverse.")
    selected = sorted(set(sources))
    if not selected:
        raise PacusageError(f"Sample {sample_id}: no evidence source was requested.")
    unsupported = sorted(set(selected).difference(EVIDENCE_SOURCES))
    if unsupported:
        raise PacusageError(f"Sample {sample_id}: unsupported evidence source {unsupported[0]!r}.")
    if layout == "SE" and "fragment_3p" in selected:
        raise PacusageError(f"Sample {sample_id}: fragment_3p requires paired-end data.")

    accumulators = [_SourceAccumulator(source) for source in selected]
    splice = _SpliceAccumulator() if splice_continuations else None
    alignment = Path(alignment)
    mode = "rc" if alignment.suffix.lower() == ".cram" else "rb"
    threads = max(1, threads)
    aliases = contig_aliases or {}
    excluded = set(excluded_contigs)

    def excluded_ids(names: tuple[str, ...]) -> frozenset[int]:
        # An excluded contig may be named as in the alignment or by its alias.
        return frozenset(
            index
            for index, name in enumerate(names)
            if name in excluded or aliases.get(name, name) in excluded
        )

    if layout == "SE":
        with pysam.AlignmentFile(
            str(alignment), mode, reference_filename=str(reference), threads=threads
        ) as handle:
            names = handle.references
            filters = (min_mapq, require_unique, exclude_duplicates, excluded_ids(names))
            _scan_single_end(handle, strandedness, filters, accumulators, splice)
    else:
        with tempfile.TemporaryDirectory(prefix="pacusage-collate-") as temporary:
            collated = Path(temporary) / "collated.bam"
            try:
                pysam.collate(
                    "-@",
                    str(threads),
                    "--reference",
                    str(reference),
                    "-o",
                    str(collated),
                    str(alignment),
                    catch_stdout=False,
                )
            except Exception as error:
                raise PacusageError(f"Could not name-collate {alignment}.") from error
            with pysam.AlignmentFile(str(collated), "rb", threads=threads) as handle:
                names = handle.references
                filters = (min_mapq, require_unique, exclude_duplicates, excluded_ids(names))
                _scan_paired_end(
                    handle, strandedness, filters, require_proper_pair, accumulators, splice
                )

    evidence = {
        accumulator.source: _finish_source(
            accumulator, sample_id, layout, strandedness, names, aliases
        )
        for accumulator in accumulators
    }
    finished_splice = (
        _finish_splice(splice, sample_id, layout, strandedness, names, aliases)
        if splice is not None
        else None
    )
    return AlignmentScan(sample_id, layout, strandedness, evidence, finished_splice)


def _scan_single_end(
    handle: pysam.AlignmentFile,
    strandedness: str,
    filters: tuple[int, bool, bool, set[str]],
    accumulators: list[_SourceAccumulator],
    splice: _SpliceAccumulator | None,
) -> None:
    min_mapq, require_unique, exclude_duplicates, excluded = filters
    for record in handle.fetch(until_eof=True):
        for accumulator in accumulators:
            accumulator.filtering["records_examined"] += 1
        if splice is not None:
            splice.filtering["splice_records_examined"] += 1
        reason = record_filter_reason(
            record, min_mapq, require_unique, exclude_duplicates, excluded
        )
        if reason:
            for accumulator in accumulators:
                accumulator.filtering[reason] += 1
            if splice is not None:
                splice.filtering[f"splice_{reason}"] += 1
            continue
        strand = transcript_strand(record, strandedness)
        key = (record.reference_id, strand)
        if splice is not None:
            splice.add(key, direct_splice_continuations(record, strand))
        clip = terminal_soft_clip(record, strand)
        poly_a_like = is_poly_a_like(clip)
        for accumulator in accumulators:
            if accumulator.source == "polyA_junction" and not poly_a_like:
                accumulator.filtering["no_poly_a_clip"] += 1
                continue
            coordinate = read_boundary(record, strand, accumulator.source)
            accumulator.add(key, coordinate, poly_a_like)


def _scan_paired_end(
    handle: pysam.AlignmentFile,
    strandedness: str,
    filters: tuple[int, bool, bool, set[str]],
    require_proper_pair: bool,
    accumulators: list[_SourceAccumulator],
    splice: _SpliceAccumulator | None,
) -> None:
    min_mapq, require_unique, exclude_duplicates, excluded = filters
    for group in query_name_groups(handle.fetch(until_eof=True)):
        for accumulator in accumulators:
            accumulator.filtering["query_groups_examined"] += 1
        if splice is not None:
            splice.filtering["splice_query_groups_examined"] += 1
        primary = [
            record
            for record in group
            if not record.is_secondary and not record.is_supplementary
        ]
        if len(primary) != 2:
            for accumulator in accumulators:
                accumulator.filtering["orphan_or_multiple_primary"] += 1
            if splice is not None:
                splice.filtering["splice_orphan_or_multiple_primary"] += 1
            continue
        first, second = primary
        reason = pair_filter_reason(
            first,
            second,
            min_mapq,
            require_unique,
            require_proper_pair,
            exclude_duplicates,
            excluded,
        )
        if reason:
            for accumulator in accumulators:
                accumulator.filtering[reason] += 1
            if splice is not None:
                splice.filtering[f"splice_{reason}"] += 1
            continue
        read1 = first if first.is_read1 else second
        strand = transcript_strand(read1, strandedness)
        key = (read1.reference_id, strand)
        if splice is not None:
            splice.add(
                key,
                {
                    edge
                    for record in (first, second)
                    for edge in direct_splice_continuations(record, strand)
                },
            )
        clip = terminal_soft_clip(read1, strand)
        poly_a_like = is_poly_a_like(clip)
        for accumulator in accumulators:
            if accumulator.source == "polyA_junction" and not poly_a_like:
                accumulator.filtering["no_poly_a_clip"] += 1
                continue
            coordinate = (
                fragment_boundary(first, second, strand)
                if accumulator.source == "fragment_3p"
                else read_boundary(read1, strand, accumulator.source)
            )
            accumulator.add(key, coordinate, poly_a_like)


def _finish_source(
    accumulator: _SourceAccumulator,
    sample_id: str,
    layout: str,
    strandedness: str,
    names: tuple[str, ...],
    aliases: dict[str, str],
) -> SourceEvidence:
    groups: dict[tuple[str, str], tuple[dict[int, int], dict[int, int]]] = {}
    for (reference_id, strand), counts in accumulator.counts.items():
        name = names[reference_id]
        group = (aliases.get(name, name), strand)
        poly_a = accumulator.poly_a.get((reference_id, strand), {})
        if group not in groups:
            groups[group] = (counts, poly_a)
            continue
        # Two alignment contigs alias to one name: merge them.
        merged_counts, merged_poly_a = groups[group]
        for coordinate, value in counts.items():
            merged_counts[coordinate] = merged_counts.get(coordinate, 0) + value
        for coordinate, value in poly_a.items():
            merged_poly_a[coordinate] = merged_poly_a.get(coordinate, 0) + value
    filtering = accumulator.filtering
    filtering["unique_observations"] = sum(len(counts) for counts, _ in groups.values())
    filtering["sample_id"] = sample_id
    filtering["layout"] = layout
    filtering["strandedness"] = strandedness
    filtering["evidence_source"] = accumulator.source
    return SourceEvidence(sample_id, accumulator.source, groups, dict(filtering))


def _finish_splice(
    accumulator: _SpliceAccumulator,
    sample_id: str,
    layout: str,
    strandedness: str,
    names: tuple[str, ...],
    aliases: dict[str, str],
) -> tuple[list[SpliceContinuation], dict[str, Any]]:
    merged: dict[tuple[str, str, int, int, int, int], int] = {}
    for (reference_id, strand, *edge), count in accumulator.counts.items():
        name = names[reference_id]
        key = (aliases.get(name, name), strand, *edge)
        merged[key] = merged.get(key, 0) + count
    continuations = [
        SpliceContinuation(
            sample_id=sample_id,
            contig=contig,
            strand=strand,
            upstream_start=upstream_start,
            upstream_end=upstream_end,
            downstream_start=downstream_start,
            downstream_end=downstream_end,
            count=count,
        )
        for (
            contig,
            strand,
            upstream_start,
            upstream_end,
            downstream_start,
            downstream_end,
        ), count in sorted(merged.items())
    ]
    filtering = accumulator.filtering
    filtering["splice_unique_continuations"] = len(continuations)
    filtering["sample_id"] = sample_id
    filtering["splice_layout"] = layout
    filtering["splice_strandedness"] = strandedness
    return continuations, dict(filtering)


def write_scan(scan: AlignmentScan, prefix: str | Path) -> Path:
    """Write the scan next to ``prefix`` and return its manifest path.

    Counters are stored as ordered key-value pairs because ``write_json``
    sorts mapping keys.
    """
    prefix = Path(prefix)
    sources: dict[str, dict[str, Any]] = {}
    for source, evidence in scan.evidence.items():
        parquet = prefix.with_name(f"{prefix.name}.{source}.parquet")
        write_evidence_parquet(evidence.observations(), parquet)
        sources[source] = {
            "parquet": parquet.name,
            "filtering": [[key, value] for key, value in evidence.filtering.items()],
        }
    splice: dict[str, Any] | None = None
    if scan.splice is not None:
        continuations, filtering = scan.splice
        path = prefix.with_name(f"{prefix.name}.splice_continuations.tsv.gz")
        write_splice_continuations(continuations, path)
        splice = {
            "continuations": path.name,
            "filtering": [[key, value] for key, value in filtering.items()],
        }
    manifest = prefix.with_name(f"{prefix.name}.filtering.json")
    write_json(
        {
            "sample_id": scan.sample_id,
            "layout": scan.layout,
            "strandedness": scan.strandedness,
            "sources": sources,
            "splice": splice,
        },
        manifest,
    )
    return manifest


def read_scan_manifest(path: str | Path) -> ScanManifest:
    path = Path(path)
    try:
        payload = json.loads(path.read_text())
    except (OSError, ValueError) as error:
        raise PacusageError(f"Cannot read alignment scan manifest {path}.") from error
    sources = {
        source: (path.parent / value["parquet"], dict(value["filtering"]))
        for source, value in payload["sources"].items()
    }
    splice_payload = payload.get("splice")
    splice = (
        None
        if splice_payload is None
        else (path.parent / splice_payload["continuations"], dict(splice_payload["filtering"]))
    )
    named = [parquet for parquet, _ in sources.values()] + ([splice[0]] if splice else [])
    missing = [item.name for item in named if not item.is_file()]
    if missing:
        raise PacusageError(
            f"Alignment scan {path} names missing files: {', '.join(missing)}. "
            "Rerun SCAN_ALIGNMENT."
        )
    return ScanManifest(
        sample_id=payload["sample_id"],
        layout=payload["layout"],
        strandedness=payload["strandedness"],
        sources=sources,
        splice=splice,
    )
