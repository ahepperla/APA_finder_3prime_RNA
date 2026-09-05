"""Assign evidence to a frozen PAC atlas and compute raw counts and PAU."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict
from collections.abc import Iterable

import numpy as np
import pandas as pd

from .models import EvidenceObservation


def quantify_exact(
    observations: Iterable[EvidenceObservation],
    atlas: list[dict[str, object]],
    assignment_radius: int,
) -> tuple[list[dict[str, object]], dict[str, int]]:
    index = _atlas_index(atlas)
    counts: defaultdict[str, int] = defaultdict(int)
    assigned = unassigned = 0
    for observation in observations:
        indexed = index.get((observation.contig, observation.strand))
        if not indexed:
            unassigned += observation.count
            continue
        coordinates, rows = indexed
        lower = bisect_left(coordinates, observation.coordinate - assignment_radius)
        upper = bisect_right(coordinates, observation.coordinate + assignment_radius)
        candidates = rows[lower:upper]
        if not candidates:
            unassigned += observation.count
            continue
        selected = min(
            candidates,
            key=lambda row: (
                abs(int(row["coordinate"]) - observation.coordinate),
                int(row["coordinate"]),
            ),
        )
        counts[str(selected["pac_id"])] += observation.count
        assigned += observation.count
    rows = [{"pac_id": str(row["pac_id"]), "count": counts[str(row["pac_id"])]} for row in atlas]
    return rows, {
        "accepted_fragments": assigned + unassigned,
        "assigned_fragments": assigned,
        "unassigned_fragments": unassigned,
        "ambiguous_fragments": 0,
    }


def quantify_proximal(
    observations: Iterable[EvidenceObservation],
    atlas: list[dict[str, object]],
    kernel: np.ndarray,
    kernel_minimum_offset: int,
    likelihood_ratio: float,
) -> tuple[list[dict[str, object]], dict[str, int]]:
    index = _oriented_atlas_index(atlas)
    nonzero = np.flatnonzero(kernel > 0)
    if not len(nonzero):
        rows = [{"pac_id": str(row["pac_id"]), "count": 0} for row in atlas]
        total = sum(observation.count for observation in observations)
        return rows, {
            "accepted_fragments": total,
            "assigned_fragments": 0,
            "unassigned_fragments": total,
            "ambiguous_fragments": 0,
        }
    minimum_offset = kernel_minimum_offset + int(nonzero[0])
    maximum_offset = kernel_minimum_offset + int(nonzero[-1])
    counts: defaultdict[str, int] = defaultdict(int)
    assigned = unassigned = ambiguous = 0
    for observation in observations:
        indexed = index.get((observation.contig, observation.strand))
        if not indexed:
            unassigned += observation.count
            continue
        coordinates, rows = indexed
        endpoint = (
            observation.coordinate if observation.strand == "+" else -observation.coordinate
        )
        lower = bisect_left(coordinates, endpoint + minimum_offset)
        upper = bisect_right(coordinates, endpoint + maximum_offset)
        likelihoods: list[tuple[float, int, str]] = []
        for coordinate, row in zip(
            coordinates[lower:upper],
            rows[lower:upper],
            strict=True,
        ):
            offset = coordinate - endpoint
            kernel_index = offset - kernel_minimum_offset
            if 0 <= kernel_index < len(kernel) and kernel[kernel_index] > 0:
                likelihoods.append(
                    (
                        float(kernel[kernel_index]),
                        int(row["coordinate"]),
                        str(row["pac_id"]),
                    )
                )
        if not likelihoods:
            unassigned += observation.count
            continue
        likelihoods.sort(key=lambda item: (-item[0], item[1]))
        best = likelihoods[0]
        second = likelihoods[1][0] if len(likelihoods) > 1 else 0.0
        if second > 0 and best[0] < likelihood_ratio * second:
            ambiguous += observation.count
            unassigned += observation.count
            continue
        counts[best[2]] += observation.count
        assigned += observation.count
    rows = [{"pac_id": str(row["pac_id"]), "count": counts[str(row["pac_id"])]} for row in atlas]
    return rows, {
        "accepted_fragments": assigned + unassigned,
        "assigned_fragments": assigned,
        "unassigned_fragments": unassigned,
        "ambiguous_fragments": ambiguous,
    }


def _atlas_index(
    atlas: Iterable[dict[str, object]],
) -> dict[tuple[str, str], tuple[list[int], list[dict[str, object]]]]:
    result: defaultdict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in atlas:
        result[(str(row["contig"]), str(row["strand"]))].append(row)
    indexed = {}
    for key in result:
        result[key].sort(key=lambda row: int(row["coordinate"]))
        indexed[key] = ([int(row["coordinate"]) for row in result[key]], result[key])
    return indexed


def _oriented_atlas_index(
    atlas: Iterable[dict[str, object]],
) -> dict[tuple[str, str], tuple[list[int], list[dict[str, object]]]]:
    result: defaultdict[tuple[str, str], list[tuple[int, dict[str, object]]]] = defaultdict(list)
    for row in atlas:
        strand = str(row["strand"])
        coordinate = int(row["coordinate"])
        oriented = coordinate if strand == "+" else -coordinate
        result[(str(row["contig"]), strand)].append((oriented, row))
    indexed = {}
    for key, values in result.items():
        values.sort(key=lambda item: item[0])
        indexed[key] = (
            [coordinate for coordinate, _ in values],
            [row for _, row in values],
        )
    return indexed


def build_count_outputs(
    per_sample_counts: dict[str, Iterable[dict[str, object]] | pd.Series],
    atlas: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample_ids = sorted(per_sample_counts)
    atlas_rows = sorted(atlas, key=lambda row: str(row["pac_id"]))
    wide = pd.DataFrame(
        {
            "gene_id": [row.get("gene_id", "") for row in atlas_rows],
            "pac_id": [str(row["pac_id"]) for row in atlas_rows],
        }
    )
    for sample_id in sample_ids:
        rows = per_sample_counts[sample_id]
        counts = (
            rows
            if isinstance(rows, pd.Series)
            else pd.Series(
                {str(row["pac_id"]): int(row["count"]) for row in rows},
                dtype="int64",
            )
        )
        counts.index = counts.index.astype(str)
        wide[sample_id] = wide["pac_id"].map(counts).fillna(0).astype("int64")
    long = wide.melt(
        id_vars=["gene_id", "pac_id"],
        value_vars=sample_ids,
        var_name="sample_id",
        value_name="count",
    )
    eligible = long[
        (long["gene_id"].astype(str) != "")
        & ~long["gene_id"].astype(str).str.contains(",", regex=False)
    ].copy()
    totals = (
        eligible.groupby(["gene_id", "sample_id"], as_index=False)["count"]
        .sum()
        .rename(columns={"count": "gene_total"})
    )
    pau = eligible.merge(totals, on=["gene_id", "sample_id"], how="left")
    pau["pau"] = pau["count"] / pau["gene_total"].replace(0, np.nan)
    verify_count_invariants(long, totals, pau)
    return wide, long, totals, pau


def verify_count_invariants(
    long_counts: pd.DataFrame, gene_totals: pd.DataFrame, pau: pd.DataFrame
) -> None:
    expected = (
        long_counts[
            (long_counts["gene_id"].astype(str) != "")
            & ~long_counts["gene_id"].astype(str).str.contains(",", regex=False)
        ]
        .groupby(["gene_id", "sample_id"])["count"]
        .sum()
        .sort_index()
    )
    observed = gene_totals.set_index(["gene_id", "sample_id"])["gene_total"].sort_index()
    if not expected.equals(observed):
        raise ValueError("PAC counts do not sum to reported gene totals.")
    positive = pau[pau["gene_total"] > 0]
    sums = positive.groupby(["gene_id", "sample_id"])["pau"].sum()
    if not np.allclose(sums.to_numpy(), 1.0):
        raise ValueError("PAU does not sum to one for one or more positive gene totals.")
