"""Condition-blind exact-boundary and proximal-tag PAC discovery."""

from __future__ import annotations

import random
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import asdict

import numpy as np

from .calibration import minimum_resolvable_separation
from .models import EvidenceObservation, PacCandidate


def cluster_exact_boundaries(
    observations: Iterable[EvidenceObservation],
    seed_radius: int,
    cluster_radius: int,
    minimum_total: int,
    minimum_sample_count: int,
    minimum_supporting_samples: int,
    known_sites: set[tuple[str, str, int]] | None = None,
    known_rescue_total: int = 5,
    known_match_radius: int = 12,
) -> tuple[list[PacCandidate], list[PacCandidate]]:
    observations = list(observations)
    known_sites = known_sites or set()
    grouped = _group_counts(observations)
    clip_counts = _group_clip_counts(observations)
    accepted: list[PacCandidate] = []
    rejected: list[PacCandidate] = []
    for (contig, strand), coordinate_counts in sorted(grouped.items()):
        coordinates = sorted(coordinate_counts)
        ranked = sorted(
            (
                _seed_metrics(coordinate, coordinates, coordinate_counts, seed_radius)
                for coordinate in coordinates
            ),
            key=lambda value: (
                -value["supporting_samples"],
                -value["capped_support"],
                -value["raw_support"],
                value["coordinate"],
            ),
        )
        assigned: set[int] = set()
        for metrics in ranked:
            coordinate = metrics["coordinate"]
            if coordinate in assigned:
                continue
            members = tuple(
                value
                for value in coordinates
                if value not in assigned and abs(value - coordinate) <= cluster_radius
            )
            assigned.update(members)
            per_sample: defaultdict[str, int] = defaultdict(int)
            for member in members:
                for sample_id, count in coordinate_counts[member].items():
                    per_sample[sample_id] += count
            total = sum(per_sample.values())
            supporting = sum(count >= minimum_sample_count for count in per_sample.values())
            capped = sum(min(count, 3) for count in per_sample.values())
            position_counts = {
                member: sum(coordinate_counts[member].values()) for member in members
            }
            fraction_within_2nt = (
                sum(
                    count
                    for member, count in position_counts.items()
                    if abs(member - coordinate) <= 2
                )
                / total
                if total
                else 0.0
            )
            width_90 = _weighted_interval_width(position_counts)
            clipped = sum(clip_counts[(contig, strand)].get(member, 0) for member in members)
            flank_count = sum(
                sum(sample_counts.values())
                for other, sample_counts in coordinate_counts.items()
                if cluster_radius < abs(other - coordinate) <= cluster_radius * 2
            )
            known = _matches_known(contig, strand, coordinate, known_sites, known_match_radius)
            primary_pass = total >= minimum_total and supporting >= minimum_supporting_samples
            rescue_pass = (
                known and total >= known_rescue_total and supporting >= minimum_supporting_samples
            )
            if primary_pass:
                status = "primary"
                reason = ""
            elif rescue_pass:
                status = "known_rescue_only"
                reason = ""
            else:
                failures = []
                if total < minimum_total:
                    failures.append(f"total_count<{minimum_total}")
                if supporting < minimum_supporting_samples:
                    failures.append(f"supporting_samples<{minimum_supporting_samples}")
                status = "rejected"
                reason = ";".join(failures)
            candidate = PacCandidate(
                contig=contig,
                strand=strand,
                coordinate=coordinate,
                total_count=total,
                supporting_samples=supporting,
                capped_support=capped,
                member_coordinates=members,
                status=status,
                rejection_reason=reason,
                fraction_within_2nt=fraction_within_2nt,
                width_90=width_90,
                poly_a_clip_fraction=clipped / total if total else 0.0,
                local_enrichment=total / max(flank_count, 1),
            )
            (accepted if status != "rejected" else rejected).append(candidate)
    return accepted, rejected


def _group_counts(
    observations: Iterable[EvidenceObservation],
) -> dict[tuple[str, str], dict[int, dict[str, int]]]:
    result: dict[tuple[str, str], dict[int, dict[str, int]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(int))
    )
    for observation in observations:
        result[(observation.contig, observation.strand)][observation.coordinate][
            observation.sample_id
        ] += observation.count
    return result


def _group_clip_counts(
    observations: Iterable[EvidenceObservation],
) -> dict[tuple[str, str], dict[int, int]]:
    result: dict[tuple[str, str], dict[int, int]] = defaultdict(lambda: defaultdict(int))
    for observation in observations:
        result[(observation.contig, observation.strand)][observation.coordinate] += (
            observation.poly_a_clip_count
        )
    return result


def _weighted_interval_width(position_counts: dict[int, int]) -> int:
    expanded = np.asarray(
        [coordinate for coordinate, count in sorted(position_counts.items()) for _ in range(count)],
        dtype=int,
    )
    if not len(expanded):
        return 0
    return int(
        np.quantile(expanded, 0.95, method="nearest")
        - np.quantile(expanded, 0.05, method="nearest")
    )


def _seed_metrics(
    coordinate: int,
    coordinates: list[int],
    coordinate_counts: dict[int, dict[str, int]],
    radius: int,
) -> dict[str, int]:
    local: defaultdict[str, int] = defaultdict(int)
    for other in coordinates:
        if abs(other - coordinate) <= radius:
            for sample_id, count in coordinate_counts[other].items():
                local[sample_id] += count
    return {
        "coordinate": coordinate,
        "supporting_samples": sum(value > 0 for value in local.values()),
        "capped_support": sum(min(value, 3) for value in local.values()),
        "raw_support": sum(local.values()),
    }


def _matches_known(
    contig: str,
    strand: str,
    coordinate: int,
    known_sites: set[tuple[str, str, int]],
    radius: int,
) -> bool:
    return any(
        known_contig == contig
        and known_strand == strand
        and abs(known_coordinate - coordinate) <= radius
        for known_contig, known_strand, known_coordinate in known_sites
    )


def discover_proximal_pacs(
    observations: Iterable[EvidenceObservation],
    kernel: np.ndarray,
    kernel_minimum_offset: int,
    overlap_threshold: float,
    minimum_total: int,
    minimum_sample_count: int,
    minimum_supporting_samples: int,
    bootstrap_replicates: int = 200,
    random_seed: int = 1729,
) -> tuple[list[PacCandidate], list[PacCandidate], int]:
    grouped = _group_counts(observations)
    resolution = minimum_resolvable_separation(kernel, overlap_threshold)
    accepted: list[PacCandidate] = []
    rejected: list[PacCandidate] = []
    offsets = np.flatnonzero(kernel > 0) + kernel_minimum_offset
    for (contig, strand), coordinate_counts in sorted(grouped.items()):
        candidate_coordinates: set[int] = set()
        for endpoint in coordinate_counts:
            for offset in offsets:
                candidate_coordinates.add(
                    endpoint - int(offset) if strand == "+" else endpoint + int(offset)
                )
        scores: dict[int, tuple[float, int, int, int]] = {}
        for coordinate in candidate_coordinates:
            sample_scores: defaultdict[str, float] = defaultdict(float)
            raw_support = 0
            for endpoint, sample_counts in coordinate_counts.items():
                delta = endpoint - coordinate if strand == "+" else coordinate - endpoint
                kernel_index = delta - kernel_minimum_offset
                if kernel_index < 0 or kernel_index >= len(kernel):
                    continue
                weight = float(kernel[kernel_index])
                if weight <= 0:
                    continue
                for sample_id, count in sample_counts.items():
                    sample_scores[sample_id] += min(count, 3) * weight
                    raw_support += count
            supporting = sum(score > 0 for score in sample_scores.values())
            capped = int(round(sum(sample_scores.values()) * 1000000))
            scores[coordinate] = (sum(sample_scores.values()), supporting, capped, raw_support)
        local_maxima = [
            coordinate
            for coordinate in sorted(scores)
            if scores[coordinate][0]
            >= max(
                scores.get(coordinate - 1, (-1, 0, 0, 0))[0],
                scores.get(coordinate + 1, (-1, 0, 0, 0))[0],
            )
        ]
        ranked = sorted(
            local_maxima,
            key=lambda coordinate: (
                -scores[coordinate][1],
                -scores[coordinate][2],
                -scores[coordinate][0],
                -scores[coordinate][3],
                coordinate,
            ),
        )
        assigned_maxima: set[int] = set()
        for coordinate in ranked:
            if coordinate in assigned_maxima:
                continue
            merged = tuple(
                value
                for value in local_maxima
                if value not in assigned_maxima and abs(value - coordinate) < resolution
            )
            assigned_maxima.update(merged)
            per_sample: defaultdict[str, int] = defaultdict(int)
            for endpoint, sample_counts in coordinate_counts.items():
                delta = endpoint - coordinate if strand == "+" else coordinate - endpoint
                kernel_index = delta - kernel_minimum_offset
                if 0 <= kernel_index < len(kernel) and kernel[kernel_index] > 0:
                    for sample_id, count in sample_counts.items():
                        per_sample[sample_id] += count
            total = sum(per_sample.values())
            supporting = sum(count >= minimum_sample_count for count in per_sample.values())
            status = (
                "primary"
                if total >= minimum_total and supporting >= minimum_supporting_samples
                else "rejected"
            )
            failures = []
            if total < minimum_total:
                failures.append(f"total_count<{minimum_total}")
            if supporting < minimum_supporting_samples:
                failures.append(f"supporting_samples<{minimum_supporting_samples}")
            sample_observations: dict[str, list[EvidenceObservation]] = defaultdict(list)
            for endpoint, sample_counts in coordinate_counts.items():
                for sample_id, count in sample_counts.items():
                    sample_observations[sample_id].append(
                        EvidenceObservation(
                            sample_id=sample_id,
                            contig=contig,
                            strand=strand,
                            coordinate=endpoint,
                            count=count,
                        )
                    )
            interval_low, interval_high, successes = bootstrap_proximal_coordinate(
                sample_observations,
                kernel,
                kernel_minimum_offset,
                strand,
                range(coordinate - resolution, coordinate + resolution + 1),
                bootstrap_replicates,
                random_seed + coordinate,
            )
            candidate = PacCandidate(
                contig=contig,
                strand=strand,
                coordinate=coordinate,
                total_count=total,
                supporting_samples=supporting,
                capped_support=sum(min(value, 3) for value in per_sample.values()),
                member_coordinates=merged,
                status=status,
                rejection_reason=";".join(failures) if status == "rejected" else "",
                coordinate_interval_low=interval_low,
                coordinate_interval_high=interval_high,
                coordinate_bootstrap_successes=successes,
            )
            (accepted if status == "primary" else rejected).append(candidate)
    return accepted, rejected, resolution


def bootstrap_proximal_coordinate(
    sample_observations: dict[str, list[EvidenceObservation]],
    kernel: np.ndarray,
    kernel_minimum_offset: int,
    strand: str,
    search_coordinates: Iterable[int],
    replicates: int,
    random_seed: int,
) -> tuple[int | None, int | None, int]:
    sample_ids = sorted(sample_observations)
    coordinates = sorted(set(search_coordinates))
    if not sample_ids or not coordinates:
        return None, None, 0
    rng = random.Random(random_seed)
    maxima: list[int] = []
    for _ in range(replicates):
        sampled = [rng.choice(sample_ids) for _ in sample_ids]
        scores: dict[int, float] = defaultdict(float)
        for sample_id in sampled:
            for observation in sample_observations[sample_id]:
                for coordinate in coordinates:
                    delta = (
                        observation.coordinate - coordinate
                        if strand == "+"
                        else coordinate - observation.coordinate
                    )
                    index = delta - kernel_minimum_offset
                    if 0 <= index < len(kernel):
                        scores[coordinate] += min(observation.count, 3) * kernel[index]
        if scores:
            maxima.append(min(scores, key=lambda value: (-scores[value], value)))
    if len(maxima) < max(1, int(replicates * 0.8)):
        return None, None, len(maxima)
    return (
        int(np.quantile(maxima, 0.025, method="nearest")),
        int(np.quantile(maxima, 0.975, method="nearest")),
        len(maxima),
    )


def candidates_as_rows(candidates: Iterable[PacCandidate]) -> list[dict[str, object]]:
    rows = []
    for candidate in candidates:
        row = asdict(candidate)
        row["member_coordinates"] = ",".join(map(str, candidate.member_coordinates))
        rows.append(row)
    return rows
