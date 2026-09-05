"""Condition-blind exact-boundary and proximal-tag PAC discovery."""

from __future__ import annotations

import math
import struct
import tempfile
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import asdict
from itertools import groupby

import numpy as np
from scipy.signal import find_peaks, oaconvolve

from .calibration import minimum_resolvable_separation
from .errors import PacusageError
from .models import EvidenceObservation, PacCandidate

_SPOOLED_OBSERVATION = struct.Struct("<qIq")


def cluster_exact_boundaries(
    observations: Iterable[EvidenceObservation],
    seed_radius: int,
    cluster_radius: int,
    minimum_total: int,
    minimum_sample_count: int,
    minimum_supporting_samples: float,
    known_sites: set[tuple[str, str, int]] | None = None,
    known_rescue_total: int = 5,
    known_match_radius: int = 12,
    observations_sorted: bool = False,
    sample_conditions: Mapping[str, str] | None = None,
) -> tuple[list[PacCandidate], list[PacCandidate]]:
    minimum_supporting_samples = _validate_support_threshold(
        minimum_supporting_samples
    )
    known_sites = known_sites or set()
    ordered: Iterable[EvidenceObservation]
    if observations_sorted:
        ordered = observations
    else:
        ordered = sorted(
            observations,
            key=lambda item: (item.contig, item.strand, item.coordinate, item.sample_id),
        )
    accepted: list[PacCandidate] = []
    rejected: list[PacCandidate] = []
    for (contig, strand), group in groupby(
        ordered,
        key=lambda item: (item.contig, item.strand),
    ):
        coordinate_counts, clip_counts = _group_exact_observations(group)
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
            lower = bisect_left(coordinates, coordinate - cluster_radius)
            upper = bisect_right(coordinates, coordinate + cluster_radius)
            members = tuple(
                value
                for value in coordinates[lower:upper]
                if value not in assigned
            )
            assigned.update(members)
            per_sample: defaultdict[str, int] = defaultdict(int)
            for member in members:
                for sample_id, count in coordinate_counts[member].items():
                    per_sample[sample_id] += count
            total = sum(per_sample.values())
            total_supporting, supporting, supporting_condition, support_pass = (
                _condition_support(
                    per_sample,
                    minimum_sample_count,
                    sample_conditions,
                    minimum_supporting_samples,
                )
            )
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
            clipped = sum(clip_counts.get(member, 0) for member in members)
            flank_count = _flank_count(
                coordinate,
                coordinates,
                coordinate_counts,
                cluster_radius,
            )
            known = _matches_known(contig, strand, coordinate, known_sites, known_match_radius)
            primary_pass = total >= minimum_total and support_pass
            rescue_pass = (
                known and total >= known_rescue_total and support_pass
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
                if not support_pass:
                    failures.append(
                        _support_failure_reason(minimum_supporting_samples)
                    )
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
                region_start=min(members),
                region_end=max(members) + 1,
                resolution_nt=1,
                total_supporting_samples=total_supporting,
                supporting_condition=supporting_condition,
            )
            (accepted if status != "rejected" else rejected).append(candidate)
    return accepted, rejected


def _group_exact_observations(
    observations: Iterable[EvidenceObservation],
) -> tuple[
    dict[int, dict[str, int]],
    dict[int, int],
]:
    counts: dict[int, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    clips: dict[int, int] = defaultdict(int)
    for observation in observations:
        counts[observation.coordinate][observation.sample_id] += observation.count
        clips[observation.coordinate] += observation.poly_a_clip_count
    return counts, clips


def _weighted_interval_width(position_counts: dict[int, int]) -> int:
    total = sum(position_counts.values())
    if not total:
        return 0
    ordered = sorted(position_counts.items())
    return _weighted_coordinate(ordered, total, 0.95) - _weighted_coordinate(
        ordered, total, 0.05
    )


def _weighted_coordinate(
    ordered_counts: list[tuple[int, int]],
    total: int,
    quantile: float,
) -> int:
    target = int(np.rint(quantile * (total - 1)))
    cumulative = 0
    for coordinate, count in ordered_counts:
        cumulative += count
        if target < cumulative:
            return coordinate
    raise ValueError("Weighted coordinate rank exceeds the available observations.")


def _seed_metrics(
    coordinate: int,
    coordinates: list[int],
    coordinate_counts: dict[int, dict[str, int]],
    radius: int,
) -> dict[str, int]:
    local: defaultdict[str, int] = defaultdict(int)
    lower = bisect_left(coordinates, coordinate - radius)
    upper = bisect_right(coordinates, coordinate + radius)
    for other in coordinates[lower:upper]:
        for sample_id, count in coordinate_counts[other].items():
            local[sample_id] += count
    return {
        "coordinate": coordinate,
        "supporting_samples": sum(value > 0 for value in local.values()),
        "capped_support": sum(min(value, 3) for value in local.values()),
        "raw_support": sum(local.values()),
    }


def _flank_count(
    coordinate: int,
    coordinates: list[int],
    coordinate_counts: dict[int, dict[str, int]],
    radius: int,
) -> int:
    if radius <= 0:
        return 0
    left = coordinates[
        bisect_left(coordinates, coordinate - 2 * radius) :
        bisect_left(coordinates, coordinate - radius)
    ]
    right = coordinates[
        bisect_right(coordinates, coordinate + radius) :
        bisect_right(coordinates, coordinate + 2 * radius)
    ]
    return sum(
        sum(coordinate_counts[value].values())
        for value in (*left, *right)
    )


def _condition_support(
    per_sample: Mapping[str, int],
    minimum_sample_count: int,
    sample_conditions: Mapping[str, str] | None,
    minimum_supporting_samples: float,
) -> tuple[int, int, str, bool]:
    qualifying = [
        sample_id
        for sample_id, count in per_sample.items()
        if count >= minimum_sample_count
    ]
    total_supporting = len(qualifying)
    if not sample_conditions:
        if minimum_supporting_samples < 1:
            raise PacusageError(
                "Fractional pac_min_supporting_samples requires sample conditions."
            )
        return (
            total_supporting,
            total_supporting,
            "",
            total_supporting >= int(minimum_supporting_samples),
        )
    missing = sorted(set(per_sample).difference(sample_conditions))
    if missing:
        raise PacusageError(
            "Evidence contains sample IDs absent from the normalized sample sheet: "
            + ", ".join(missing[:10])
        )
    by_condition = Counter(sample_conditions[sample_id] for sample_id in qualifying)
    if not by_condition:
        return 0, 0, "", False
    condition_sizes = Counter(sample_conditions.values())
    fractional = minimum_supporting_samples < 1
    supporting_condition = ""
    supporting = 0
    support_pass = False
    best_score = -1.0
    for condition in sorted(condition_sizes):
        condition_support = by_condition[condition]
        required = _required_supporting_samples(
            minimum_supporting_samples,
            condition_sizes[condition],
        )
        condition_pass = condition_support >= required
        score = (
            condition_support / condition_sizes[condition]
            if fractional
            else float(condition_support)
        )
        ranking = (condition_pass, score, condition_support)
        best_ranking = (support_pass, best_score, supporting)
        if condition_support > 0 and ranking > best_ranking:
            supporting_condition = condition
            supporting = condition_support
            support_pass = condition_pass
            best_score = score
    return total_supporting, supporting, supporting_condition, support_pass


def _validate_support_threshold(value: float) -> float:
    threshold = float(value)
    if not math.isfinite(threshold):
        raise PacusageError("pac_min_supporting_samples must be a finite number.")
    if threshold <= 0:
        raise PacusageError("pac_min_supporting_samples must be greater than zero.")
    if threshold >= 1 and not threshold.is_integer():
        raise PacusageError(
            "pac_min_supporting_samples must be a fraction in (0, 1) "
            "or a whole-number sample count."
        )
    return threshold


def _required_supporting_samples(threshold: float, condition_size: int) -> int:
    if threshold < 1:
        return math.ceil(threshold * condition_size)
    return int(threshold)


def _support_failure_reason(threshold: float) -> str:
    if threshold < 1:
        return f"supporting_sample_fraction<{threshold:g}"
    return f"supporting_samples<{int(threshold)}"


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
    minimum_supporting_samples: float,
    bin_size: int = 25,
    assignment_likelihood_ratio: float = 3.0,
    observations_sorted: bool = False,
    sample_conditions: Mapping[str, str] | None = None,
) -> tuple[list[PacCandidate], list[PacCandidate], int]:
    minimum_supporting_samples = _validate_support_threshold(
        minimum_supporting_samples
    )
    if bin_size < 1:
        raise ValueError("Proximal discovery bin size must be at least one nucleotide.")
    nonzero = np.flatnonzero(kernel > 0)
    if not len(nonzero):
        raise PacusageError("The calibration kernel contains no positive weights.")
    resolution = minimum_resolvable_separation(kernel, overlap_threshold)
    effective_bin_size = min(bin_size, max(1, resolution))
    binned_kernel, kernel_minimum_bin = _bin_kernel(
        kernel,
        kernel_minimum_offset,
        effective_bin_size,
    )
    kernel_minimum_nonzero = kernel_minimum_offset + int(nonzero[0])
    kernel_maximum_nonzero = kernel_minimum_offset + int(nonzero[-1])
    ordered: Iterable[EvidenceObservation]
    if observations_sorted:
        ordered = observations
    else:
        ordered = sorted(
            observations,
            key=lambda item: (item.contig, item.strand, item.coordinate, item.sample_id),
        )
    accepted: list[PacCandidate] = []
    rejected: list[PacCandidate] = []
    for (contig, strand), group in groupby(
        ordered,
        key=lambda item: (item.contig, item.strand),
    ):
        group_accepted, group_rejected = _discover_proximal_group(
            contig,
            strand,
            group,
            kernel,
            kernel_minimum_offset,
            kernel_minimum_nonzero,
            kernel_maximum_nonzero,
            binned_kernel,
            kernel_minimum_bin,
            resolution,
            effective_bin_size,
            assignment_likelihood_ratio,
            minimum_total,
            minimum_sample_count,
            minimum_supporting_samples,
            sample_conditions,
        )
        accepted.extend(group_accepted)
        rejected.extend(group_rejected)
    accepted.sort(key=lambda item: (item.contig, item.coordinate, item.strand))
    rejected.sort(key=lambda item: (item.contig, item.coordinate, item.strand))
    return accepted, rejected, resolution


def _discover_proximal_group(
    contig: str,
    strand: str,
    observations: Iterable[EvidenceObservation],
    kernel: np.ndarray,
    kernel_minimum_offset: int,
    kernel_minimum_nonzero: int,
    kernel_maximum_nonzero: int,
    binned_kernel: np.ndarray,
    kernel_minimum_bin: int,
    resolution: int,
    bin_size: int,
    assignment_likelihood_ratio: float,
    minimum_total: int,
    minimum_sample_count: int,
    minimum_supporting_samples: float,
    sample_conditions: Mapping[str, str] | None,
) -> tuple[list[PacCandidate], list[PacCandidate]]:
    binned_counts: defaultdict[int, int] = defaultdict(int)
    sample_indexes: dict[str, int] = {}
    with tempfile.TemporaryFile() as spool:
        for observation in observations:
            oriented = _oriented_coordinate(observation.coordinate, strand)
            bin_index = _nearest_bin(oriented, bin_size)
            binned_counts[bin_index] += min(observation.count, 3)
            sample_index = sample_indexes.setdefault(
                observation.sample_id,
                len(sample_indexes),
            )
            spool.write(
                _SPOOLED_OBSERVATION.pack(
                    oriented,
                    sample_index,
                    observation.count,
                )
            )

        peaks = _regional_peaks(
            binned_counts,
            binned_kernel,
            kernel_minimum_bin,
            resolution,
            bin_size,
        )
        peaks = [
            peak for peak in peaks if _genomic_coordinate(peak[0], strand) >= 0
        ]
        if not peaks:
            return [], []
        peak_coordinates = [value[0] for value in peaks]
        counts = np.zeros((len(peaks), len(sample_indexes)), dtype=np.int64)
        spool.seek(0)
        while block := spool.read(_SPOOLED_OBSERVATION.size * 65536):
            if len(block) % _SPOOLED_OBSERVATION.size:
                raise PacusageError("Temporary proximal-discovery data are truncated.")
            for endpoint, sample_index, count in _SPOOLED_OBSERVATION.iter_unpack(block):
                selected = _best_proximal_peak(
                    endpoint,
                    peak_coordinates,
                    kernel,
                    kernel_minimum_offset,
                    kernel_minimum_nonzero,
                    kernel_maximum_nonzero,
                    assignment_likelihood_ratio,
                )
                if selected is not None:
                    counts[selected, sample_index] += count

    total_supporting, condition_support, supporting_conditions, support_passes = (
        _condition_support_arrays(
            counts,
            sample_indexes,
            minimum_sample_count,
            sample_conditions,
            minimum_supporting_samples,
        )
    )
    accepted: list[PacCandidate] = []
    rejected: list[PacCandidate] = []
    left_width = resolution // 2
    right_width = resolution - left_width
    for index, (oriented_coordinate, _, members) in enumerate(peaks):
        coordinate = _genomic_coordinate(oriented_coordinate, strand)
        per_sample = counts[index]
        total = int(per_sample.sum())
        supporting = int(condition_support[index])
        capped = int(np.minimum(per_sample, 3).sum())
        failures = []
        if total < minimum_total:
            failures.append(f"total_count<{minimum_total}")
        if not bool(support_passes[index]):
            failures.append(_support_failure_reason(minimum_supporting_samples))
        status = "primary" if not failures else "rejected"
        genomic_members = tuple(
            sorted(_genomic_coordinate(value, strand) for value in members)
        )
        candidate = PacCandidate(
            contig=contig,
            strand=strand,
            coordinate=coordinate,
            total_count=total,
            supporting_samples=supporting,
            capped_support=capped,
            member_coordinates=genomic_members,
            status=status,
            rejection_reason=";".join(failures),
            region_start=max(0, coordinate - left_width),
            region_end=coordinate + right_width,
            resolution_nt=resolution,
            total_supporting_samples=int(total_supporting[index]),
            supporting_condition=str(supporting_conditions[index]),
        )
        (accepted if status == "primary" else rejected).append(candidate)
    return accepted, rejected


def _condition_support_arrays(
    counts: np.ndarray,
    sample_indexes: Mapping[str, int],
    minimum_sample_count: int,
    sample_conditions: Mapping[str, str] | None,
    minimum_supporting_samples: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    qualifying = counts >= minimum_sample_count
    total_supporting = qualifying.sum(axis=1)
    if not sample_conditions:
        if minimum_supporting_samples < 1:
            raise PacusageError(
                "Fractional pac_min_supporting_samples requires sample conditions."
            )
        return (
            total_supporting,
            total_supporting.copy(),
            np.full(len(counts), "", dtype=object),
            total_supporting >= int(minimum_supporting_samples),
        )
    missing = sorted(set(sample_indexes).difference(sample_conditions))
    if missing:
        raise PacusageError(
            "Evidence contains sample IDs absent from the normalized sample sheet: "
            + ", ".join(missing[:10])
        )
    condition_sizes = Counter(sample_conditions.values())
    condition_samples: defaultdict[str, list[int]] = defaultdict(list)
    for sample_id, sample_index in sample_indexes.items():
        condition_samples[sample_conditions[sample_id]].append(sample_index)
    maximum = np.zeros(len(counts), dtype=np.int64)
    winners = np.full(len(counts), "", dtype=object)
    passes = np.zeros(len(counts), dtype=bool)
    best_scores = np.full(len(counts), -1.0)
    fractional = minimum_supporting_samples < 1
    for condition in sorted(condition_sizes):
        indexes = condition_samples[condition]
        support = (
            qualifying[:, indexes].sum(axis=1)
            if indexes
            else np.zeros(len(counts), dtype=np.int64)
        )
        required = _required_supporting_samples(
            minimum_supporting_samples,
            condition_sizes[condition],
        )
        condition_passes = support >= required
        scores = (
            support / condition_sizes[condition]
            if fractional
            else support.astype(float)
        )
        same_pass_state = condition_passes == passes
        replace = (support > 0) & (
            (condition_passes & ~passes)
            | (
                same_pass_state
                & (
                    (scores > best_scores)
                    | ((scores == best_scores) & (support > maximum))
                )
            )
        )
        maximum[replace] = support[replace]
        winners[replace] = condition
        passes[replace] = condition_passes[replace]
        best_scores[replace] = scores[replace]
    return total_supporting, maximum, winners, passes


def _regional_peaks(
    binned_counts: dict[int, int],
    binned_kernel: np.ndarray,
    kernel_minimum_bin: int,
    resolution: int,
    bin_size: int,
) -> list[tuple[int, float, tuple[int, ...]]]:
    if not binned_counts:
        return []
    kernel_span = len(binned_kernel) - 1
    merge_distance_bins = max(0, (resolution - 1) // bin_size)
    selected: list[tuple[int, float, tuple[int, ...]]] = []
    for block in _active_bin_blocks(binned_counts, kernel_span):
        signal_start = block[0][0]
        signal_end = block[-1][0]
        signal = np.zeros(signal_end - signal_start + 1, dtype=float)
        for bin_index, count in block:
            signal[bin_index - signal_start] = count
        score = oaconvolve(signal, binned_kernel, mode="full")
        tolerance = max(float(score.max()) * 1e-12, np.finfo(float).eps)
        score[score < tolerance] = 0
        peak_indexes, _ = find_peaks(score, height=tolerance)
        boundary_peaks = []
        if len(score) == 1 and score[0] > 0:
            boundary_peaks.append(0)
        elif len(score) > 1:
            if score[0] > 0 and score[0] >= score[1]:
                boundary_peaks.append(0)
            if score[-1] > 0 and score[-1] >= score[-2]:
                boundary_peaks.append(len(score) - 1)
        if boundary_peaks:
            peak_indexes = np.unique(
                np.concatenate((peak_indexes, np.asarray(boundary_peaks, dtype=int)))
            )
        if not len(peak_indexes) and score.max() > 0:
            peak_indexes = np.asarray([int(np.argmax(score))])
        ranked = sorted(
            (int(index) for index in peak_indexes),
            key=lambda index: (
                -float(score[index]),
                signal_start + kernel_minimum_bin + index,
            ),
        )
        suppressed = np.zeros(len(score), dtype=bool)
        retained = []
        for index in ranked:
            if suppressed[index]:
                continue
            retained.append(index)
            lower = max(0, index - merge_distance_bins)
            upper = min(len(score), index + merge_distance_bins + 1)
            suppressed[lower:upper] = True
        retained.sort()
        retained_bins = [signal_start + kernel_minimum_bin + index for index in retained]
        member_map: dict[int, list[int]] = {value: [] for value in retained_bins}
        for index in peak_indexes:
            peak_bin = signal_start + kernel_minimum_bin + int(index)
            nearest = _nearest_value(peak_bin, retained_bins)
            if abs(peak_bin - nearest) * bin_size < resolution:
                member_map[nearest].append(peak_bin * bin_size)
        for index, peak_bin in zip(retained, retained_bins, strict=True):
            coordinate = peak_bin * bin_size
            selected.append(
                (
                    coordinate,
                    float(score[index]),
                    tuple(sorted(set(member_map[peak_bin] or [coordinate]))),
                )
            )
    selected.sort(key=lambda value: value[0])
    return selected


def _active_bin_blocks(
    counts: dict[int, int],
    maximum_gap: int,
) -> Iterator[list[tuple[int, int]]]:
    block: list[tuple[int, int]] = []
    previous: int | None = None
    for item in sorted(counts.items()):
        if previous is not None and item[0] - previous > maximum_gap:
            yield block
            block = []
        block.append(item)
        previous = item[0]
    if block:
        yield block


def _bin_kernel(
    kernel: np.ndarray,
    kernel_minimum_offset: int,
    bin_size: int,
) -> tuple[np.ndarray, int]:
    weights: defaultdict[int, float] = defaultdict(float)
    for index, weight in enumerate(kernel):
        if weight <= 0:
            continue
        offset = kernel_minimum_offset + index
        weights[_nearest_bin(offset, bin_size)] += float(weight)
    if not weights:
        raise PacusageError("The calibration kernel contains no positive weights.")
    minimum = min(weights)
    maximum = max(weights)
    result = np.zeros(maximum - minimum + 1, dtype=float)
    for bin_index, weight in weights.items():
        result[bin_index - minimum] = weight
    return result / result.sum(), minimum


def _best_proximal_peak(
    endpoint: int,
    peak_coordinates: list[int],
    kernel: np.ndarray,
    kernel_minimum_offset: int,
    kernel_minimum_nonzero: int,
    kernel_maximum_nonzero: int,
    likelihood_ratio: float,
) -> int | None:
    lower = bisect_left(peak_coordinates, endpoint + kernel_minimum_nonzero)
    upper = bisect_right(peak_coordinates, endpoint + kernel_maximum_nonzero)
    likelihoods = []
    for index in range(lower, upper):
        offset = peak_coordinates[index] - endpoint
        weight = float(kernel[offset - kernel_minimum_offset])
        if weight > 0:
            likelihoods.append((weight, index))
    if not likelihoods:
        return None
    likelihoods.sort(key=lambda value: (-value[0], peak_coordinates[value[1]]))
    best = likelihoods[0]
    second = likelihoods[1][0] if len(likelihoods) > 1 else 0.0
    if second > 0 and best[0] < likelihood_ratio * second:
        return None
    return best[1]


def _nearest_value(value: int, ordered: list[int]) -> int:
    insertion = bisect_left(ordered, value)
    candidates = ordered[max(0, insertion - 1) : min(len(ordered), insertion + 1)]
    return min(candidates, key=lambda candidate: (abs(candidate - value), candidate))


def _nearest_bin(value: int, bin_size: int) -> int:
    magnitude = (abs(value) + bin_size // 2) // bin_size
    return int(magnitude if value >= 0 else -magnitude)


def _oriented_coordinate(coordinate: int, strand: str) -> int:
    return coordinate if strand == "+" else -coordinate


def _genomic_coordinate(coordinate: int, strand: str) -> int:
    return coordinate if strand == "+" else -coordinate


def candidates_as_rows(candidates: Iterable[PacCandidate]) -> list[dict[str, object]]:
    rows = []
    for candidate in candidates:
        row = asdict(candidate)
        row["member_coordinates"] = ",".join(map(str, candidate.member_coordinates))
        rows.append(row)
    return rows
