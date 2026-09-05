"""Calibration metrics and deterministic proximal-tag kernels."""

from __future__ import annotations

from bisect import bisect_left
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import numpy as np

from .errors import PacusageError
from .models import EvidenceObservation


@dataclass(frozen=True)
class CalibrationMetrics:
    sample_id: str
    evidence_source: str
    calibration_genes: int
    observations: int
    median_offset: float
    modal_offset: int
    central_low: float
    central_high: float
    adjacent_boundary_fraction: float
    poly_a_clip_fraction: float
    reproducibility: float = float("nan")
    classification: str = "ambiguous"
    reason: str = ""


def observation_offsets(
    observations: Iterable[EvidenceObservation],
    ends: dict[tuple[str, str], list[tuple[int, str]]],
    maximum_distance: int,
    minimum_count_per_gene: int,
) -> tuple[Counter[int], int, int]:
    offsets: Counter[int] = Counter()
    genes: Counter[str] = Counter()
    matched: dict[str, Counter[int]] = {}
    clipped_by_gene: Counter[str] = Counter()
    end_index = {
        key: ([coordinate for coordinate, _ in values], values)
        for key, values in ends.items()
    }
    for observation in observations:
        indexed = end_index.get((observation.contig, observation.strand))
        if not indexed:
            continue
        coordinates, values = indexed
        insertion = bisect_left(coordinates, observation.coordinate)
        candidate_indexes = []
        if insertion > 0:
            previous_coordinate = coordinates[insertion - 1]
            candidate_indexes.append(bisect_left(coordinates, previous_coordinate))
        if insertion < len(values):
            candidate_indexes.append(insertion)
        candidates = [values[index] for index in candidate_indexes]
        nearest_coordinate, gene_id = min(
            candidates, key=lambda item: abs(item[0] - observation.coordinate)
        )
        distance = abs(nearest_coordinate - observation.coordinate)
        if distance > maximum_distance:
            continue
        offset = (
            nearest_coordinate - observation.coordinate
            if observation.strand == "+"
            else observation.coordinate - nearest_coordinate
        )
        matched.setdefault(gene_id, Counter())[offset] += observation.count
        clipped_by_gene[gene_id] += min(observation.poly_a_clip_count, observation.count)
        genes[gene_id] += observation.count
    eligible = {gene for gene, count in genes.items() if count >= minimum_count_per_gene}
    for gene_id in eligible:
        offsets.update(matched[gene_id])
    clips = sum(clipped_by_gene[gene_id] for gene_id in eligible)
    return offsets, len(eligible), clips


def calculate_metrics(
    sample_id: str,
    evidence_source: str,
    offsets: Mapping[int, int] | Iterable[int],
    calibration_genes: int,
    poly_a_clips: int,
    quantile_low: float,
    quantile_high: float,
) -> CalibrationMetrics:
    counts = _offset_counts(offsets)
    total = sum(counts.values())
    if not total:
        return CalibrationMetrics(
            sample_id=sample_id,
            evidence_source=evidence_source,
            calibration_genes=calibration_genes,
            observations=0,
            median_offset=float("nan"),
            modal_offset=0,
            central_low=float("nan"),
            central_high=float("nan"),
            adjacent_boundary_fraction=0.0,
            poly_a_clip_fraction=0.0,
            reason="no observations near calibration transcript ends",
        )
    mode = min(
        (offset for offset, count in counts.items() if count == max(counts.values())),
        default=0,
    )
    return CalibrationMetrics(
        sample_id=sample_id,
        evidence_source=evidence_source,
        calibration_genes=calibration_genes,
        observations=total,
        median_offset=_weighted_quantile(counts, 0.5),
        modal_offset=mode,
        central_low=_weighted_quantile(counts, quantile_low),
        central_high=_weighted_quantile(counts, quantile_high),
        adjacent_boundary_fraction=(
            sum(count for offset, count in counts.items() if abs(offset) <= 1) / total
        ),
        poly_a_clip_fraction=poly_a_clips / total,
    )


def classify_metrics(metrics: CalibrationMetrics, params: dict) -> CalibrationMetrics:
    values = metrics.__dict__.copy()
    if metrics.calibration_genes < int(params["calibration_min_genes"]):
        values.update(
            classification="ambiguous",
            reason=(
                f"{metrics.calibration_genes} calibration genes; "
                f"need {params['calibration_min_genes']}"
            ),
        )
    elif not np.isfinite(metrics.reproducibility) or metrics.reproducibility < float(
        params["calibration_min_kernel_correlation"]
    ):
        values.update(classification="ambiguous", reason="kernel reproducibility below threshold")
    elif (
        abs(metrics.median_offset) <= float(params["exact_max_median_abs_offset"])
        and metrics.central_high - metrics.central_low <= float(params["exact_max_central_width"])
        and metrics.adjacent_boundary_fraction >= float(params["exact_min_boundary_fraction"])
    ):
        values.update(classification="exact", reason="exact-boundary criteria passed")
    elif (
        metrics.median_offset >= float(params["proximal_min_median_upstream_offset"])
        and metrics.central_low >= -float(params["calibration_max_distance"])
        and metrics.central_high <= float(params["calibration_max_distance"])
    ):
        values.update(classification="proximal", reason="proximal-tag criteria passed")
    else:
        values.update(classification="ambiguous", reason="neither endpoint model passed")
    return CalibrationMetrics(**values)


def empirical_kernel(
    offsets: Mapping[int, int] | Iterable[int], minimum: int, maximum: int
) -> np.ndarray:
    values = np.zeros(maximum - minimum + 1, dtype=float)
    for offset, count in _offset_counts(offsets).items():
        if minimum <= offset <= maximum:
            values[offset - minimum] += count
    if not values.any():
        return values
    smoothed = np.convolve(values, np.asarray([1, 2, 3, 2, 1], dtype=float), mode="same")
    return smoothed / smoothed.sum()


def _offset_counts(offsets: Mapping[int, int] | Iterable[int]) -> Counter[int]:
    if isinstance(offsets, Mapping):
        return Counter(
            {int(offset): int(count) for offset, count in offsets.items() if int(count) > 0}
        )
    return Counter(int(offset) for offset in offsets)


def _weighted_quantile(counts: Mapping[int, int], quantile: float) -> float:
    total = sum(counts.values())
    position = quantile * (total - 1)
    lower_rank = int(np.floor(position))
    upper_rank = int(np.ceil(position))

    def value_at(rank: int) -> int:
        cumulative = 0
        for value, count in sorted(counts.items()):
            cumulative += count
            if rank < cumulative:
                return value
        raise ValueError("Weighted quantile rank exceeds the available observations.")

    lower = value_at(lower_rank)
    upper = value_at(upper_rank)
    return float(lower + (upper - lower) * (position - lower_rank))


def kernel_correlations(kernels: list[np.ndarray]) -> list[float]:
    if len(kernels) < 2:
        return [1.0] * len(kernels)
    correlations: list[float] = []
    for index, kernel in enumerate(kernels):
        others = [item for other_index, item in enumerate(kernels) if other_index != index]
        pooled = np.mean(others, axis=0)
        if np.std(kernel) == 0 or np.std(pooled) == 0:
            correlations.append(0.0)
        else:
            correlations.append(float(np.corrcoef(kernel, pooled)[0, 1]))
    return correlations


def pooled_kernel(kernels: list[np.ndarray]) -> np.ndarray:
    if not kernels:
        raise PacusageError("Cannot pool an empty collection of calibration kernels.")
    pooled = np.mean(kernels, axis=0)
    return pooled / pooled.sum() if pooled.sum() else pooled


def kernel_overlap(kernel: np.ndarray, separation: int) -> float:
    if separation <= 0:
        return 1.0
    shifted = np.zeros_like(kernel)
    shifted[separation:] = kernel[:-separation] if separation < len(kernel) else 0
    return float(np.minimum(kernel, shifted).sum())


def minimum_resolvable_separation(kernel: np.ndarray, threshold: float) -> int:
    for separation in range(1, len(kernel)):
        if kernel_overlap(kernel, separation) <= threshold:
            return separation
    return len(kernel)
