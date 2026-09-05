import numpy as np

from pacusage.calibration import (
    CalibrationMetrics,
    calculate_metrics,
    classify_metrics,
    empirical_kernel,
    kernel_correlations,
    kernel_overlap,
    observation_offsets,
)
from pacusage.models import EvidenceObservation
from pacusage.parameters import normalize_parameters


def params() -> dict:
    return normalize_parameters({"input": "x", "assembly": "x", "fasta": "x", "gtf": "x"})


def test_empirical_kernel_is_normalized_and_deterministic() -> None:
    first = empirical_kernel([0, 0, 1, -1], -5, 5)
    second = empirical_kernel([0, 0, 1, -1], -5, 5)
    assert np.allclose(first, second)
    assert np.isclose(first.sum(), 1)
    assert kernel_overlap(first, 0) == 1


def test_exact_and_proximal_classification() -> None:
    values = params()
    values["calibration_min_genes"] = 1
    exact = CalibrationMetrics("s", "read_3p", 2, 20, 0, 0, -1, 1, 0.9, 0, 0.95)
    proximal = CalibrationMetrics("s", "read_3p", 2, 20, 50, 50, 30, 70, 0, 0, 0.95)
    assert classify_metrics(exact, values).classification == "exact"
    assert classify_metrics(proximal, values).classification == "proximal"


def test_leave_one_out_kernel_correlations() -> None:
    kernel = empirical_kernel([0, 0, 1], -3, 3)
    values = kernel_correlations([kernel, kernel.copy(), kernel.copy()])
    assert np.allclose(values, 1)


def test_calibration_offsets_remain_aggregated_for_high_read_counts() -> None:
    observations = [
        EvidenceObservation("sample", "chr1", "+", 100, 1_000_000, 250_000)
    ]
    offsets, genes, clips = observation_offsets(
        observations,
        {("chr1", "+"): [(100, "gene1")]},
        maximum_distance=10,
        minimum_count_per_gene=1,
    )
    assert offsets == {0: 1_000_000}
    assert genes == 1
    assert clips == 250_000
    metrics = calculate_metrics("sample", "read_3p", offsets, genes, clips, 0.05, 0.95)
    assert metrics.observations == 1_000_000
    assert metrics.median_offset == 0
    assert metrics.poly_a_clip_fraction == 0.25
