import numpy as np

from pacusage.calibration import (
    CalibrationMetrics,
    classify_metrics,
    empirical_kernel,
    kernel_correlations,
    kernel_overlap,
)
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
