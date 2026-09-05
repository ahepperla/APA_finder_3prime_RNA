import numpy as np
import pytest

from pacusage.calibration import minimum_resolvable_separation
from pacusage.clustering import (
    _regional_peaks,
    cluster_exact_boundaries,
    discover_proximal_pacs,
)
from pacusage.errors import PacusageError
from pacusage.models import EvidenceObservation


def observation(
    sample: str,
    coordinate: int,
    count: int,
    strand: str = "+",
) -> EvidenceObservation:
    return EvidenceObservation(sample, "chr1", strand, coordinate, count)


def test_exact_seed_ranking_is_deterministic() -> None:
    observations = [
        observation("a", 100, 3),
        observation("b", 101, 3),
        observation("a", 120, 9),
    ]
    accepted, rejected = cluster_exact_boundaries(
        observations,
        seed_radius=2,
        cluster_radius=5,
        minimum_total=2,
        minimum_sample_count=2,
        minimum_supporting_samples=1,
    )
    assert [candidate.coordinate for candidate in accepted] == [100, 120]
    assert accepted[0].member_coordinates == (100, 101)
    assert accepted[0].fraction_within_2nt == 1
    assert accepted[0].width_90 == 1
    assert not rejected


def test_radius_boundary_is_inclusive() -> None:
    observations = [observation("a", 100, 2), observation("b", 112, 2)]
    accepted, _ = cluster_exact_boundaries(observations, 0, 12, 1, 1, 1)
    assert len(accepted) == 1
    assert accepted[0].member_coordinates == (100, 112)


def test_proximal_resolution_merges_unresolved_maxima() -> None:
    kernel = np.array([0.1, 0.2, 0.4, 0.2, 0.1])
    resolution = minimum_resolvable_separation(kernel, 0.5)
    observations = [
        observation("a", 100, 3),
        observation("b", 101, 3),
    ]
    accepted, _, observed_resolution = discover_proximal_pacs(
        observations, kernel, -2, 0.5, 1, 1, 1
    )
    assert observed_resolution == resolution
    assert accepted
    assert len({value for candidate in accepted for value in candidate.member_coordinates}) >= 1


@pytest.mark.parametrize(
    ("strand", "endpoint", "expected_coordinate"),
    [("+", 90, 100), ("-", 110, 100)],
)
def test_proximal_positive_offset_points_toward_the_pac(
    strand: str,
    endpoint: int,
    expected_coordinate: int,
) -> None:
    observations = [
        observation("a", endpoint, 3, strand),
        observation("b", endpoint, 3, strand),
    ]
    accepted, rejected, resolution = discover_proximal_pacs(
        observations,
        np.asarray([1.0]),
        kernel_minimum_offset=10,
        overlap_threshold=0.5,
        minimum_total=1,
        minimum_sample_count=1,
        minimum_supporting_samples=1,
        bin_size=1,
    )
    assert not rejected
    assert resolution == 1
    assert [candidate.coordinate for candidate in accepted] == [expected_coordinate]
    assert accepted[0].region_start == expected_coordinate
    assert accepted[0].region_end == expected_coordinate + 1


def test_regional_peaks_merge_only_inside_calibrated_resolution() -> None:
    kernel = np.zeros(41)
    kernel[20] = 1
    peaks = _regional_peaks(
        {0: 3, 10: 3, 40: 3},
        kernel,
        kernel_minimum_bin=-20,
        resolution=315,
        bin_size=25,
    )
    assert [peak[0] for peak in peaks] == [0, 1000]
    assert peaks[0][2] == (0, 250)


def test_proximal_discovery_rejects_an_empty_kernel() -> None:
    with pytest.raises(PacusageError, match="no positive weights"):
        discover_proximal_pacs([], np.zeros(3), -1, 0.5, 1, 1, 1)
