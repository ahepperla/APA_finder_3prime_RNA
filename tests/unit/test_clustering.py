import numpy as np
import pytest

from pacusage.calibration import minimum_resolvable_separation
from pacusage.clustering import (
    _regional_peaks,
    cluster_exact_boundaries,
    discover_proximal_pacs,
    filter_constitutive_readthrough,
)
from pacusage.errors import PacusageError
from pacusage.models import EvidenceObservation, PacCandidate, SpliceContinuation


def observation(
    sample: str,
    coordinate: int,
    count: int,
    strand: str = "+",
) -> EvidenceObservation:
    return EvidenceObservation(sample, "chr1", strand, coordinate, count)


def candidate(coordinate: int) -> PacCandidate:
    return PacCandidate("chr1", "+", coordinate, 10, 2, 6, (coordinate,))


def continuation(sample: str, count: int = 2) -> SpliceContinuation:
    return SpliceContinuation(sample, "chr1", "+", 100, 200, 250, 350, count)


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


def test_exact_support_must_come_from_one_condition() -> None:
    observations = [
        observation("control_1", 100, 3),
        observation("treatment_1", 100, 3),
    ]
    accepted, rejected = cluster_exact_boundaries(
        observations,
        0,
        12,
        1,
        1,
        2,
        sample_conditions={
            "control_1": "control",
            "treatment_1": "treatment",
        },
    )
    assert not accepted
    assert len(rejected) == 1
    assert rejected[0].supporting_samples == 1
    assert rejected[0].total_supporting_samples == 2


def test_exact_fractional_support_uses_each_condition_size() -> None:
    observations = [
        observation("large_1", 100, 1),
        observation("large_2", 100, 1),
        observation("large_3", 100, 1),
        observation("small_1", 100, 1),
        observation("small_2", 100, 1),
    ]
    sample_conditions = {
        **{f"large_{index}": "large" for index in range(1, 11)},
        "small_1": "small",
        "small_2": "small",
    }
    accepted, rejected = cluster_exact_boundaries(
        observations,
        seed_radius=0,
        cluster_radius=12,
        minimum_total=1,
        minimum_sample_count=1,
        minimum_supporting_samples=0.75,
        sample_conditions=sample_conditions,
    )
    assert not rejected
    assert len(accepted) == 1
    assert accepted[0].supporting_condition == "small"
    assert accepted[0].supporting_samples == 2
    assert accepted[0].total_supporting_samples == 5


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


def test_proximal_support_must_come_from_one_condition() -> None:
    observations = [
        observation("control_1", 90, 3),
        observation("treatment_1", 90, 3),
    ]
    accepted, rejected, _ = discover_proximal_pacs(
        observations,
        np.asarray([1.0]),
        kernel_minimum_offset=10,
        overlap_threshold=0.5,
        minimum_total=1,
        minimum_sample_count=1,
        minimum_supporting_samples=2,
        bin_size=1,
        sample_conditions={
            "control_1": "control",
            "treatment_1": "treatment",
        },
    )
    assert not accepted
    assert len(rejected) == 1
    assert rejected[0].supporting_samples == 1
    assert rejected[0].total_supporting_samples == 2


def test_proximal_fractional_support_rounds_up() -> None:
    sample_conditions = {
        f"sample_{index}": "condition" for index in range(1, 6)
    }
    rejected_observations = [
        observation(f"sample_{index}", 90, 1) for index in range(1, 3)
    ]
    accepted, rejected, _ = discover_proximal_pacs(
        rejected_observations,
        np.asarray([1.0]),
        kernel_minimum_offset=10,
        overlap_threshold=0.5,
        minimum_total=1,
        minimum_sample_count=1,
        minimum_supporting_samples=0.5,
        bin_size=1,
        sample_conditions=sample_conditions,
    )
    assert not accepted
    assert rejected[0].rejection_reason == "supporting_sample_fraction<0.5"

    accepted_observations = [
        observation(f"sample_{index}", 90, 1) for index in range(1, 4)
    ]
    accepted, rejected, _ = discover_proximal_pacs(
        accepted_observations,
        np.asarray([1.0]),
        kernel_minimum_offset=10,
        overlap_threshold=0.5,
        minimum_total=1,
        minimum_sample_count=1,
        minimum_supporting_samples=0.5,
        bin_size=1,
        sample_conditions=sample_conditions,
    )
    assert not rejected
    assert accepted[0].supporting_samples == 3


def test_constitutive_readthrough_rejects_only_upstream_blocks_with_universal_support() -> None:
    conditions = {
        "control_1": "control",
        "control_2": "control",
        "treatment_1": "treatment",
        "treatment_2": "treatment",
    }
    accepted, rejected = filter_constitutive_readthrough(
        [candidate(150), candidate(350)],
        [continuation(sample) for sample in conditions],
        conditions,
        minimum_junction_count=2,
        minimum_replicate_support="all",
    )
    assert [item.coordinate for item in accepted] == [350]
    assert [item.coordinate for item in rejected] == [150]
    assert rejected[0].rejection_reason == "constitutive_readthrough"


def test_constitutive_readthrough_keeps_a_condition_specific_terminal_candidate() -> None:
    conditions = {
        "control_1": "control",
        "control_2": "control",
        "treatment_1": "treatment",
        "treatment_2": "treatment",
    }
    accepted, rejected = filter_constitutive_readthrough(
        [candidate(150)],
        [continuation("control_1"), continuation("control_2")],
        conditions,
        minimum_junction_count=2,
        minimum_replicate_support="all",
    )
    assert [item.coordinate for item in accepted] == [150]
    assert not rejected


def test_constitutive_readthrough_support_can_use_a_per_condition_fraction() -> None:
    conditions = {
        **{f"control_{index}": "control" for index in range(1, 4)},
        **{f"treatment_{index}": "treatment" for index in range(1, 4)},
    }
    continuations = [
        continuation("control_1"),
        continuation("control_2"),
        continuation("treatment_1"),
        continuation("treatment_2"),
    ]
    accepted, rejected = filter_constitutive_readthrough(
        [candidate(150)],
        continuations,
        conditions,
        minimum_junction_count=2,
        minimum_replicate_support=0.5,
    )
    assert not accepted
    assert [item.coordinate for item in rejected] == [150]
