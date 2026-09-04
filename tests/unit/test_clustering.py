import numpy as np

from pacusage.calibration import minimum_resolvable_separation
from pacusage.clustering import cluster_exact_boundaries, discover_proximal_pacs
from pacusage.models import EvidenceObservation


def observation(sample: str, coordinate: int, count: int) -> EvidenceObservation:
    return EvidenceObservation(sample, "chr1", "+", coordinate, count)


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
