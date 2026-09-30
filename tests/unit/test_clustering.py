import numpy as np
import pytest

from pacusage.calibration import minimum_resolvable_separation
from pacusage.clustering import (
    _regional_peaks,
    candidates_as_rows,
    cluster_exact_boundaries,
    discover_proximal_pacs,
    filter_constitutive_readthrough,
    filter_internal_exon_ends,
    proximal_pile_shift,
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


def candidate(
    coordinate: int, strand: str = "+", conditions: tuple[str, ...] = ()
) -> PacCandidate:
    return PacCandidate(
        "chr1", strand, coordinate, 10, 2, 6, (coordinate,), supporting_conditions=conditions
    )


# Samples "a" and "b" as replicates of one condition.
ONE_CONDITION = {"a": "C", "b": "C"}
# Two replicates each of a control and a treatment.
TWO_CONDITIONS = {
    "control_1": "control",
    "control_2": "control",
    "treatment_1": "treatment",
    "treatment_2": "treatment",
}
BOTH = ("control", "treatment")


def continuation(sample: str, count: int = 2) -> SpliceContinuation:
    """Reads splicing from the plus-strand exon block [100, 200) to [250, 350)."""
    return SpliceContinuation(sample, "chr1", "+", 100, 200, 250, 350, count)


def readthrough(
    candidates: list[PacCandidate],
    continuations: list[SpliceContinuation],
    minimum: int = 2,
    shift_bins: int = 0,
) -> tuple[list[int], list[int]]:
    accepted, rejected = filter_constitutive_readthrough(
        candidates, continuations, TWO_CONDITIONS, minimum, shift_bins, bin_size=25
    )
    assert all(item.rejection_reason == "constitutive_readthrough" for item in rejected)
    return [item.coordinate for item in accepted], [item.coordinate for item in rejected]


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
        sample_conditions=ONE_CONDITION,
    )
    assert [candidate.coordinate for candidate in accepted] == [100, 120]
    assert accepted[0].member_coordinates == (100, 101)
    assert accepted[0].fraction_within_2nt == 1
    assert accepted[0].width_90 == 1
    assert not rejected


def test_radius_boundary_is_inclusive() -> None:
    observations = [observation("a", 100, 2), observation("b", 112, 2)]
    accepted, _ = cluster_exact_boundaries(
        observations, 0, 12, 1, 1, 1, sample_conditions=ONE_CONDITION
    )
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
        observations, kernel, -2, 0.5, 1, 1, 1, sample_conditions=ONE_CONDITION
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
        sample_conditions=ONE_CONDITION,
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
        discover_proximal_pacs([], np.zeros(3), -1, 0.5, 1, 1, 1, sample_conditions={})


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


def test_constitutive_readthrough_pools_each_supporting_condition() -> None:
    # One continuation read per sample: two per condition, pooled.
    pooled = [continuation(sample, 1) for sample in TWO_CONDITIONS]
    assert readthrough([candidate(150, conditions=BOTH)], pooled) == ([], [150])
    assert readthrough([candidate(150, conditions=BOTH)], pooled, minimum=3) == ([150], [])
    # The next block starts at 250, and 350 is more than a bin past this one.
    assert readthrough([candidate(350, conditions=BOTH)], pooled) == ([350], [])


def test_constitutive_readthrough_asks_only_the_conditions_that_support_a_candidate() -> None:
    treated = [continuation("treatment_1"), continuation("treatment_2")]
    # The control has no continuation reads, but it does not support the
    # candidate, so its low coverage cannot protect it.
    assert readthrough([candidate(150, conditions=("treatment",))], treated) == ([], [150])
    # A supporting condition without continuation keeps the candidate: the
    # exon may end there in that condition.
    assert readthrough([candidate(150, conditions=BOTH)], treated) == ([150], [])
    # A candidate with no supporting condition is never rejected.
    assert readthrough([candidate(150)], treated) == ([150], [])


def test_constitutive_readthrough_looks_at_the_read_pile_behind_a_candidate() -> None:
    reads = [continuation(sample) for sample in TWO_CONDITIONS]
    # 250 is two bins past the block's end, but its reads pile up two bins
    # upstream of it, inside the block.
    assert readthrough([candidate(250, conditions=BOTH)], reads, shift_bins=2) == ([], [250])
    assert readthrough([candidate(250, conditions=BOTH)], reads, shift_bins=0) == ([250], [])


def test_constitutive_readthrough_treats_the_donor_alike_on_both_strands() -> None:
    plus = [continuation(sample) for sample in TWO_CONDITIONS]
    # The block's 3' end is 200 on the plus strand; up to a bin past it counts.
    assert readthrough(
        [candidate(value, conditions=BOTH) for value in (200, 225, 250)], plus
    ) == ([250], [200, 225])
    minus = [
        SpliceContinuation(sample, "chr1", "-", 300, 400, 100, 200, 2)
        for sample in TWO_CONDITIONS
    ]
    # On the minus strand the block's 3' end is its start, 300.
    assert readthrough(
        [candidate(value, "-", BOTH) for value in (250, 275, 300)], minus
    ) == ([250], [275, 300])


def test_constitutive_readthrough_rejects_unknown_samples() -> None:
    with pytest.raises(PacusageError, match="absent from the normalized sample sheet: other"):
        readthrough([candidate(150, conditions=BOTH)], [continuation("other")])


def test_internal_exon_end_filter_expects_the_pile_downstream_of_the_donor() -> None:
    donors = {("chr1", "+"): [1000]}
    candidates = [candidate(value) for value in (950, 975, 1000, 1025, 1050, 1150)]
    accepted, rejected = filter_internal_exon_ends(candidates, donors, 0, 25)
    assert [item.coordinate for item in rejected] == [975, 1000, 1025]
    assert [item.coordinate for item in accepted] == [950, 1050, 1150]
    assert {item.rejection_reason for item in rejected} == {"internal_exon_end"}
    assert {item.status for item in rejected} == {"rejected"}
    # A kernel that peaks six bins downstream moves the expected candidate.
    accepted, rejected = filter_internal_exon_ends(candidates, donors, 6, 25)
    assert [item.coordinate for item in rejected] == [1150]


def test_internal_exon_end_filter_on_the_minus_strand() -> None:
    # A minus-strand donor is an exon start; its pile's peak lies at smaller
    # coordinates.
    donors = {("chr1", "-"): [1000]}
    candidates = [candidate(value, "-") for value in (800, 825, 850, 875, 900, 1000)]
    accepted, rejected = filter_internal_exon_ends(candidates, donors, 6, 25)
    assert [item.coordinate for item in rejected] == [825, 850, 875]
    assert [item.coordinate for item in accepted] == [800, 900, 1000]
    # Donors on the other strand never match.
    plus_only = {("chr1", "+"): [850]}
    assert filter_internal_exon_ends(candidates, plus_only, 0, 25)[1] == []


def test_a_read_pile_peaks_where_the_pile_shift_says() -> None:
    # A kernel peaking 150 nt downstream of the read ends, as in the
    # Plasmidsaurus fixture.
    offsets = np.arange(20, 281)
    weights = np.minimum(offsets - 19, 281 - offsets).astype(float)
    kernel = weights / weights.sum()
    bin_size, shift_bins = proximal_pile_shift(kernel, 20, 0.5, 25)
    assert (bin_size, shift_bins) == (25, 6)
    reads = [observation(sample, 1000, 5) for sample in ("a", "b")]
    reads += [observation(sample, 1000, 5, "-") for sample in ("a", "b")]
    accepted, _, _ = discover_proximal_pacs(
        reads, kernel, 20, 0.5, 1, 1, 1, bin_size=25, sample_conditions=ONE_CONDITION
    )
    assert [(item.strand, item.coordinate) for item in accepted] == [
        ("-", 1000 - shift_bins * bin_size),
        ("+", 1000 + shift_bins * bin_size),
    ]


def test_discovery_records_every_supporting_condition() -> None:
    conditions = {**ONE_CONDITION, "c": "D", "d": "D", "e": "E", "f": "E"}
    reads = [observation(sample, 100, 3) for sample in ("a", "b", "c", "d", "e")]
    accepted, _, _ = discover_proximal_pacs(
        reads, np.asarray([1.0]), 0, 0.5, 1, 2, 2, bin_size=1, sample_conditions=conditions
    )
    # Condition E has one qualifying replicate of two, short of the two needed.
    assert accepted[0].supporting_conditions == ("C", "D")
    exact, _ = cluster_exact_boundaries(
        reads, 2, 2, 1, 2, 2, sample_conditions=conditions
    )
    assert exact[0].supporting_conditions == ("C", "D")
    row = next(candidates_as_rows(accepted))
    assert row["supporting_conditions"] == "C;D"
