import numpy as np

from pacusage.models import EvidenceObservation
from pacusage.quantification import build_count_outputs, quantify_exact, quantify_proximal

ATLAS = [
    {"pac_id": "p1", "contig": "chr1", "strand": "+", "coordinate": 100, "gene_id": "g1"},
    {"pac_id": "p2", "contig": "chr1", "strand": "+", "coordinate": 120, "gene_id": "g1"},
]


def test_exact_assignment_uses_nearest_and_conserves_counts() -> None:
    observations = [
        EvidenceObservation("s", "chr1", "+", 101, 4),
        EvidenceObservation("s", "chr1", "+", 111, 2),
        EvidenceObservation("s", "chr1", "+", 200, 1),
    ]
    rows, qc = quantify_exact(observations, ATLAS, 12)
    assert {row["pac_id"]: row["count"] for row in rows} == {"p1": 4, "p2": 2}
    assert qc["accepted_fragments"] == qc["assigned_fragments"] + qc["unassigned_fragments"]


def test_proximal_ambiguous_assignment_stays_unassigned() -> None:
    kernel = np.zeros(21)
    kernel[0] = 0.5
    kernel[20] = 0.5
    observations = [EvidenceObservation("s", "chr1", "+", 110, 3)]
    _, qc = quantify_proximal(observations, ATLAS, kernel, -10, 3)
    assert qc["ambiguous_fragments"] == 3
    assert qc["assigned_fragments"] == 0


def test_proximal_assignment_uses_transcript_oriented_offsets() -> None:
    atlas = [
        {"pac_id": "plus", "contig": "chr1", "strand": "+", "coordinate": 100},
        {"pac_id": "minus", "contig": "chr1", "strand": "-", "coordinate": 200},
    ]
    observations = [
        EvidenceObservation("s", "chr1", "+", 90, 4),
        EvidenceObservation("s", "chr1", "-", 210, 5),
    ]
    rows, qc = quantify_proximal(observations, atlas, np.asarray([1.0]), 10, 3)
    assert {row["pac_id"]: row["count"] for row in rows} == {
        "plus": 4,
        "minus": 5,
    }
    assert qc["assigned_fragments"] == 9


def test_pau_sums_to_one_without_pseudocounts() -> None:
    per_sample = {
        "a": [{"pac_id": "p1", "count": 3}, {"pac_id": "p2", "count": 1}],
        "b": [{"pac_id": "p1", "count": 0}, {"pac_id": "p2", "count": 0}],
    }
    _, _, totals, pau = build_count_outputs(per_sample, ATLAS)
    assert totals.loc[totals["sample_id"] == "a", "gene_total"].iloc[0] == 4
    assert np.isclose(pau.loc[pau["sample_id"] == "a", "pau"].sum(), 1)
    assert pau.loc[pau["sample_id"] == "b", "pau"].isna().all()
