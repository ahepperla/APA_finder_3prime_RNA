import numpy as np
import pytest

from pacusage.models import EvidenceObservation
from pacusage.quantification import build_count_outputs, quantify_exact, quantify_proximal


def atlas_row(pac_id: str, strand: str, coordinate: int, gene_id: str, gene_name: str) -> dict:
    return {
        "pac_id": pac_id, "gene_id": gene_id, "gene_name": gene_name, "chrom": "chr1",
        "start": coordinate, "end": coordinate + 1, "strand": strand,
        "locus": f"chr1:{coordinate + 1}-{coordinate + 1}", "coordinate": coordinate,
    }


ATLAS = [
    atlas_row("p1", "+", 100, "g1", "GeneOne"),
    atlas_row("p2", "+", 120, "g1", "GeneOne"),
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
    atlas = [atlas_row("plus", "+", 100, "", ""), atlas_row("minus", "-", 200, "", "")]
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


def test_count_tables_carry_gene_names_after_gene_ids() -> None:
    atlas = [*ATLAS, atlas_row("p3", "+", 900, "", "")]
    counts = {"p1": 3, "p2": 1, "p3": 2}
    per_sample = {"a": [{"pac_id": pac_id, "count": count} for pac_id, count in counts.items()]}
    wide, long, totals, pau = build_count_outputs(per_sample, atlas)
    identity = ["pac_id", "gene_id", "gene_name", "chrom", "start", "end", "strand", "locus"]
    assert list(wide.columns) == [*identity, "a"]
    assert list(long.columns) == ["pac_id", "gene_id", "gene_name", "sample_id", "count"]
    assert list(totals.columns) == ["gene_id", "gene_name", "sample_id", "gene_total"]
    assert list(pau.columns) == [
        "pac_id", "gene_id", "gene_name", "sample_id", "count", "gene_total", "pau",
    ]
    # Rows keep the atlas's order, not pac_id string order.
    assert wide["pac_id"].tolist() == ["p1", "p2", "p3"]
    assert wide["locus"].tolist() == ["chr1:101-101", "chr1:121-121", "chr1:901-901"]
    assert wide["gene_name"].tolist() == ["GeneOne", "GeneOne", ""]
    assert totals.to_dict("records") == [
        {"gene_id": "g1", "gene_name": "GeneOne", "sample_id": "a", "gene_total": 4}
    ]
    assert set(pau["gene_name"]) == {"GeneOne"}


def test_a_gene_with_two_names_is_rejected() -> None:
    atlas = [ATLAS[0], {**ATLAS[1], "gene_name": "OtherName"}]
    per_sample = {"a": [{"pac_id": "p1", "count": 3}, {"pac_id": "p2", "count": 1}]}
    with pytest.raises(ValueError, match="Gene g1 has more than one gene_name"):
        build_count_outputs(per_sample, atlas)
