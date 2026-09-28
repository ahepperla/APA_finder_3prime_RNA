import numpy as np
import pandas as pd
import pytest

from pacusage.statistics import (
    add_bh_fdr,
    cmh_kmer_test,
    motif_usage_scores,
)


def test_bh_fdr_is_deterministic() -> None:
    adjusted = add_bh_fdr([0.01, 0.04, 0.03])
    assert np.allclose(adjusted, [0.03, 0.04, 0.04])


def test_motif_scores_give_each_gene_equal_weight() -> None:
    atlas = pd.DataFrame(
        {
            "pac_id": ["a1", "a2", "b1", "b2"],
            "primary_pas_motif": ["AATAAA", "", "AATAAA", ""],
            "primary_motif_class": ["canonical", "none", "canonical", "none"],
            "known_rescue_only": [False] * 4,
            "primary_pas_motif_rna": ["AAUAAA", "none", "AAUAAA", "none"],
        }
    )
    pau = pd.DataFrame(
        {
            "gene_id": ["a", "a", "b", "b"],
            "pac_id": ["a1", "a2", "b1", "b2"],
            "sample_id": ["s"] * 4,
            "count": [90, 10, 1, 9],
            "gene_total": [100, 100, 10, 10],
            "pau": [0.9, 0.1, 0.1, 0.9],
        }
    )
    scores = motif_usage_scores(pau, atlas, minimum_gene_total=5)
    canonical = scores.loc[scores["primary_motif_class"] == "canonical", "motif_usage"].iloc[0]
    assert np.isclose(canonical, 0.5)
    assert set(scores["primary_pas_motif_rna"]) == {"AAUAAA", "none"}


def test_motif_scores_keep_exact_variants_separate() -> None:
    atlas = pd.DataFrame(
        {
            "pac_id": ["p1", "p2"],
            "primary_pas_motif": ["TATAAA", "AGTAAA"],
            "primary_motif_class": ["other_variant", "other_variant"],
            "known_rescue_only": [False, False],
            "primary_pas_motif_rna": ["UAUAAA", "AGUAAA"],
        }
    )
    pau = pd.DataFrame(
        {
            "gene_id": ["g1", "g1"],
            "pac_id": ["p1", "p2"],
            "sample_id": ["s1", "s1"],
            "count": [3, 1],
            "gene_total": [4, 4],
            "pau": [0.75, 0.25],
        }
    )

    scores = motif_usage_scores(pau, atlas, minimum_gene_total=1)

    assert set(scores["primary_pas_motif_rna"]) == {"UAUAAA", "AGUAAA"}
    assert len(scores) == 2


def test_cmh_kmer_test_against_r_mantelhaen_test() -> None:
    # Test case (a): Two strata as described in the task
    # stratum 1: a=3, b=1, c=2, d=4
    # stratum 2: a=2, b=2, c=1, d=5
    event_presence = {
        "gene1": [True, True, True, False],
        "gene2": [True, True, False, False],
    }
    background_presence = {
        "gene1": [True, True, False, False, False, False],
        "gene2": [True, False, False, False, False, False],
    }

    result = cmh_kmer_test(event_presence, background_presence)

    assert np.isclose(result["common_odds_ratio"], 5.5)
    assert np.isclose(result["ci_low"], 0.7254731846, rtol=1e-6)
    assert np.isclose(result["ci_high"], 41.69692366, rtol=1e-6)
    assert np.isclose(result["p_value"], 0.1041180298, rtol=1e-6)
    assert result["informative_genes"] == 2


def test_cmh_kmer_test_with_twenty_identical_strata() -> None:
    # Test case (b): 20 identical strata with a=3, b=1, c=2, d=4
    event_presence = {f"gene{i}": [True, True, True, False] for i in range(20)}
    background_presence = {f"gene{i}": [True, True, False, False, False, False] for i in range(20)}

    result = cmh_kmer_test(event_presence, background_presence)

    assert np.isclose(result["common_odds_ratio"], 6.0)
    assert np.isclose(result["ci_low"], 3.187330757, rtol=1e-6)
    assert np.isclose(result["ci_high"], 11.29471735, rtol=1e-6)
    assert np.isclose(result["p_value"], 4.320463058e-08, rtol=1e-6)


def test_motif_scores_use_only_genes_covered_in_every_sample() -> None:
    # g2 has 100 reads in s1 but only 5 in s2, so it leaves both samples.
    atlas = pd.DataFrame(
        {
            "pac_id": ["p1", "p2", "p3", "p4"],
            "primary_motif_class": ["canonical", "no_recognized_motif"] * 2,
            "known_rescue_only": [False] * 4,
            "primary_pas_motif_rna": ["AAUAAA", "none"] * 2,
        }
    )
    counts = {("g1", "s1"): (60, 40), ("g1", "s2"): (30, 70), ("g2", "s1"): (90, 10),
              ("g2", "s2"): (1, 4)}
    rows = []
    for (gene_id, sample_id), values in counts.items():
        pac_ids = ("p1", "p2") if gene_id == "g1" else ("p3", "p4")
        for pac_id, count in zip(pac_ids, values, strict=True):
            rows.append({"gene_id": gene_id, "pac_id": pac_id, "sample_id": sample_id,
                         "count": count, "gene_total": sum(values), "pau": count / sum(values)})
    pau = pd.DataFrame(rows)

    scores = motif_usage_scores(pau, atlas, minimum_gene_total=10)

    usage = {(row.sample_id, row.primary_pas_motif_rna): row.motif_usage
             for row in scores.itertuples()}
    assert usage == pytest.approx({("s1", "AAUAAA"): 0.6, ("s1", "none"): 0.4,
                                   ("s2", "AAUAAA"): 0.3, ("s2", "none"): 0.7})
    assert set(scores["informative_genes"]) == {1}


def test_motif_scores_class_level_pools_hexamers_within_each_gene() -> None:
    # g1 has two other_variant hexamers and a canonical PAC; g2 has one
    # other_variant PAC and one without a motif; g3 has canonical and none.
    rows = [
        ("g1", "p1", "TATAAA", "other_variant", 3),
        ("g1", "p2", "AGTAAA", "other_variant", 1),
        ("g1", "p3", "AATAAA", "canonical", 4),
        ("g2", "p4", "TATAAA", "other_variant", 1),
        ("g2", "p5", "", "no_recognized_motif", 3),
        ("g3", "p6", "AATAAA", "canonical", 2),
        ("g3", "p7", "", "no_recognized_motif", 2),
    ]
    atlas = pd.DataFrame(
        {
            "pac_id": [row[1] for row in rows],
            "primary_motif_class": [row[3] for row in rows],
            "known_rescue_only": [False] * len(rows),
            "primary_pas_motif_rna": [row[2].replace("T", "U") for row in rows],
        }
    )
    totals = {"g1": 8, "g2": 4, "g3": 4}
    pau = pd.DataFrame(
        {
            "gene_id": [row[0] for row in rows],
            "pac_id": [row[1] for row in rows],
            "sample_id": ["s1"] * len(rows),
            "count": [row[4] for row in rows],
            "gene_total": [totals[row[0]] for row in rows],
        }
    )
    pau["pau"] = pau["count"] / pau["gene_total"]

    classes = motif_usage_scores(pau, atlas, minimum_gene_total=1, level="class")
    assert list(classes.columns) == [
        "sample_id", "primary_motif_class", "motif_usage", "transformed_motif_usage",
        "informative_genes",
    ]
    by_class = classes.set_index("primary_motif_class")
    # other_variant: g1 (3 + 1) / 8 = 0.5 and g2 1 / 4 = 0.25, so the mean is 0.375.
    # canonical: g1 4 / 8 and g3 2 / 4, both 0.5. none: g2 0.75 and g3 0.5.
    assert by_class["motif_usage"].to_dict() == pytest.approx(
        {"other_variant": 0.375, "canonical": 0.5, "no_recognized_motif": 0.625}
    )
    assert by_class["informative_genes"].to_dict() == {
        "other_variant": 2, "canonical": 2, "no_recognized_motif": 2,
    }
    assert np.allclose(
        classes["transformed_motif_usage"], np.arcsin(np.sqrt(classes["motif_usage"]))
    )

    # The same data at motif level keep the two other_variant hexamers apart.
    motifs = motif_usage_scores(pau, atlas, minimum_gene_total=1)
    by_motif = motifs.set_index("primary_pas_motif_rna")
    # UAUAAA: g1 3 / 8 and g2 1 / 4; AGUAAA: g1 1 / 8 only.
    assert by_motif.loc["UAUAAA", "motif_usage"] == pytest.approx(0.3125)
    assert by_motif.loc["AGUAAA", "motif_usage"] == pytest.approx(0.125)
    assert by_motif.loc["AGUAAA", "informative_genes"] == 1


def test_motif_scores_rejects_unknown_level() -> None:
    atlas = pd.DataFrame(
        {
            "pac_id": ["p1"],
            "primary_pas_motif": ["TATAAA"],
            "primary_motif_class": ["other_variant"],
            "known_rescue_only": [False],
            "primary_pas_motif_rna": ["UAUAAA"],
        }
    )
    pau = pd.DataFrame(
        {
            "gene_id": ["g1"],
            "pac_id": ["p1"],
            "sample_id": ["s1"],
            "count": [1],
            "gene_total": [1],
            "pau": [1.0],
        }
    )

    with pytest.raises(ValueError, match="level must be 'motif' or 'class'"):
        motif_usage_scores(pau, atlas, minimum_gene_total=1, level="unknown")
