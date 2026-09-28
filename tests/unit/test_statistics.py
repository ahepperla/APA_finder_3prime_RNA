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
