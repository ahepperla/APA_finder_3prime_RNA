import numpy as np
import pandas as pd
import pytest

from pacusage.errors import PacusageError
from pacusage.statistics import (
    add_bh_fdr,
    deterministic_seed,
    motif_usage_scores,
    validate_design_matrix,
)


def test_design_detects_condition_covariate_confounding() -> None:
    samples = pd.DataFrame({"condition": ["C", "C", "T", "T"], "batch": ["a", "a", "b", "b"]})
    with pytest.raises(PacusageError, match="not full rank"):
        validate_design_matrix(samples, ["batch"])


def test_seed_and_bh_are_deterministic() -> None:
    assert deterministic_seed(1729, "atlas", "g1") == deterministic_seed(1729, "atlas", "g1")
    adjusted = add_bh_fdr([0.01, 0.04, 0.03])
    assert np.allclose(adjusted, [0.03, 0.04, 0.04])


def test_motif_scores_give_each_gene_equal_weight() -> None:
    atlas = pd.DataFrame(
        {
            "pac_id": ["a1", "a2", "b1", "b2"],
            "primary_pas_motif": ["AATAAA", "", "AATAAA", ""],
            "primary_motif_class": ["canonical", "none", "canonical", "none"],
            "known_rescue_only": [False] * 4,
            "candidate_status": ["primary"] * 4,
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
    scores = motif_usage_scores(pau, atlas)
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
            "candidate_status": ["primary", "primary"],
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

    scores = motif_usage_scores(pau, atlas)

    assert set(scores["primary_pas_motif_rna"]) == {"UAUAAA", "AGUAAA"}
    assert len(scores) == 2
