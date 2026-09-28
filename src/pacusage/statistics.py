"""Motif usage scores and the exploratory k-mer test."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd
from scipy import stats


def motif_usage_scores(
    pau: pd.DataFrame,
    atlas: pd.DataFrame,
    minimum_gene_total: int,
    excluded_rescue_sites: bool = True,
) -> pd.DataFrame:
    columns = ["pac_id", "primary_motif_class", "known_rescue_only", "primary_pas_motif_rna"]
    merged = pau.merge(atlas[columns], on="pac_id", how="inner")
    merged["primary_pas_motif_rna"] = (
        merged["primary_pas_motif_rna"].fillna("").astype(str).replace("", "none")
    )
    # A fixed, condition-blind gene set: genes with at least minimum_gene_total
    # reads in every sample.
    covered = merged.groupby("gene_id")["gene_total"].min() >= minimum_gene_total
    merged = merged[merged["gene_id"].isin(covered.index[covered])]
    if excluded_rescue_sites:
        merged = merged[~merged["known_rescue_only"].astype(bool)]
    retained_total = merged.groupby(["gene_id", "sample_id"])["count"].transform("sum")
    merged = merged[retained_total > 0].copy()
    merged["renormalized_pau"] = merged["count"] / retained_total[retained_total > 0]
    site_counts = merged.groupby(["gene_id", "sample_id"])["pac_id"].transform("nunique")
    merged = merged[site_counts >= 2]
    motif_keys = ["primary_pas_motif_rna", "primary_motif_class"]
    by_gene = merged.groupby(
        ["sample_id", "gene_id", *motif_keys], as_index=False
    )["renormalized_pau"].sum()
    scores = (
        by_gene.groupby(["sample_id", *motif_keys], as_index=False)["renormalized_pau"]
        .mean()
        .rename(columns={"renormalized_pau": "motif_usage"})
    )
    scores["transformed_motif_usage"] = np.arcsin(np.sqrt(scores["motif_usage"]))
    genes = by_gene.groupby(["sample_id", *motif_keys])["gene_id"].nunique()
    scores["informative_genes"] = [
        int(genes.loc[(row.sample_id, row.primary_pas_motif_rna, row.primary_motif_class)])
        for row in scores.itertuples()
    ]
    return scores


def cmh_kmer_test(
    event_presence: dict[str, list[bool]],
    background_presence: dict[str, list[bool]],
) -> dict[str, float | int]:
    tables = []
    for gene_id in sorted(set(event_presence).intersection(background_presence)):
        event = event_presence[gene_id]
        background = background_presence[gene_id]
        table = np.array(
            [
                [sum(event), len(event) - sum(event)],
                [sum(background), len(background) - sum(background)],
            ],
            dtype=float,
        )
        if table.sum(axis=0).min() == 0 or table.sum(axis=1).min() == 0:
            continue
        tables.append(table)
    if not tables:
        return {
            "common_odds_ratio": float("nan"),
            "ci_low": float("nan"),
            "ci_high": float("nan"),
            "p_value": float("nan"),
            "informative_genes": 0,
        }
    # Per stratum: a, b are event PACs with and without the k-mer; c, d are
    # background PACs with and without it.
    a, b, c, d = (np.array([table.flat[cell] for table in tables]) for cell in range(4))
    n = a + b + c + d
    r, s = a * d / n, b * c / n
    odds_ratio = r.sum() / s.sum() if s.sum() else float("inf")
    score = sum(table[0, 0] - table[0].sum() * table[:, 0].sum() / table.sum() for table in tables)
    variance = sum(
        np.prod(table.sum(axis=0))
        * np.prod(table.sum(axis=1))
        / (table.sum() ** 2 * (table.sum() - 1))
        for table in tables
        if table.sum() > 1
    )
    chi_square = score**2 / variance if variance else float("nan")
    p_value = float(stats.chi2.sf(chi_square, 1)) if np.isfinite(chi_square) else float("nan")
    ci_low = ci_high = float("nan")
    if r.sum() > 0 and s.sum() > 0:
        # Robins-Breslow-Greenland variance of log(OR_MH), as in R's
        # mantelhaen.test.
        p, q = (a + d) / n, (b + c) / n
        log_variance = (
            (p * r).sum() / (2 * r.sum() ** 2)
            + (p * s + q * r).sum() / (2 * r.sum() * s.sum())
            + (q * s).sum() / (2 * s.sum() ** 2)
        )
        half_width = stats.norm.ppf(0.975) * np.sqrt(log_variance)
        ci_low = float(np.exp(np.log(odds_ratio) - half_width))
        ci_high = float(np.exp(np.log(odds_ratio) + half_width))
    return {
        "common_odds_ratio": float(odds_ratio),
        "ci_low": ci_low,
        "ci_high": ci_high,
        "p_value": p_value,
        "informative_genes": len(tables),
    }


def add_bh_fdr(values: Iterable[float]) -> np.ndarray:
    pvalues = np.asarray(list(values), dtype=float)
    result = np.full(len(pvalues), np.nan)
    finite = np.flatnonzero(np.isfinite(pvalues))
    if not len(finite):
        return result
    order = finite[np.argsort(pvalues[finite])]
    adjusted = pvalues[order] * len(order) / np.arange(1, len(order) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result[order] = np.minimum(adjusted, 1.0)
    return result
