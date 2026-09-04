"""Statistical input validation and researcher-facing summaries."""

from __future__ import annotations

import hashlib
from collections.abc import Iterable

import numpy as np
import pandas as pd
from scipy import stats

from .errors import PacusageError


def deterministic_seed(base_seed: int, *parts: str) -> int:
    payload = "\0".join([str(base_seed), *parts]).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "big")


def validate_design_matrix(
    samples: pd.DataFrame, covariates: list[str], include_condition: bool = True
) -> pd.DataFrame:
    missing = [name for name in covariates if name not in samples.columns]
    if missing:
        raise PacusageError("Model covariates absent from sample data: " + ", ".join(missing))
    columns = []
    names = ["intercept"]
    columns.append(np.ones(len(samples)))
    for covariate in covariates:
        if samples[covariate].isna().any() or (samples[covariate].astype(str) == "").any():
            raise PacusageError(f"Covariate {covariate!r} is incomplete in this comparison family.")
        encoded = pd.get_dummies(samples[covariate].astype(str), prefix=covariate, drop_first=True)
        for name in encoded:
            columns.append(encoded[name].to_numpy(dtype=float))
            names.append(name)
    if include_condition:
        encoded = pd.get_dummies(
            samples["condition"].astype(str), prefix="condition", drop_first=True
        )
        for name in encoded:
            columns.append(encoded[name].to_numpy(dtype=float))
            names.append(name)
    matrix = np.column_stack(columns)
    rank = np.linalg.matrix_rank(matrix)
    if rank < matrix.shape[1]:
        raise PacusageError(
            "The model design is not full rank. Condition is confounded with one or more "
            f"covariates ({', '.join(covariates) or 'none'})."
        )
    return pd.DataFrame(matrix, columns=names, index=samples.index)


def filter_testable_features(
    counts: pd.DataFrame,
    minimum_gene_total: int,
    minimum_site_count: int,
    minimum_site_usage: float,
    minimum_supporting_samples: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_columns = [name for name in counts if name not in {"gene_id", "pac_id"}]
    keep_rows = []
    reasons = []
    for gene_id, group in counts.groupby("gene_id", sort=False):
        gene_total = int(group[sample_columns].to_numpy().sum())
        eligible_sites = []
        for _, row in group.iterrows():
            total = int(row[sample_columns].sum())
            supporting = int((row[sample_columns] > 0).sum())
            usage = total / gene_total if gene_total else 0.0
            failures = []
            if total < minimum_site_count:
                failures.append(f"site_count<{minimum_site_count}")
            if supporting < minimum_supporting_samples:
                failures.append(f"supporting_samples<{minimum_supporting_samples}")
            if usage < minimum_site_usage:
                failures.append(f"site_usage<{minimum_site_usage}")
            if not failures:
                eligible_sites.append(str(row["pac_id"]))
            else:
                reasons.append(
                    {
                        "gene_id": gene_id,
                        "pac_id": row["pac_id"],
                        "tested": False,
                        "reason": ";".join(failures),
                    }
                )
        gene_failures = []
        if gene_total < minimum_gene_total:
            gene_failures.append(f"gene_total<{minimum_gene_total}")
        if len(eligible_sites) < 2:
            gene_failures.append("fewer_than_2_testable_pacs")
        if gene_failures:
            reasons.extend(
                {
                    "gene_id": gene_id,
                    "pac_id": pac_id,
                    "tested": False,
                    "reason": ";".join(gene_failures),
                }
                for pac_id in eligible_sites
            )
        else:
            keep_rows.extend(eligible_sites)
    kept = counts[counts["pac_id"].astype(str).isin(keep_rows)].copy()
    return kept, pd.DataFrame(reasons, columns=["gene_id", "pac_id", "tested", "reason"])


def classify_event(row: pd.Series, params: dict) -> str:
    control_detected = int(row["control_supporting_samples"]) >= int(
        params["event_min_supporting_samples"]
    )
    treatment_detected = int(row["treatment_supporting_samples"]) >= int(
        params["event_min_supporting_samples"]
    )
    positive = float(row["delta_pau"]) >= float(params["min_abs_delta_pau"])
    negative = float(row["delta_pau"]) <= -float(params["min_abs_delta_pau"])
    significant = _finite_and_at_most(
        row.get("gene_fdr"), float(params["gene_fdr"])
    ) and _finite_and_at_most(row.get("pac_fdr"), float(params["site_fdr"]))
    stable = not bool(row.get("zero_boundary_unstable", False))
    confident = str(row.get("confidence", "")) != "low" and not bool(
        row.get("internal_priming_flag", False)
    )
    gained_detection = (
        not control_detected
        and treatment_detected
        and float(row["fitted_control_pau"]) <= float(params["event_max_control_pau"])
        and float(row["fitted_treatment_pau"]) >= float(params["event_min_treatment_pau"])
        and positive
    )
    lost_detection = (
        control_detected
        and not treatment_detected
        and float(row["fitted_treatment_pau"]) <= float(params["event_max_control_pau"])
        and float(row["fitted_control_pau"]) >= float(params["event_min_treatment_pau"])
        and negative
    )
    if gained_detection:
        return "gained" if significant and stable and confident else "gained_candidate"
    if lost_detection:
        return "lost" if significant and stable and confident else "lost_candidate"
    if control_detected and treatment_detected and significant and positive:
        return "increased_usage"
    if control_detected and treatment_detected and significant and negative:
        return "decreased_usage"
    return "none"


def _finite_and_at_most(value: object, threshold: float) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return np.isfinite(number) and number <= threshold


def motif_usage_scores(
    pau: pd.DataFrame,
    atlas: pd.DataFrame,
    excluded_rescue_sites: bool = True,
) -> pd.DataFrame:
    columns = [
        "pac_id",
        "primary_pas_motif",
        "primary_motif_class",
        "known_rescue_only",
        "candidate_status",
    ]
    if "primary_pas_motif_rna" in atlas:
        columns.append("primary_pas_motif_rna")
    merged = pau.merge(atlas[columns], on="pac_id", how="inner")
    if "primary_pas_motif_rna" not in merged:
        merged["primary_pas_motif_rna"] = (
            merged["primary_pas_motif"].fillna("").astype(str).str.upper().str.replace("T", "U")
        )
    merged["primary_pas_motif_rna"] = (
        merged["primary_pas_motif_rna"].fillna("").astype(str).replace("", "none")
    )
    if excluded_rescue_sites:
        merged = merged[
            ~merged["known_rescue_only"].astype(bool)
            & (merged["candidate_status"].astype(str) != "motif_assisted_rescue")
        ]
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
    numerator = sum(table[0, 0] * table[1, 1] / table.sum() for table in tables)
    denominator = sum(table[0, 1] * table[1, 0] / table.sum() for table in tables)
    odds_ratio = numerator / denominator if denominator else float("inf")
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
    se = np.sqrt(sum(1 / cell for table in tables for cell in table.flat if cell > 0))
    log_or = np.log(odds_ratio) if 0 < odds_ratio < np.inf else float("nan")
    return {
        "common_odds_ratio": float(odds_ratio),
        "ci_low": float(np.exp(log_or - 1.96 * se)) if np.isfinite(log_or) else float("nan"),
        "ci_high": float(np.exp(log_or + 1.96 * se)) if np.isfinite(log_or) else float("nan"),
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
