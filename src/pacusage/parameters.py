"""Parameter defaults, normalization, and lightweight schema validation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml
from jsonschema import Draft7Validator

from .errors import PacusageError

DEFAULTS: dict[str, Any] = {
    "outdir": "results",
    "known_pacs": None,
    "chromosome_aliases": None,
    "layout": "auto",
    "strandedness": "auto",
    "library_profile": "generic_3prime",
    "endpoint_model": "auto",
    "evidence_source": "auto",
    "min_mapq": 20,
    "require_unique": True,
    "require_proper_pair": True,
    "exclude_duplicates": True,
    "excluded_contigs": [],
    "min_replicates_per_condition": 2,
    "insufficient_replicates_policy": "error",
    "strand_min_informative_fragments": 10000,
    "strand_max_sampled_fragments": 200000,
    "strand_decision_fraction": 0.80,
    "calibration_min_genes": 100,
    "calibration_max_distance": 1000,
    "calibration_quantile_low": 0.05,
    "calibration_quantile_high": 0.95,
    "calibration_min_kernel_correlation": 0.80,
    "calibration_min_model_margin": 0.10,
    "exact_max_median_abs_offset": 2,
    "exact_max_central_width": 12,
    "exact_min_boundary_fraction": 0.50,
    "proximal_min_median_upstream_offset": 3,
    "pac_seed_radius": 2,
    "pac_cluster_radius": 12,
    "pac_min_total_count": 10,
    "pac_min_sample_count": 2,
    "pac_min_supporting_samples": 2,
    "known_pac_rescue_total": 5,
    "known_pac_match_radius": 12,
    "proximal_kernel_overlap_threshold": 0.50,
    "proximal_assignment_likelihood_ratio": 3.0,
    "pac_coordinate_bootstrap_replicates": 200,
    "max_downstream_distance": 5000,
    "pas_scan_upstream_far": 50,
    "pas_scan_upstream_near": 5,
    "pas_core_upstream_far": 35,
    "pas_core_upstream_near": 10,
    "pas_motif_catalog": None,
    "internal_priming_window": 20,
    "internal_priming_max_a_run": 6,
    "internal_priming_max_a_fraction": 0.60,
    "min_gene_total": 20,
    "min_site_count": 5,
    "min_site_usage": 0.01,
    "min_test_supporting_samples": 2,
    "model_covariates": [],
    "gene_fdr": 0.05,
    "site_fdr": 0.05,
    "min_abs_delta_pau": 0.10,
    "dm_bootstrap_replicates": 200,
    "dm_bootstrap_min_success_fraction": 0.80,
    "dm_zero_sensitivity_repeats": 5,
    "dm_zero_max_delta_pau_spread": 0.02,
    "random_seed": 1729,
    "event_min_treatment_pau": 0.05,
    "event_max_control_pau": 0.01,
    "event_min_supporting_samples": 2,
    "motif_preference_min_genes": 50,
    "motif_kmer_length": 6,
    "run_kmer_enrichment": True,
    "save_prepared_reference": False,
    "save_prepared_alignments": False,
    "save_intermediates": False,
}

REQUIRED = ("input", "assembly", "fasta", "gtf")
CHOICES = {
    "layout": {"auto", "SE", "PE"},
    "strandedness": {"auto", "forward", "reverse"},
    "library_profile": {"generic_3prime", "plasmidsaurus_3prime", "exact_boundary"},
    "endpoint_model": {"auto", "exact_boundary", "proximal_tag"},
    "evidence_source": {"auto", "read_3p", "read_5p", "fragment_3p", "polyA_junction"},
    "insufficient_replicates_policy": {"error", "warn"},
}


def load_parameters(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    with path.open() as handle:
        if path.suffix.lower() == ".json":
            supplied = json.load(handle)
        else:
            supplied = yaml.safe_load(handle) or {}
    if not isinstance(supplied, dict):
        raise PacusageError(f"Parameter file {path} must contain a mapping.")
    params = normalize_parameters(supplied)
    validate_against_schema(params)
    return params


def normalize_parameters(supplied: dict[str, Any]) -> dict[str, Any]:
    params = {**DEFAULTS, **supplied}
    missing = [key for key in REQUIRED if params.get(key) in (None, "")]
    if missing:
        raise PacusageError(
            "Missing required parameters: "
            + ", ".join(missing)
            + ". Supply them in analysis.yaml or on the Nextflow command line."
        )
    for key, choices in CHOICES.items():
        if params[key] not in choices:
            expected = ", ".join(sorted(choices))
            raise PacusageError(f"Invalid {key}={params[key]!r}; expected one of: {expected}.")
    if not 0 < float(params["strand_decision_fraction"]) <= 1:
        raise PacusageError("strand_decision_fraction must be in (0, 1].")
    if (
        not 0
        <= float(params["calibration_quantile_low"])
        < float(params["calibration_quantile_high"])
        <= 1
    ):
        raise PacusageError("Calibration quantiles must satisfy 0 <= low < high <= 1.")
    if int(params["min_replicates_per_condition"]) < 1:
        raise PacusageError("min_replicates_per_condition must be at least 1.")
    if not isinstance(params["model_covariates"], list):
        raise PacusageError("model_covariates must be a YAML list.")
    return params


def write_resolved_parameters(params: dict[str, Any], path: str | Path) -> None:
    with Path(path).open("w") as handle:
        yaml.safe_dump(params, handle, sort_keys=True)


def validate_against_schema(params: dict[str, Any]) -> None:
    schema_path = Path(__file__).resolve().parents[2] / "nextflow_schema.json"
    if not schema_path.is_file():
        return
    with schema_path.open() as handle:
        schema = json.load(handle)
    errors = sorted(Draft7Validator(schema).iter_errors(params), key=lambda error: list(error.path))
    if errors:
        details = "; ".join(
            f"{'.'.join(map(str, error.path)) or 'parameters'}: {error.message}"
            for error in errors[:10]
        )
        raise PacusageError(f"Parameter schema validation failed: {details}")
