import json
import re
from pathlib import Path

import pytest

from pacusage.errors import PacusageError
from pacusage.parameters import (
    coerce_command_line_types,
    read_supplied_parameters,
    resolve_parameters,
)

REQUIRED = {"input": "samples.tsv", "assembly": "test", "fasta": "genome.fa", "gtf": "genes.gtf"}


def resolve_encoded(tmp_path: Path, schema: dict, **values: object) -> dict:
    """Resolve parameters as VALIDATE_INPUTS does, from Nextflow's JSON encoding.

    Nextflow 25.10 and later encode every command-line value as a string.
    """
    path = tmp_path / "input_params.json"
    path.write_text(json.dumps({**REQUIRED, **values}))
    return resolve_parameters(read_supplied_parameters(path), schema)


def test_required_parameters(schema) -> None:
    with pytest.raises(PacusageError, match="Missing required parameters"):
        resolve_parameters({}, schema)


def test_supplied_values_override_schema_defaults(resolved_params) -> None:
    params = resolved_params(min_mapq=30)
    assert params["min_mapq"] == 30
    assert params["pac_cluster_radius"] == 12
    assert params["proximal_bin_size"] == 25
    assert params["internal_exon_end_filter"] is True
    assert params["potential_internal_priming_withheld_calls"] is True
    assert params["dm_bootstrap_include_candidates"] is True
    assert params["statistics_bootstrap_batch_size"] == 500
    assert params["bind_paths"] == []


def test_invalid_values_name_the_parameter_and_its_meaning(resolved_params) -> None:
    message = r"proximal_bin_size: 0 is less than the minimum of 1 \(Bin width"
    with pytest.raises(PacusageError, match=message):
        resolved_params(proximal_bin_size=0)


def test_fractional_pac_support_is_accepted(resolved_params) -> None:
    assert resolved_params(pac_min_supporting_samples=0.5)["pac_min_supporting_samples"] == 0.5


@pytest.mark.parametrize("value", [0, -0.5, 1.5, float("nan"), True, "invalid"])
def test_invalid_pac_support_threshold_is_rejected(resolved_params, value: object) -> None:
    with pytest.raises(PacusageError, match="pac_min_supporting_samples"):
        resolved_params(pac_min_supporting_samples=value)


@pytest.mark.parametrize("value", [0, 1.5, float("nan"), True, "invalid"])
def test_invalid_constitutive_readthrough_junction_count_is_rejected(
    resolved_params, value: object
) -> None:
    with pytest.raises(PacusageError, match="constitutive_readthrough_min_junction_count"):
        resolved_params(constitutive_readthrough_min_junction_count=value)


def test_integral_floats_become_integers_for_integer_parameters(resolved_params) -> None:
    params = resolved_params(min_mapq=30.0, pac_min_supporting_samples=2.0, gene_fdr=0.05)
    assert params["min_mapq"] == 30 and type(params["min_mapq"]) is int
    assert type(params["pac_min_supporting_samples"]) is int
    assert type(params["gene_fdr"]) is float


def test_calibration_quantiles_must_be_ordered(resolved_params) -> None:
    with pytest.raises(PacusageError, match="calibration_quantile_low must be below"):
        resolved_params(calibration_quantile_low=0.9, calibration_quantile_high=0.1)


def test_command_line_integers_become_integers(tmp_path: Path, schema) -> None:
    params = resolve_encoded(
        tmp_path, schema, min_mapq="30", statistics_bootstrap_cpus="4", dm_bootstrap_replicates="0"
    )
    keys = ("min_mapq", "statistics_bootstrap_cpus", "dm_bootstrap_replicates")
    values = [params[key] for key in keys]
    assert values == [30, 4, 0]
    assert [type(value) for value in values] == [int, int, int]


@pytest.mark.parametrize("value", ["3.5", "30x", "1e3", ""])
def test_command_line_integers_reject_other_values(tmp_path: Path, schema, value: str) -> None:
    message = re.escape(f"min_mapq: '{value}' is not of type 'integer'")
    with pytest.raises(PacusageError, match=message):
        resolve_encoded(tmp_path, schema, min_mapq=value)


def test_command_line_numbers_become_numbers(tmp_path: Path, schema) -> None:
    params = resolve_encoded(
        tmp_path,
        schema,
        gene_fdr="0.01",
        calibration_min_model_margin="1e-1",
        exact_max_central_width="12",
    )
    assert params["gene_fdr"] == 0.01 and type(params["gene_fdr"]) is float
    assert params["calibration_min_model_margin"] == 0.1
    # An integer string stays an integer, as the same value would in YAML.
    assert params["exact_max_central_width"] == 12
    assert type(params["exact_max_central_width"]) is int


@pytest.mark.parametrize("value", ["nan", "inf", "0.01x", "1_000"])
def test_command_line_numbers_reject_other_values(tmp_path: Path, schema, value: str) -> None:
    message = re.escape(f"gene_fdr: '{value}' is not of type 'number'")
    with pytest.raises(PacusageError, match=message):
        resolve_encoded(tmp_path, schema, gene_fdr=value)


def test_command_line_booleans_become_booleans(tmp_path: Path, schema) -> None:
    params = resolve_encoded(
        tmp_path,
        schema,
        save_prepared_alignments="true",
        require_unique="FALSE",
        run_kmer_enrichment="False",
    )
    assert params["save_prepared_alignments"] is True
    assert params["require_unique"] is False
    assert params["run_kmer_enrichment"] is False


@pytest.mark.parametrize("value", ["yes", "1", "0", "on"])
def test_command_line_booleans_reject_other_words(tmp_path: Path, schema, value: str) -> None:
    message = re.escape(f"require_unique: '{value}' is not of type 'boolean'")
    with pytest.raises(PacusageError, match=message):
        resolve_encoded(tmp_path, schema, require_unique=value)


def test_command_line_support_threshold_accepts_counts_and_fractions(
    tmp_path: Path, schema
) -> None:
    # The support threshold accepts a whole-number sample count or a fraction.
    resolved = resolve_encoded(tmp_path, schema, pac_min_supporting_samples="2")
    assert resolved["pac_min_supporting_samples"] == 2
    fraction = resolve_encoded(tmp_path, schema, pac_min_supporting_samples="0.5")
    assert fraction["pac_min_supporting_samples"] == 0.5


def test_command_line_strings_are_not_split_into_lists(tmp_path: Path, schema) -> None:
    # A list comes from a params file; a command-line string is not split.
    message = re.escape("excluded_contigs: 'chrM' is not of type 'array'")
    with pytest.raises(PacusageError, match=message):
        resolve_encoded(tmp_path, schema, excluded_contigs="chrM")
    listed = resolve_encoded(tmp_path, schema, excluded_contigs=["chrM"])
    assert listed["excluded_contigs"] == ["chrM"]


def test_typed_and_string_parameters_are_unchanged(tmp_path: Path, schema) -> None:
    params = resolve_encoded(
        tmp_path, schema, min_mapq=30, require_unique=False, gene_fdr=0.01, outdir="30"
    )
    assert [params[key] for key in ("min_mapq", "require_unique", "gene_fdr", "outdir")] == [
        30,
        False,
        0.01,
        "30",
    ]
    unknown = {"not_a_parameter": "30"}
    assert coerce_command_line_types(unknown, schema) == unknown


@pytest.mark.parametrize(
    "key", ["potential_internal_priming_withheld_calls", "internal_exon_end_filter"]
)
def test_switches_are_booleans_from_the_command_line(tmp_path: Path, schema, key: str) -> None:
    assert resolve_encoded(tmp_path, schema, **{key: "false"})[key] is False
    assert resolve_encoded(tmp_path, schema, **{key: "true"})[key] is True
    with pytest.raises(PacusageError, match=key):
        resolve_encoded(tmp_path, schema, **{key: "sometimes"})


def test_excluded_contigs_default_includes_mitochondrial_contigs(resolved_params) -> None:
    params = resolved_params()
    assert params["excluded_contigs"] == ["chrM", "MT", "chrMT"]


def test_excluded_contigs_empty_list_overrides_default(resolved_params) -> None:
    params = resolved_params(excluded_contigs=[])
    assert params["excluded_contigs"] == []


def test_excluded_contigs_custom_list_replaces_default(resolved_params) -> None:
    params = resolved_params(excluded_contigs=["chrY"])
    assert params["excluded_contigs"] == ["chrY"]
