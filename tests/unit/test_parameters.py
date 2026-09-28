import json
import re
from pathlib import Path

import pytest

from pacusage import parameters
from pacusage.errors import PacusageError
from pacusage.parameters import coerce_command_line_types, load_parameters, normalize_parameters

REQUIRED = {"input": "samples.tsv", "assembly": "test", "fasta": "genome.fa", "gtf": "genes.gtf"}


def load_encoded(tmp_path: Path, **values: object) -> dict:
    """Load parameters as VALIDATE_INPUTS does, from Nextflow's JSON encoding.

    Nextflow 25.10 and later encode every command-line value as a string.
    """
    path = tmp_path / "input_params.json"
    path.write_text(json.dumps({**REQUIRED, **values}))
    return load_parameters(path)


def test_required_parameters() -> None:
    with pytest.raises(PacusageError, match="Missing required parameters"):
        normalize_parameters({})


def test_command_values_override_defaults() -> None:
    params = normalize_parameters(
        {
            "input": "samples.tsv",
            "assembly": "test",
            "fasta": "genome.fa",
            "gtf": "genes.gtf",
            "min_mapq": 30,
        }
    )
    assert params["min_mapq"] == 30
    assert params["pac_cluster_radius"] == 12
    assert params["proximal_bin_size"] == 25
    assert params["constitutive_readthrough_min_replicate_support"] == "all"
    assert params["dm_bootstrap_include_candidates"] is True
    assert params["statistics_cpus"] == 8
    assert params["statistics_bootstrap_cpus"] == 4
    assert params["statistics_bootstrap_batch_size"] == 500
    assert params["statistics_bootstrap_max_forks"] == 8
    assert params["bind_paths"] == []


def test_proximal_bin_size_must_be_positive() -> None:
    with pytest.raises(PacusageError, match="proximal_bin_size"):
        normalize_parameters(
            {
                "input": "samples.tsv",
                "assembly": "test",
                "fasta": "genome.fa",
                "gtf": "genes.gtf",
                "proximal_bin_size": 0,
            }
        )


def test_fractional_pac_support_is_accepted() -> None:
    params = normalize_parameters(
        {
            "input": "samples.tsv",
            "assembly": "test",
            "fasta": "genome.fa",
            "gtf": "genes.gtf",
            "pac_min_supporting_samples": 0.5,
        }
    )
    assert params["pac_min_supporting_samples"] == 0.5


@pytest.mark.parametrize("value", [0, -0.5, 1.5, float("nan"), True, "invalid"])
def test_invalid_pac_support_threshold_is_rejected(value: object) -> None:
    with pytest.raises(PacusageError, match="pac_min_supporting_samples"):
        normalize_parameters(
            {
                "input": "samples.tsv",
                "assembly": "test",
                "fasta": "genome.fa",
                "gtf": "genes.gtf",
                "pac_min_supporting_samples": value,
            }
        )


@pytest.mark.parametrize("value", [0, 1.5, float("nan"), True, "invalid"])
def test_invalid_constitutive_readthrough_replicate_support_is_rejected(value: object) -> None:
    with pytest.raises(PacusageError, match="constitutive_readthrough_min_replicate_support"):
        normalize_parameters(
            {
                "input": "samples.tsv",
                "assembly": "test",
                "fasta": "genome.fa",
                "gtf": "genes.gtf",
                "constitutive_readthrough_min_replicate_support": value,
            }
        )


def test_constitutive_readthrough_replicate_support_uses_unambiguous_all_or_counts() -> None:
    base = {
        "input": "samples.tsv",
        "assembly": "test",
        "fasta": "genome.fa",
        "gtf": "genes.gtf",
    }
    assert normalize_parameters(
        {**base, "constitutive_readthrough_min_replicate_support": "all"}
    )["constitutive_readthrough_min_replicate_support"] == "all"
    assert normalize_parameters(
        {**base, "constitutive_readthrough_min_replicate_support": 1}
    )["constitutive_readthrough_min_replicate_support"] == 1


@pytest.mark.parametrize("value", [0, 1.5, float("nan"), True, "invalid"])
def test_invalid_constitutive_readthrough_junction_count_is_rejected(value: object) -> None:
    with pytest.raises(PacusageError, match="constitutive_readthrough_min_junction_count"):
        normalize_parameters(
            {
                "input": "samples.tsv",
                "assembly": "test",
                "fasta": "genome.fa",
                "gtf": "genes.gtf",
                "constitutive_readthrough_min_junction_count": value,
            }
        )


def test_command_line_integers_become_integers(tmp_path: Path) -> None:
    params = load_encoded(
        tmp_path, min_mapq="30", statistics_bootstrap_cpus="4", dm_bootstrap_replicates="0"
    )
    keys = ("min_mapq", "statistics_bootstrap_cpus", "dm_bootstrap_replicates")
    values = [params[key] for key in keys]
    assert values == [30, 4, 0]
    assert [type(value) for value in values] == [int, int, int]


@pytest.mark.parametrize("value", ["3.5", "30x", "1e3", ""])
def test_command_line_integers_reject_other_values(tmp_path: Path, value: str) -> None:
    message = re.escape(f"min_mapq: '{value}' is not of type 'integer'")
    with pytest.raises(PacusageError, match=message):
        load_encoded(tmp_path, min_mapq=value)


def test_command_line_numbers_become_numbers(tmp_path: Path) -> None:
    params = load_encoded(
        tmp_path, gene_fdr="0.01", calibration_min_model_margin="1e-1", exact_max_central_width="12"
    )
    assert params["gene_fdr"] == 0.01 and type(params["gene_fdr"]) is float
    assert params["calibration_min_model_margin"] == 0.1
    # An integer string stays an integer, as the same value would in YAML.
    assert params["exact_max_central_width"] == 12
    assert type(params["exact_max_central_width"]) is int


@pytest.mark.parametrize("value", ["nan", "inf", "0.01x", "1_000"])
def test_command_line_numbers_reject_other_values(tmp_path: Path, value: str) -> None:
    message = re.escape(f"gene_fdr: '{value}' is not of type 'number'")
    with pytest.raises(PacusageError, match=message):
        load_encoded(tmp_path, gene_fdr=value)


def test_command_line_booleans_become_booleans(tmp_path: Path) -> None:
    params = load_encoded(
        tmp_path,
        save_prepared_alignments="true",
        require_unique="FALSE",
        run_kmer_enrichment="False",
    )
    assert params["save_prepared_alignments"] is True
    assert params["require_unique"] is False
    assert params["run_kmer_enrichment"] is False


@pytest.mark.parametrize("value", ["yes", "1", "0", "on"])
def test_command_line_booleans_reject_other_words(tmp_path: Path, value: str) -> None:
    message = re.escape(f"require_unique: '{value}' is not of type 'boolean'")
    with pytest.raises(PacusageError, match=message):
        load_encoded(tmp_path, require_unique=value)


def test_command_line_support_thresholds_accept_counts_and_fractions(tmp_path: Path) -> None:
    # Both accept a whole-number sample count or a fraction; the replicate
    # support also accepts "all".
    assert load_encoded(tmp_path, pac_min_supporting_samples="2")["pac_min_supporting_samples"] == 2
    fraction = load_encoded(tmp_path, pac_min_supporting_samples="0.5")
    assert fraction["pac_min_supporting_samples"] == 0.5
    support = "constitutive_readthrough_min_replicate_support"
    assert load_encoded(tmp_path, **{support: "all"})[support] == "all"
    assert load_encoded(tmp_path, **{support: "0.5"})[support] == 0.5
    assert load_encoded(tmp_path, **{support: "3"})[support] == 3


def test_command_line_strings_are_not_split_into_lists(tmp_path: Path) -> None:
    # A list comes from a params file; a command-line string is not split.
    message = re.escape("excluded_contigs: 'chrM' is not of type 'array'")
    with pytest.raises(PacusageError, match=message):
        load_encoded(tmp_path, excluded_contigs="chrM")
    assert load_encoded(tmp_path, excluded_contigs=["chrM"])["excluded_contigs"] == ["chrM"]


def test_typed_and_string_parameters_are_unchanged(tmp_path: Path) -> None:
    params = load_encoded(tmp_path, min_mapq=30, require_unique=False, gene_fdr=0.01, outdir="30")
    assert [params[key] for key in ("min_mapq", "require_unique", "gene_fdr", "outdir")] == [
        30,
        False,
        0.01,
        "30",
    ]
    assert coerce_command_line_types({"not_a_parameter": "30"}, None) == {"not_a_parameter": "30"}


def test_types_come_from_the_defaults_without_a_schema(tmp_path: Path, monkeypatch) -> None:
    # The container installs the package without nextflow_schema.json, and
    # bool("false") would be True downstream.
    monkeypatch.setattr(parameters, "SCHEMA_PATH", tmp_path / "missing_schema.json")
    params = load_encoded(tmp_path, min_mapq="30", require_unique="false", gene_fdr="0.01")
    assert params["min_mapq"] == 30 and type(params["min_mapq"]) is int
    assert params["require_unique"] is False
    assert params["gene_fdr"] == 0.01
