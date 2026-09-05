import pytest

from pacusage.errors import PacusageError
from pacusage.parameters import normalize_parameters


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
