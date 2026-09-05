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
