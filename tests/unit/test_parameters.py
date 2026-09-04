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
