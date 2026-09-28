"""Fixtures shared by the unit and integration tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from pacusage.parameters import load_schema, resolve_parameters, write_resolved_parameters

SCHEMA = Path(__file__).resolve().parents[1] / "nextflow_schema.json"
PLACEHOLDERS = {
    "input": "samples.tsv",
    "assembly": "test",
    "fasta": "genome.fa",
    "gtf": "genes.gtf",
}


@pytest.fixture(scope="session")
def schema() -> dict:
    return load_schema(SCHEMA)


@pytest.fixture(scope="session")
def resolved_params(schema):
    """Resolve parameters as VALIDATE_INPUTS does.

    Call it with parameter overrides; with ``path``, it also writes the
    resolved parameters there, as the resolved_params.yaml later steps read.
    """

    def resolve(path: Path | None = None, **supplied) -> dict:
        params = resolve_parameters({**PLACEHOLDERS, **supplied}, schema)
        if path is not None:
            write_resolved_parameters(params, path)
        return params

    return resolve
