"""Parameter resolution against nextflow_schema.json.

VALIDATE_INPUTS resolves the parameters once: schema defaults, overridden by
what Nextflow passes (analysis.yaml, then the command line). It writes
resolved_params.yaml, which every later step reads as it is.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

import yaml
from jsonschema import Draft7Validator

from .errors import PacusageError

REQUIRED = ("input", "assembly", "fasta", "gtf")
_INTEGER = re.compile(r"[+-]?\d+")
_NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")


def load_schema(path: str | Path) -> dict[str, Any]:
    with Path(path).open() as handle:
        return json.load(handle)


def read_supplied_parameters(path: str | Path) -> dict[str, Any]:
    """The parameters Nextflow passes, as JSON or YAML."""
    path = Path(path)
    with path.open() as handle:
        supplied = json.load(handle) if path.suffix.lower() == ".json" else yaml.safe_load(handle)
    if not isinstance(supplied, dict):
        raise PacusageError(f"Parameter file {path} must contain a mapping.")
    return supplied


def resolve_parameters(supplied: dict[str, Any], schema: dict[str, Any]) -> dict[str, Any]:
    """Schema defaults overridden by ``supplied``, validated against the schema."""
    properties = _schema_properties(schema)
    defaults = {key: spec["default"] for key, spec in properties.items() if "default" in spec}
    params = {**defaults, **coerce_command_line_types(supplied, schema)}
    missing = [key for key in REQUIRED if params.get(key) in (None, "")]
    if missing:
        raise PacusageError(
            "Missing required parameters: "
            + ", ".join(missing)
            + ". Supply them in analysis.yaml or on the Nextflow command line."
        )
    # JSON Schema range checks pass NaN, since every comparison with it is false.
    for key, value in params.items():
        if isinstance(value, float) and not math.isfinite(value):
            raise PacusageError(f"{key} must be a finite number, not {value}.")
    errors = sorted(
        Draft7Validator(schema).iter_errors(params), key=lambda error: list(error.path)
    )
    if errors:
        details = []
        for error in errors[:10]:
            key = ".".join(map(str, error.path)) or "parameters"
            name = str(error.path[0]) if error.path else ""
            description = properties.get(name, {}).get("description", "")
            details.append(f"{key}: {error.message}" + (f" ({description})" if description else ""))
        raise PacusageError("Invalid parameters: " + "; ".join(details))
    if not params["calibration_quantile_low"] < params["calibration_quantile_high"]:
        raise PacusageError("calibration_quantile_low must be below calibration_quantile_high.")
    # The allowance applies to the default count only, so an explicit count
    # beside it would silently ignore it.
    if params.get("min_site_usage_samples") is not None and params.get("min_site_usage_dropouts"):
        raise PacusageError(
            "Set min_site_usage_samples or min_site_usage_dropouts, not both: the dropout "
            "allowance lowers the default count, the size of the family's smallest condition."
        )
    # YAML writes 2.0 for an integer parameter given as 2.0; keep it an integer.
    for key, types in _declared_types(schema).items():
        value = params.get(key)
        if "integer" in types and isinstance(value, float) and value.is_integer():
            params[key] = int(value)
    return params


def coerce_command_line_types(
    supplied: dict[str, Any], schema: dict[str, Any]
) -> dict[str, Any]:
    """Convert command-line strings to the types their parameters declare.

    Nextflow 25.10 and later pass every command-line value as a string, so
    `--min_mapq 30` arrives as "30" and a bare `--save_prepared_alignments` as
    "true". A params file keeps its YAML types. Only an exact integer, a finite
    number, or true/false (in any case) is converted, and only for a parameter
    that accepts that type; a list parameter is never split. Anything else is
    left for validation to report.
    """
    declared = _declared_types(schema)
    coerced = dict(supplied)
    for key, value in supplied.items():
        accepted = declared.get(key, set())
        if not isinstance(value, str) or not accepted or "array" in accepted:
            continue
        text = value.strip()
        if "boolean" in accepted and text.lower() in {"true", "false"}:
            coerced[key] = text.lower() == "true"
        elif accepted & {"integer", "number"} and _INTEGER.fullmatch(text):
            coerced[key] = int(text)
        elif "number" in accepted and _NUMBER.fullmatch(text) and math.isfinite(float(text)):
            coerced[key] = float(text)
    return coerced


def read_resolved_parameters(path: str | Path) -> dict[str, Any]:
    with Path(path).open() as handle:
        return yaml.safe_load(handle)


def write_resolved_parameters(params: dict[str, Any], path: str | Path) -> None:
    with Path(path).open("w") as handle:
        yaml.safe_dump(params, handle, sort_keys=True)


def _schema_properties(schema: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        key: specification
        for group in schema.get("definitions", {}).values()
        for key, specification in group.get("properties", {}).items()
    }


def _declared_types(schema: dict[str, Any]) -> dict[str, set[str]]:
    """The JSON types each parameter accepts."""
    declared: dict[str, set[str]] = {}
    for key, specification in _schema_properties(schema).items():
        types: set[str] = set()
        for branch in specification.get("anyOf", [specification]):
            value = branch.get("type")
            types.update(value if isinstance(value, list) else [value] if value else [])
        declared[key] = types
    return declared
