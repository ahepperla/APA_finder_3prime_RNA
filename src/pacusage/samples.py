"""Sample-sheet parsing and condition/control validation."""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from .errors import PacusageError
from .models import Sample

REQUIRED_COLUMNS = ("sample_id", "alignment", "condition", "control")
OPTIONAL_COLUMNS = {
    "replicate",
    "batch",
    "donor",
    "layout",
    "strandedness",
    "library_profile",
    "evidence_source",
}
SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def read_and_validate_samples(
    path: str | Path, params: dict[str, Any], check_files: bool = True
) -> tuple[list[Sample], list[dict[str, str]]]:
    sheet_path = Path(path).expanduser().resolve()
    if not sheet_path.is_file():
        raise PacusageError(f"Sample sheet does not exist or is not readable: {sheet_path}")

    with sheet_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PacusageError(f"Sample sheet is empty: {sheet_path}")
        missing = [name for name in REQUIRED_COLUMNS if name not in reader.fieldnames]
        if missing:
            raise PacusageError(
                f"Sample sheet {sheet_path} is missing columns: {', '.join(missing)}."
            )
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]

    if not rows:
        raise PacusageError(f"Sample sheet has a header but no samples: {sheet_path}")

    ids = [row["sample_id"] for row in rows]
    duplicates = sorted(sample_id for sample_id, count in Counter(ids).items() if count > 1)
    if duplicates:
        raise PacusageError(f"Duplicate sample_id values: {', '.join(duplicates)}.")
    unsafe = [sample_id for sample_id in ids if not SAFE_SAMPLE_ID.fullmatch(sample_id)]
    if unsafe:
        raise PacusageError(
            "sample_id must contain only letters, digits, dot, underscore, and hyphen; "
            f"invalid: {', '.join(unsafe)}."
        )

    conditions: dict[str, set[str]] = defaultdict(set)
    counts: Counter[str] = Counter()
    for row in rows:
        if not row["condition"]:
            raise PacusageError(f"Sample {row['sample_id']} has an empty condition.")
        conditions[row["condition"]].add(row["control"])
        counts[row["condition"]] += 1
    inconsistent = {
        condition: values for condition, values in conditions.items() if len(values) > 1
    }
    if inconsistent:
        detail = "; ".join(f"{key}: {sorted(value)}" for key, value in sorted(inconsistent.items()))
        raise PacusageError(f"Each condition must have one control value; found {detail}.")

    mapping = {condition: next(iter(values)) for condition, values in conditions.items()}
    control_conditions = {condition for condition, control in mapping.items() if control == ""}
    for condition, control in mapping.items():
        if not control:
            continue
        if control == condition:
            raise PacusageError(f"Condition {condition!r} cannot name itself as its control.")
        if control not in mapping:
            raise PacusageError(
                f"Condition {condition!r} references absent control condition {control!r}."
            )
        if control not in control_conditions:
            raise PacusageError(
                f"Condition {condition!r} references {control!r}, but {control!r} is not "
                "a control condition with blank control fields."
            )
    unreferenced = sorted(
        control for control in control_conditions if control not in set(mapping.values())
    )
    if unreferenced:
        raise PacusageError(
            "Every control condition must be referenced by a treatment; unreferenced: "
            + ", ".join(unreferenced)
            + "."
        )

    minimum = int(params["min_replicates_per_condition"])
    insufficient = {condition: count for condition, count in counts.items() if count < minimum}
    if insufficient and params["insufficient_replicates_policy"] == "error":
        detail = ", ".join(
            f"{condition}={count}" for condition, count in sorted(insufficient.items())
        )
        raise PacusageError(
            f"Conditions below min_replicates_per_condition={minimum}: {detail}. "
            "Add biological replicates or set insufficient_replicates_policy: warn "
            "for an explicitly exploratory run."
        )

    model_covariates = set(params["model_covariates"])
    missing_covariates = sorted(model_covariates.difference(reader.fieldnames or []))
    if missing_covariates:
        raise PacusageError(
            "model_covariates names missing from sample sheet: " + ", ".join(missing_covariates)
        )

    samples: list[Sample] = []
    checks: list[dict[str, str]] = []
    for row in rows:
        alignment = (sheet_path.parent / row["alignment"]).resolve()
        if check_files and not alignment.is_file():
            raise PacusageError(
                f"Sample {row['sample_id']} alignment does not exist: {alignment}. "
                "Paths are resolved relative to the sample sheet."
            )
        if alignment.suffix.lower() not in {".bam", ".cram"}:
            raise PacusageError(
                f"Sample {row['sample_id']} alignment must end in .bam or .cram: {alignment}"
            )
        for covariate in model_covariates:
            if not row.get(covariate, ""):
                raise PacusageError(
                    f"Sample {row['sample_id']} is missing model covariate {covariate!r}."
                )
        layout = row.get("layout") or params["layout"]
        strandedness = row.get("strandedness") or params["strandedness"]
        profile = row.get("library_profile") or params["library_profile"]
        evidence_source = row.get("evidence_source") or params["evidence_source"]
        _validate_sample_choice(row["sample_id"], "layout", layout, {"auto", "SE", "PE"})
        _validate_sample_choice(
            row["sample_id"], "strandedness", strandedness, {"auto", "forward", "reverse"}
        )
        _validate_sample_choice(
            row["sample_id"],
            "library_profile",
            profile,
            {"generic_3prime", "plasmidsaurus_3prime", "exact_boundary"},
        )
        _validate_sample_choice(
            row["sample_id"],
            "evidence_source",
            evidence_source,
            {"auto", "read_3p", "read_5p", "fragment_3p", "polyA_junction"},
        )
        custom = {
            key: value
            for key, value in row.items()
            if key not in REQUIRED_COLUMNS and key not in OPTIONAL_COLUMNS
        }
        sample = Sample(
            sample_id=row["sample_id"],
            alignment=str(alignment),
            condition=row["condition"],
            control=row["control"],
            control_condition=row["control"] or row["condition"],
            replicate=row.get("replicate", ""),
            batch=row.get("batch", ""),
            donor=row.get("donor", ""),
            layout=layout,
            strandedness=strandedness,
            library_profile=profile,
            evidence_source=evidence_source,
            covariates=custom,
        )
        samples.append(sample)
        checks.append(
            {
                "sample_id": sample.sample_id,
                "check": "sample_sheet",
                "status": "PASS",
                "detail": "validated",
                "exploratory_insufficient_replicates": str(
                    sample.condition in insufficient
                ).lower(),
            }
        )
    return samples, checks


def _validate_sample_choice(sample_id: str, field: str, value: str, choices: set[str]) -> None:
    if value not in choices:
        raise PacusageError(
            f"Sample {sample_id} has invalid {field}={value!r}; "
            f"expected one of: {', '.join(sorted(choices))}."
        )


def write_normalized_samples(samples: Iterable[Sample], path: str | Path) -> None:
    rows = [sample.as_dict() for sample in samples]
    base = [
        "sample_id",
        "alignment",
        "condition",
        "control",
        "control_condition",
        "replicate",
        "batch",
        "donor",
        "layout",
        "strandedness",
        "library_profile",
        "evidence_source",
    ]
    extra = sorted({key for row in rows for key in row}.difference(base))
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=base + extra)
        writer.writeheader()
        writer.writerows(rows)


def control_mapping_rows(samples: Iterable[Sample]) -> list[dict[str, str]]:
    mapping = {sample.condition: sample.control_condition for sample in samples}
    return [
        {
            "condition": condition,
            "control_condition": control,
            "role": "control" if condition == control else "treatment",
        }
        for condition, control in sorted(mapping.items())
    ]
