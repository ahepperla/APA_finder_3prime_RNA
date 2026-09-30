"""Check that docs/output_columns.md documents every column of every table a
run publishes, and nothing else.

The guide gives each table a "### " section headed by its path pattern in
backticks, with one table row per column. The identity columns are documented
once, under "### Identity columns". A section whose body has a line starting
"Same columns as `PATTERN`" takes that pattern's columns. Rows after a line
starting "Columns that only some runs have" are optional: a table may lack
them, but every column it has must be documented.
"""

from __future__ import annotations

import gzip
import re
from collections import Counter
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from pacusage.models import IDENTITY_COLUMNS

GUIDE = Path(__file__).resolve().parents[2] / "docs" / "output_columns.md"
IDENTITY_HEADING = "Identity columns"
OPTIONAL_MARKER = "Columns that only some runs have"
TABLE_SUFFIXES = (".tsv", ".tsv.gz", ".parquet")
# Nextflow's own execution reports.
UNPUBLISHED_DIRECTORIES = {"pipeline_info"}
SAMPLE, COMPARISON, FAMILY = "SAMPLE", "CONDITION_vs_CONTROL", "FAMILY"


Section = tuple[list[str], list[str]]


def read_guide(text: str) -> tuple[list[str], dict[str, Section]]:
    """The identity columns, and each table pattern's documented columns:
    those every run has, then the optional ones."""
    identity: list[str] = []
    sections: dict[str, Section] = {}
    inherited: dict[str, str] = {}
    current: list[str] | None = None
    section: Section | None = None
    patterns: list[str] = []
    for line in text.splitlines():
        if line.startswith("#"):
            current, section, patterns = None, None, []
            if line.startswith("### "):
                heading = line[4:].strip()
                if heading == IDENTITY_HEADING:
                    current = identity
                else:
                    patterns = [name for name in re.findall(r"`([^`]+)`", heading) if "/" in name]
                    if patterns:
                        section = ([], [])
                        current = section[0]
                    for pattern in patterns:
                        assert pattern not in sections, f"guide: {pattern} has two sections"
                        sections[pattern] = section
        elif current is not None and line.startswith("| `"):
            current.extend(re.findall(r"`([^`]+)`", line.split("|")[1]))
        elif section is not None and line.startswith(OPTIONAL_MARKER):
            current = section[1]
        elif patterns and (match := re.match(r"Same columns as `([^`]+)`", line)):
            for pattern in patterns:
                inherited[pattern] = match.group(1)
    for pattern, source in inherited.items():
        assert source in sections, f"guide: {pattern} takes the columns of unknown {source}"
        required, optional = sections[pattern]
        sections[pattern] = (required + sections[source][0], optional + sections[source][1])
    return identity, sections


def table_header(path: Path) -> list[str]:
    if path.name.endswith(".parquet"):
        return list(pq.read_schema(path).names)
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as handle:
        line = handle.readline().rstrip("\n")
    return line.split("\t") if line else []


def run_names(root: Path) -> dict[str, dict[str, str]]:
    """The run's sample, comparison, and family names, each kind mapped to
    its placeholder."""
    samples = pd.read_csv(root / "manifest" / "normalized_samples.tsv", sep="\t", dtype=str)
    mapping = pd.read_csv(root / "qc" / "control_mapping.tsv", sep="\t", dtype=str)
    treatments = mapping[mapping["condition"] != mapping["control_condition"]]
    return {
        SAMPLE: {sample_id: SAMPLE for sample_id in samples["sample_id"]},
        COMPARISON: {
            f"{row.condition}_vs_{row.control_condition}": COMPARISON
            for row in treatments.itertuples()
        },
        FAMILY: {family: FAMILY for family in treatments["control_condition"]},
    }


def table_pattern(relative: str, names: dict[str, dict[str, str]]) -> str:
    """A table's path with its leading name replaced by the placeholder:
    comparison and family names in statistics/, comparison names in motifs/
    and figures/, and sample names elsewhere."""
    directory, _, filename = relative.rpartition("/")
    kinds = {"statistics": (COMPARISON, FAMILY), "motifs": (COMPARISON,), "figures": (COMPARISON,)}
    candidates = {}
    for kind in kinds.get(directory, (SAMPLE,)):
        candidates |= names[kind]
    for name in sorted(candidates, key=len, reverse=True):
        if filename.startswith(f"{name}."):
            filename = candidates[name] + filename[len(name):]
            break
    return f"{directory}/{filename}" if directory else filename


def check_column_guide(root: Path, guide: Path = GUIDE) -> int:
    """Assert that the guide documents exactly the columns of every table
    under root; returns the number of tables checked."""
    identity, sections = read_guide(guide.read_text())
    assert identity == list(IDENTITY_COLUMNS), f"guide: identity columns are {identity}"
    names = run_names(root)
    # Per-sample and per-comparison columns, as in counts/pac_counts.tsv.gz.
    column_names = names[SAMPLE] | names[COMPARISON]
    problems = []
    for pattern, (required, optional) in sections.items():
        repeated = sorted(
            name for name, count in Counter(required + optional).items() if count > 1
        )
        if repeated:
            problems.append(f"{pattern}: documented twice: {', '.join(repeated)}")
    tables = sorted(
        path
        for path in root.rglob("*")
        if path.is_file()
        and path.name.endswith(TABLE_SUFFIXES)
        and not UNPUBLISHED_DIRECTORIES & set(path.relative_to(root).parts)
    )
    assert tables, f"no tables under {root}"
    for path in tables:
        relative = path.relative_to(root).as_posix()
        pattern = table_pattern(relative, names)
        if pattern not in sections:
            problems.append(f"{relative}: no section for {pattern}")
            continue
        header = [column_names.get(column, column) for column in table_header(path)]
        published = set(header) - set(identity)
        required, optional = (set(columns) for columns in sections[pattern])
        if published - required - optional:
            missing = ", ".join(sorted(published - required - optional))
            problems.append(f"{relative}: undocumented columns: {missing}")
        if required - published:
            absent = ", ".join(sorted(required - published))
            problems.append(f"{relative}: documented columns it lacks: {absent}")
    assert not problems, "docs/output_columns.md is out of date:\n" + "\n".join(problems)
    return len(tables)
