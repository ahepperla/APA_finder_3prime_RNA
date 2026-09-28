"""The Nextflow configuration agrees with the modules and the parameter schema."""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MAX_HOURS = 72


def label_blocks(config_text: str) -> dict[str, str]:
    """The body of each ``withLabel: name { ... }`` block (no nested braces)."""
    return dict(re.findall(r"withLabel:\s*(\w+)\s*\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}", config_text))


def test_every_module_label_is_configured_and_every_label_is_used() -> None:
    used = {
        match
        for path in (ROOT / "modules" / "local").glob("*.nf")
        for match in re.findall(r"label '(\w+)'", path.read_text())
    }
    configured = set(label_blocks((ROOT / "conf" / "base.config").read_text()))
    assert used == configured


def test_memory_and_time_scale_with_the_attempt_within_the_ceiling() -> None:
    config = (ROOT / "conf" / "base.config").read_text()
    attempts = int(re.search(r"maxRetries\s*=\s*(\d+)", config).group(1)) + 1
    for label, body in label_blocks(config).items():
        assert re.search(r"memory\s*=\s*\{[^}]*\*\s*task\.attempt\s*\}", body), label
        time = re.search(r"time\s*=\s*\{\s*(\d+)\.h\s*\*\s*task\.attempt\s*\}", body)
        assert time, f"{label}: time must be {{ N.h * task.attempt }}"
        assert int(time.group(1)) * attempts <= MAX_HOURS, label


def test_nextflow_config_defaults_match_the_schema() -> None:
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    expected = {
        key: specification.get("default")
        for group in schema["definitions"].values()
        for key, specification in group["properties"].items()
    }
    config = (ROOT / "nextflow.config").read_text()
    block = re.search(r"^params \{\n(.*?)^\}", config, re.MULTILINE | re.DOTALL).group(1)
    observed = {}
    for key, value in re.findall(r"^\s+(\w+)\s*=\s*(.+?)\s*$", block, re.MULTILINE):
        # The container default is relative to the pipeline directory.
        value = value.replace('"${projectDir}/', '"')
        observed[key] = json.loads(value.replace("'", '"'))
    observed.pop("help")
    assert observed == expected
