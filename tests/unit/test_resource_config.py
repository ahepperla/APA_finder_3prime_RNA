import re
from pathlib import Path


def resource_problems(config_text: str, max_hours: float = 72.0) -> list[str]:
    """
    Find resource configuration problems in Nextflow config text.

    Checks each withLabel block for:
    - memory and time fields that don't scale with task.attempt
    - time values that cannot be parsed
    - time values that exceed max_hours after retries
    """
    problems: list[str] = []

    # Parse maxRetries from the config
    max_retries_match = re.search(r'maxRetries\s*=\s*(\d+)', config_text)
    max_retries = int(max_retries_match.group(1)) if max_retries_match else 0
    total_attempts = max_retries + 1

    # Find all withLabel blocks with proper brace matching
    start = 0
    while True:
        start = config_text.find('withLabel:', start)
        if start == -1:
            break

        # Extract label name
        label_match = re.match(r'withLabel:\s*(\w+)\s*\{', config_text[start:])
        if not label_match:
            start += 1
            continue

        label_name = label_match.group(1)
        brace_start = start + label_match.end() - 1  # Position of opening brace

        # Find matching closing brace
        brace_count = 1
        pos = brace_start + 1
        while pos < len(config_text) and brace_count > 0:
            if config_text[pos] == '{':
                brace_count += 1
            elif config_text[pos] == '}':
                brace_count -= 1
            pos += 1

        block_content = config_text[brace_start + 1 : pos - 1]
        start = pos

        # Extract memory and time lines - match key = value anywhere in block
        memory_match = re.search(r'memory\s*=\s*(.+?)(?:\n|$)', block_content)
        time_match = re.search(r'time\s*=\s*(.+?)(?:\n|$)', block_content)

        # Check memory scaling
        if not memory_match:
            problems.append(f"{label_name}: memory field missing")
        else:
            memory_value = memory_match.group(1).strip()
            if 'task.attempt' not in memory_value:
                problems.append(f"{label_name}: memory does not scale with task.attempt")

        # Check time scaling
        if not time_match:
            problems.append(f"{label_name}: time field missing")
        else:
            time_value = time_match.group(1).strip()
            if 'task.attempt' not in time_value:
                problems.append(f"{label_name}: time does not scale with task.attempt")
            else:
                # Try to parse the base time value
                # Expected format: { <number>.<unit> * task.attempt }
                base_match = re.search(r'{\s*(\d+(?:\.\d+)?)\s*\.\s*([mhd])\s*\*', time_value)
                if not base_match:
                    problems.append(f"{label_name}: cannot parse time base value")
                else:
                    try:
                        value = float(base_match.group(1))
                        unit = base_match.group(2)

                        # Convert to hours
                        if unit == 'm':
                            hours = value / 60.0
                        elif unit == 'h':
                            hours = value
                        elif unit == 'd':
                            hours = value * 24.0
                        else:
                            problems.append(f"{label_name}: unknown time unit {unit}")
                            continue

                        # Check if total time (base * attempts) exceeds max_hours
                        total_hours = hours * total_attempts
                        if total_hours > max_hours:
                            problems.append(
                                f"{label_name}: time {hours}h * {total_attempts} attempts = "
                                f"{total_hours}h exceeds {max_hours}h"
                            )
                    except ValueError:
                        problems.append(f"{label_name}: cannot parse time value as number")

    return problems


def test_base_config_resources_scale_with_attempt() -> None:
    """Test that base.config has all labels with proper resource scaling."""
    config_path = Path(__file__).resolve().parents[2] / "conf" / "base.config"
    config_text = config_path.read_text()

    problems = resource_problems(config_text)
    assert problems == [], f"Resource problems found: {problems}"

    # Verify all expected labels are present
    found_labels = set(re.findall(r'withLabel:\s*(\w+)\s*\{', config_text))
    expected_labels = {
        "low",
        "medium",
        "high",
        "statistics_high",
        "statistics_bootstrap",
        "serial_medium",
        "serial_high",
    }
    assert found_labels == expected_labels, (
        f"Label mismatch: found {found_labels}, expected {expected_labels}"
    )


def test_fixed_walltime_is_reported() -> None:
    """Test that fixed (non-scaling) walltime is detected."""
    config_text = """
process {
  maxRetries = 2
  withLabel: statistics_high {
    cpus = 8
    memory = { 32.GB * task.attempt }
    time = 3.d
  }
}
"""
    problems = resource_problems(config_text)
    assert len(problems) == 1
    assert "statistics_high" in problems[0]


def test_scaled_walltime_over_ceiling_is_reported() -> None:
    """Test that scaled walltime exceeding ceiling is detected."""
    config_text = """
process {
  maxRetries = 2
  withLabel: test_label {
    cpus = 8
    memory = { 32.GB * task.attempt }
    time = { 30.h * task.attempt }
  }
}
"""
    problems = resource_problems(config_text)
    # 30h * 3 attempts = 90h > 72h
    assert len(problems) == 1
    assert "test_label" in problems[0]
