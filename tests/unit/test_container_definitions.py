from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_container_builds_use_miniforge_conda() -> None:
    definitions = [
        (ROOT / "Dockerfile").read_text(),
        (ROOT / "containers" / "Apptainer.def").read_text(),
    ]

    for definition in definitions:
        assert "micromamba" not in definition
        assert "/opt/conda/bin/conda env create" in definition
        assert "--no-build-isolation" in definition
