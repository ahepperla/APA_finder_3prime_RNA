from pathlib import Path

import pytest

from pacusage.errors import PacusageError
from pacusage.parameters import normalize_parameters
from pacusage.samples import read_and_validate_samples


def parameters(tmp_path: Path) -> dict:
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nAAAA\n")
    gtf.write_text('chr1\ttest\tgene\t1\t4\t.\t+\t.\tgene_id "g1";\n')
    return normalize_parameters(
        {
            "input": str(tmp_path / "samples.tsv"),
            "assembly": "test",
            "fasta": str(fasta),
            "gtf": str(gtf),
        }
    )


def touch_alignments(tmp_path: Path, names: list[str]) -> None:
    for name in names:
        (tmp_path / name).touch()


def test_control_mapping_and_relative_paths(tmp_path: Path) -> None:
    names = ["c1.bam", "c2.bam", "t1.bam", "t2.bam"]
    touch_alignments(tmp_path, names)
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\talignment\tcondition\tcontrol\n"
        "C1\tc1.bam\tControl\t\n"
        "C2\tc2.bam\tControl\t\n"
        "T1\tt1.bam\tTreatment\tControl\n"
        "T2\tt2.bam\tTreatment\tControl\n"
    )
    samples, _ = read_and_validate_samples(sheet, parameters(tmp_path))
    assert samples[0].control_condition == "Control"
    assert samples[2].control_condition == "Control"
    assert samples[0].alignment == str((tmp_path / "c1.bam").resolve())


def test_treatment_cannot_reference_treatment(tmp_path: Path) -> None:
    names = ["a1.bam", "a2.bam", "b1.bam", "b2.bam", "c1.bam", "c2.bam"]
    touch_alignments(tmp_path, names)
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\talignment\tcondition\tcontrol\n"
        "A1\ta1.bam\tA\t\nA2\ta2.bam\tA\t\n"
        "B1\tb1.bam\tB\tA\nB2\tb2.bam\tB\tA\n"
        "C1\tc1.bam\tC\tB\nC2\tc2.bam\tC\tB\n"
    )
    with pytest.raises(PacusageError, match="not a control condition"):
        read_and_validate_samples(sheet, parameters(tmp_path))


def test_warn_policy_marks_low_replication(tmp_path: Path) -> None:
    touch_alignments(tmp_path, ["c.bam", "t.bam"])
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\talignment\tcondition\tcontrol\n"
        "C\tc.bam\tControl\t\nT\tt.bam\tTreatment\tControl\n"
    )
    params = parameters(tmp_path)
    params["insufficient_replicates_policy"] = "warn"
    _, checks = read_and_validate_samples(sheet, params)
    assert all(row["exploratory_insufficient_replicates"] == "true" for row in checks)
