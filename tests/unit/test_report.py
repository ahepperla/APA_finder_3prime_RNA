from pathlib import Path

import numpy as np

from pacusage.report import _histogram_svg, build_report
from pacusage.tableio import write_tsv


def test_diagnostic_histogram_is_labeled_and_responsive() -> None:
    chart = _histogram_svg(np.array([0.67, 1.0, 1.0]))

    assert "PAC-level p-value distribution" in chart
    assert "PAC-level p-value" in chart
    assert "PAC count" in chart
    assert "class='histogram-layout'" in chart
    assert "3 finite tests" in chart
    assert chart.count("class='histogram-bin'") == 3


def test_report_tolerates_unavailable_model_statistics(tmp_path: Path) -> None:
    statistics = tmp_path / "statistics"
    statistics.mkdir()
    write_tsv(
        [
            {
                "gene_id": "gene-1",
                "feature_id": "pac-1",
                "pvalue_pac": None,
                "model_status": "drimseq_add_uniform",
                "zero_boundary_unstable": None,
                "bootstrap_successes": None,
                "delta_pau_ci_low": None,
                "delta_pau_ci_high": None,
            }
        ],
        statistics / "treatment_vs_control.pacs.tsv.gz",
    )

    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output)

    report = output.read_text()
    assert "Model diagnostics" in report
    assert "No finite PAC-level p-values were available." in report
    assert "gene-1" in report


def test_bootstrap_interval_metric_counts_finite_intervals(tmp_path: Path) -> None:
    statistics = tmp_path / "statistics"
    statistics.mkdir()
    write_tsv(
        [
            {
                "gene_id": "gene-1",
                "feature_id": "pac-1",
                "pvalue_pac": 0.01,
                "model_status": "drimseq",
                "zero_boundary_unstable": False,
                "bootstrap_successes": 40,
                "delta_pau_ci_low": 0.10,
                "delta_pau_ci_high": 0.30,
            },
            {
                "gene_id": "gene-2",
                "feature_id": "pac-2",
                "pvalue_pac": 0.02,
                "model_status": "drimseq_add_uniform",
                "zero_boundary_unstable": False,
                "bootstrap_successes": 14,
                "delta_pau_ci_low": None,
                "delta_pau_ci_high": None,
            },
            {
                "gene_id": "gene-3",
                "feature_id": "pac-3",
                "pvalue_pac": 0.50,
                "model_status": "drimseq",
                "zero_boundary_unstable": False,
                "bootstrap_successes": 0,
                "delta_pau_ci_low": None,
                "delta_pau_ci_high": None,
            },
        ],
        statistics / "treatment_vs_control.pacs.tsv.gz",
    )

    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output)

    report = output.read_text()
    assert "Model diagnostics" in report
    assert "<strong>1</strong>Bootstrap intervals" in report
    assert "<strong>2</strong>Bootstrap intervals" not in report


def _kernel_diagnostics(status: str, reason: str) -> dict[str, object]:
    return {
        "endpoint_model": "proximal_tag",
        "evidence_source": "read_3p",
        "kernel_modes": 3 if status == "warning" else 1,
        "minimum_resolvable_separation": 2 if status == "warning" else 77,
        "median_central_interval_width": 100 if status == "warning" else 177,
        "spread_to_resolution_ratio": 50.0 if status == "warning" else 2.299,
        "status": status,
        "reason": reason,
    }


def test_calibration_warning_leads_the_report(tmp_path: Path) -> None:
    qc = tmp_path / "qc"
    qc.mkdir()
    write_tsv(
        [_kernel_diagnostics("warning", "the pooled kernel has 3 <separated> modes")],
        qc / "calibration_kernel_diagnostics.tsv",
    )
    write_tsv([{"check": "sample_sheet", "status": "passed"}], qc / "input_validation.tsv")
    write_tsv([{"sample_id": "S1", "classification": "exact"}], qc / "library_calibration.tsv")
    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output)

    report = output.read_text()
    assert "<p>the pooled kernel has 3 &lt;separated&gt; modes</p>" in report
    assert "<separated>" not in report
    headings = [
        "<section class='warning'><h2>Calibration warning</h2>",
        "<h2>Input validation</h2>",
        "<h2>Library calibration</h2>",
        "<h2>Calibration kernel</h2>",
    ]
    positions = [report.index(heading) for heading in headings]
    assert positions == sorted(positions)


def test_ok_calibration_kernel_is_tabulated_without_a_warning(tmp_path: Path) -> None:
    qc = tmp_path / "qc"
    qc.mkdir()
    write_tsv([_kernel_diagnostics("ok", "")], qc / "calibration_kernel_diagnostics.tsv")
    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output)

    report = output.read_text()
    assert "Calibration warning" not in report
    assert "<h2>Calibration kernel</h2>" in report


def test_report_tolerates_missing_calibration_diagnostics(tmp_path: Path) -> None:
    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output)

    report = output.read_text()
    assert "Calibration warning" not in report
    assert "<h2>Calibration kernel</h2>" not in report
