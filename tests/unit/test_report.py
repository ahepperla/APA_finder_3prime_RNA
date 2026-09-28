from pathlib import Path

import numpy as np
import pandas as pd

from pacusage.report import _histogram_svg, build_report
from pacusage.tableio import gzip_compression, write_tsv


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
    build_report(tmp_path, output, 1)

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
    build_report(tmp_path, output, 1)

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
    build_report(tmp_path, output, 1)

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
    build_report(tmp_path, output, 1)

    report = output.read_text()
    assert "Calibration warning" not in report
    assert "<h2>Calibration kernel</h2>" in report


def test_report_tolerates_missing_calibration_diagnostics(tmp_path: Path) -> None:
    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 1)

    report = output.read_text()
    assert "Calibration warning" not in report
    assert "<h2>Calibration kernel</h2>" not in report


def test_pau_qc_uses_only_genes_covered_in_every_sample(tmp_path: Path) -> None:
    # Gene A has at least 10 reads everywhere. Gene B has 5 reads in S2, so
    # the whole gene stays out of the QC although its PAU is defined; gene C
    # has none in S2, where its PAU is undefined.
    gene_a = {"S1": [0.5, 0.3, 0.2], "S2": [0.6, 0.3, 0.1], "S3": [0.2, 0.3, 0.5]}
    gene_b = {"S1": [0.9, 0.1], "S2": [0.8, 0.2], "S3": [0.1, 0.9]}
    gene_c = {"S1": [0.7, 0.3], "S2": [None, None], "S3": [0.4, 0.6]}
    rows = []
    for gene_id, values, totals in (
        ("geneA", gene_a, {"S1": 100, "S2": 80, "S3": 90}),
        ("geneB", gene_b, {"S1": 50, "S2": 5, "S3": 60}),
        ("geneC", gene_c, {"S1": 40, "S2": 0, "S3": 30}),
    ):
        for sample_id, usage in values.items():
            for index, pau in enumerate(usage):
                rows.append(
                    {
                        "gene_id": gene_id,
                        "pac_id": f"{gene_id}.{index}",
                        "sample_id": sample_id,
                        "gene_total": totals[sample_id],
                        "pau": pau,
                    }
                )
    path = tmp_path / "observed_pau.tsv.gz"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression=gzip_compression(path))

    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 10)

    report = output.read_text()
    assert "Genes with at least 10 reads in every sample (1 genes, 3 PACs)." in report
    expected = round(float(np.corrcoef(gene_a["S1"], gene_a["S3"])[0, 1]), 3)
    assert expected == -0.929
    assert f"<td>{expected}</td>" in report


def test_statistical_filtering_summarizes_each_family(tmp_path: Path) -> None:
    def family_table(family: str, tested: list[bool]) -> None:
        rows = [
            {
                "family": family,
                "gene_id": f"g{index}",
                "pac_id": f"{family}.p{index}",
                "tested": "TRUE" if value else "FALSE",
                "reason": "" if value else f"site_count<5 ({family}.p{index})",
            }
            for index, value in enumerate(tested)
        ]
        path = tmp_path / f"{family}.statistical_filtering.tsv.gz"
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression=gzip_compression(path))

    family_table("DMSO", [True, False, False])
    family_table("Vehicle", [True, True, False])
    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 1)

    report = output.read_text()
    assert "<p>DMSO: 3 PACs, 1 tested, 2 not tested.</p>" in report
    assert "<p>Vehicle: 3 PACs, 2 tested, 1 not tested.</p>" in report
    for pac_id in ("DMSO.p1", "DMSO.p2", "Vehicle.p2"):
        assert f"site_count&lt;5 ({pac_id})" in report
    assert "DMSO.p0" not in report and "Vehicle.p0" not in report


def test_statistical_filtering_table_is_capped(tmp_path: Path) -> None:
    rows = [
        {"family": "DMSO", "gene_id": "g", "pac_id": f"p{index}", "tested": False,
         "reason": "site_count<5"}
        for index in range(250)
    ]
    path = tmp_path / "DMSO.statistical_filtering.tsv.gz"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression=gzip_compression(path))
    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 1)

    report = output.read_text()
    assert "Showing the first 200 of 250 untested PACs." in report
    assert report.count("<td>site_count&lt;5</td>") == 200
