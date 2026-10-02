import base64
import re
from pathlib import Path

import numpy as np
import pandas as pd

from pacusage.report import _figure_sections, _histogram_svg, build_report
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
                "pac_id": "pac-1",
                "pac_pvalue": None,
                "model_status": "fitted_with_zero_count_stabilization",
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
                "pac_id": "pac-1",
                "pac_pvalue": 0.01,
                "model_status": "fitted",
                "zero_boundary_unstable": False,
                "bootstrap_successes": 40,
                "delta_pau_ci_low": 0.10,
                "delta_pau_ci_high": 0.30,
            },
            {
                "gene_id": "gene-2",
                "pac_id": "pac-2",
                "pac_pvalue": 0.02,
                "model_status": "fitted_with_zero_count_stabilization",
                "zero_boundary_unstable": False,
                "bootstrap_successes": 14,
                "delta_pau_ci_low": None,
                "delta_pau_ci_high": None,
            },
            {
                "gene_id": "gene-3",
                "pac_id": "pac-3",
                "pac_pvalue": 0.50,
                "model_status": "fitted",
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
    assert "<strong>3</strong>Finite PAC p-values" in report
    assert "<strong>1</strong>Tests with zero-count stabilization" in report


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
        "<h2>Input checks</h2>",
        "<h2>Where reads end, relative to known transcript ends</h2>",
        "<h2>Read-end offset profile</h2>",
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
    assert "<h2>Read-end offset profile</h2>" in report


def test_report_tolerates_missing_calibration_diagnostics(tmp_path: Path) -> None:
    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 1)

    report = output.read_text()
    assert "Calibration warning" not in report
    assert "<h2>Read-end offset profile</h2>" not in report


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
    # The PCA is a figure from the figures step, no longer a table.
    assert "PAU principal components" not in report
    assert report.count("<table") == 1


def test_pau_pca_figure_follows_the_correlation_table(tmp_path: Path) -> None:
    rows = [
        {"gene_id": "g", "pac_id": f"g.{index}", "sample_id": sample, "gene_total": 50,
         "pau": pau}
        for sample, usage in (("S1", [0.2, 0.8]), ("S2", [0.4, 0.6]), ("S3", [0.9, 0.1]))
        for index, pau in enumerate(usage)
    ]
    path = tmp_path / "observed_pau.tsv.gz"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression=gzip_compression(path))
    write_figures(tmp_path, ["pau_pca", "event_counts"])
    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 10)
    report = output.read_text()
    start = report.index("<h2>PAU principal components</h2>")
    end = report.index("</section>", start)
    assert report.index("<h2>PAU sample correlation</h2>") < start
    assert start < report.index("<h2>Treatment-control figures</h2>")
    pca = report[start:end]
    assert figure_sources(pca) == [encoded("pau_pca")]
    assert figure_links(pca) == ["../figures/pau_pca.pdf"]
    assert "<table" not in pca
    # The treatment-control figures leave the PCA out.
    assert report.count(encoded("pau_pca")) == 1


def test_pau_pca_figure_is_shown_without_observed_pau(tmp_path: Path) -> None:
    write_figures(tmp_path, ["pau_pca"])
    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 10)
    report = output.read_text()
    assert "<h2>PAU sample correlation</h2>" not in report
    assert figure_sources(report) == [encoded("pau_pca")]
    assert "Treatment-control figures" not in report


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


def test_gene_plot_title_shows_name_and_id(tmp_path: Path) -> None:
    """Gene plot title should show both gene name and gene ID."""
    counts = tmp_path / "counts"
    counts.mkdir()
    rows = [
        {
            "gene_id": "ENSG00000000001",
            "gene_name": "BRCA1",
            "pac_id": "pac1",
            "sample_id": "S1",
            "count": 100,
            "gene_total": 200,
            "pau": 0.5,
        }
    ]
    path = counts / "observed_pau.tsv.gz"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression=gzip_compression(path))

    write_tsv(
        [{"pac_id": "pac1", "coordinate": 1000, "strand": "+"}],
        tmp_path / "pacs.v1.metadata.tsv.gz",
    )
    write_tsv(
        [{"sample_id": "S1", "condition": "ctrl"}],
        tmp_path / "normalized_samples.tsv",
    )

    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 1)

    report = output.read_text()
    # Title should show: <strong>BRCA1</strong> (ENSG00000000001)
    assert "<strong>BRCA1</strong> (ENSG00000000001)" in report


def test_gene_plot_title_omits_id_when_name_equals_id(tmp_path: Path) -> None:
    """When gene_name equals gene_id, show only the ID once."""
    counts = tmp_path / "counts"
    counts.mkdir()
    rows = [
        {
            "gene_id": "GENE001",
            "gene_name": "GENE001",
            "pac_id": "pac1",
            "sample_id": "S1",
            "count": 100,
            "gene_total": 200,
            "pau": 0.5,
        }
    ]
    path = counts / "observed_pau.tsv.gz"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression=gzip_compression(path))

    write_tsv(
        [{"pac_id": "pac1", "coordinate": 1000, "strand": "+"}],
        tmp_path / "pacs.v1.metadata.tsv.gz",
    )
    write_tsv(
        [{"sample_id": "S1", "condition": "ctrl"}],
        tmp_path / "normalized_samples.tsv",
    )

    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 1)

    report = output.read_text()
    # Should show: <strong>GENE001</strong> (not with ID repeated)
    assert "<strong>GENE001</strong>" in report
    # Make sure it's not showing the (ID) part
    assert "<strong>GENE001</strong> (GENE001)" not in report


def test_gene_plot_title_html_escapes_name(tmp_path: Path) -> None:
    """Gene names with HTML characters should be escaped."""
    counts = tmp_path / "counts"
    counts.mkdir()
    rows = [
        {
            "gene_id": "GENE<001>",
            "gene_name": "BR<CA1>",
            "pac_id": "pac1",
            "sample_id": "S1",
            "count": 100,
            "gene_total": 200,
            "pau": 0.5,
        }
    ]
    path = counts / "observed_pau.tsv.gz"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression=gzip_compression(path))

    write_tsv(
        [{"pac_id": "pac1", "coordinate": 1000, "strand": "+"}],
        tmp_path / "pacs.v1.metadata.tsv.gz",
    )
    write_tsv(
        [{"sample_id": "S1", "condition": "ctrl"}],
        tmp_path / "normalized_samples.tsv",
    )

    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 1)

    report = output.read_text()
    # Should escape the < and > in both name and ID
    assert "&lt;CA1&gt;" in report
    assert "&lt;001&gt;" in report
    # Make sure raw HTML tags are not present
    assert "<CA1>" not in report
    assert "<001>" not in report


def test_motif_class_preference_tables_are_rendered(tmp_path: Path) -> None:
    """*.preference_class.tsv.gz tables are rendered under 'PolyA-signal preference by class'."""
    motifs = tmp_path / "motifs"
    motifs.mkdir()

    # Create a preference_class table
    rows = [
        {"motif_id": "m1", "class": "C1", "preference": 0.5},
        {"motif_id": "m2", "class": "C2", "preference": 0.7},
    ]
    path = motifs / "treatment_vs_control.preference_class.tsv.gz"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression=gzip_compression(path))

    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 1)

    report = output.read_text()
    assert "<h2>PolyA-signal preference by class</h2>" in report
    assert "<h2>treatment vs control: polyA-signal preference by class</h2>" in report
    assert "<td>C1</td>" in report and "<td>C2</td>" in report


PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def write_figures(root: Path, stems: list[str]) -> None:
    # Each file's bytes name the file, so a misplaced image is caught.
    for stem in stems:
        (root / f"{stem}.png").write_bytes(PNG_SIGNATURE + stem.encode())


def figure_sources(section: str) -> list[str]:
    return re.findall(r"src='data:image/png;base64,([^']*)'", section)


def figure_links(section: str) -> list[str]:
    return re.findall(r"href='([^']*)'", section)


def encoded(stem: str) -> str:
    return base64.b64encode(PNG_SIGNATURE + stem.encode()).decode("ascii")


def test_figure_section_orders_comparisons_then_summaries(tmp_path: Path) -> None:
    # Written out of order, to show that the section sorts comparisons.
    write_figures(
        tmp_path,
        [
            "B_vs_A.shifts_by_gene_region", "B_vs_A.volcano", "B_vs_A.distal_usage",
            "A_vs_B.distal_usage", "A_vs_B.shifts_by_gene_region", "A_vs_B.volcano",
            "effect_vs_coverage", "concordance", "pau_pca", "concordance_matrix",
            "event_counts", "apa_pattern_grid",
        ],
    )
    section = _figure_sections(tmp_path)
    expected = [
        "A_vs_B.volcano", "A_vs_B.distal_usage", "A_vs_B.shifts_by_gene_region",
        "B_vs_A.volcano", "B_vs_A.distal_usage", "B_vs_A.shifts_by_gene_region",
        "event_counts", "apa_pattern_grid", "concordance_matrix", "concordance",
        "effect_vs_coverage",
    ]
    assert section.startswith("<section><h2>Treatment-control figures</h2>")
    assert re.findall(r"<h3>(.*?)</h3>", section) == ["A vs B", "B vs A", "All comparisons"]
    assert figure_sources(section) == [encoded(stem) for stem in expected]
    assert figure_links(section) == [f"../figures/{stem}.pdf" for stem in expected]
    assert re.findall(r"alt='([^']*)'", section)[:3] == [
        "Volcano plot for A vs B",
        "Distal PAC usage for A vs B",
        "Shifts by gene region for A vs B",
    ]


def test_figure_section_skips_a_missing_figure(tmp_path: Path) -> None:
    write_figures(tmp_path, ["T_vs_C.volcano", "T_vs_C.distal_usage"])
    section = _figure_sections(tmp_path)
    assert figure_sources(section) == [encoded("T_vs_C.volcano"), encoded("T_vs_C.distal_usage")]
    assert section.count("<figure>") == 2
    assert "<h3>All comparisons</h3>" not in section


def test_figure_links_are_percent_encoded(tmp_path: Path) -> None:
    write_figures(tmp_path, ["Drug 10%_vs_DMSO.volcano"])
    section = _figure_sections(tmp_path)
    assert figure_links(section) == ["../figures/Drug%2010%25_vs_DMSO.volcano.pdf"]
    assert "<h3>Drug 10% vs DMSO</h3>" in section


def test_report_has_no_figure_section_without_figures(tmp_path: Path) -> None:
    assert _figure_sections(tmp_path) == ""
    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 1)
    assert "Treatment-control figures" not in output.read_text()


def test_figure_section_comes_before_the_calls_tables(tmp_path: Path) -> None:
    write_figures(tmp_path, ["A_vs_B.volcano"])
    write_tsv(
        [{"pac_id": "p1", "gene_id": "g1", "event_type": "increased_usage"}],
        tmp_path / "A_vs_B.calls.tsv.gz",
    )
    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 1)
    report = output.read_text()
    assert report.index("<h2>Treatment-control figures</h2>") < report.index(
        "<h2>A vs B: PACs with a call</h2>"
    )


def test_statistics_tables_are_titled_by_comparison_and_described(tmp_path: Path) -> None:
    # Underscores inside a condition's name survive; only _vs_ is replaced.
    for suffix in ("calls", "pacs"):
        write_tsv(
            [{"pac_id": "p1", "gene_id": "g1", "event_type": "increased_usage"}],
            tmp_path / f"Drug_high_vs_DMSO.{suffix}.tsv.gz",
        )
    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 1)
    report = output.read_text()
    calls = report.index(
        "<h2>Drug_high vs DMSO: PACs with a call</h2><p class='chart-note'>PACs that gained"
    )
    tested = report.index(
        "<h2>Drug_high vs DMSO: every tested PAC</h2><p class='chart-note'>Every tested PAC,"
    )
    assert calls < tested
    assert "Drug high" not in report


def test_report_opens_with_the_glossary_and_the_column_guide_link(tmp_path: Path) -> None:
    output = tmp_path / "report" / "index.html"
    build_report(tmp_path, output, 1)
    report = output.read_text()
    main = report[report.index("<main>"):]
    assert main.startswith("<main><p class='chart-note'>A PAC is a polyadenylation site")
    assert (
        "<a href='https://github.com/ahepperla/APA_finder_3prime_RNA/blob/main/docs/"
        "output_columns.md'>docs/output_columns.md</a>"
    ) in main


def test_report_with_figures_is_byte_identical_across_runs(tmp_path: Path) -> None:
    write_figures(tmp_path, ["A_vs_B.volcano", "A_vs_B.distal_usage", "event_counts"])
    first = tmp_path / "first" / "index.html"
    second = tmp_path / "second" / "index.html"
    build_report(tmp_path, first, 1)
    build_report(tmp_path, second, 1)
    assert first.read_bytes() == second.read_bytes()
