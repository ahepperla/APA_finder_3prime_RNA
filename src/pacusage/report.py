"""Build a portable, self-contained HTML analysis report."""

from __future__ import annotations

import base64
import html
import urllib.parse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

# Figures drawn by scripts/plot_usage_figures.R, in report order, as
# (kind, label, caption). The numbers are drawn in the figures themselves.
COMPARISON_FIGURES = (
    (
        "volcano",
        "Volcano plot",
        "Each point is a tested PAC: its change in fitted PAU against its PAC-level "
        "p-value. Colored points have a confirmed call and open circles a candidate "
        "call; grey PACs have no PAC call. Triangles at the top have p-values of 0. "
        "Up to 20 genes in each direction are labeled with their names.",
    ),
    (
        "distal_usage",
        "Distal PAC usage",
        "Each point is a tested gene: the fitted usage of its most 3' tested PAC in the "
        "last exon or downstream of the gene, control against treatment. Colors show the gene's "
        "APA pattern from the .genes table; filled points have a confirmed call at the "
        "distal PAC, and diamonds mark a pattern that only PACs flagged for possible "
        "internal priming support. Up to 20 genes whose distal PAC rose and 20 whose "
        "distal PAC fell, the largest changes, are labeled.",
    ),
    (
        "calls_by_gene_region",
        "PAC calls by gene region",
        "Confirmed PAC calls by the part of the gene the PAC lies in (gene_region). Left "
        "of zero, PACs that lost usage (lost or decreased); right of zero, PACs that "
        "gained usage (gained or increased). Each label gives the number of PACs tested "
        "in that region.",
    ),
)
SUMMARY_FIGURES = (
    (
        "event_counts",
        "Events per comparison",
        "PAC events count PACs, with candidates in lighter shades. Gene events count "
        "genes that pass the gene-level screen, from the dominant_switch and "
        "active_pacs_change columns of the .genes tables. The lower panel counts genes "
        "per APA pattern (apa_pattern); a gene with two patterns counts in both. Patterns "
        "that only PACs flagged for possible internal priming support, the "
        "_potential_internal_priming patterns, are counted apart.",
    ),
    (
        "apa_pattern_grid",
        "APA patterns shared between comparisons",
        "Genes with an APA pattern in two or more comparisons, up to 50: those shared by "
        "the most comparisons first, then by best gene FDR. Each cell shows the gene's "
        "first pattern in that comparison: a + marks two or more, and a * a pattern that "
        "only PACs flagged for possible internal priming support. Grey marks a gene not "
        "tested there. figures/apa_patterns_by_comparison.tsv.gz lists every tested gene.",
    ),
    (
        "concordance_matrix",
        "Correlation between comparisons",
        "Pearson correlation of the change in PAU at the PACs that two comparisons both "
        "tested. Comparisons against the same control share its estimate, so they "
        "correlate positively without any shared biology. Chained comparisons, where one's "
        "treatment is the other's control, estimate that condition from the same samples "
        "with opposite signs, so they correlate negatively.",
    ),
    (
        "concordance",
        "Concordance between comparisons",
        "One panel per pair of comparisons: each PAC tested in both, its change in PAU in "
        "one against the other. Colored PACs have a confirmed call in both comparisons or "
        "in one. The solid line is y = x and the dashed line y = -x. With more than six "
        "comparisons, only pairs that share a control or a condition are drawn.",
    ),
    (
        "effect_vs_coverage",
        "Change in usage against gene coverage",
        "Change in fitted PAU against the reads at the gene in the less-covered group, "
        "one panel per comparison. The dashed line is min_gene_total, the coverage a "
        "gained or lost call requires.",
    ),
)
COLUMN_GUIDE_URL = (
    "https://github.com/ahepperla/APA_finder_3prime_RNA/blob/main/docs/output_columns.md"
)
GLOSSARY = (
    "A PAC is a polyadenylation site cluster: one place where transcripts end. PAU is a "
    "PAC's share of its gene's reads; fitted PAU comes from the model, observed PAU from "
    "the counts. Gene FDR says whether a gene's usage changed, and PAC FDR which of its "
    "PACs changed. The column guide explains every column: "
)
PCA_CAPTION = (
    "Each point is a sample, placed by its observed PAU at the PACs of genes with at "
    "least min_gene_total reads in every sample, each PAC centered across samples. "
    "Replicates should sit together; samples that separate by condition differ in PAC "
    "usage. figures/pau_pca.tsv has the coordinates."
)


def build_report(
    results_root: str | Path, output_html: str | Path, minimum_gene_total: int
) -> None:
    root = Path(results_root)
    sections = [
        _calibration_warning_section(_locate(root, "calibration_kernel_diagnostics.tsv")),
        _table_section(
            "Input checks", _locate(root, "input_validation.tsv"),
            description="Each sample and file check, with PASS or the problem found.",
        ),
        _table_section(
            "Reference genome preparation", _locate(root, "reference_preparation.tsv"),
            description="How the reference FASTA was linked and indexed.",
        ),
        _combined_table_section(
            "Alignment preparation", list(root.rglob("*.alignment_preparation.tsv")),
            description="How each sample's alignment was linked, sorted, or indexed.",
        ),
        _table_section(
            "Which control each condition is compared with", _locate(root, "control_mapping.tsv"),
            description="Each condition's direct control, and whether it is a treatment, a "
            "control, or both.",
        ),
        _table_section(
            "Where reads end, relative to known transcript ends",
            _locate(root, "library_calibration.tsv"),
            description="Per sample: how far read ends fall from annotated transcript ends, "
            "which decides how PACs are found.",
        ),
        _table_section(
            "Read-end offset profile", _locate(root, "calibration_kernel_diagnostics.tsv"),
            description="The pooled profile of how far reads end from their PAC, which places "
            "PACs, with the resolution it allows and any warning.",
        ),
        _combined_table_section(
            "Library strandedness", list(root.rglob("*.strandedness.tsv")),
            description="The read orientation found for each sample.",
        ),
        _combined_table_section(
            "Reads kept and removed", list(root.rglob("*.fragment_filtering.tsv")),
            description="Per sample: alignment records examined, fragments kept, and the "
            "spliced reads recorded for the readthrough filter.",
        ),
        _table_section(
            "PAC discovery", _locate(root, "pac_discovery.tsv"),
            description="Candidate PACs accepted and rejected, by filter. "
            "atlas/rejected_candidates.tsv.gz gives each rejection's reason.",
        ),
        _statistical_filtering_section(list(root.rglob("*.statistical_filtering.tsv.gz"))),
        _combined_table_section(
            "Reads assigned to PACs", list(root.rglob("*.quantification.tsv")),
            description="Per sample: reads assigned to a PAC, left unassigned, or ambiguous.",
        ),
        _atlas_summary(_locate(root, "pacs.v1.metadata.tsv.gz")),
        _pau_qc_sections(root, minimum_gene_total),
        _model_diagnostic_sections(root),
        _figure_sections(root),
        _statistics_sections(root),
        _top_genes_section(root),
        _gene_plot_section(root),
        _table_section(
            "PolyA-signal usage by sample", _locate(root, "motif_scores.tsv"),
            description="Each sample's usage of PACs grouped by the polyA-signal hexamer "
            "upstream of them.",
        ),
        _motif_sections(root),
    ]
    body = "\n".join(section for section in sections if section)
    glossary = (
        f"<p class='chart-note'>{html.escape(GLOSSARY)}"
        f"<a href='{COLUMN_GUIDE_URL}'>docs/output_columns.md</a>.</p>"
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>PACusage report</title>
<style>
:root {{ color-scheme: light; --ink:#1d252c; --muted:#5e6b73; --line:#d8dee2;
  --accent:#176b63; --warn:#9b4d16; --paper:#ffffff; --wash:#f5f7f6; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; color:var(--ink); background:var(--wash);
  font:14px/1.45 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; }}
header {{ background:#14342f; color:white; padding:28px max(24px,calc((100vw - 1240px)/2)); }}
h1 {{ margin:0 0 4px; font-size:30px; letter-spacing:0; }}
header p {{ margin:0; color:#d4e4df; }}
main {{ max-width:1240px; margin:0 auto; padding:24px; }}
section {{ margin:0 0 26px; }}
h2 {{ font-size:19px; margin:0 0 10px; }}
h3 {{ font-size:14px; margin:0; }}
.table-wrap {{ overflow:auto; border:1px solid var(--line); background:var(--paper); }}
table {{ border-collapse:collapse; width:100%; font-size:12px; }}
th,td {{ border-bottom:1px solid var(--line); padding:7px 9px; text-align:left;
  white-space:nowrap; }}
th {{ position:sticky; top:0; background:#e9efed; z-index:1; }}
tr:last-child td {{ border-bottom:0; }}
.empty {{ color:var(--muted); }}
.metrics {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr)); gap:8px; }}
.metric {{ border-left:4px solid var(--accent); background:var(--paper); padding:12px; }}
.metric strong {{ display:block; font-size:23px; }}
.warning {{ border-left:4px solid var(--warn); background:#fbf1e8; padding:12px 14px; }}
.chart-panel {{ margin-top:10px; border:1px solid var(--line); background:var(--paper);
  padding:14px; }}
.chart-note {{ color:var(--muted); margin:2px 0 10px; }}
.histogram-layout {{ display:grid; grid-template-columns:20px 34px minmax(0,1fr);
  grid-template-rows:220px 24px 24px; width:min(760px,100%); }}
.histogram-y-label {{ grid-column:1; grid-row:1; align-self:center; justify-self:center;
  color:var(--ink); font-size:11px; writing-mode:vertical-rl; transform:rotate(180deg); }}
.histogram-y-axis {{ grid-column:2; grid-row:1; position:relative;
  border-right:1px solid var(--muted); }}
.histogram-y-tick {{ position:absolute; right:8px; color:var(--muted); font-size:11px;
  line-height:1; transform:translateY(50%); }}
.histogram-plot {{ grid-column:3; grid-row:1; position:relative;
  border-bottom:1px solid var(--muted); overflow:hidden; }}
.histogram-grid-line {{ position:absolute; left:0; right:0; border-top:1px solid #e4e9e7; }}
.histogram-bars {{ position:absolute; inset:0; display:grid; align-items:end;
  gap:clamp(4px,1.5vw,12px); padding:0 clamp(4px,1vw,10px); }}
.histogram-bin {{ height:100%; min-width:0; display:flex; align-items:flex-end; }}
.histogram-bar {{ width:100%; min-height:0; background:var(--accent); color:white;
  display:flex; justify-content:center; align-items:flex-start; padding-top:4px;
  font-size:11px; font-weight:600; }}
.histogram-x-ticks {{ grid-column:3; grid-row:2; display:flex;
  justify-content:space-between; color:var(--muted); font-size:11px; padding-top:5px; }}
.histogram-x-label {{ grid-column:3; grid-row:3; text-align:center; font-size:12px; }}
input[type=search] {{ width:min(420px,100%); border:1px solid #aeb9bd; padding:8px 10px;
  margin:0 0 8px; background:white; }}
.figure-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr));
  gap:14px; align-items:start; }}
figure {{ margin:0; }}
.single-figure figure {{ max-width:640px; }}
figure img {{ display:block; width:100%; height:auto; border:1px solid var(--line);
  background:white; }}
figcaption {{ color:var(--muted); font-size:12px; margin-top:6px; }}
</style>
</head>
<body>
<header><h1>PACusage</h1><p>Polyadenylation site discovery and usage analysis</p></header>
<main>{glossary}
{body}</main>
<script>
for (const input of document.querySelectorAll('[data-table-filter]')) {{
  input.addEventListener('input', () => {{
    const table = document.getElementById(input.dataset.tableFilter);
    const query = input.value.toLowerCase();
    for (const row of table.tBodies[0].rows) {{
      row.hidden = !row.textContent.toLowerCase().includes(query);
    }}
  }});
}}
</script>
</body></html>
"""
    output = Path(output_html)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document)


def _locate(root: Path, name: str) -> Path:
    direct = list(root.rglob(name))
    return direct[0] if direct else root / name


def _read_table(path: Path, **options: object) -> pd.DataFrame | None:
    if not path.is_file():
        return None
    try:
        return pd.read_csv(path, sep="\t", **options)
    except (pd.errors.EmptyDataError, OSError):
        return None


def _calibration_warning_section(path: Path) -> str:
    frame = _read_table(path)
    if frame is None:
        return ""
    reasons = frame.loc[frame["status"] == "warning", "reason"].fillna("")
    if reasons.empty:
        return ""
    paragraphs = "".join(f"<p>{html.escape(str(reason))}</p>" for reason in reasons)
    return (
        f"<section class='warning'><h2>Calibration warning</h2>{paragraphs}"
        "<p>See qc/calibration_kernel_diagnostics.tsv. The run continued; check that the "
        "calibration kernel suits this library before trusting proximal-tag PACs.</p></section>"
    )


def _table_section(title: str, path: Path, limit: int = 200, description: str = "") -> str:
    frame = _read_table(path, nrows=limit + 1)
    if frame is None:
        return ""
    truncated = len(frame) > limit
    table_id = "table-" + "".join(character for character in title.lower() if character.isalnum())
    table = frame.head(limit).fillna("").to_html(index=False, escape=True, table_id=table_id)
    note = (
        f"<p class='empty'>Showing the first {limit:,} rows.</p>"
        if truncated
        else ""
    )
    return (
        f"<section><h2>{html.escape(title)}</h2>{_description(description)}"
        f"<input type='search' placeholder='Filter rows' data-table-filter='{table_id}'>"
        f"<div class='table-wrap'>{table}</div>{note}</section>"
    )


def _description(text: str) -> str:
    return f"<p class='chart-note'>{html.escape(text)}</p>" if text else ""


def _combined_table_section(
    title: str, paths: list[Path], limit: int = 200, description: str = ""
) -> str:
    frames = [frame for path in paths if (frame := _read_table(path)) is not None]
    if not frames:
        return ""
    frame = pd.concat(frames, ignore_index=True, sort=False)
    table_id = "table-" + "".join(character for character in title.lower() if character.isalnum())
    table = frame.head(limit).fillna("").to_html(index=False, escape=True, table_id=table_id)
    return (
        f"<section><h2>{html.escape(title)}</h2>{_description(description)}"
        f"<input type='search' placeholder='Filter rows' data-table-filter='{table_id}'>"
        f"<div class='table-wrap'>{table}</div></section>"
    )


def _atlas_summary(path: Path) -> str:
    columns = {
        "gene_id",
        "known_pac",
        "internal_priming_flag",
        "gene_region",
        "confidence",
    }
    frame = _read_table(path, usecols=lambda name: name in columns)
    if frame is None:
        return ""
    regions = Counter(frame.get("gene_region", pd.Series(dtype=str)).fillna("unassigned"))
    confidence = Counter(frame.get("confidence", pd.Series(dtype=str)).fillna("unknown"))
    metrics = {
        "PACs in atlas": len(frame),
        "Genes represented": (
            frame.get("gene_id", pd.Series(dtype=str)).replace("", pd.NA).nunique()
        ),
        "Known PAC matches": int(frame.get("known_pac", pd.Series(dtype=bool)).fillna(False).sum()),
        "Internal priming flags": int(
            frame.get("internal_priming_flag", pd.Series(dtype=bool)).fillna(False).sum()
        ),
    }
    boxes = "".join(
        f"<div class='metric'><strong>{value:,}</strong>{html.escape(label)}</div>"
        for label, value in metrics.items()
    )
    details = "".join(
        f"<p class='empty'>{html.escape(label)}: "
        + html.escape(", ".join(f"{name} {count:,}" for name, count in sorted(counts.items())))
        + "</p>"
        for label, counts in (("PACs by gene region", regions), ("PACs by confidence", confidence))
    )
    return (
        "<section><h2>PAC atlas</h2>"
        + _description(
            "The PACs found, pooled across every sample, before any comparison was tested."
        )
        + f"<div class='metrics'>{boxes}</div>{details}</section>"
    )


def _statistical_filtering_section(paths: list[Path], limit: int = 200) -> str:
    """Per-family counts of tested PACs, then the first untested PACs and why."""
    frames = [frame for path in sorted(paths) if (frame := _read_table(path)) is not None]
    if not frames:
        return ""
    frame = pd.concat(frames, ignore_index=True)
    tested = frame["tested"].astype(str).str.lower() == "true"
    summary = "".join(
        f"<p>{html.escape(str(family))}: {len(rows):,} PACs, {int(rows.sum()):,} tested, "
        f"{int((~rows).sum()):,} not tested.</p>"
        for family, rows in tested.groupby(frame["family"], sort=True)
    )
    untested = frame.loc[~tested]
    table_id = "table-statisticalfiltering"
    table = untested.head(limit).fillna("").to_html(index=False, escape=True, table_id=table_id)
    note = (
        f"<p class='empty'>Showing the first {limit:,} of {len(untested):,} untested PACs.</p>"
        if len(untested) > limit
        else ""
    )
    return (
        "<section><h2>PACs left out of testing, and why</h2>"
        + _description(
            "Per comparison family: PACs tested and not tested. The table lists untested "
            "PACs with the filter that removed each."
        )
        + f"{summary}"
        f"<input type='search' placeholder='Filter rows' data-table-filter='{table_id}'>"
        f"<div class='table-wrap'>{table}</div>{note}</section>"
    )


def _pau_qc_sections(root: Path, minimum_gene_total: int) -> str:
    """Sample PAU correlation on genes covered in every sample, and the PCA.

    A gene without reads in a sample has no PAU there, so only genes with at
    least ``minimum_gene_total`` reads in every sample take part. The PCA is
    plot_usage_figures.R's pau_pca figure, drawn by the same rule.
    """
    figure = _locate(root, "pau_pca.png")
    pca = (
        "<section><h2>PAU principal components</h2><div class='chart-panel single-figure'>"
        + _figure_html(figure, "PAU principal components", PCA_CAPTION)
        + "</div></section>"
        if figure.is_file()
        else ""
    )
    columns = {"gene_id", "pac_id", "sample_id", "pau", "gene_total"}
    frame = _read_table(
        _locate(root, "observed_pau.tsv.gz"),
        usecols=lambda name: name in columns,
    )
    if frame is None or frame.empty:
        return pca
    covered = frame.groupby("gene_id")["gene_total"].min() >= minimum_gene_total
    frame = frame[frame["gene_id"].isin(covered.index[covered])]
    matrix = frame.pivot_table(
        index=["gene_id", "pac_id"], columns="sample_id", values="pau"
    ).dropna()
    if matrix.empty or matrix.shape[1] < 2:
        return pca
    correlation = matrix.corr().round(3)
    correlation.insert(0, "sample_id", correlation.index)
    correlation.index = range(len(correlation))

    genes = matrix.index.get_level_values("gene_id").nunique()
    note = (
        f"Genes with at least {minimum_gene_total} reads in every sample "
        f"({genes:,} genes, {len(matrix):,} PACs)."
    )
    return _frame_section("PAU sample correlation", correlation, note=note) + pca


def _model_diagnostic_sections(root: Path) -> str:
    pac_files = sorted(root.rglob("*.pacs.tsv.gz"))
    if not pac_files:
        return ""
    columns = {
        "pac_pvalue",
        "model_status",
        "zero_boundary_unstable",
        "delta_pau_ci_low",
        "delta_pau_ci_high",
    }
    tests = finite_pvalues = stabilized = unstable = intervals = 0
    pvalue_parts: list[np.ndarray] = []
    for path in pac_files:
        try:
            chunks = pd.read_csv(
                path,
                sep="\t",
                usecols=lambda name: name in columns,
                chunksize=100_000,
            )
            for frame in chunks:
                tests += len(frame)
                pvalues = pd.to_numeric(
                    frame.get("pac_pvalue", pd.Series(dtype=float)), errors="coerce"
                ).dropna()
                finite_pvalues += len(pvalues)
                pvalue_parts.append(pvalues.to_numpy(dtype=float))
                stabilized += int(
                    frame.get("model_status", pd.Series(dtype=str))
                    .astype(str)
                    .str.contains("zero_count_stabilization")
                    .sum()
                )
                unstable += int(
                    frame.get("zero_boundary_unstable", pd.Series(dtype=bool))
                    .astype(str)
                    .str.lower()
                    .isin(["true", "t", "1"])
                    .sum()
                )
                low = pd.to_numeric(
                    frame.get("delta_pau_ci_low", pd.Series(np.nan, index=frame.index)),
                    errors="coerce",
                )
                high = pd.to_numeric(
                    frame.get("delta_pau_ci_high", pd.Series(np.nan, index=frame.index)),
                    errors="coerce",
                )
                intervals += int((np.isfinite(low) & np.isfinite(high)).sum())
        except (pd.errors.EmptyDataError, OSError):
            continue
    if not tests:
        return ""
    metrics = {
        "PAC tests": tests,
        "Finite PAC p-values": finite_pvalues,
        "Tests with zero-count stabilization": stabilized,
        "Unstable at a zero boundary": unstable,
        "Bootstrap intervals": intervals,
    }
    boxes = "".join(
        f"<div class='metric'><strong>{value:,}</strong>{html.escape(label)}</div>"
        for label, value in metrics.items()
    )
    pvalues = np.concatenate(pvalue_parts) if pvalue_parts else np.asarray([], dtype=float)
    histogram = _histogram_svg(pvalues)
    return (
        "<section><h2>Model diagnostics</h2>"
        + _description(
            "PAC tests across comparisons, how many needed stabilizing, and the spread of "
            "their p-values."
        )
        + f"<div class='metrics'>{boxes}</div>"
        f"{histogram}</section>"
    )


def _histogram_svg(values: np.ndarray, bins: int = 10) -> str:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    finite = finite[(finite >= 0) & (finite <= 1)]
    if not len(finite):
        return "<p class='empty'>No finite PAC-level p-values were available.</p>"

    bin_count = min(bins, max(3, int(np.ceil(np.sqrt(len(finite))))))
    counts, edges = np.histogram(finite, bins=bin_count, range=(0, 1))
    maximum = max(int(counts.max()), 1)
    if maximum <= 4:
        y_ticks = list(range(maximum + 1))
    else:
        tick_step = max(1, int(np.ceil(maximum / 4)))
        y_ticks = list(range(0, maximum + 1, tick_step))
        if y_ticks[-1] != maximum:
            y_ticks.append(maximum)

    y_axis = []
    grid = []
    for tick in y_ticks:
        position = tick / maximum * 100
        y_axis.append(
            f"<span class='histogram-y-tick' style='bottom:{position:.2f}%'>{tick}</span>"
        )
        grid.append(
            f"<span class='histogram-grid-line' style='bottom:{position:.2f}%'></span>"
        )

    bars = []
    for index, count in enumerate(counts):
        bar_height = 100 * int(count) / maximum
        range_label = f"{edges[index]:.2f}-{edges[index + 1]:.2f}"
        count_label = str(int(count)) if count else ""
        empty_style = "background:transparent;padding-top:0;" if not count else ""
        bars.append(
            "<div class='histogram-bin'>"
            f"<div class='histogram-bar' style='height:{bar_height:.2f}%;{empty_style}' "
            f"title='p-value {range_label}: {int(count)} PACs' "
            f"aria-label='p-value {range_label}: {int(count)} PACs'>{count_label}</div></div>"
        )

    return (
        "<div class='chart-panel'><h3>PAC-level p-value distribution</h3>"
        f"<p class='chart-note'>{len(finite):,} finite tests across all comparisons. "
        "Counts are printed inside non-empty bins.</p>"
        "<div class='histogram-layout' role='img' "
        "aria-label='Histogram of finite PAC-level model p-values between zero and one'>"
        "<div class='histogram-y-label'>PAC count</div>"
        f"<div class='histogram-y-axis'>{''.join(y_axis)}</div>"
        f"<div class='histogram-plot'>{''.join(grid)}"
        f"<div class='histogram-bars' style='grid-template-columns:repeat({bin_count},"
        f"minmax(0,1fr))'>{''.join(bars)}</div></div>"
        "<div class='histogram-x-ticks'><span>0.00</span><span>0.25</span>"
        "<span>0.50</span><span>0.75</span><span>1.00</span></div>"
        "<div class='histogram-x-label'>PAC-level p-value</div></div></div>"
    )


def _top_genes_section(root: Path) -> str:
    frames = []
    for path in sorted(root.rglob("*.genes.tsv.gz")):
        frame = _read_table(path)
        if frame is None:
            continue
        frame["comparison"] = path.name.removesuffix(".genes.tsv.gz")
        frames.append(frame)
    if not frames:
        return ""
    values = pd.concat(frames, ignore_index=True, sort=False)
    values["gene_fdr"] = pd.to_numeric(values["gene_fdr"], errors="coerce")
    return _frame_section(
        "Genes with the lowest gene FDR",
        values.sort_values("gene_fdr").head(50),
        note="Across comparisons, the 50 genes with the strongest evidence that their PAC "
        "usage changed, from the .genes tables.",
    )


def _gene_plot_section(root: Path) -> str:
    pau_columns = {"gene_id", "gene_name", "pac_id", "sample_id", "count", "pau"}
    pau = _read_table(
        _locate(root, "observed_pau.tsv.gz"),
        usecols=lambda name: name in pau_columns,
    )
    atlas = _read_table(
        _locate(root, "pacs.v1.metadata.tsv.gz"),
        usecols=lambda name: name in {"pac_id", "coordinate", "strand"},
    )
    samples = _read_table(
        _locate(root, "normalized_samples.tsv"),
        usecols=lambda name: name in {"sample_id", "condition"},
    )
    if pau is None or atlas is None or samples is None or pau.empty:
        return ""
    merged = pau.merge(atlas[["pac_id", "coordinate", "strand"]], on="pac_id", how="left").merge(
        samples[["sample_id", "condition"]], on="sample_id", how="left"
    )
    gene_totals = merged.groupby("gene_id")["count"].sum().sort_values(ascending=False)
    plots = [
        _gene_svg(gene_id, merged[merged["gene_id"] == gene_id])
        for gene_id in gene_totals.head(8).index
    ]
    return "<section><h2>PAC usage in the most-read genes</h2>" + "".join(plots) + "</section>"


def _gene_svg(gene_id: str, values: pd.DataFrame) -> str:
    gene_name = str(values["gene_name"].iloc[0])
    summary = (
        values.groupby(["condition", "pac_id", "coordinate"], as_index=False)["pau"]
        .mean()
        .fillna(0)
    )
    conditions = sorted(summary["condition"].astype(str).unique())
    pacs = sorted(summary["pac_id"].astype(str).unique())
    colors = ["#176b63", "#b85c38", "#466995", "#7a6c5d", "#6b5b95"]
    width, height = 620, 170
    group_width = width / max(len(pacs), 1)
    bar_width = max(5, (group_width - 14) / max(len(conditions), 1))
    bars = []
    for pac_index, pac_id in enumerate(pacs):
        for condition_index, condition in enumerate(conditions):
            subset = summary[
                (summary["pac_id"].astype(str) == pac_id)
                & (summary["condition"].astype(str) == condition)
            ]
            usage = float(subset["pau"].iloc[0]) if len(subset) else 0.0
            x = pac_index * group_width + 8 + condition_index * bar_width
            bar_height = usage * 110
            bars.append(
                f"<rect x='{x:.1f}' y='{130 - bar_height:.1f}' width='{bar_width - 2:.1f}' "
                f"height='{bar_height:.1f}' fill='{colors[condition_index % len(colors)]}'>"
                f"<title>{html.escape(condition)}: {usage:.3f}</title></rect>"
            )
    legend = " | ".join(conditions)
    # A gene the annotation leaves unnamed carries its ID as its name.
    if gene_name == gene_id:
        title = f"<strong>{html.escape(gene_id)}</strong>"
        label = gene_id
    else:
        title = f"<strong>{html.escape(gene_name)}</strong> ({html.escape(gene_id)})"
        label = f"{gene_name} ({gene_id})"
    return (
        f"<div>{title}"
        f"<span class='empty'> {html.escape(legend)}</span>"
        f"<svg viewBox='0 0 {width} {height}' role='img' "
        f"aria-label='PAC usage for {html.escape(label)}' "
        "style='display:block;max-width:620px;background:white;border:1px solid #d8dee2'>"
        + "".join(bars)
        + "</svg></div>"
    )


def _frame_section(title: str, frame: pd.DataFrame, limit: int = 200, note: str = "") -> str:
    if frame.empty:
        return ""
    table_id = "table-" + "".join(character for character in title.lower() if character.isalnum())
    table = frame.head(limit).fillna("").to_html(index=False, escape=True, table_id=table_id)
    note_html = f"<p class='empty'>{html.escape(note)}</p>" if note else ""
    return (
        f"<section><h2>{html.escape(title)}</h2>{note_html}"
        f"<input type='search' placeholder='Filter rows' data-table-filter='{table_id}'>"
        f"<div class='table-wrap'>{table}</div></section>"
    )


def _figure_sections(root: Path) -> str:
    images: dict[str, Path] = {}
    for path in sorted(root.rglob("*.png")):
        images.setdefault(path.name, path)
    comparisons = sorted(
        name.removesuffix(".volcano.png") for name in images if name.endswith(".volcano.png")
    )
    panels = []
    for comparison in comparisons:
        title = _comparison_title(comparison)
        figures = [
            _figure_html(images[f"{comparison}.{kind}.png"], f"{label} for {title}", caption)
            for kind, label, caption in COMPARISON_FIGURES
            if f"{comparison}.{kind}.png" in images
        ]
        panels.append(_figure_panel(title, figures))
    summary = [
        _figure_html(images[f"{kind}.png"], label, caption)
        for kind, label, caption in SUMMARY_FIGURES
        if f"{kind}.png" in images
    ]
    if summary:
        panels.append(_figure_panel("All comparisons", summary))
    if not panels:
        return ""
    return (
        "<section><h2>Treatment-control figures</h2>"
        "<p class='chart-note'>Drawn from each comparison's .pacs and .genes tables. "
        "The same figures are in figures/ as vector PDFs.</p>" + "".join(panels) + "</section>"
    )


def _figure_panel(title: str, figures: list[str]) -> str:
    return (
        f"<div class='chart-panel'><h3>{html.escape(title)}</h3>"
        f"<div class='figure-grid'>{''.join(figures)}</div></div>"
    )


def _figure_html(path: Path, alt: str, caption: str) -> str:
    # The report sits in report/, beside the figures/ folder that holds the PDFs.
    pdf = urllib.parse.quote(path.name.removesuffix(".png") + ".pdf")
    image = base64.b64encode(path.read_bytes()).decode("ascii")
    return (
        f"<figure><img src='data:image/png;base64,{image}' alt='{html.escape(alt)}'>"
        f"<figcaption>{html.escape(caption)} <a href='../figures/{pdf}'>PDF</a>"
        "</figcaption></figure>"
    )


def _comparison_title(stem: str) -> str:
    """A comparison's file stem, CONDITION_vs_CONTROL, as a title."""
    return stem.replace("_vs_", " vs ")


STATISTICS_TABLES = (
    (
        ".calls.tsv.gz",
        "PACs with a call",
        "PACs that gained or lost usage, or rose or fell significantly, and the candidates. "
        "Columns are described in docs/output_columns.md.",
    ),
    (
        ".pacs.tsv.gz",
        "every tested PAC",
        "Every tested PAC, with its change in usage, its tests, and its annotation.",
    ),
)


def _statistics_sections(root: Path) -> str:
    return "\n".join(
        _table_section(
            f"{_comparison_title(path.name.removesuffix(suffix))}: {label}", path,
            description=description,
        )
        for suffix, label, description in STATISTICS_TABLES
        for path in sorted(root.rglob(f"*{suffix}"))
    )


def _motif_sections(root: Path) -> str:
    def tables(suffix: str, label: str, description: str) -> list[str]:
        return [
            _table_section(
                f"{_comparison_title(path.name.removesuffix(suffix))}: {label}", path,
                description=description,
            )
            for path in sorted(root.rglob(f"*{suffix}"))
        ]

    sections = tables(
        ".preference.tsv.gz", "polyA-signal preference",
        "Per polyA-signal hexamer: how the usage of PACs with that signal upstream changed, "
        "treatment against control, with its FDR.",
    )
    class_tables = tables(
        ".preference_class.tsv.gz", "polyA-signal preference by class",
        "The same, with the hexamers grouped into classes: canonical, common variant, other "
        "variant, and no recognized signal.",
    )
    if class_tables:
        sections.append(
            "<section><h2>PolyA-signal preference by class</h2>"
            + "\n".join(class_tables)
            + "</section>"
        )
    return "\n".join(sections)
