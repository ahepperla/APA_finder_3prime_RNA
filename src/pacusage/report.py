"""Build a portable, self-contained HTML analysis report."""

from __future__ import annotations

import html
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


def build_report(results_root: str | Path, output_html: str | Path) -> None:
    root = Path(results_root)
    sections = [
        _table_section("Input validation", _locate(root, "input_validation.tsv")),
        _table_section("Reference preparation", _locate(root, "reference_preparation.tsv")),
        _combined_table_section(
            "Alignment preparation", list(root.rglob("*.alignment_preparation.tsv"))
        ),
        _table_section("Condition to control mapping", _locate(root, "control_mapping.tsv")),
        _table_section("Library calibration", _locate(root, "library_calibration.tsv")),
        _combined_table_section("Strandedness", list(root.rglob("*.strandedness.tsv"))),
        _combined_table_section("Fragment filtering", list(root.rglob("*.fragment_filtering.tsv"))),
        _table_section("PAC discovery", _locate(root, "pac_discovery.tsv")),
        _table_section("Statistical filtering", _locate(root, "statistical_filtering.tsv")),
        _combined_table_section("Quantification", list(root.rglob("*.quantification.tsv"))),
        _atlas_summary(_locate(root, "pacs.v1.metadata.tsv.gz")),
        _pau_qc_sections(root),
        _model_diagnostic_sections(root),
        _statistics_sections(root / "statistics"),
        _top_genes_section(root),
        _gene_plot_section(root),
        _table_section("Motif usage by sample", _locate(root, "motif_scores.tsv")),
        _motif_sections(root / "motifs"),
    ]
    body = "\n".join(section for section in sections if section)
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
</style>
</head>
<body>
<header><h1>PACusage</h1><p>Polyadenylation site discovery and usage analysis</p></header>
<main>{body}</main>
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


def _read_table(path: Path) -> pd.DataFrame | None:
    if not path.is_file():
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except (pd.errors.EmptyDataError, OSError):
        return None


def _table_section(title: str, path: Path, limit: int = 200) -> str:
    frame = _read_table(path)
    if frame is None:
        return ""
    table_id = "table-" + "".join(character for character in title.lower() if character.isalnum())
    table = frame.head(limit).fillna("").to_html(index=False, escape=True, table_id=table_id)
    note = (
        f"<p class='empty'>Showing the first {limit:,} of {len(frame):,} rows.</p>"
        if len(frame) > limit
        else ""
    )
    return (
        f"<section><h2>{html.escape(title)}</h2>"
        f"<input type='search' placeholder='Filter rows' data-table-filter='{table_id}'>"
        f"<div class='table-wrap'>{table}</div>{note}</section>"
    )


def _combined_table_section(title: str, paths: list[Path], limit: int = 200) -> str:
    frames = [frame for path in paths if (frame := _read_table(path)) is not None]
    if not frames:
        return ""
    frame = pd.concat(frames, ignore_index=True, sort=False)
    table_id = "table-" + "".join(character for character in title.lower() if character.isalnum())
    table = frame.head(limit).fillna("").to_html(index=False, escape=True, table_id=table_id)
    return (
        f"<section><h2>{html.escape(title)}</h2>"
        f"<input type='search' placeholder='Filter rows' data-table-filter='{table_id}'>"
        f"<div class='table-wrap'>{table}</div></section>"
    )


def _atlas_summary(path: Path) -> str:
    frame = _read_table(path)
    if frame is None:
        return ""
    assignment = Counter(frame.get("assignment_class", pd.Series(dtype=str)).fillna("unassigned"))
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
    details = html.escape(
        json.dumps(
            {"assignment_class": dict(assignment), "confidence": dict(confidence)},
            sort_keys=True,
        )
    )
    return (
        "<section><h2>PAC atlas</h2>"
        f"<div class='metrics'>{boxes}</div><p class='empty'>{details}</p></section>"
    )


def _pau_qc_sections(root: Path) -> str:
    frame = _read_table(_locate(root, "observed_pau.tsv.gz"))
    if frame is None or frame.empty:
        return ""
    matrix = frame.pivot_table(
        index=["gene_id", "pac_id"], columns="sample_id", values="pau", fill_value=0
    )
    correlation = matrix.corr().round(3)
    correlation.insert(0, "sample_id", correlation.index)
    correlation.index = range(len(correlation))
    centered = matrix.T.to_numpy(dtype=float, copy=True)
    centered -= centered.mean(axis=0, keepdims=True)
    if centered.shape[0] >= 2 and centered.shape[1] >= 1:
        u, singular, _ = np.linalg.svd(centered, full_matrices=False)
        scores = u[:, :2] * singular[:2]
        if scores.shape[1] == 1:
            scores = np.column_stack([scores[:, 0], np.zeros(scores.shape[0])])
        pca = pd.DataFrame(
            {
                "sample_id": matrix.columns,
                "PC1": scores[:, 0].round(4),
                "PC2": scores[:, 1].round(4),
            }
        )
    else:
        pca = pd.DataFrame(columns=["sample_id", "PC1", "PC2"])
    return _frame_section("PAU sample correlation", correlation) + _frame_section(
        "PAU principal components", pca
    )


def _model_diagnostic_sections(root: Path) -> str:
    pac_files = sorted(root.rglob("*.pacs.tsv.gz"))
    frames = [frame for path in pac_files if (frame := _read_table(path)) is not None]
    if not frames:
        return ""
    values = pd.concat(frames, ignore_index=True, sort=False)
    metrics = {
        "PAC tests": len(values),
        "Finite PAC p-values": int(
            pd.to_numeric(values.get("pvalue_pac"), errors="coerce").notna().sum()
        ),
        "Stabilized fits": int(
            values.get("model_status", pd.Series(dtype=str))
            .astype(str)
            .str.contains("add_uniform")
            .sum()
        ),
        "Unstable boundaries": int(
            values.get("zero_boundary_unstable", pd.Series(dtype=bool))
            .astype(str)
            .str.lower()
            .isin(["true", "t", "1"])
            .sum()
        ),
        "Bootstrap intervals": int(
            pd.to_numeric(values.get("bootstrap_successes"), errors="coerce").fillna(0).gt(0).sum()
        ),
    }
    boxes = "".join(
        f"<div class='metric'><strong>{value:,}</strong>{html.escape(label)}</div>"
        for label, value in metrics.items()
    )
    pvalues = pd.to_numeric(values.get("pvalue_pac"), errors="coerce").dropna().to_numpy()
    histogram = _histogram_svg(pvalues)
    return (
        f"<section><h2>Model diagnostics</h2><div class='metrics'>{boxes}</div>"
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
    sort_column = "gene_fdr" if "gene_fdr" in values else "pvalue"
    values[sort_column] = pd.to_numeric(values[sort_column], errors="coerce")
    return _frame_section("Top genes", values.sort_values(sort_column).head(50))


def _gene_plot_section(root: Path) -> str:
    pau = _read_table(_locate(root, "observed_pau.tsv.gz"))
    atlas = _read_table(_locate(root, "pacs.v1.metadata.tsv.gz"))
    samples = _read_table(_locate(root, "normalized_samples.tsv"))
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
    return "<section><h2>Top gene usage profiles</h2>" + "".join(plots) + "</section>"


def _gene_svg(gene_id: str, values: pd.DataFrame) -> str:
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
    return (
        f"<div><strong>{html.escape(str(gene_id))}</strong>"
        f"<span class='empty'> {html.escape(legend)}</span>"
        f"<svg viewBox='0 0 {width} {height}' role='img' "
        f"aria-label='PAC usage for {html.escape(str(gene_id))}' "
        "style='display:block;max-width:620px;background:white;border:1px solid #d8dee2'>"
        + "".join(bars)
        + "</svg></div>"
    )


def _frame_section(title: str, frame: pd.DataFrame, limit: int = 200) -> str:
    if frame.empty:
        return ""
    table_id = "table-" + "".join(character for character in title.lower() if character.isalnum())
    table = frame.head(limit).fillna("").to_html(index=False, escape=True, table_id=table_id)
    return (
        f"<section><h2>{html.escape(title)}</h2>"
        f"<input type='search' placeholder='Filter rows' data-table-filter='{table_id}'>"
        f"<div class='table-wrap'>{table}</div></section>"
    )


def _statistics_sections(directory: Path) -> str:
    root = directory if directory.is_dir() else directory.parent
    return "\n".join(
        _table_section(path.name.replace(".tsv.gz", "").replace("_", " "), path)
        for pattern in ("*.events.tsv.gz", "*.pacs.tsv.gz")
        for path in sorted(root.rglob(pattern))
    )


def _motif_sections(directory: Path) -> str:
    root = directory if directory.is_dir() else directory.parent
    return "\n".join(
        _table_section(path.name.replace(".tsv.gz", "").replace("_", " "), path)
        for path in sorted(root.rglob("*.preference.tsv.gz"))
    )
