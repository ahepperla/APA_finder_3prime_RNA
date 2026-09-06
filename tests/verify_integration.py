#!/usr/bin/env python3
"""Assertions for stable biological behavior of the generated test run."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "results-test"


def main() -> None:
    counts = pd.read_csv(ROOT / "counts" / "pac_counts.tsv.gz", sep="\t")
    pau = pd.read_csv(ROOT / "counts" / "observed_pau.tsv.gz", sep="\t")
    treatment = pd.read_csv(ROOT / "statistics" / "TreatmentA_vs_DMSO.pacs.tsv.gz", sep="\t")
    second_family = pd.read_csv(ROOT / "statistics" / "TreatmentB_vs_Vehicle.pacs.tsv.gz", sep="\t")
    nested_family = pd.read_csv(
        ROOT / "statistics" / "Rescue_vs_TreatmentA.pacs.tsv.gz", sep="\t"
    )
    gained = treatment[treatment["pac_id"].astype(str).str.endswith(".350")].iloc[0]
    assert gained["delta_pau"] > 0.2
    assert gained["event_type"] in {"gained", "gained_candidate"}
    assert gained["bootstrap_successes"] > 0
    assert len(second_family) >= 2
    assert len(nested_family) >= 2
    assert "primary_pas_motif_rna" in treatment
    assert "AAUAAA" in set(treatment["primary_pas_motif_rna"].dropna())

    positive = pau[pau["gene_total"] > 0]
    sums = positive.groupby(["gene_id", "sample_id"])["pau"].sum().to_numpy()
    assert np.allclose(sums, 1)

    sample_columns = [column for column in counts if column not in {"gene_id", "pac_id"}]
    scaled = counts.copy()
    scaled[sample_columns] *= 3
    original_usage = counts[sample_columns].div(
        counts.groupby("gene_id")[sample_columns].transform("sum")
    )
    scaled_usage = scaled[sample_columns].div(
        scaled.groupby("gene_id")[sample_columns].transform("sum")
    )
    assert np.allclose(original_usage.fillna(0), scaled_usage.fillna(0))

    report = ROOT / "report" / "index.html"
    assert report.stat().st_size > 10000
    report_text = report.read_text()
    assert "PAU sample correlation" in report_text
    assert "PAC-level p-value distribution" in report_text
    assert "primary_pas_motif_rna" in report_text
    assert "AAUAAA" in report_text

    traces = list(ROOT.parent.glob("trace-*.txt"))
    assert traces
    latest_trace = max(traces, key=lambda path: path.stat().st_mtime)
    trace_text = latest_trace.read_text()
    assert trace_text.count("PACUSAGE:STATISTICS:FIT_USAGE_MODEL") == 3
    assert trace_text.count("PACUSAGE:STATISTICS:MERGE_USAGE_MODELS") == 1


if __name__ == "__main__":
    main()
