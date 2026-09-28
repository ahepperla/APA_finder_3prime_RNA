#!/usr/bin/env python3
"""Assertions for stable biological behavior of the generated test run."""

from __future__ import annotations

import hashlib
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

REPOSITORY = Path(__file__).resolve().parents[1]
ROOT = REPOSITORY / "results-test"
FIXTURES = REPOSITORY / "tests" / "fixtures"
COMPARISONS = ("TreatmentA_vs_DMSO", "TreatmentB_vs_Vehicle", "Rescue_vs_TreatmentA")
FAMILIES = ("DMSO", "Vehicle", "TreatmentA")


def pac(contig: str, strand: str, coordinate: int) -> str:
    return f"PACv1.synthetic.{contig}.{strand}.{coordinate}"


def statistics_table(name: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / "statistics" / name, sep="\t")


def one_row(table: pd.DataFrame, pac_id: str) -> pd.Series:
    rows = table[table["pac_id"].astype(str) == pac_id]
    assert len(rows) == 1, f"expected one row for {pac_id}, found {len(rows)}"
    return rows.iloc[0]


def number(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def check_event(table: pd.DataFrame, pac_id: str, expected: str, label: str) -> pd.Series:
    row = one_row(table, pac_id)
    observed = row["event_type"]
    assert observed == expected, f"{label}: {pac_id} is {observed}, not {expected}"
    return row


def check_confirmed(row: pd.Series, params: dict, label: str) -> None:
    assert number(row["pac_fdr"]) <= params["site_fdr"], f"{label}: pac_fdr {row['pac_fdr']}"
    assert number(row["gene_fdr"]) <= params["gene_fdr"], f"{label}: gene_fdr {row['gene_fdr']}"


def check_interval(row: pd.Series, params: dict, label: str) -> None:
    low = number(row["delta_pau_ci_low"])
    high = number(row["delta_pau_ci_high"])
    delta = number(row["delta_pau"])
    fraction = params["dm_bootstrap_min_success_fraction"]
    minimum = math.ceil(fraction * params["dm_bootstrap_replicates"])
    assert row["bootstrap_status"] == "ok", f"{label}: bootstrap_status {row['bootstrap_status']}"
    assert np.isfinite(low) and np.isfinite(high), f"{label}: no interval"
    assert low <= delta <= high, f"{label}: interval [{low}, {high}] excludes {delta}"
    successes = number(row["bootstrap_successes"])
    assert successes >= minimum, f"{label}: {successes} bootstrap successes"


def task_directories(trace_text: str, process: str) -> dict[str, Path]:
    """Work directories of one process's tasks, keyed by task tag."""
    lines = trace_text.splitlines()
    header = lines[0].split("\t")
    directories = {}
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t"), strict=True))
        match = re.fullmatch(rf"PACUSAGE:\S*{process} \((.+)\)", row["name"])
        if match:
            candidates = list((REPOSITORY / "work").glob(f"{row['hash']}*"))
            assert len(candidates) == 1, f"no single work directory for {row['name']}"
            directories[match.group(1)] = candidates[0]
    return directories


def check_alignment_handling(trace_text: str) -> None:
    # One scan per sample replaces per-source calibration passes.
    assert trace_text.count("PACUSAGE:DISCOVERY:SCAN_ALIGNMENT") == 10
    assert "CALIBRATE_SAMPLE" not in trace_text
    assert trace_text.count("PACUSAGE:RECORD_INPUT_CHECKSUMS") == 1
    samples = pd.read_csv(FIXTURES / "samples.tsv", sep="\t", dtype=str)
    sources = dict(zip(samples["sample_id"], samples["alignment"], strict=True))
    prepared = task_directories(trace_text, "PREPARE_ALIGNMENT")
    assert set(prepared) == set(sources)
    for sample_id, source in sources.items():
        output = prepared[sample_id] / f"{sample_id}{Path(source).suffix}"
        if sample_id == "DMSO_1":
            # The only unsorted fixture is sorted into the work directory.
            assert output.is_file() and not output.is_symlink(), output
        else:
            assert output.is_symlink(), f"{output} is a copy, not a link"
            assert output.resolve() == (FIXTURES / source).resolve(), output
    for process in ("INFER_STRANDEDNESS", "SCAN_ALIGNMENT"):
        for sample_id, directory in task_directories(trace_text, process).items():
            copies = [
                path.name
                for path in directory.iterdir()
                if path.suffix in {".bam", ".cram"} and not path.is_symlink()
            ]
            assert copies == [], f"{process} ({sample_id}) holds alignment copies: {copies}"


def check_input_checksums() -> None:
    rows = pd.read_csv(ROOT / "manifest" / "input_checksums.tsv", sep="\t", dtype=str)
    samples = pd.read_csv(FIXTURES / "samples.tsv", sep="\t", dtype=str)
    assert list(rows["role"]) == ["sample_sheet", "fasta", "annotation"] + ["alignment"] * len(
        samples
    )
    alignments = rows[rows["role"] == "alignment"]
    expected = [str((FIXTURES / name).resolve()) for name in samples["alignment"]]
    assert list(alignments["path"]) == expected
    for path, digest in zip(alignments["path"], alignments["sha256"], strict=True):
        assert digest == hashlib.sha256(Path(path).read_bytes()).hexdigest(), path


def main() -> None:
    params = yaml.safe_load((ROOT / "manifest" / "resolved_params.yaml").read_text())
    tables = {name: statistics_table(f"{name}.pacs.tsv.gz") for name in COMPARISONS}
    events = {name: statistics_table(f"{name}.events.tsv.gz") for name in COMPARISONS}
    treatment = tables["TreatmentA_vs_DMSO"]

    # chr1 PAC 350 has no reads in DMSO and is used in TreatmentA.
    gained = check_event(treatment, pac("chr1", "+", 350), "gained", "gene_plus")
    check_confirmed(gained, params, "gene_plus")
    check_interval(gained, params, "gene_plus")
    assert number(gained["delta_pau"]) >= 0.3
    assert gained["delta_pau_ci_low"] > 0
    assert gained["raw_control_counts"] == "DMSO_1=0,DMSO_2=0"
    assert gained["model_status"] == "drimseq_add_uniform"
    gained_events = events["TreatmentA_vs_DMSO"]
    is_gained = (gained_events["pac_id"] == pac("chr1", "+", 350)) & (
        gained_events["event_type"] == "gained"
    )
    assert int(is_gained.sum()) == 1

    # bg01's middle PAC is silent only in TreatmentA: lost against DMSO and
    # gained again in the nested Rescue comparison.
    lost = check_event(treatment, pac("chr2", "+", 1350), "lost", "bg01")
    check_confirmed(lost, params, "bg01 lost")
    regained = check_event(tables["Rescue_vs_TreatmentA"], pac("chr2", "+", 1350), "gained", "bg01")
    check_confirmed(regained, params, "bg01 regained")
    # bg02's middle PAC is used only in TreatmentB.
    check_confirmed(
        check_event(tables["TreatmentB_vs_Vehicle"], pac("chr2", "-", 2050), "gained", "bg02"),
        params,
        "bg02",
    )
    # bg09 gains an internal-priming-like PAC: detected but not confirmed.
    primed = check_event(treatment, pac("chr2", "+", 9350), "gained_candidate", "bg09")
    assert primed["confidence"] == "low"
    assert str(primed["internal_priming_flag"]).lower() == "true"
    # bg03-bg06 shift usage from the first toward the last PAC in TreatmentA.
    for index in range(3, 7):
        origin = 1000 * index
        strand = "+" if index % 2 else "-"
        first, last = (origin + 300, origin + 400) if strand == "+" else (origin, origin + 100)
        label = f"bg{index:02d}"
        for coordinate, expected in ((first, "decreased_usage"), (last, "increased_usage")):
            row = check_event(treatment, pac("chr2", strand, coordinate), expected, label)
            check_confirmed(row, params, label)

    # Exact nulls: bg07 (TreatmentB is three times Vehicle) and bg08 (3:2 in
    # every sample) show no change.
    null_checks = [("TreatmentB_vs_Vehicle", "bg07")] + [(name, "bg08") for name in COMPARISONS]
    for name, gene_id in null_checks:
        rows = tables[name][tables[name]["gene_id"] == gene_id]
        assert len(rows) >= 2, f"{name}: {gene_id} is missing"
        assert set(rows["event_type"]) == {"none"}, f"{name}: {gene_id} has events"
        smallest = rows["pvalue_gene"].min()
        assert (rows["pvalue_gene"] > 0.5).all(), f"{name}: {gene_id} gene p is {smallest}"
        assert (rows["delta_pau"].abs() < 0.01).all(), f"{name}: {gene_id} delta is not zero"

    # PAC-level FDRs exist wherever stageR confirms genes.
    designed = {
        "TreatmentA_vs_DMSO": ["gene_plus", "bg01", "bg03", "bg04", "bg05", "bg06", "bg09"],
        "TreatmentB_vs_Vehicle": ["bg02"],
        "Rescue_vs_TreatmentA": ["bg01"],
    }
    for name, table in tables.items():
        rows = table[table["gene_id"].isin(designed[name])]
        assert set(rows["gene_id"]) == set(designed[name]), f"{name}: designed genes are missing"
        assert rows["pac_fdr"].notna().all(), f"{name}: designed genes lack pac_fdr"
        screened = table[table["gene_fdr"] <= params["site_fdr"] / 2]
        assert screened["pac_fdr"].notna().all(), f"{name}: screened genes lack pac_fdr"
        assert table["pvalue_pac"].notna().mean() >= 0.95, f"{name}: too many PAC tests are missing"
        assert "gene_minus" not in set(table["gene_id"]), f"{name}: gene_minus was tested"
        assert len(table) >= 2
    assert "primary_pas_motif_rna" in treatment
    assert "AAUAAA" in set(treatment["primary_pas_motif_rna"].dropna())

    precision = statistics_table("gene_precision.tsv.gz")
    for family in FAMILIES:
        values = precision.loc[precision["family"] == family, "precision"]
        assert len(values) >= 25, f"{family}: only {len(values)} genes have a precision"
        assert ((values > 0) & np.isfinite(values)).mean() >= 0.95, f"{family}: precision missing"

    # Raw counts equal the reads written into the fixture, PAC by PAC.
    counts = pd.read_csv(ROOT / "counts" / "pac_counts.tsv.gz", sep="\t")
    expected = pd.read_csv(FIXTURES / "expected_pac_counts.tsv", sep="\t")
    samples = [column for column in expected if column not in {"gene_id", "pac_id"}]
    observed = counts.set_index("pac_id").sort_index()
    wanted = expected.set_index("pac_id").sort_index()
    assert list(observed.index) == list(wanted.index), "atlas PACs differ from the fixture"
    assert (observed[samples].astype(int) == wanted[samples].astype(int)).all().all()
    assert (observed["gene_id"] == wanted["gene_id"]).all()

    pau = pd.read_csv(ROOT / "counts" / "observed_pau.tsv.gz", sep="\t")
    positive = pau[pau["gene_total"] > 0]
    sums = positive.groupby(["gene_id", "sample_id"])["pau"].sum().to_numpy()
    assert np.allclose(sums, 1)

    report = ROOT / "report" / "index.html"
    assert report.stat().st_size > 10000
    report_text = report.read_text()
    assert "PAU sample correlation" in report_text
    assert "PAC-level p-value distribution" in report_text
    assert "primary_pas_motif_rna" in report_text
    assert "AAUAAA" in report_text
    finite_intervals = sum(
        int((table["delta_pau_ci_low"].notna() & table["delta_pau_ci_high"].notna()).sum())
        for table in tables.values()
    )
    assert f"<strong>{finite_intervals:,}</strong>Bootstrap intervals" in report_text

    traces = list(REPOSITORY.glob("trace-*.txt"))
    assert traces
    trace_text = max(traces, key=lambda path: path.stat().st_mtime).read_text()
    check_alignment_handling(trace_text)
    check_input_checksums()
    assert trace_text.count("PACUSAGE:STATISTICS:FIT_USAGE_MODEL") == 3
    assert trace_text.count("PACUSAGE:STATISTICS:FINALIZE_USAGE_MODEL") == 3
    assert trace_text.count("PACUSAGE:STATISTICS:MERGE_USAGE_MODELS") == 1
    tags = re.findall(r"BOOTSTRAP_USAGE_INTERVALS \(([^:)]+):", trace_text)
    for family in FAMILIES:
        assert tags.count(family) >= 1, f"{family} has no bootstrap task"
    assert tags.count("DMSO") >= 2, "the DMSO family was not scattered into several batches"


if __name__ == "__main__":
    main()
