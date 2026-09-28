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

from pacusage.cli import _kmer_rows
from pacusage.models import IDENTITY_COLUMNS

REPOSITORY = Path(__file__).resolve().parents[2]
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
    # One scan per sample reads each alignment once for every evidence source.
    assert trace_text.count("PACUSAGE:DISCOVERY:SCAN_ALIGNMENT") == 10
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


def expected_gene_name(gene_id: str) -> str:
    """The fixture annotation's names, given on its gene lines only."""
    return {"gene_plus": "GenePlus", "gene_minus": "GeneMinus"}.get(gene_id, gene_id.upper())


def is_true(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "t", "1"})


def check_output_layout(tables: dict[str, pd.DataFrame], params: dict) -> None:
    """Every PAC table starts with one identity block that locates the PAC as
    the atlas BED does; gene names come from the annotation's gene lines."""
    def read(path: Path) -> pd.DataFrame:
        # Blank names stay blank rather than becoming NaN.
        return pd.read_csv(path, sep="\t", keep_default_na=False)

    atlas = read(ROOT / "atlas" / "pacs.v1.metadata.tsv.gz")
    bed = pd.read_csv(
        ROOT / "atlas" / "pacs.v1.bed.gz",
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "pac_id", "score", "strand"],
    ).set_index("pac_id")
    assert (bed["score"] <= 1000).all(), "BED scores exceed 1000"
    statistics = ROOT / "statistics"
    located = {
        "atlas": atlas,
        "pac_motifs": read(ROOT / "motifs" / "pac_motifs.tsv.gz"),
        "pac_counts": read(ROOT / "counts" / "pac_counts.tsv.gz"),
        "fitted_pau": read(statistics / "fitted_pau.tsv.gz"),
        **{f"{name}.pacs": read(statistics / f"{name}.pacs.tsv.gz") for name in tables},
        **{
            f"{family}.filtering": read(statistics / f"{family}.statistical_filtering.tsv.gz")
            for family in FAMILIES
        },
    }
    for name, table in located.items():
        assert list(table.columns[: len(IDENTITY_COLUMNS)]) == IDENTITY_COLUMNS, name
        rows = table.drop_duplicates("pac_id").set_index("pac_id")
        interval = bed.loc[rows.index]
        assert (rows["chrom"] == interval["chrom"]).all(), f"{name}: chrom differs from the BED"
        assert (rows["start"] == interval["start"]).all(), f"{name}: start differs from the BED"
        assert (rows["end"] == interval["end"]).all(), f"{name}: end differs from the BED"
        loci = rows["chrom"] + ":" + (rows["start"] + 1).astype(str) + "-" + rows["end"].astype(str)
        assert (rows["locus"] == loci).all(), f"{name}: locus is not the 1-based interval"
        named = rows[rows["gene_id"].astype(str).str.fullmatch(r"[^,]+")]
        assert (named["gene_name"] == named["gene_id"].map(expected_gene_name)).all(), name
    removed = {"feature_id", "site_class", "family", "contig", "region_start", "region_end"}
    for name, table in tables.items():
        assert not removed & set(table.columns), f"{name}: {sorted(removed & set(table.columns))}"
    assert not {"contig", "region_start", "region_end", "resolution_group"} & set(atlas.columns)
    genes = statistics_table("TreatmentA_vs_DMSO.genes.tsv.gz")
    assert list(genes.columns[:2]) == ["gene_id", "gene_name"]
    assert (genes["gene_name"] == genes["gene_id"].map(expected_gene_name)).all()

    # Descriptive labels are given only in genes that pass the gene screen.
    descriptive = {"dominant_switch", "complexity_gain", "complexity_loss"}
    for name, table in tables.items():
        labelled = table[table["event_type"].isin(descriptive)]
        assert (labelled["gene_fdr"] <= params["gene_fdr"]).all(), f"{name}: unscreened labels"

    # The k-mer background holds only PACs tested in the comparison: the
    # published tables equal a rerun on the atlas restricted that way, and the
    # restriction removes at least one untested atlas PAC from an event gene.
    usable = atlas[
        ~is_true(atlas["known_rescue_only"])
        & ~is_true(atlas["ambiguous_gene_assignment"])
        & (atlas["gene_id"].astype(str) != "")
    ]
    restricted_somewhere = False
    for name, table in tables.items():
        selected = set(table.loc[table["event_type"].isin(["gained", "increased_usage"]), "pac_id"])
        tested = usable[usable["pac_id"].isin(set(table["pac_id"]))]
        event_genes = set(tested.loc[tested["pac_id"].isin(selected), "gene_id"])
        in_event_genes = usable["gene_id"].isin(event_genes)
        untested = usable[in_event_genes & ~usable["pac_id"].isin(tested["pac_id"])]
        restricted_somewhere |= len(untested) > 0
        expected = {
            row["kmer"]: (row["event_pacs"], row["background_pacs"], row["informative_genes"])
            for row in _kmer_rows(tested, selected, int(params["motif_kmer_length"]))
        }
        published = pd.read_csv(ROOT / "motifs" / f"{name}.kmer_enrichment.tsv.gz", sep="\t")
        observed = {
            row.kmer: (row.event_pacs, row.background_pacs, row.informative_genes)
            for row in published.itertuples()
        }
        assert observed == expected, f"{name}: k-mer background differs from the tested PACs"
    assert restricted_somewhere, "no comparison exercises the tested-PAC background"

    # Motif preference is also reported per motif class.
    class_tables = sorted((ROOT / "motifs").glob("*.preference_class.tsv.gz"))
    assert len(class_tables) == len(COMPARISONS), [path.name for path in class_tables]
    for path in class_tables:
        columns = list(pd.read_csv(path, sep="\t", nrows=0).columns)
        assert columns[:2] == ["primary_motif_class", "condition"], (path.name, columns)
        assert "primary_pas_motif_rna" not in columns, path.name
    assert params["excluded_contigs"] == ["chrM", "MT", "chrMT"], params["excluded_contigs"]


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

    # Every PAC gets one filtering row per family, tested exactly when it
    # appears in that family's comparisons.
    assert not (ROOT / "qc" / "statistical_filtering.tsv").exists()
    for family in FAMILIES:
        filtering = statistics_table(f"{family}.statistical_filtering.tsv.gz")
        assert sorted(filtering["pac_id"]) == sorted(counts["pac_id"]), family
        tested = set(filtering.loc[filtering["tested"], "pac_id"])
        for name in COMPARISONS:
            if name.endswith(f"_vs_{family}"):
                assert set(tables[name]["pac_id"]) == tested, name
        assert (filtering.loc[~filtering["tested"], "reason"].str.len() > 0).all(), family
    assert not list((ROOT / "statistics").glob("*.preference*")), "preference tables in statistics/"
    assert len(list((ROOT / "motifs").glob("*.preference.tsv.gz"))) == len(COMPARISONS)

    versions = pd.read_csv(ROOT / "manifest" / "software_versions.tsv", sep="\t")
    assert {"pacusage", "pysam", "R", "DRIMSeq", "stageR", "limma"} <= set(versions["software"])

    pau = pd.read_csv(ROOT / "counts" / "observed_pau.tsv.gz", sep="\t")
    positive = pau[pau["gene_total"] > 0]
    sums = positive.groupby(["gene_id", "sample_id"])["pau"].sum().to_numpy()
    assert np.allclose(sums, 1)

    # Exact-boundary discovery does not assign reads through the kernel, so the
    # kernel is described but not judged.
    kernel = pd.read_csv(ROOT / "qc" / "calibration_kernel_diagnostics.tsv", sep="\t")
    assert list(kernel["endpoint_model"]) == ["exact_boundary"], list(kernel["endpoint_model"])
    assert list(kernel["status"]) == ["not_applicable"], list(kernel["status"])

    report = ROOT / "report" / "index.html"
    assert report.stat().st_size > 10000
    report_text = report.read_text()
    assert "<h2>Calibration kernel</h2>" in report_text
    assert "Calibration warning" not in report_text
    assert "PAU sample correlation" in report_text
    assert "<h2>Top gene usage profiles</h2>" in report_text
    assert "<h2>Statistical filtering</h2>" in report_text
    assert "PAC-level p-value distribution" in report_text
    assert "primary_pas_motif_rna" in report_text
    assert "AAUAAA" in report_text
    finite_intervals = sum(
        int((table["delta_pau_ci_low"].notna() & table["delta_pau_ci_high"].notna()).sum())
        for table in tables.values()
    )
    assert f"<strong>{finite_intervals:,}</strong>Bootstrap intervals" in report_text

    check_output_layout(tables, params)

    trace_text = (ROOT / "pipeline_info" / "execution_trace.txt").read_text()
    check_alignment_handling(trace_text)
    check_input_checksums()
    assert trace_text.count("PACUSAGE:STATISTICS:FIT_USAGE_MODEL") == 3
    assert trace_text.count("PACUSAGE:STATISTICS:FINALIZE_USAGE_MODEL") == 3
    assert trace_text.count("PACUSAGE:STATISTICS:MERGE_USAGE_MODELS") == 1
    tags = re.findall(r"BOOTSTRAP_USAGE_INTERVALS \(([^:)]+):", trace_text)
    for family in FAMILIES:
        assert tags.count(family) >= 1, f"{family} has no bootstrap task"
    assert tags.count("DMSO") >= 2, "the DMSO family was not scattered into several batches"
    # Each batch file becomes one bootstrap task, including a family's only batch.
    fits = task_directories(trace_text, "FIT_USAGE_MODEL")
    batches = [len(list(path.glob("family-*/bootstrap-batches/*.rds"))) for path in fits.values()]
    assert sum(batches) == len(tags), f"{sum(batches)} batch files, {len(tags)} bootstrap tasks"
    assert 1 in batches, "no family wrote a single bootstrap batch"


if __name__ == "__main__":
    main()
