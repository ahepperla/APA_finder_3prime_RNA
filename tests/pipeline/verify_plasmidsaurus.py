#!/usr/bin/env python3
"""Assertions for the Plasmidsaurus-like fixture run, which uses proximal tags.

Usage: verify_plasmidsaurus.py RESULTS_DIRECTORY
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from column_guide import check_column_guide

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "plasmidsaurus"
SAMPLES = ["DMSO_1", "DMSO_2", "TRA_1", "TRA_2"]


def main(root: Path) -> None:
    params = yaml.safe_load((root / "manifest" / "resolved_params.yaml").read_text())
    resolution = json.loads((root / "manifest" / "library_resolution.json").read_text())
    selected = (resolution["library_profile"], resolution["evidence_source"])
    assert selected == ("plasmidsaurus_3prime", "read_3p"), selected
    assert resolution["endpoint_model"] == "proximal_tag", resolution["endpoint_model"]
    calibration = pd.read_csv(root / "qc" / "library_calibration.tsv", sep="\t")
    assert list(calibration["classification"]) == ["proximal"] * len(SAMPLES)
    # Only the 24 single-PAC genes have one annotated end.
    assert list(calibration["calibration_genes"]) == [24] * len(SAMPLES)

    expected = pd.read_csv(FIXTURE / "expected_pacs.tsv", sep="\t")
    discovery = pd.read_csv(root / "qc" / "pac_discovery.tsv", sep="\t").iloc[0]
    assert discovery["endpoint_model"] == "proximal_tag"
    assert int(discovery["accepted_pacs"]) == len(expected), discovery["accepted_pacs"]
    # Only the designed artifacts are rejected: reads at an internal exon's
    # donor, and reads inside an exon that spliced reads continue from.
    wanted = pd.read_csv(FIXTURE / "expected_rejected.tsv", sep="\t")
    rejected = pd.read_csv(root / "atlas" / "rejected_candidates.tsv.gz", sep="\t")
    observed = list(rejected[["strand", "coordinate", "rejection_reason"]].itertuples(index=False))
    assert observed == list(wanted.itertuples(index=False)), observed
    assert int(discovery["rejected_candidates"]) == len(wanted), discovery["rejected_candidates"]
    assert int(discovery["internal_exon_end_rejected"]) == 1
    assert int(discovery["constitutive_readthrough_rejected"]) == 1
    assert set(rejected["supporting_conditions"]) == {"DMSO;TreatmentA"}
    resolution_nt = int(discovery["minimum_resolvable_separation"])
    assert resolution_nt > 1, "a proximal-tag atlas cannot resolve single nucleotides"
    # Read ends scatter around one site each, so the kernel has one mode and
    # calibration does not warn. Discovery used the same kernel.
    kernel = pd.read_csv(root / "qc" / "calibration_kernel_diagnostics.tsv", sep="\t").iloc[0]
    assert kernel["status"] == "ok", kernel["reason"]
    assert int(kernel["kernel_modes"]) == 1, kernel["kernel_modes"]
    assert int(kernel["minimum_resolvable_separation"]) == resolution_nt

    # Designed sites sit on the discovery bins, so each representative is its
    # site. A region spans the assay resolution and holds one designed PAC.
    atlas = pd.read_csv(root / "atlas" / "pacs.v1.metadata.tsv.gz", sep="\t")
    pacs = expected.merge(
        atlas, on=["strand", "coordinate"], suffixes=("", "_atlas"), validate="one_to_one"
    )
    assert len(pacs) == len(atlas) == len(expected)
    assert (pacs["gene_id"] == pacs["gene_id_atlas"]).all(), "PACs assigned to the wrong gene"
    assert (pacs["resolution_nt"] == resolution_nt).all()
    assert (pacs["gene_name"] == pacs["gene_id"].str.upper()).all(), "gene names are missing"
    # A proximal-tag PAC's start and end are its resolution region.
    assert ((atlas["end"] - atlas["start"]) == resolution_nt).all()
    loci = atlas["chrom"] + ":" + (atlas["start"] + 1).astype(str) + "-" + atlas["end"].astype(str)
    assert (atlas["locus"] == loci).all(), "locus is not the 1-based region"
    for row in atlas.itertuples():
        inside = expected[
            (expected["strand"] == row.strand)
            & (expected["coordinate"] >= row.start)
            & (expected["coordinate"] < row.end)
        ]
        assert len(inside) == 1, f"{row.pac_id}: region holds {len(inside)} designed PACs"

    # Proximal-tag BED records span the resolution region, not one nucleotide.
    bed = pd.read_csv(
        root / "atlas" / "pacs.v1.bed.gz",
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "pac_id", "score", "strand"],
    )
    regions = bed.merge(atlas, on="pac_id", suffixes=("", "_atlas"), validate="one_to_one")
    assert len(regions) == len(atlas)
    assert (regions["start"] == regions["start_atlas"]).all()
    assert (regions["end"] == regions["end_atlas"]).all()

    # PACs are 400 nt apart, beyond the kernel's reach, so every read is
    # assigned to its own PAC and raw counts equal the reads written.
    counts = pd.read_csv(root / "counts" / "pac_counts.tsv.gz", sep="\t")
    table = pacs[["pac_id", *SAMPLES]].merge(
        counts, on="pac_id", suffixes=("_written", ""), validate="one_to_one"
    )
    for sample_id in SAMPLES:
        assert (table[sample_id] == table[f"{sample_id}_written"]).all(), sample_id

    tests = pd.read_csv(root / "statistics" / "TreatmentA_vs_DMSO.pacs.tsv.gz", sep="\t")
    # Every multi-PAC gene is tested. Single-PAC genes, the calibration genes
    # and the artifact hosts, cannot be.
    tested = set(pacs.loc[~pacs["design"].isin(["calibration", "artifact_host"]), "pac_id"])
    assert set(tests["pac_id"]) == tested, sorted(set(tests["pac_id"]) ^ tested)
    tests = tests.merge(
        pacs[["pac_id", "design", "rank_from_distal"]], on="pac_id", validate="one_to_one"
    )
    fraction = params["dm_bootstrap_min_success_fraction"]
    minimum_successes = int(np.ceil(fraction * params["dm_bootstrap_replicates"]))
    # Shifted genes move from 0.75 to 0.25 distal usage in TreatmentA.
    shifted = tests[tests["design"] == "shifted"]
    assert shifted["gene_id"].nunique() == 3, sorted(shifted["gene_id"].unique())
    for row in shifted.itertuples():
        label = f"{row.gene_id} {row.pac_id}"
        assert row.gene_fdr <= params["gene_fdr"], f"{label}: gene_fdr {row.gene_fdr}"
        assert row.pac_fdr <= params["site_fdr"], f"{label}: pac_fdr {row.pac_fdr}"
        proximal = row.rank_from_distal == 1
        expected_event = "increased_usage" if proximal else "decreased_usage"
        assert row.event_type == expected_event, f"{label}: {row.event_type}"
        change = row.delta_pau if proximal else -row.delta_pau
        assert 0.4 <= change <= 0.6, f"{label}: delta_pau {row.delta_pau}"
        assert row.bootstrap_status == "ok", f"{label}: bootstrap {row.bootstrap_status}"
        assert row.delta_pau_ci_low <= row.delta_pau <= row.delta_pau_ci_high, label
        assert row.bootstrap_successes >= minimum_successes, label
    # Exact nulls use their PACs 3:2 in every sample.
    nulls = tests[tests["design"] == "exact_null"]
    assert nulls["gene_id"].nunique() == 2, sorted(nulls["gene_id"].unique())
    assert set(nulls["event_type"]) == {"none"}, set(nulls["event_type"])
    assert (nulls["delta_pau"].abs() < 0.01).all(), list(nulls["delta_pau"])
    assert (nulls["gene_pvalue"] > 0.5).all(), list(nulls["gene_pvalue"])

    # Shifted genes lose usage at their distal PAC, their 3' UTRs shorten, and
    # their most-used PAC switches from distal to proximal.
    distal = pd.read_csv(root / "figures" / "TreatmentA_vs_DMSO.distal_usage.tsv.gz", sep="\t")
    distal = distal.set_index("gene_id")
    shifted_genes = sorted(shifted["gene_id"].unique())
    expected_distal = sorted(shifted.loc[shifted["rank_from_distal"] == 0, "pac_id"])
    assert sorted(distal.loc[shifted_genes, "distal_pac_id"]) == expected_distal
    assert set(distal.loc[shifted_genes, "direction"]) == {"distal_down"}
    assert set(distal.loc[shifted_genes, "apa_pattern"]) == {"utr_shortening"}
    null_genes = sorted(nulls["gene_id"].unique())
    assert set(distal.loc[null_genes, "direction"]) == {"none"}
    assert set(distal.loc[null_genes, "apa_pattern"]) == {"none"}
    genes = pd.read_csv(root / "statistics" / "TreatmentA_vs_DMSO.genes.tsv.gz", sep="\t")
    switched = genes.loc[genes["dominant_switch"].astype(str).str.lower() == "true", "gene_id"]
    assert set(shifted_genes) <= set(switched), sorted(switched)
    figures = {path.name for path in (root / "figures").iterdir()}
    stems = [f"TreatmentA_vs_DMSO.{kind}" for kind in ("volcano", "distal_usage", "site_classes")]
    stems += [
        "event_counts", "apa_pattern_grid", "concordance_matrix", "concordance",
        "effect_vs_coverage", "pau_pca",
    ]
    expected_figures = {f"{stem}.{suffix}" for stem in stems for suffix in ("pdf", "png")}
    expected_figures |= {
        "TreatmentA_vs_DMSO.distal_usage.tsv.gz", "apa_patterns_by_comparison.tsv.gz",
        "concordance.tsv.gz", "pau_pca.tsv",
    }
    assert figures == expected_figures, sorted(figures ^ expected_figures)
    # One comparison: no pairs, and no gene shared between comparisons.
    assert pd.read_csv(root / "figures" / "concordance.tsv.gz", sep="\t").empty
    grid = pd.read_csv(root / "figures" / "apa_patterns_by_comparison.tsv.gz", sep="\t")
    assert list(grid.columns) == [
        "gene_id", "gene_name", "patterned_comparisons", "TreatmentA_vs_DMSO"
    ]
    assert set(grid.loc[grid["gene_id"].isin(shifted_genes), "patterned_comparisons"]) == {1}
    pca = pd.read_csv(root / "figures" / "pau_pca.tsv", sep="\t")
    assert len(pca) == len(pd.read_csv(root / "manifest" / "normalized_samples.tsv", sep="\t"))

    report = root / "report" / "index.html"
    assert report.stat().st_size > 10000
    report_text = report.read_text()
    assert "<h2>Read-end offset profile</h2>" in report_text
    assert "Calibration warning" not in report_text
    assert report_text.count("<img src='data:image/png;base64,") == len(stems)
    # run_nextflow.sh publishes the per-sample counts, so the column guide's
    # section for them is checked too.
    assert str(params["save_intermediates"]).lower() == "true", params["save_intermediates"]
    per_sample = sorted((root / "counts" / "per_sample").glob("*.pac_counts.tsv.gz"))
    assert [path.name for path in per_sample] == [f"{s}.pac_counts.tsv.gz" for s in SAMPLES]
    tables = check_column_guide(root)
    print(f"Plasmidsaurus-like run verified: {len(atlas)} proximal-tag PACs, {tables} tables.")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
