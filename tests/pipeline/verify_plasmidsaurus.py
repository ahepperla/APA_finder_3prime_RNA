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
    assert int(discovery["rejected_candidates"]) == 0, discovery["rejected_candidates"]
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
    # Every multi-PAC gene is tested, and single-PAC genes cannot be.
    tested = set(pacs.loc[pacs["design"] != "calibration", "pac_id"])
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
    assert (nulls["pvalue_gene"] > 0.5).all(), list(nulls["pvalue_gene"])

    report = root / "report" / "index.html"
    assert report.stat().st_size > 10000
    report_text = report.read_text()
    assert "<h2>Calibration kernel</h2>" in report_text
    assert "Calibration warning" not in report_text
    print(f"Plasmidsaurus-like run verified: {len(atlas)} proximal-tag PACs.")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
