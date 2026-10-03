#!/usr/bin/env python3
"""Assertions for stable biological behavior of the generated test run."""

from __future__ import annotations

import base64
import hashlib
import itertools
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from column_guide import check_column_guide

from pacusage.cli import _kmer_rows
from pacusage.models import IDENTITY_COLUMNS

REPOSITORY = Path(__file__).resolve().parents[2]
ROOT = REPOSITORY / "results-test"
FIXTURES = REPOSITORY / "tests" / "fixtures"
COMPARISONS = ("TreatmentA_vs_DMSO", "TreatmentB_vs_Vehicle", "Rescue_vs_TreatmentA")
FAMILIES = ("DMSO", "Vehicle", "TreatmentA")
COMPARISON_FIGURES = ("volcano", "distal_usage", "shifts_by_gene_region")
SUMMARY_FIGURES = (
    "event_counts", "apa_pattern_grid", "concordance_matrix", "concordance",
    "effect_vs_coverage", "pau_pca",
)
FIGURE_TABLES = ("apa_patterns_by_comparison.tsv.gz", "concordance.tsv.gz", "pau_pca.tsv")
# Figure sizes in pixels at 200 dpi, from scripts/plot_usage_figures.R; the
# summaries grow with the three comparisons. The concordance matrix has two
# columns and two rows for their three pairs. The grid's height also grows
# with its genes, so grid_pixels() computes it. The shifts figure grows with
# the region pairs that shifts take in any comparison (three in the fixture).
# The PCA widens for its legend column, from the longest condition name
# (TreatmentA).
FIGURE_PIXELS = {
    "volcano": (1500, 1100),
    "distal_usage": (1400, 1150),
    "shifts_by_gene_region": (1300, 600),
    "event_counts": (1500, 1800),
    "effect_vs_coverage": (1500, 800),
    "concordance": (1100, 1050),
    "concordance_matrix": (1075, 1075),
    "pau_pca": (1350, 1100),
}
CONFIRMED_EVENTS = {"gained", "increased_usage", "lost", "decreased_usage"}
DISTAL_COLUMNS = [
    "gene_id", "gene_name", "condition", "control_condition", "direction", "apa_pattern",
    "distal_pac_id", "distal_locus", "distal_gene_region", "fitted_control_distal_pau",
    "fitted_treatment_distal_pau", "delta_distal_pau", "distal_event_type", "gene_fdr",
    "distal_pac_fdr", "tested_pacs",
]
DISTAL_DIRECTIONS = {
    "gained": "distal_up",
    "increased_usage": "distal_up",
    "lost": "distal_down",
    "decreased_usage": "distal_down",
    "gained_candidate": "distal_up_candidate",
    "lost_candidate": "distal_down_candidate",
}
# The genes the fixture designs; the others are seeded nulls, whose chance
# results are checked only against an independent recomputation.
DESIGNED_GENES = {"gene_plus", *(f"bg{index:02d}" for index in range(1, 10)), "ipa01", "ale01"}
# The designed genes whose distal PAC has a call; the other designed genes
# have none. bg01 loses its middle PAC in TreatmentA, which raises the distal
# share.
EXPECTED_DISTAL_CALLS = {
    "TreatmentA_vs_DMSO": {
        "bg01": "distal_up", "bg03": "distal_up", "bg05": "distal_up", "ale01": "distal_up",
        "bg04": "distal_down", "bg06": "distal_down", "gene_plus": "distal_down",
        "ipa01": "distal_down",
    },
    "TreatmentB_vs_Vehicle": {"bg02": "distal_down", "gene_plus": "distal_down"},
    "Rescue_vs_TreatmentA": {
        "bg04": "distal_up", "bg06": "distal_up", "gene_plus": "distal_up", "ipa01": "distal_up",
        "bg01": "distal_down", "bg03": "distal_down", "bg05": "distal_down",
        "ale01": "distal_down",
    },
}
# Designed genes with a dominant switch, more active PACs, and fewer active
# PACs; the other designed genes have none.
EXPECTED_GENE_EVENTS = {
    # ipa01's intronic PAC has fitted PAU 0.09 in DMSO and 0.57 in TreatmentA,
    # so at the default active_pac_min_pau of 0.10 it is active only in the
    # treatment.
    "TreatmentA_vs_DMSO": (
        {"bg01", "bg03", "bg04", "bg05", "bg06", "gene_plus", "ipa01", "ale01"},
        {"bg09", "gene_plus", "ipa01"},
        {"bg01"},
    ),
    "TreatmentB_vs_Vehicle": ({"gene_plus"}, {"bg02"}, set()),
    "Rescue_vs_TreatmentA": (
        {"bg01", "bg03", "bg04", "bg05", "bg06", "ipa01", "ale01"}, {"bg01"}, {"bg09"}
    ),
}
# Each designed gene's APA pattern; any designed gene not listed is "none".
EXPECTED_APA_PATTERNS = {
    "TreatmentA_vs_DMSO": {
        "gene_plus": "utr_shortening", "bg04": "utr_shortening", "bg06": "utr_shortening",
        "bg01": "utr_lengthening", "bg03": "utr_lengthening", "bg05": "utr_lengthening",
        "bg09": "unclassified_change", "ipa01": "intronic_gain", "ale01": "alternative_last_exon",
    },
    "TreatmentB_vs_Vehicle": {"gene_plus": "utr_shortening", "bg02": "utr_shortening"},
    "Rescue_vs_TreatmentA": {
        "gene_plus": "utr_lengthening", "bg04": "utr_lengthening", "bg06": "utr_lengthening",
        "bg01": "utr_shortening", "bg03": "utr_shortening", "bg05": "utr_shortening",
        "bg09": "unclassified_change", "ipa01": "intronic_loss", "ale01": "alternative_last_exon",
    },
}
# Each designed gene's shift (DESIGNED_SHIFT_COLUMNS): its direction, then the
# region and call of the PAC its usage moved from and of the one it moved to.
# bg09 has a confirmed call on one side only. Designed genes not listed have
# no shift.
EXPECTED_SHIFTS = {
    "TreatmentA_vs_DMSO": {
        "ale01": ("distal", "last_exon", "decreased_usage", "last_exon", "increased_usage"),
        "bg01": ("distal", "last_exon", "lost", "last_exon", "increased_usage"),
        "bg03": ("distal", "last_exon", "decreased_usage", "last_exon", "increased_usage"),
        "bg04": ("proximal", "last_exon", "decreased_usage", "last_exon", "increased_usage"),
        "bg05": ("distal", "last_exon", "decreased_usage", "last_exon", "increased_usage"),
        "bg06": ("proximal", "last_exon", "decreased_usage", "last_exon", "increased_usage"),
        "bg09": ("distal", "last_exon", "decreased_usage", "last_exon", "gained_candidate"),
        "gene_plus": ("proximal", "last_exon", "decreased_usage", "last_exon", "gained"),
        "ipa01": ("proximal", "last_exon", "decreased_usage", "intron", "increased_usage"),
    },
    "TreatmentB_vs_Vehicle": {
        "bg02": ("proximal", "last_exon", "decreased_usage", "last_exon", "gained"),
        "gene_plus": ("proximal", "last_exon", "decreased_usage", "last_exon", "increased_usage"),
    },
    "Rescue_vs_TreatmentA": {
        "ale01": ("proximal", "last_exon", "decreased_usage", "last_exon", "increased_usage"),
        "bg01": ("proximal", "last_exon", "decreased_usage", "last_exon", "gained"),
        "bg03": ("proximal", "last_exon", "decreased_usage", "last_exon", "increased_usage"),
        "bg04": ("distal", "last_exon", "decreased_usage", "last_exon", "increased_usage"),
        "bg05": ("proximal", "last_exon", "decreased_usage", "last_exon", "increased_usage"),
        "bg06": ("distal", "last_exon", "decreased_usage", "last_exon", "increased_usage"),
        "bg09": ("proximal", "last_exon", "lost_candidate", "last_exon", "increased_usage"),
        "gene_plus": ("distal", "last_exon", "decreased_usage", "last_exon", "increased_usage"),
        "ipa01": ("distal", "intron", "decreased_usage", "last_exon", "increased_usage"),
    },
}
PAC_EVENT_TYPES = {
    "gained", "lost", "increased_usage", "decreased_usage", "gained_candidate", "lost_candidate",
    "none",
}
GENE_REGIONS = ("last_exon", "internal_exon", "intron", "downstream_of_gene", "intergenic")
APA_PATTERN_CLASSES = (
    "intronic_gain", "intronic_loss", "alternative_last_exon", "utr_shortening",
    "utr_lengthening",
)
POTENTIAL_INTERNAL_PRIMING = "_potential_internal_priming"
APA_METRIC_COLUMNS = ("delta_intronic_share", "delta_utr_distal_share", "last_exon_switch")
SHIFT_COLUMNS = (
    "shift_direction", "shift_from_pac_id", "shift_from_gene_region", "shift_from_event_type",
    "shift_to_pac_id", "shift_to_gene_region", "shift_to_event_type",
)
# Each side of a shift names a PAC, its gene region, and its call.
SHIFT_FIELDS = ("pac_id", "gene_region", "event_type")
DESIGNED_SHIFT_COLUMNS = (
    "shift_direction", "shift_from_gene_region", "shift_from_event_type", "shift_to_gene_region",
    "shift_to_event_type",
)
APA_TOLERANCE = 1e-9

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


def png_size(path: Path) -> tuple[int, int]:
    data = path.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n" and data[12:16] == b"IHDR", path.name
    return int.from_bytes(data[16:20], "big"), int.from_bytes(data[20:24], "big")


def pick_pac(rows: pd.DataFrame, sign: int) -> pd.Series:
    """The row whose delta_pau goes furthest in the direction of sign, ties
    going to the lower pac_fdr (missing last), then the PAC ID."""
    rows = rows.assign(change=-sign * rows["delta_pau"].astype(float))
    return rows.sort_values(["change", "pac_fdr", "pac_id"], na_position="last").iloc[0]


def expected_shifts(pacs: pd.DataFrame, coordinates: dict) -> dict[str, dict]:
    """Each gene's shift, recomputed from its PAC rows: its confirmed gain and
    loss with the largest changes or, on a side without a confirmed call, its
    other PAC with the largest fitted change that way."""
    shifts = {}
    for gene_id, rows in pacs.groupby("gene_id"):
        up = rows["event_type"].isin({"gained", "increased_usage"})
        down = rows["event_type"].isin({"lost", "decreased_usage"})
        if up.any():
            to = pick_pac(rows[up], 1)
            others = rows[down] if down.any() else rows[rows["pac_id"] != to["pac_id"]]
            source = pick_pac(others, -1)
        elif down.any():
            source = pick_pac(rows[down], -1)
            to = pick_pac(rows[rows["pac_id"] != source["pac_id"]], 1)
        else:
            continue
        # Transcript orientation: larger is more 3'.
        position = {
            side: coordinates[row["pac_id"]] * (1 if row["strand"] == "+" else -1)
            for side, row in (("from", source), ("to", to))
        }
        shifts[gene_id] = {
            "shift_direction": "distal" if position["to"] > position["from"] else "proximal",
            **{f"shift_from_{field}": source[field] for field in SHIFT_FIELDS},
            **{f"shift_to_{field}": to[field] for field in SHIFT_FIELDS},
        }
    return shifts


def check_gene_shifts(tables: dict[str, pd.DataFrame]) -> None:
    """Exactly the genes with a confirmed call have a shift, each one follows
    from the gene's PAC rows, and the designed genes shift as designed."""
    atlas = pd.read_csv(
        ROOT / "atlas" / "pacs.v1.metadata.tsv.gz", sep="\t", usecols=["pac_id", "coordinate"]
    )
    coordinates = dict(zip(atlas["pac_id"], atlas["coordinate"], strict=True))
    for name, pacs in tables.items():
        genes = statistics_table(f"{name}.genes.tsv.gz").set_index("gene_id")
        expected = expected_shifts(pacs, coordinates)
        shifted = genes.index[genes["shift_direction"].notna()]
        assert sorted(shifted) == sorted(expected), f"{name}: shifted genes {sorted(shifted)}"
        for gene_id, fields in expected.items():
            observed = genes.loc[gene_id, list(SHIFT_COLUMNS)].to_dict()
            assert observed == fields, f"{name}: {gene_id} shift {observed}, not {fields}"
        unshifted = genes.loc[genes["shift_direction"].isna(), list(SHIFT_COLUMNS)]
        assert unshifted.isna().all().all(), f"{name}: a gene without a shift has shift fields"
        designed = {
            gene_id: tuple(genes.loc[gene_id, column] for column in DESIGNED_SHIFT_COLUMNS)
            for gene_id in sorted(DESIGNED_GENES & set(shifted))
        }
        assert designed == EXPECTED_SHIFTS[name], f"{name}: designed shifts {designed}"


def check_gene_events(tables: dict[str, pd.DataFrame], params: dict) -> None:
    """The genes tables flag every gene with a switch or a change in its
    number of active PACs, and the PAC tables carry only PAC calls."""
    for name, (switched, more, fewer) in EXPECTED_GENE_EVENTS.items():
        genes = statistics_table(f"{name}.genes.tsv.gz")
        assert list(genes.columns[:14]) == [
            "gene_id", "gene_name", "condition", "control_condition", "dominant_switch",
            "active_pacs_change", "apa_pattern",
        ] + list(SHIFT_COLUMNS), name
        assert set(genes["active_pacs_change"]) <= {"more", "fewer", "none"}, name
        flagged = {
            "switch": set(genes.loc[is_true(genes["dominant_switch"]), "gene_id"]),
            "more": set(genes.loc[genes["active_pacs_change"] == "more", "gene_id"]),
            "fewer": set(genes.loc[genes["active_pacs_change"] == "fewer", "gene_id"]),
        }
        kinds = zip(("switch", "more", "fewer"), (switched, more, fewer), strict=True)
        for kind, expected in kinds:
            observed = flagged[kind] & DESIGNED_GENES
            assert observed == expected, f"{name}: designed {kind} genes {sorted(observed)}"
        # Every gene's flags follow from its PAC rows, recomputed here.
        pacs = tables[name]
        first = pacs.drop_duplicates("gene_id").set_index("gene_id")
        fitted = pacs.groupby("gene_id")[["fitted_control_pau", "fitted_treatment_pau"]].apply(
            lambda rows: bool(np.isfinite(rows.to_numpy(dtype=float)).all())
        )
        screened = first["gene_fdr"] <= params["gene_fdr"]
        switch = screened & (first["dominant_pac_control"] != first["dominant_pac_treatment"])
        switch &= first["dominant_pac_control"].notna() & first["dominant_pac_treatment"].notna()
        # Active PACs have at least active_pac_min_pau of the gene's fitted
        # usage in a group.
        minimum = params["active_pac_min_pau"]
        for group in ("control", "treatment"):
            active = (pacs[f"fitted_{group}_pau"] >= minimum).groupby(pacs["gene_id"]).sum()
            recorded = first[f"{group}_active_pacs"]
            assert (recorded == active.reindex(first.index)).all(), f"{name}: {group} active"
        change = first["treatment_active_pacs"] - first["control_active_pacs"]
        change = change.where(screened & fitted.reindex(first.index), 0)
        assert flagged["switch"] == set(first.index[switch]), f"{name}: dominant_switch"
        assert flagged["more"] == set(first.index[change > 0]), f"{name}: more active PACs"
        assert flagged["fewer"] == set(first.index[change < 0]), f"{name}: fewer active PACs"
        assert set(pacs["event_type"]) <= PAC_EVENT_TYPES, f"{name}: {set(pacs['event_type'])}"


def truthy(values: pd.Series) -> np.ndarray:
    return values.astype(str).str.lower().isin(["true", "t", "1"]).to_numpy()


def withheld_calls(rows: pd.DataFrame, params: dict) -> tuple[np.ndarray, np.ndarray]:
    """Gains and losses on flagged PACs that the flag alone withheld: a
    candidate that is significant, stable, covered, and not exploratory."""
    none = np.zeros(len(rows), dtype=bool)
    if not params["potential_internal_priming_withheld_calls"]:
        return none, none
    events = rows["event_type"].to_numpy()
    gene_fdr = rows["gene_fdr"].to_numpy(dtype=float)
    pac_fdr = rows["pac_fdr"].to_numpy(dtype=float)
    trusted = (
        truthy(rows["internal_priming_flag"])
        & (gene_fdr <= params["gene_fdr"])
        & (pac_fdr <= params["site_fdr"])
        & ~truthy(rows["zero_boundary_unstable"])
        & ~truthy(rows["exploratory_insufficient_replicates"])
    )
    covered = params["min_gene_total"]
    gains = trusted & (events == "gained_candidate")
    gains &= rows["control_gene_total"].to_numpy(dtype=float) >= covered
    losses = trusted & (events == "lost_candidate")
    losses &= rows["treatment_gene_total"].to_numpy(dtype=float) >= covered
    return gains, losses


def expected_apa_patterns(pacs: pd.DataFrame, coordinates: dict, params: dict) -> dict:
    """APA patterns recomputed from a PAC table, independently of the R code:
    region shares of the fitted usage, gated by confirmed calls on PACs that
    are not low confidence. A pattern that holds only once flagged PACs'
    calls, withheld ones included, count takes the potential suffix."""
    threshold = params["apa_pattern_min_change"] - APA_TOLERANCE
    minimum = params["event_min_treatment_pau"] - APA_TOLERANCE
    results = {}
    for gene_id, rows in pacs.groupby("gene_id", sort=False):
        metrics = dict.fromkeys(APA_METRIC_COLUMNS, math.nan)
        control = rows["fitted_control_pau"].to_numpy(dtype=float)
        treatment = rows["fitted_treatment_pau"].to_numpy(dtype=float)
        if not (np.isfinite(control).all() and np.isfinite(treatment).all()):
            results[gene_id] = ("none", metrics)
            continue
        coordinate = rows["pac_id"].map(coordinates).to_numpy(dtype=float)
        oriented = -coordinate if rows["strand"].iloc[0] == "-" else coordinate
        regions = rows["gene_region"].to_numpy()
        upstream = np.isin(regions, ["intron", "internal_exon"])
        three_prime = np.flatnonzero(np.isin(regions, ["last_exon", "downstream_of_gene"]))
        events = rows["event_type"].to_numpy()
        confident = rows["confidence"].to_numpy() != "low"
        up = np.isin(events, ["gained", "increased_usage"])
        down = np.isin(events, ["lost", "decreased_usage"])
        flagged = truthy(rows["internal_priming_flag"])
        withheld_gains, withheld_losses = withheld_calls(rows, params)
        flagged_up = flagged & (up | withheld_gains)
        flagged_down = flagged & (down | withheld_losses)
        called = up | down | withheld_gains | withheld_losses
        delta_intronic = treatment[upstream].sum() - control[upstream].sum()
        metrics["delta_intronic_share"] = delta_intronic
        last_exon = rows["last_exon_locus"].fillna("").to_numpy()[three_prime]
        exons = [three_prime[last_exon == name] for name in dict.fromkeys(last_exon)]
        exon_control = [control[index].sum() for index in exons]
        exon_treatment = [treatment[index].sum() for index in exons]
        exon_delta = [
            after - before for before, after in zip(exon_control, exon_treatment, strict=True)
        ]
        if len(exons) >= 2:
            metrics["last_exon_switch"] = min(max(exon_delta), -min(exon_delta))
        main: np.ndarray = np.asarray([], dtype=int)
        distal = -1
        if exons:
            main_exon = min(
                range(len(exons)),
                key=lambda k: (-(exon_control[k] + exon_treatment[k]), -oriented[exons[k]].max()),
            )
            main = exons[main_exon]
            if (
                len(main) >= 2
                and exon_control[main_exon] >= minimum
                and exon_treatment[main_exon] >= minimum
            ):
                distal = int(main[np.argmax(oriented[main])])
                metrics["delta_utr_distal_share"] = (
                    treatment[distal] / exon_treatment[main_exon]
                    - control[distal] / exon_control[main_exon]
                )
        gene_fdr = rows["gene_fdr"].iloc[0]
        if not (np.isfinite(gene_fdr) and gene_fdr <= params["gene_fdr"]) or not called.any():
            results[gene_id] = ("none", metrics)
            continue
        gene = {
            "delta_intronic": delta_intronic, "upstream": upstream, "exons": exons,
            "exon_delta": exon_delta, "switch": metrics["last_exon_switch"],
            "utr": metrics["delta_utr_distal_share"], "control": control,
            "treatment": treatment, "main": main, "distal": distal,
        }
        found = supported_patterns(gene, up & confident, down & confident, threshold)
        wider = supported_patterns(
            gene, (up & confident) | flagged_up, (down & confident) | flagged_down, threshold
        )
        found += [pattern + POTENTIAL_INTERNAL_PRIMING for pattern in wider if pattern not in found]
        results[gene_id] = (";".join(found) if found else "unclassified_change", metrics)
    return results


def supported_patterns(
    gene: dict, gate_up: np.ndarray, gate_down: np.ndarray, threshold: float
) -> list[str]:
    """The patterns that the calls gate_up and gate_down support, in order."""
    found = []
    delta_intronic, upstream = gene["delta_intronic"], gene["upstream"]
    if delta_intronic >= threshold and (gate_up & upstream).any():
        found.append("intronic_gain")
    if delta_intronic <= -threshold and (gate_down & upstream).any():
        found.append("intronic_loss")
    switch = gene["switch"]
    if np.isfinite(switch) and switch >= threshold:
        gaining = gene["exons"][int(np.argmax(gene["exon_delta"]))]
        losing = gene["exons"][int(np.argmin(gene["exon_delta"]))]
        if gate_up[gaining].any() and gate_down[losing].any():
            found.append("alternative_last_exon")
    utr = gene["utr"]
    if np.isfinite(utr) and abs(utr) >= threshold:
        shifted = abs(delta_intronic) >= threshold or (np.isfinite(switch) and switch >= threshold)
        main, distal = gene["main"], gene["distal"]
        main_change = gene["treatment"][main].sum() - gene["control"][main].sum()
        use_up = not shifted or main_change < 0
        use_down = not shifted or main_change > 0
        proximal = np.setdiff1d(main, [distal])
        if utr < 0:
            if (use_up and gate_up[proximal].any()) or (use_down and gate_down[distal]):
                found.append("utr_shortening")
        elif (use_up and gate_up[distal]) or (use_down and gate_down[proximal].any()):
            found.append("utr_lengthening")
    return found


def check_apa_patterns(tables: dict[str, pd.DataFrame], params: dict) -> None:
    """Every gene's APA pattern and numbers match a recomputation, and the
    designed genes get the patterns their design implies. The atlas records
    the chr3 genes' last exons."""
    atlas = pd.read_csv(
        ROOT / "atlas" / "pacs.v1.metadata.tsv.gz", sep="\t", keep_default_na=False
    )
    last_exons = dict(zip(atlas["pac_id"], atlas["last_exon_locus"], strict=True))
    regions = dict(zip(atlas["pac_id"], atlas["gene_region"], strict=True))
    assert (regions[pac("chr3", "+", 1150)], last_exons[pac("chr3", "+", 1150)]) == ("intron", "")
    assert last_exons[pac("chr3", "+", 2400)] == "chr3:2201-2400"
    assert last_exons[pac("chr3", "+", 5000)] == "chr3:4801-5000"
    assert last_exons[pac("chr3", "+", 5600)] == "chr3:5401-5600"
    assert last_exons[pac("chr1", "+", 350)] == "chr1:101-400"
    in_last_exon = atlas["gene_region"] == "last_exon"
    assert (atlas.loc[in_last_exon, "last_exon_locus"] != "").all(), "a last-exon PAC lacks a locus"
    assert set(atlas["gene_region"]) <= set(GENE_REGIONS), set(atlas["gene_region"])
    coordinates = dict(zip(atlas["pac_id"], atlas["coordinate"].astype(float), strict=True))
    for name, table in tables.items():
        genes = statistics_table(f"{name}.genes.tsv.gz").set_index("gene_id")
        expected = expected_apa_patterns(table, coordinates, params)
        assert set(expected) == set(genes.index), f"{name}: genes differ from the PAC table"
        for gene_id, (pattern, metrics) in expected.items():
            row = genes.loc[gene_id]
            label = f"{name}: {gene_id}"
            assert row["apa_pattern"] == pattern, f"{label} is {row['apa_pattern']}, not {pattern}"
            for column, value in metrics.items():
                observed = number(row[column])
                same = (math.isnan(value) and math.isnan(observed)) or abs(observed - value) < 1e-9
                assert same, f"{name}: {gene_id} {column} is {observed}, not {value}"
        designed = {
            gene_id: genes.loc[gene_id, "apa_pattern"]
            for gene_id in sorted(DESIGNED_GENES & set(genes.index))
        }
        wanted = {gene_id: EXPECTED_APA_PATTERNS[name].get(gene_id, "none") for gene_id in designed}
        assert designed == wanted, f"{name}: designed patterns {designed}"
        valid = {
            *APA_PATTERN_CLASSES,
            *(pattern + POTENTIAL_INTERNAL_PRIMING for pattern in APA_PATTERN_CLASSES),
            "unclassified_change",
            "none",
        }
        for pattern in genes["apa_pattern"]:
            assert set(pattern.split(";")) <= valid, f"{name}: invalid pattern {pattern}"


def grid_pixels(figures: Path) -> tuple[int, int]:
    """The shared-genes grid: 3 + 0.75 inches per comparison wide, and 2.375
    inches, 0.875 for the longest label line ("vs TreatmentA"), and 0.125 per
    drawn gene high."""
    grid = pd.read_csv(figures / "apa_patterns_by_comparison.tsv.gz", sep="\t")
    rows = min(50, int((grid["patterned_comparisons"] >= 2).sum()))
    return 200 * 3 + 150 * len(COMPARISONS), round(200 * (2.375 + 0.875 + 0.125 * rows))


def check_figures(tables: dict[str, pd.DataFrame], trace_text: str, report_text: str) -> None:
    """Figures exist for every comparison, carry no dates, and the distal-usage
    table agrees with the PAC tables and the fixture's designed genes."""
    figures = ROOT / "figures"
    stems = [f"{name}.{kind}" for name in COMPARISONS for kind in COMPARISON_FIGURES]
    stems += list(SUMMARY_FIGURES)
    expected = {f"{stem}.{suffix}" for stem in stems for suffix in ("pdf", "png")}
    expected |= {f"{name}.distal_usage.tsv.gz" for name in COMPARISONS}
    expected |= set(FIGURE_TABLES)
    observed = {path.name for path in figures.iterdir()}
    assert observed == expected, sorted(observed ^ expected)
    assert trace_text.count("PACUSAGE:PLOT_FIGURES") == 1
    pixels = {**FIGURE_PIXELS, "apa_pattern_grid": grid_pixels(figures)}
    for stem in stems:
        pdf = (figures / f"{stem}.pdf").read_bytes()
        assert pdf.startswith(b"%PDF-"), stem
        for field in (b"/CreationDate", b"/ModDate", b"/Producer"):
            assert field not in pdf, f"{stem}.pdf records {field.decode()}"
        kind = stem.rsplit(".", 1)[-1]
        assert png_size(figures / f"{stem}.png") == pixels[kind], stem
        image = base64.b64encode((figures / f"{stem}.png").read_bytes()).decode("ascii")
        assert f"src='data:image/png;base64,{image}'" in report_text, f"{stem} is not in the report"
    assert report_text.count("<img src='data:image/png;base64,") == len(stems)

    atlas = pd.read_csv(
        ROOT / "atlas" / "pacs.v1.metadata.tsv.gz", sep="\t", usecols=["pac_id", "coordinate"]
    )
    coordinates = dict(zip(atlas["pac_id"], atlas["coordinate"], strict=True))
    for name, table in tables.items():
        distal = pd.read_csv(figures / f"{name}.distal_usage.tsv.gz", sep="\t")
        assert list(distal.columns) == DISTAL_COLUMNS, name
        ordered = distal.sort_values(["gene_fdr", "gene_id"], na_position="last", kind="stable")
        assert list(distal["gene_id"]) == list(ordered["gene_id"]), f"{name}: rows are not sorted"
        # The distal PAC, recomputed: the most 3' tested PAC in the last exon
        # or downstream of the gene.
        eligible = table[table["gene_region"].isin(["last_exon", "downstream_of_gene"])].copy()
        coordinate = eligible["pac_id"].map(coordinates)
        eligible["position"] = np.where(eligible["strand"] == "+", coordinate, -coordinate)
        chosen = eligible.sort_values(["gene_id", "position"]).groupby("gene_id").tail(1)
        chosen = chosen.set_index("gene_id")
        rows = distal.set_index("gene_id")
        assert sorted(rows.index) == sorted(chosen.index), name
        chosen = chosen.loc[rows.index]
        assert list(rows["distal_pac_id"]) == list(chosen["pac_id"]), name
        assert np.allclose(rows["delta_distal_pau"], chosen["delta_pau"]), name
        assert list(rows["tested_pacs"]) == list(table.groupby("gene_id").size().loc[rows.index])
        directions = chosen["event_type"].map(DISTAL_DIRECTIONS).fillna("none")
        assert list(rows["direction"]) == list(directions), name
        calls = {
            gene_id: direction
            for gene_id, direction in rows.loc[rows["direction"] != "none", "direction"].items()
            if gene_id in DESIGNED_GENES
        }
        assert calls == EXPECTED_DISTAL_CALLS[name], f"{name}: distal calls {calls}"
        genes = statistics_table(f"{name}.genes.tsv.gz").set_index("gene_id")
        assert list(rows["apa_pattern"]) == list(genes.loc[rows.index, "apa_pattern"]), name


def check_pattern_grid() -> None:
    """The grid table lists every tested gene's pattern in each comparison,
    recomputed from the genes tables, most shared first, then by best FDR."""
    stems = sorted(COMPARISONS)
    grid = pd.read_csv(
        ROOT / "figures" / "apa_patterns_by_comparison.tsv.gz", sep="\t", dtype=str,
        keep_default_na=False,
    )
    assert list(grid.columns) == ["gene_id", "gene_name", "patterned_comparisons", *stems]
    patterns: dict[str, dict[str, str]] = {}
    names: dict[str, str] = {}
    fdrs: dict[str, list[float]] = {}
    for stem in stems:
        genes = pd.read_csv(
            ROOT / "statistics" / f"{stem}.genes.tsv.gz", sep="\t", keep_default_na=False
        )
        patterns[stem] = dict(zip(genes["gene_id"], genes["apa_pattern"], strict=True))
        fdr = pd.to_numeric(genes["gene_fdr"].replace("", np.nan))
        columns = (genes["gene_id"], genes["gene_name"], fdr)
        for gene_id, gene_name, value in zip(*columns, strict=True):
            names.setdefault(gene_id, gene_name)
            fdrs.setdefault(gene_id, [])
            if not math.isnan(value):
                fdrs[gene_id].append(float(value))
    rows = []
    for gene_id, values in fdrs.items():
        cells = [patterns[stem].get(gene_id, "") for stem in stems]
        patterned = sum(cell not in ("", "none") for cell in cells)
        best = min(values) if values else math.inf
        key = (-patterned, not values, best, gene_id)
        rows.append((key, [gene_id, names[gene_id], str(patterned), *cells]))
    expected = [row for _, row in sorted(rows)]
    assert grid.values.tolist() == expected
    # The fixture's reversal genes change in TreatmentA and back in Rescue.
    shared = set(grid.loc[grid["patterned_comparisons"].astype(int) >= 2, "gene_id"])
    assert {"gene_plus", "bg01", "ipa01", "ale01"} <= shared, sorted(shared)


def check_concordance(tables: dict[str, pd.DataFrame]) -> None:
    """Each pair of comparisons' shared PACs, calls in both, and correlation,
    recomputed from the PAC tables. Rescue_vs_TreatmentA reverses
    TreatmentA_vs_DMSO and shares TreatmentA's estimate with the opposite
    sign, so that chained pair correlates strongly negatively."""
    summary = pd.read_csv(ROOT / "figures" / "concordance.tsv.gz", sep="\t")
    assert list(summary.columns) == [
        "comparison_a", "comparison_b", "relation", "shared_pacs", "called_in_both", "pearson_r"
    ]
    expected = []
    for first, second in itertools.combinations(sorted(COMPARISONS), 2):
        treatment_a, control_a = first.split("_vs_")
        treatment_b, control_b = second.split("_vs_")
        if control_a == control_b:
            relation = "shared_control"
        elif treatment_a == control_b or treatment_b == control_a:
            relation = "chained"
        else:
            relation = "unrelated"
        left = tables[first].set_index("pac_id")
        right = tables[second].set_index("pac_id")
        shared = left.index.intersection(right.index)
        x = left.loc[shared, "delta_pau"].to_numpy(dtype=float)
        y = right.loc[shared, "delta_pau"].to_numpy(dtype=float)
        kept = np.isfinite(x) & np.isfinite(y)
        both = (
            left.loc[shared, "event_type"].isin(CONFIRMED_EVENTS).to_numpy()
            & right.loc[shared, "event_type"].isin(CONFIRMED_EVENTS).to_numpy()
            & kept
        )
        r = float(np.corrcoef(x[kept], y[kept])[0, 1]) if kept.sum() >= 3 else math.nan
        expected.append((first, second, relation, int(kept.sum()), int(both.sum()), r))
    observed = list(summary.itertuples(index=False, name=None))
    assert [row[:5] for row in observed] == [row[:5] for row in expected], observed
    assert np.allclose([row[5] for row in observed], [row[5] for row in expected], rtol=1e-9)
    relations = {(row[0], row[1]): (row[2], row[5]) for row in observed}
    relation, r = relations[("Rescue_vs_TreatmentA", "TreatmentA_vs_DMSO")]
    assert relation == "chained" and r < -0.5, (relation, r)


def check_pau_pca(params: dict, report_text: str) -> None:
    """The PAU PCA, recomputed by the rule the report used to tabulate it,
    matches up to each component's sign, and the sign follows the figure
    script's rule. The report shows the figure where the table was."""
    pau = pd.read_csv(ROOT / "counts" / "observed_pau.tsv.gz", sep="\t")
    covered = pau.groupby("gene_id")["gene_total"].min() >= params["min_gene_total"]
    frame = pau[pau["gene_id"].isin(covered.index[covered])]
    matrix = frame.pivot_table(index="pac_id", columns="sample_id", values="pau").dropna()
    centered = matrix.T.to_numpy(dtype=float)
    centered = centered - centered.mean(axis=0, keepdims=True)
    u, singular, _ = np.linalg.svd(centered, full_matrices=False)
    pca = pd.read_csv(ROOT / "figures" / "pau_pca.tsv", sep="\t")
    samples = pd.read_csv(ROOT / "manifest" / "normalized_samples.tsv", sep="\t")
    assert list(pca["sample_id"]) == sorted(samples["sample_id"]) == list(matrix.columns)
    conditions = dict(zip(samples["sample_id"], samples["condition"], strict=True))
    assert list(pca["condition"]) == [conditions[sample] for sample in pca["sample_id"]]
    variance = singular**2 / np.sum(singular**2)
    for index, column in enumerate(("PC1", "PC2")):
        scores = pca[column].to_numpy(dtype=float)
        direct = u[:, index] * singular[index]
        sign = 1.0 if np.dot(scores, direct) >= 0 else -1.0
        assert np.allclose(scores, sign * direct, atol=1e-9), column
        # The first PAC by ID within a relative 1e-6 of the largest loading
        # has a positive loading.
        loading = centered.T @ scores
        largest = np.flatnonzero(np.abs(loading) >= np.abs(loading).max() * (1 - 1e-6))[0]
        assert loading[largest] > 0, column
        assert np.allclose(pca[f"pc{index + 1}_variance_fraction"], variance[index]), column
    # TreatmentA's designed changes dominate the first component.
    treated = pca["condition"] == "TreatmentA"
    assert (np.sign(pca.loc[treated, "PC1"]) != np.sign(pca.loc[~treated, "PC1"].mean())).all()
    start = report_text.index("<h2>PAU principal components</h2>")
    section = report_text[start : report_text.index("</section>", start)]
    assert report_text.index("<h2>PAU sample correlation</h2>") < start
    assert start < report_text.index("<h2>Treatment-control figures</h2>")
    image = base64.b64encode((ROOT / "figures" / "pau_pca.png").read_bytes()).decode("ascii")
    assert f"src='data:image/png;base64,{image}'" in section
    assert "<table" not in section


# Genes each comparison leaves untested for depth: off01 is silent in
# TreatmentA, and low01 has 2 reads per PAC there.
EXPECTED_WITHOUT_DEPTH = {
    "TreatmentA_vs_DMSO": {"off01": "turned_off", "low01": "too_low_in_treatment"},
    "TreatmentB_vs_Vehicle": {},
    "Rescue_vs_TreatmentA": {"off01": "turned_on", "low01": "too_low_in_control"},
}


def check_genes_without_depth(tables: dict[str, pd.DataFrame], counts: pd.DataFrame) -> None:
    """The comparisons skip the genes without depth, and the tables' CPM uses
    each sample's assigned reads, the column totals of the count table."""
    sheet = pd.read_csv(ROOT / "manifest" / "normalized_samples.tsv", sep="\t")
    library = counts[list(sheet["sample_id"])].sum()
    for sample_id in sheet["sample_id"]:
        quantification = pd.read_csv(ROOT / "qc" / f"{sample_id}.quantification.tsv", sep="\t")
        assigned = int(quantification["assigned_fragments"].iloc[0])
        assert assigned == int(library[sample_id]), sample_id
    for name, expected in EXPECTED_WITHOUT_DEPTH.items():
        skipped = statistics_table(f"{name}.genes_without_depth.tsv.gz")
        observed = dict(zip(skipped["gene_id"], skipped["depth_status"], strict=True))
        assert observed == expected, f"{name}: {observed}"
        assert not set(skipped["gene_id"]) & set(tables[name]["gene_id"]), name
        treatment, control = name.split("_vs_")
        for _, row in skipped.iterrows():
            gene = counts[counts["gene_id"] == row["gene_id"]]
            for group, condition in (("control", control), ("treatment", treatment)):
                ids = list(sheet.loc[sheet["condition"] == condition, "sample_id"])
                reads = gene[ids].sum()
                assert int(row[f"{group}_gene_total"]) == int(reads.sum()), (name, row["gene_id"])
                cpm = (reads / library[ids] * 1e6).mean()
                assert np.isclose(row[f"{group}_mean_cpm"], cpm), (name, row["gene_id"], group)


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
    # Names that version 0.4.0 replaced, by the tables that had them.
    atlas_names = {"assignment_class", "last_exon", "supporting_condition"}
    gene_names = {"pvalue", "lr", "df", "complexity_change"}
    renamed = {
        "atlas/pacs.v1.metadata.tsv.gz": atlas_names,
        "atlas/rejected_candidates.tsv.gz": atlas_names,
        "statistics/*.genes.tsv.gz": gene_names,
        "statistics/*.gene_omnibus.tsv.gz": gene_names,
        "motifs/*.tsv.gz": {"p_value"},
    }
    removed = {"feature_id", "site_class", "family", "contig", "region_start", "region_end"}
    removed |= atlas_names | gene_names | {
        "pvalue_pac", "pvalue_gene", "control_detected_complexity", "treatment_detected_complexity",
    }
    for name, table in tables.items():
        assert not removed & set(table.columns), f"{name}: {sorted(removed & set(table.columns))}"
    for pattern, names in renamed.items():
        paths = sorted(ROOT.glob(pattern))
        assert paths, pattern
        for path in paths:
            columns = set(pd.read_csv(path, sep="\t", nrows=0).columns)
            assert not names & columns, f"{path.name}: {sorted(names & columns)}"
    assert not list((ROOT / "statistics").glob("*.events.tsv.gz")), "an events table remains"
    assert not {"contig", "region_start", "region_end", "resolution_group"} & set(atlas.columns)
    genes = statistics_table("TreatmentA_vs_DMSO.genes.tsv.gz")
    assert list(genes.columns[:2]) == ["gene_id", "gene_name"]
    assert (genes["gene_name"] == genes["gene_id"].map(expected_gene_name)).all()

    # Gene-level events are given only in genes that pass the gene screen.
    for name in tables:
        genes = statistics_table(f"{name}.genes.tsv.gz")
        changed = genes["active_pacs_change"] != "none"
        described = genes[is_true(genes["dominant_switch"]) | changed]
        assert (described["gene_fdr"] <= params["gene_fdr"]).all(), f"{name}: unscreened events"

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
    calls = {name: statistics_table(f"{name}.calls.tsv.gz") for name in COMPARISONS}
    treatment = tables["TreatmentA_vs_DMSO"]

    # chr1 PAC 350 has no reads in DMSO and is used in TreatmentA.
    gained = check_event(treatment, pac("chr1", "+", 350), "gained", "gene_plus")
    check_confirmed(gained, params, "gene_plus")
    check_interval(gained, params, "gene_plus")
    assert number(gained["delta_pau"]) >= 0.3
    assert gained["delta_pau_ci_low"] > 0
    assert gained["raw_control_counts"] == "DMSO_1=0,DMSO_2=0"
    assert gained["model_status"] == "fitted_with_zero_count_stabilization"
    gained_calls = calls["TreatmentA_vs_DMSO"]
    is_gained = (gained_calls["pac_id"] == pac("chr1", "+", 350)) & (
        gained_calls["event_type"] == "gained"
    )
    assert int(is_gained.sum()) == 1
    # A calls table is its comparison's PAC rows with a call.
    for name, table in tables.items():
        called = table[table["event_type"] != "none"].reset_index(drop=True)
        pd.testing.assert_frame_equal(calls[name], called, obj=f"{name}.calls")

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
        smallest = rows["gene_pvalue"].min()
        assert (rows["gene_pvalue"] > 0.5).all(), f"{name}: {gene_id} gene p is {smallest}"
        assert (rows["delta_pau"].abs() < 0.01).all(), f"{name}: {gene_id} delta is not zero"

    # PAC-level FDRs exist wherever stageR confirms genes.
    designed = {
        "TreatmentA_vs_DMSO": [
            "gene_plus", "bg01", "bg03", "bg04", "bg05", "bg06", "bg09", "ipa01", "ale01"
        ],
        "TreatmentB_vs_Vehicle": ["bg02"],
        "Rescue_vs_TreatmentA": ["bg01", "ipa01", "ale01"],
    }
    for name, table in tables.items():
        rows = table[table["gene_id"].isin(designed[name])]
        assert set(rows["gene_id"]) == set(designed[name]), f"{name}: designed genes are missing"
        assert rows["pac_fdr"].notna().all(), f"{name}: designed genes lack pac_fdr"
        screened = table[table["gene_fdr"] <= params["site_fdr"] / 2]
        assert screened["pac_fdr"].notna().all(), f"{name}: screened genes lack pac_fdr"
        assert table["pac_pvalue"].notna().mean() >= 0.95, f"{name}: too many PAC tests are missing"
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
                # A comparison tests the family's PACs except those of genes
                # it lists as without depth.
                skipped = statistics_table(f"{name}.genes_without_depth.tsv.gz")
                listed = filtering["gene_id"].isin(skipped["gene_id"])
                without = set(filtering.loc[listed, "pac_id"])
                assert without <= tested, name
                assert set(tables[name]["pac_id"]) == tested - without, name
        assert (filtering.loc[~filtering["tested"], "reason"].str.len() > 0).all(), family
    assert not list((ROOT / "statistics").glob("*.preference*")), "preference tables in statistics/"
    assert len(list((ROOT / "motifs").glob("*.preference.tsv.gz"))) == len(COMPARISONS)

    versions = pd.read_csv(ROOT / "manifest" / "software_versions.tsv", sep="\t")
    assert {"pacusage", "pysam", "R", "DRIMSeq", "stageR", "limma", "ggplot2", "ggrepel"} <= set(
        versions["software"]
    )

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
    assert "<h2>Read-end offset profile</h2>" in report_text
    assert "Calibration warning" not in report_text
    assert "PAU sample correlation" in report_text
    assert "<h2>PAC usage in the most-read genes</h2>" in report_text
    assert "<h2>PACs left out of testing, and why</h2>" in report_text
    assert "<h2>Genes with the lowest gene FDR</h2>" in report_text
    for name in COMPARISONS:
        title = name.replace("_vs_", " vs ")
        assert f"<h2>{title}: PACs with a call</h2>" in report_text, name
        assert f"<h2>{title}: every tested PAC</h2>" in report_text, name
        assert f"<h2>{title}: genes not tested for depth</h2>" in report_text, name
    assert "docs/output_columns.md</a>" in report_text
    assert "PAC-level p-value distribution" in report_text
    assert "primary_pas_motif_rna" in report_text
    assert "AAUAAA" in report_text
    finite_intervals = sum(
        int((table["delta_pau_ci_low"].notna() & table["delta_pau_ci_high"].notna()).sum())
        for table in tables.values()
    )
    assert f"<strong>{finite_intervals:,}</strong>Bootstrap intervals" in report_text

    check_output_layout(tables, params)
    check_genes_without_depth(tables, counts)
    check_column_guide(ROOT)

    trace_text = (ROOT / "pipeline_info" / "execution_trace.txt").read_text()
    check_alignment_handling(trace_text)
    check_input_checksums()
    check_gene_events(tables, params)
    check_gene_shifts(tables)
    check_apa_patterns(tables, params)
    check_figures(tables, trace_text, report_text)
    check_pattern_grid()
    check_concordance(tables)
    check_pau_pca(params, report_text)
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
