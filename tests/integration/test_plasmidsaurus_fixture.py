"""The Plasmidsaurus-like fixture must calibrate as proximal and yield its PACs.

The fixture's sample sheet once reused the exact-boundary alignments, whose
reads end on each PAC. Calibrated against the annotated distal ends of genes
with PACs 50 nt apart, that kernel had spikes at 0, 50, and 100 nt. Every
endpoint then matched several peaks equally well, and proximal discovery
accepted no PAC.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pacusage.cli import main
from pacusage.evidence import (
    extract_evidence,
    extract_splice_continuations,
    write_evidence,
    write_splice_continuations,
)
from pacusage.tableio import read_tsv

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "plasmidsaurus"


@pytest.fixture(scope="module")
def discovery(tmp_path_factory: pytest.TempPathFactory) -> dict[str, object]:
    """Calibrate the fixture's samples and discover PACs as the pipeline does."""
    work = tmp_path_factory.mktemp("plasmidsaurus")
    # A link keeps htslib from writing an index beside the fixture FASTA.
    fasta = work / "genome.fa"
    fasta.symlink_to(FIXTURE / "genome.fa")
    annotation = FIXTURE / "genes.gtf"
    samples = read_tsv(FIXTURE / "samples.tsv")
    sheet = work / "normalized_samples.tsv"
    sheet.write_text(
        "sample_id\tcondition\n"
        + "".join(f"{row['sample_id']}\t{row['condition']}\n" for row in samples)
    )
    params = work / "params.json"
    params.write_text(
        json.dumps(
            {
                "input": str(FIXTURE / "samples.tsv"),
                "assembly": "synthetic",
                "fasta": str(fasta),
                "gtf": str(annotation),
                # The test profile's value for this small fixture.
                "calibration_min_genes": 20,
            }
        )
    )
    alignments = [str(FIXTURE / row["alignment"]) for row in samples]
    resolutions = []
    for row in samples:
        resolution = work / f"{row['sample_id']}.resolution.json"
        resolution.write_text(
            json.dumps(
                {
                    "sample_id": row["sample_id"],
                    "library_profile": row["library_profile"],
                    "layout": "SE",
                    "strandedness": "forward",
                    "evidence_source": "read_3p",
                    "endpoint_model": "proximal_tag",
                }
            )
        )
        resolutions.append(str(resolution))
    arguments = ["calibrate", "--samples", str(sheet), "--alignments", *alignments]
    arguments += ["--resolutions", *resolutions, "--reference", str(fasta)]
    arguments += ["--annotation", str(annotation), "--params", str(params)]
    arguments += ["--output", str(work / "library_calibration.tsv")]
    arguments += ["--kernel", str(work / "kernel.tsv"), "--resolution", str(work / "run.json")]
    assert main(arguments) == 0

    evidence, continuations = [], []
    for row, alignment in zip(samples, alignments, strict=True):
        sample_id = row["sample_id"]
        observations, _ = extract_evidence(sample_id, alignment, fasta, "SE", "forward", "read_3p")
        write_evidence(observations, work / f"{sample_id}.tsv.gz", work / f"{sample_id}.parquet")
        splices, _ = extract_splice_continuations(sample_id, alignment, fasta, "SE", "forward")
        write_splice_continuations(splices, work / f"{sample_id}.splice.tsv.gz")
        evidence.append(str(work / f"{sample_id}.parquet"))
        continuations.append(str(work / f"{sample_id}.splice.tsv.gz"))
    arguments = ["cluster", "--evidence", *evidence, "--splice-continuations", *continuations]
    arguments += ["--samples", str(sheet), "--resolution", str(work / "run.json")]
    arguments += ["--kernel", str(work / "kernel.tsv"), "--params", str(params)]
    arguments += ["--accepted", str(work / "accepted.tsv.gz")]
    arguments += ["--rejected", str(work / "rejected.tsv.gz"), "--qc", str(work / "discovery.tsv")]
    assert main(arguments) == 0
    return {
        "calibration": read_tsv(work / "library_calibration.tsv"),
        "resolution": json.loads((work / "run.json").read_text()),
        "accepted": read_tsv(work / "accepted.tsv.gz"),
        "rejected": read_tsv(work / "rejected.tsv.gz"),
        "summary": read_tsv(work / "discovery.tsv")[0],
    }


def test_fixture_calibrates_as_proximal_tag(discovery) -> None:
    assert discovery["resolution"]["endpoint_model"] == "proximal_tag"
    assert discovery["resolution"]["evidence_source"] == "read_3p"
    for metric in discovery["calibration"]:
        assert metric["classification"] == "proximal", metric["sample_id"]
        # Only the 24 single-PAC genes have one annotated end.
        assert int(metric["calibration_genes"]) == 24
        # Read ends sit 20-280 nt upstream, with the median and mode at 150.
        assert float(metric["median_offset"]) == 150.0
        assert float(metric["reproducibility"]) >= 0.99


def test_fixture_yields_exactly_its_designed_pacs(discovery) -> None:
    expected = read_tsv(FIXTURE / "expected_pacs.tsv")
    designed = [(row["strand"], int(row["coordinate"])) for row in expected]
    accepted = [
        (row["strand"], int(row["coordinate"]), int(row["region_start"]), int(row["region_end"]))
        for row in discovery["accepted"]
    ]
    assert discovery["rejected"] == []
    assert int(discovery["summary"]["accepted_pacs"]) == len(designed) == 64
    # A resolution group spans several nucleotides; it is not a cleavage site.
    assert int(discovery["summary"]["minimum_resolvable_separation"]) > 1
    # Sites sit on the discovery bins, so each representative is its site, and
    # each site lies inside its own region only.
    assert sorted((strand, coordinate) for strand, coordinate, _, _ in accepted) == sorted(designed)
    for strand, coordinate, start, end in accepted:
        inside = [site for site in designed if site[0] == strand and start <= site[1] < end]
        assert inside == [(strand, coordinate)]
