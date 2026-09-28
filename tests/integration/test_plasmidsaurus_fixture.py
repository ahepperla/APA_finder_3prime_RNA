"""The Plasmidsaurus-like fixture calibrates as proximal-tag and yields its PACs.

Its reads end 20-280 nt upstream of their PACs (see
tests/fixtures/build_plasmidsaurus_fixture.py). The last test calibrates the
exact-boundary fixture under the Plasmidsaurus profile instead: those reads end
on PACs 50 nt apart, so the kernel has spikes at 0, 50, and 100 nt and the
calibration warns.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pacusage.cli import main
from pacusage.tableio import read_tsv

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "plasmidsaurus"
EXACT_FIXTURE = FIXTURE.parent


def _calibrate(
    work: Path, reference: Path, annotation: Path, samples: list[dict], resolved_params
) -> Path:
    """Calibrate SE forward read_3p samples as the Plasmidsaurus profile would."""
    # A link keeps htslib from writing an index beside the fixture FASTA.
    fasta = work / "genome.fa"
    fasta.symlink_to(reference)
    sheet = work / "normalized_samples.tsv"
    sheet.write_text(
        "sample_id\tcondition\n"
        + "".join(f"{row['sample_id']}\t{row['condition']}\n" for row in samples)
    )
    params = work / "resolved_params.yaml"
    # The test profile's value for these small fixtures.
    resolved_params(params, calibration_min_genes=20)

    transcript_ends = work / "calibration_transcript_ends.tsv"
    assert main([
        "calibration-reference",
        "--annotation", str(annotation),
        "--output", str(transcript_ends),
    ]) == 0

    resolutions = []
    calibrations = []
    for row in samples:
        resolution = work / f"{row['sample_id']}.resolution.json"
        resolution.write_text(
            json.dumps(
                {
                    "sample_id": row["sample_id"],
                    "library_profile": "plasmidsaurus_3prime",
                    "layout": "SE",
                    "strandedness": "forward",
                    "evidence_source": "read_3p",
                    "endpoint_model": "proximal_tag",
                }
            )
        )
        resolutions.append(str(resolution))

        calibration = work / f"{row['sample_id']}.calibration.json"
        assert main([
            "scan-alignment",
            "--sample-id", row["sample_id"],
            "--alignment", row["alignment"],
            "--resolution", str(resolution),
            "--reference", str(fasta),
            "--transcript-ends", str(transcript_ends),
            "--params", str(params),
            "--calibration", str(calibration),
            "--output-prefix", str(work / f"{row['sample_id']}.scan"),
        ]) == 0
        calibrations.append(str(calibration))

    assert main([
        "aggregate-calibration",
        "--samples", str(sheet),
        "--calibrations", *calibrations,
        "--params", str(params),
        "--output", str(work / "library_calibration.tsv"),
        "--kernel", str(work / "kernel.tsv"),
        "--resolution", str(work / "run.json"),
        "--kernel-diagnostics", str(work / "kernel_diagnostics.tsv"),
    ]) == 0

    return params


@pytest.fixture(scope="module")
def discovery(tmp_path_factory: pytest.TempPathFactory, resolved_params) -> dict[str, object]:
    """Calibrate the fixture's samples and discover PACs as the pipeline does."""
    work = tmp_path_factory.mktemp("plasmidsaurus")
    samples = [
        {**row, "alignment": str(FIXTURE / row["alignment"])}
        for row in read_tsv(FIXTURE / "samples.tsv")
    ]
    params = _calibrate(
        work, FIXTURE / "genome.fa", FIXTURE / "genes.gtf", samples, resolved_params
    )
    sheet = work / "normalized_samples.tsv"

    evidence, continuations = [], []
    for row in samples:
        sample_id = row["sample_id"]
        scan_manifest = work / f"{sample_id}.scan.filtering.json"

        assert main([
            "extract-evidence",
            "--sample-id", sample_id,
            "--scan", str(scan_manifest),
            "--resolution", str(work / "run.json"),
            "--params", str(params),
            "--tsv", str(work / f"{sample_id}.tsv.gz"),
            "--parquet", str(work / f"{sample_id}.parquet"),
            "--plus-track", str(work / f"{sample_id}.plus.bedgraph"),
            "--minus-track", str(work / f"{sample_id}.minus.bedgraph"),
            "--splice-continuations", str(work / f"{sample_id}.splice.tsv.gz"),
            "--qc", str(work / f"{sample_id}.extract_qc.tsv"),
        ]) == 0

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
        "diagnostics": read_tsv(work / "kernel_diagnostics.tsv"),
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


def test_fixture_kernel_passes_the_calibration_warning(discovery) -> None:
    [row] = discovery["diagnostics"]
    assert (row["status"], row["kernel_modes"], row["reason"]) == ("ok", "1", "")
    # The diagnostics describe the kernel discovery used.
    resolution = discovery["summary"]["minimum_resolvable_separation"]
    assert row["minimum_resolvable_separation"] == resolution


def test_exact_boundary_reads_calibrated_as_plasmidsaurus_warn(
    tmp_path: Path, resolved_params
) -> None:
    # The original failure: reads that end on PACs 50 nt apart calibrate to a
    # kernel of PAC spacings, which now warns before discovery finds nothing.
    conditions = {"DMSO_1": "DMSO", "DMSO_2": "DMSO", "TRA_1": "TreatmentA", "TRA_2": "TreatmentA"}
    samples = [
        {
            "sample_id": sample_id,
            "condition": condition,
            "alignment": str(EXACT_FIXTURE / f"{sample_id}.bam"),
        }
        for sample_id, condition in conditions.items()
    ]
    _calibrate(
        tmp_path, EXACT_FIXTURE / "genome.fa", EXACT_FIXTURE / "genes.gtf", samples, resolved_params
    )
    [row] = read_tsv(tmp_path / "kernel_diagnostics.tsv")
    assert (row["status"], row["kernel_modes"], row["minimum_resolvable_separation"]) == (
        "warning",
        "3",
        "2",
    )
    assert row["reason"].startswith(
        "the pooled kernel has 3 separated modes; the samples' median central interval"
    )
