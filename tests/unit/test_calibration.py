import json
from pathlib import Path

import numpy as np

from pacusage.calibration import (
    CalibrationMetrics,
    calculate_metrics,
    classify_metrics,
    empirical_kernel,
    kernel_correlations,
    kernel_overlap,
    observation_offsets,
)
from pacusage.cli import main
from pacusage.models import EvidenceObservation
from pacusage.parameters import normalize_parameters


def params() -> dict:
    return normalize_parameters({"input": "x", "assembly": "x", "fasta": "x", "gtf": "x"})


def test_empirical_kernel_is_normalized_and_deterministic() -> None:
    first = empirical_kernel([0, 0, 1, -1], -5, 5)
    second = empirical_kernel([0, 0, 1, -1], -5, 5)
    assert np.allclose(first, second)
    assert np.isclose(first.sum(), 1)
    assert kernel_overlap(first, 0) == 1


def test_exact_and_proximal_classification() -> None:
    values = params()
    values["calibration_min_genes"] = 1
    exact = CalibrationMetrics("s", "read_3p", 2, 20, 0, 0, -1, 1, 0.9, 0, 0.95)
    proximal = CalibrationMetrics("s", "read_3p", 2, 20, 50, 50, 30, 70, 0, 0, 0.95)
    assert classify_metrics(exact, values).classification == "exact"
    assert classify_metrics(proximal, values).classification == "proximal"


def test_leave_one_out_kernel_correlations() -> None:
    kernel = empirical_kernel([0, 0, 1], -3, 3)
    values = kernel_correlations([kernel, kernel.copy(), kernel.copy()])
    assert np.allclose(values, 1)


def test_calibration_offsets_remain_aggregated_for_high_read_counts() -> None:
    observations = [
        EvidenceObservation("sample", "chr1", "+", 100, 1_000_000, 250_000)
    ]
    offsets, genes, clips = observation_offsets(
        observations,
        {("chr1", "+"): [(100, "gene1")]},
        maximum_distance=10,
        minimum_count_per_gene=1,
    )
    assert offsets == {0: 1_000_000}
    assert genes == 1
    assert clips == 250_000
    metrics = calculate_metrics("sample", "read_3p", offsets, genes, clips, 0.05, 0.95)
    assert metrics.observations == 1_000_000
    assert metrics.median_offset == 0
    assert metrics.poly_a_clip_fraction == 0.25


def test_binary_search_preserves_first_gene_at_duplicate_end_coordinate() -> None:
    observations = [
        EvidenceObservation("sample", "chr1", "+", 110, 3),
        EvidenceObservation("sample", "chr1", "+", 300, 1),
    ]
    offsets, genes, _ = observation_offsets(
        observations,
        {("chr1", "+"): [(100, "gene1"), (100, "gene2"), (300, "gene1")]},
        maximum_distance=200,
        minimum_count_per_gene=4,
    )
    assert offsets == {-10: 3, 0: 1}
    assert genes == 1


def test_split_calibration_matches_legacy_outputs(tmp_path: Path, monkeypatch) -> None:
    samples = tmp_path / "samples.tsv"
    samples.write_text("sample_id\nsample_a\nsample_b\n")
    annotation = tmp_path / "genes.gtf"
    annotation.write_text(
        'chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1";\n'
    )
    reference = tmp_path / "genome.fa"
    reference.write_text(">chr1\n" + "A" * 200 + "\n")
    params_path = tmp_path / "params.json"
    params_path.write_text(
        json.dumps(
            {
                "input": str(samples),
                "assembly": "test",
                "fasta": str(reference),
                "gtf": str(annotation),
                "calibration_min_genes": 1,
                "pac_min_sample_count": 1,
            }
        )
    )

    resolutions = []
    alignments = []
    for sample_id in ("sample_a", "sample_b"):
        resolution = tmp_path / f"{sample_id}.resolution.json"
        resolution.write_text(
            json.dumps(
                {
                    "sample_id": sample_id,
                    "library_profile": "exact_boundary",
                    "layout": "SE",
                    "strandedness": "forward",
                    "evidence_source": "read_3p",
                    "endpoint_model": "exact_boundary",
                }
            )
        )
        resolutions.append(str(resolution))
        alignments.append(str(tmp_path / f"{sample_id}.bam"))

    def fake_extract_evidence(sample_id: str, *_args, **_kwargs):
        count = 10 if sample_id == "sample_a" else 12
        return [EvidenceObservation(sample_id, "chr1", "+", 100, count)], {}

    monkeypatch.setattr("pacusage.cli.extract_evidence", fake_extract_evidence)

    legacy = tmp_path / "legacy"
    legacy.mkdir()
    assert (
        main(
            [
                "calibrate",
                "--samples",
                str(samples),
                "--alignments",
                *alignments,
                "--resolutions",
                *resolutions,
                "--reference",
                str(reference),
                "--annotation",
                str(annotation),
                "--params",
                str(params_path),
                "--output",
                str(legacy / "library_calibration.tsv"),
                "--kernel",
                str(legacy / "calibration_kernel.tsv"),
                "--resolution",
                str(legacy / "library_resolution.json"),
            ]
        )
        == 0
    )

    split = tmp_path / "split"
    split.mkdir()
    transcript_ends = split / "calibration_transcript_ends.tsv"
    assert (
        main(
            [
                "calibration-reference",
                "--annotation",
                str(annotation),
                "--output",
                str(transcript_ends),
            ]
        )
        == 0
    )
    summaries = []
    for sample_id, alignment, resolution in zip(
        ("sample_a", "sample_b"), alignments, resolutions, strict=True
    ):
        summary = split / f"{sample_id}.calibration.json"
        assert (
            main(
                [
                    "calibrate-sample",
                    "--sample-id",
                    sample_id,
                    "--alignment",
                    alignment,
                    "--resolution",
                    resolution,
                    "--reference",
                    str(reference),
                    "--transcript-ends",
                    str(transcript_ends),
                    "--params",
                    str(params_path),
                    "--output",
                    str(summary),
                ]
            )
            == 0
        )
        summaries.append(str(summary))
    assert (
        main(
            [
                "aggregate-calibration",
                "--samples",
                str(samples),
                "--calibrations",
                *summaries,
                "--params",
                str(params_path),
                "--output",
                str(split / "library_calibration.tsv"),
                "--kernel",
                str(split / "calibration_kernel.tsv"),
                "--resolution",
                str(split / "library_resolution.json"),
            ]
        )
        == 0
    )

    for name in (
        "library_calibration.tsv",
        "calibration_kernel.tsv",
        "library_resolution.json",
    ):
        assert (split / name).read_bytes() == (legacy / name).read_bytes()
