import json
from pathlib import Path

import numpy as np
import pytest
from alignment_builders import aligned_segment, write_alignment, write_reference

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


# Plus-strand genes end at 1000, 3000 and 5000; minus-strand genes at 7000 and
# 9000 (interbase). Each case gives one source a reproducible kernel and every
# other source kernels that do not overlap between the two samples, so source
# selection is unambiguous while every candidate source is still scanned.
PLUS_ENDS = (1000, 3000, 5000)
MINUS_ENDS = (7000, 9000)
CONTIGS = [("chr1", 12000)]


def calibration_annotation(path: Path) -> Path:
    lines = [
        f'chr1\ttest\texon\t{end - 999}\t{end}\t.\t+\t.\tgene_id "p{end}"; '
        f'transcript_id "tp{end}";\n'
        for end in PLUS_ENDS
    ] + [
        f'chr1\ttest\texon\t{end + 1}\t{end + 1000}\t.\t-\t.\tgene_id "m{end}"; '
        f'transcript_id "tm{end}";\n'
        for end in MINUS_ENDS
    ]
    path.write_text("".join(lines))
    return path


def single_end_calibration_reads(length: int, clip_shift: int, clip_count: int) -> list:
    """Reads ending at each gene end; clipped reads end ``clip_shift`` downstream."""
    reads = []
    for end in PLUS_ENDS:
        reads += [
            aligned_segment(f"p{end}_{index}", 0, 0, end - length, f"{length}M")
            for index in range(20)
        ]
        clipped_end = end + clip_shift
        reads += [
            aligned_segment(
                f"p{end}_clip{index}", 0, 0, clipped_end - length, f"{length}M10S",
                right="A" * 10,
            )
            for index in range(clip_count)
        ]
    for end in MINUS_ENDS:
        reads += [
            aligned_segment(f"m{end}_{index}", 16, 0, end, f"{length}M") for index in range(20)
        ]
        reads += [
            aligned_segment(
                f"m{end}_clip{index}", 16, 0, end - clip_shift, f"10S{length}M", left="T" * 10
            )
            for index in range(clip_count)
        ]
    return reads


def paired_end_calibration_reads(fragment: int) -> list:
    """Fragments whose 3' ends sit on gene ends; read 1 ends ``fragment`` upstream."""
    reads = []
    for end in PLUS_ENDS:
        for index in range(20):
            clip = {"right": "A" * 10} if index < 5 else {}
            cigar = "30M10S" if index < 5 else "30M"
            reads += [
                aligned_segment(f"p{end}_{index}", 99, 0, end - fragment, cigar, **clip),
                aligned_segment(f"p{end}_{index}", 147, 0, end - 30, "30M"),
            ]
    for end in MINUS_ENDS:
        for index in range(20):
            clip = {"left": "T" * 10} if index < 5 else {}
            cigar = "10S30M" if index < 5 else "30M"
            reads += [
                aligned_segment(f"m{end}_{index}", 83, 0, end + fragment - 30, cigar, **clip),
                aligned_segment(f"m{end}_{index}", 163, 0, end, "30M"),
            ]
    return reads


CALIBRATION_CASES = {
    # polyA_junction and read_3p agree in both samples.
    "exact_boundary_SE": (
        "exact_boundary",
        "SE",
        "exact_boundary",
        lambda sample: single_end_calibration_reads(30, 0, 10),
    ),
    # read_3p wins; read_5p and polyA_junction differ between samples.
    "generic_SE": (
        "generic_3prime",
        "SE",
        "auto",
        lambda sample: single_end_calibration_reads(
            30 if sample == "sample_a" else 60, 3 if sample == "sample_a" else -3, 5
        ),
    ),
    # fragment_3p wins; the read-level sources differ between samples.
    "generic_PE": (
        "generic_3prime",
        "PE",
        "auto",
        lambda sample: paired_end_calibration_reads(200 if sample == "sample_a" else 300),
    ),
}


@pytest.mark.parametrize("case", CALIBRATION_CASES)
def test_split_calibration_matches_legacy_outputs(tmp_path: Path, case: str) -> None:
    profile, layout, endpoint_model, reads = CALIBRATION_CASES[case]
    samples = tmp_path / "samples.tsv"
    samples.write_text("sample_id\nsample_a\nsample_b\n")
    annotation = calibration_annotation(tmp_path / "genes.gtf")
    reference = write_reference(tmp_path / "genome.fa", CONTIGS)
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
                    "library_profile": profile,
                    "layout": layout,
                    "strandedness": "forward",
                    "evidence_source": "auto",
                    "endpoint_model": endpoint_model,
                }
            )
        )
        resolutions.append(str(resolution))
        alignment = write_alignment(tmp_path / f"{sample_id}.bam", CONTIGS, reads(sample_id))
        alignments.append(str(alignment))

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
                    "scan-alignment",
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
                    "--calibration",
                    str(summary),
                    "--output-prefix",
                    str(split / f"{sample_id}.scan"),
                ]
            )
            == 0
        )
        summaries.append(str(summary))
        scanned = json.loads((split / f"{sample_id}.scan.filtering.json").read_text())
        expected_sources = {
            "exact_boundary_SE": ["polyA_junction", "read_3p"],
            "generic_SE": ["polyA_junction", "read_3p", "read_5p"],
            "generic_PE": ["fragment_3p", "polyA_junction", "read_3p", "read_5p"],
        }[case]
        assert sorted(scanned["sources"]) == expected_sources
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
    selected = json.loads((split / "library_resolution.json").read_text())
    assert (selected["evidence_source"], selected["endpoint_model"]) == {
        "exact_boundary_SE": ("polyA_junction", "exact_boundary"),
        "generic_SE": ("read_3p", "exact_boundary"),
        "generic_PE": ("fragment_3p", "exact_boundary"),
    }[case]
