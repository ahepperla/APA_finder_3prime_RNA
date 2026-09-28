import json
from dataclasses import asdict
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
    kernel_diagnostics,
    kernel_modes,
    kernel_overlap,
    minimum_resolvable_separation,
    observation_offsets,
    pooled_kernel,
)
from pacusage.cli import _read_kernel, main
from pacusage.models import EvidenceObservation
from pacusage.tableio import read_tsv


def test_empirical_kernel_is_normalized_and_deterministic() -> None:
    first = empirical_kernel([0, 0, 1, -1], -5, 5)
    second = empirical_kernel([0, 0, 1, -1], -5, 5)
    assert np.allclose(first, second)
    assert np.isclose(first.sum(), 1)
    assert kernel_overlap(first, 0) == 1


def test_exact_and_proximal_classification(resolved_params) -> None:
    values = resolved_params()
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


# The outputs of the serial calibration that the split path replaced,
# recorded while both existed. Per sample: evidence source, calibration genes,
# observations, then median_offset, central_low, central_high,
# adjacent_boundary_fraction, poly_a_clip_fraction, and reproducibility.
TRIANGLE = {-2: 1 / 9, -1: 2 / 9, 0: 3 / 9, 1: 2 / 9, 2: 1 / 9}
CALIBRATION_EXPECTED = {
    "exact_boundary_SE": {
        "resolution": ("exact_boundary", "polyA_junction", "exact_boundary", "SE"),
        "samples": {
            "sample_a": ("polyA_junction", 5, 50, (0.0, 0.0, 0.0, 1.0, 1.0, 1.0)),
            "sample_b": ("polyA_junction", 5, 50, (0.0, 0.0, 0.0, 1.0, 1.0, 1.0)),
        },
        "kernel": TRIANGLE,
        "diagnostics": (1, 2, 0.0, 0.0),
    },
    "generic_SE": {
        "resolution": ("generic_3prime", "read_3p", "exact_boundary", "SE"),
        "samples": {
            "sample_a": ("read_3p", 5, 125, (0.0, -3.0, 0.0, 0.8, 0.2, 0.9463258650628944)),
            "sample_b": ("read_3p", 5, 125, (0.0, 0.0, 3.0, 0.8, 0.2, 0.9463258650628944)),
        },
        "kernel": {
            offset: weight / 90
            for offset, weight in zip(
                range(-5, 6), (1, 2, 3, 10, 17, 24, 17, 10, 3, 2, 1), strict=True
            )
        },
        "diagnostics": (1, 3, 3.0, 1.0),
    },
    "generic_PE": {
        "resolution": ("generic_3prime", "fragment_3p", "exact_boundary", "PE"),
        "samples": {
            "sample_a": ("fragment_3p", 5, 100, (0.0, 0.0, 0.0, 1.0, 0.25, 1.0)),
            "sample_b": ("fragment_3p", 5, 100, (0.0, 0.0, 0.0, 1.0, 0.25, 1.0)),
        },
        "kernel": TRIANGLE,
        "diagnostics": (1, 2, 0.0, 0.0),
    },
}
METRIC_COLUMNS = (
    "median_offset",
    "central_low",
    "central_high",
    "adjacent_boundary_fraction",
    "poly_a_clip_fraction",
    "reproducibility",
)


@pytest.mark.parametrize("case", CALIBRATION_CASES)
def test_split_calibration_outputs(tmp_path: Path, case: str, resolved_params) -> None:
    profile, layout, endpoint_model, reads = CALIBRATION_CASES[case]
    samples = tmp_path / "samples.tsv"
    samples.write_text("sample_id\nsample_a\nsample_b\n")
    annotation = calibration_annotation(tmp_path / "genes.gtf")
    reference = write_reference(tmp_path / "genome.fa", CONTIGS)
    params_path = tmp_path / "resolved_params.yaml"
    resolved_params(params_path, calibration_min_genes=1, pac_min_sample_count=1)

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
                "--kernel-diagnostics",
                str(split / "calibration_kernel_diagnostics.tsv"),
            ]
        )
        == 0
    )

    expected = CALIBRATION_EXPECTED[case]
    selected = json.loads((split / "library_resolution.json").read_text())
    profile_name, source, model, sample_layout = expected["resolution"]
    assert (selected["library_profile"], selected["evidence_source"]) == (profile_name, source)
    assert selected["endpoint_model"] == model
    assert {value["layout"] for value in selected["samples"].values()} == {sample_layout}
    metrics = {row["sample_id"]: row for row in read_tsv(split / "library_calibration.tsv")}
    assert set(metrics) == set(expected["samples"])
    for sample_id, (row_source, genes, observations, values) in expected["samples"].items():
        row = metrics[sample_id]
        assert (row["evidence_source"], row["classification"]) == (row_source, "exact")
        assert (int(row["calibration_genes"]), int(row["observations"])) == (genes, observations)
        assert [float(row[name]) for name in METRIC_COLUMNS] == pytest.approx(values)
    kernel = {
        int(row["offset"]): float(row["weight"])
        for row in read_tsv(split / "calibration_kernel.tsv")
    }
    assert kernel == pytest.approx(expected["kernel"])
    # Every case resolves to exact boundaries, which the kernel warning skips.
    (diagnostics,) = read_tsv(split / "calibration_kernel_diagnostics.tsv")
    assert diagnostics["status"] == "not_applicable"
    assert [
        int(diagnostics["kernel_modes"]),
        int(diagnostics["minimum_resolvable_separation"]),
        float(diagnostics["median_central_interval_width"]),
        float(diagnostics["spread_to_resolution_ratio"]),
    ] == pytest.approx(expected["diagnostics"])


# The pooled kernel from calibrating exact-boundary reads, PACs 50 nt apart,
# under the Plasmidsaurus profile: offsets -2..102, zero outside the spikes.
SPACING_KERNEL = {
    -2: 0.044121495611877164, -1: 0.08824299122375433, 0: 0.13236448683563146,
    1: 0.08824299122375433, 2: 0.044121495611877164, 48: 0.031574023768461475,
    49: 0.06314804753692295, 50: 0.09472207130538442, 51: 0.06314804753692295,
    52: 0.031574023768461475, 98: 0.03541559173077248, 99: 0.07083118346154496,
    100: 0.10624677519231743, 101: 0.07083118346154496, 102: 0.03541559173077248,
}


def dense(weights: dict[int, float]) -> np.ndarray:
    kernel = np.zeros(max(weights) - min(weights) + 1)
    for offset, weight in weights.items():
        kernel[offset - min(weights)] = weight
    return kernel


def gaussian(center: float, sd: float, length: int = 1001) -> np.ndarray:
    positions = np.arange(length)
    return np.exp(-0.5 * ((positions - center) / sd) ** 2)


def test_kernel_of_pac_spacings_warns() -> None:
    row = kernel_diagnostics(dense(SPACING_KERNEL), [100.0] * 4, "proximal_tag", "read_3p", 0.5)
    assert row["status"] == "warning"
    assert (row["kernel_modes"], row["minimum_resolvable_separation"]) == (3, 2)
    assert row["spread_to_resolution_ratio"] == 50.0
    assert row["reason"] == (
        "the pooled kernel has 3 separated modes; the samples' median central interval "
        "(100 nt) is 50.0 times the kernel's minimum resolvable separation (2 nt)"
    )


def test_single_site_kernel_is_ok() -> None:
    kernel = gaussian(300, 60)
    row = kernel_diagnostics(kernel / kernel.sum(), [197.4] * 4, "proximal_tag", "read_3p", 0.5)
    assert (row["status"], row["kernel_modes"], row["reason"]) == ("ok", 1, "")
    assert row["minimum_resolvable_separation"] == 81


@pytest.mark.parametrize(
    ("width", "status", "reason"),
    [
        # Exactly 6 times the 81 nt resolution is still acceptable.
        (486.0, "ok", ""),
        (
            500.0,
            "warning",
            "the samples' median central interval (500 nt) is 6.2 times the kernel's "
            "minimum resolvable separation (81 nt)",
        ),
    ],
)
def test_spread_warns_above_six_times_the_resolution(width, status, reason) -> None:
    kernel = gaussian(300, 60)
    row = kernel_diagnostics(kernel / kernel.sum(), [width] * 4, "proximal_tag", "read_3p", 0.5)
    assert (row["status"], row["kernel_modes"], row["reason"]) == (status, 1, reason)


def test_flat_top_with_ripples_and_tail_bumps_is_one_mode() -> None:
    # Like the Plasmidsaurus-like fixture: a triangle whose flat top has
    # equal-height ripples, with narrow bumps at a fifth of the height at each
    # tail. Raw local maxima would count several modes here.
    offsets = np.arange(0, 301)
    kernel = np.minimum(offsets, 300 - offsets).astype(float)
    ripples = np.where(np.abs(offsets - 150) < 10, 1.5 * np.cos(offsets), 0)
    kernel = np.minimum(kernel, 140.0) + ripples
    for tail in (10, 290):
        kernel[tail - 2 : tail + 3] += 28.0
    kernel /= kernel.sum()
    assert kernel_modes(kernel, 1) > 1
    row = kernel_diagnostics(kernel, [220.0] * 4, "proximal_tag", "read_3p", 0.5)
    assert (row["status"], row["kernel_modes"]) == ("ok", 1)


def test_a_shoulder_above_half_its_height_is_one_mode() -> None:
    # Smoothing at the 89 nt resolution leaves a dip between the peak and its
    # shoulder 140 nt away, but the dip stays well above half the shoulder.
    kernel = gaussian(300, 40) + 0.6 * gaussian(440, 40)
    row = kernel_diagnostics(kernel / kernel.sum(), [240.0] * 4, "proximal_tag", "read_3p", 0.5)
    assert (row["status"], row["kernel_modes"]) == ("ok", 1)
    assert row["minimum_resolvable_separation"] == 89


def test_a_minor_site_is_not_a_mode_but_widens_the_spread() -> None:
    # A site at 15% of the main one's height, 300 nt away, stays below a
    # quarter of the highest mode. It still shrinks the resolution to 55 nt
    # and widens the central interval, so the spread check warns alone.
    kernel = gaussian(300, 40) + 0.15 * gaussian(600, 40)
    row = kernel_diagnostics(kernel / kernel.sum(), [375.0] * 4, "proximal_tag", "read_3p", 0.5)
    assert (row["status"], row["kernel_modes"], row["minimum_resolvable_separation"]) == (
        "warning",
        1,
        55,
    )
    assert row["reason"].startswith("the samples' median central interval (375 nt)")


def test_two_sites_with_a_valley_below_half_the_lower_peak_are_two_modes() -> None:
    # After smoothing at the 70 nt resolution, the valley between sites 165 nt
    # apart falls to 37% of the lower peak.
    kernel = gaussian(300, 40) + gaussian(465, 40)
    row = kernel_diagnostics(kernel / kernel.sum(), [267.0] * 4, "proximal_tag", "read_3p", 0.5)
    assert (row["status"], row["kernel_modes"], row["minimum_resolvable_separation"]) == (
        "warning",
        2,
        70,
    )
    assert row["reason"] == "the pooled kernel has 2 separated modes"


def test_two_separated_modes_warn_without_a_wide_spread() -> None:
    kernel = gaussian(200, 40) + 0.3 * gaussian(450, 40)
    row = kernel_diagnostics(kernel / kernel.sum(), [150.0] * 4, "proximal_tag", "read_3p", 0.5)
    assert (row["status"], row["kernel_modes"]) == ("warning", 2)
    assert row["reason"] == "the pooled kernel has 2 separated modes"


def test_exact_boundary_kernels_are_described_but_not_judged() -> None:
    row = kernel_diagnostics(dense(SPACING_KERNEL), [100.0] * 4, "exact_boundary", "read_3p", 0.5)
    assert (row["status"], row["reason"]) == ("not_applicable", "")
    assert (row["kernel_modes"], row["minimum_resolvable_separation"]) == (3, 2)


def test_diagnostics_describe_the_kernel_discovery_reads_back(
    tmp_path: Path, resolved_params
) -> None:
    # These offsets put the kernel's overlap at 3 nt exactly on the 0.5
    # threshold. Summed over the full pooled array, the tie rounds the other
    # way, giving a resolution that discovery never uses.
    offsets = {191: 3, 193: 3, 196: 1, 197: 1}
    full = pooled_kernel([empirical_kernel(offsets, -1000, 1000)] * 2)
    assert minimum_resolvable_separation(full, 0.5) == 3, "the offsets no longer tie"
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("sample_id\tcondition\nS1\tA\nS2\tA\n")
    params = tmp_path / "resolved_params.yaml"
    resolved_params(params)
    summaries = []
    for sample_id in ("S1", "S2"):
        metrics = calculate_metrics(sample_id, "read_3p", offsets, 100, 0, 0.05, 0.95)
        resolution = {
            "sample_id": sample_id,
            "library_profile": "plasmidsaurus_3prime",
            "layout": "SE",
            "strandedness": "forward",
            "evidence_source": "read_3p",
            "endpoint_model": "proximal_tag",
        }
        summary = tmp_path / f"{sample_id}.calibration.json"
        summary.write_text(
            json.dumps(
                {
                    "sample_id": sample_id,
                    "resolution": resolution,
                    "sources": {
                        "read_3p": {"metrics": asdict(metrics), "offset_counts": offsets}
                    },
                }
            )
        )
        summaries.append(str(summary))
    arguments = ["aggregate-calibration", "--samples", str(sheet), "--calibrations", *summaries]
    arguments += ["--params", str(params), "--output", str(tmp_path / "calibration.tsv")]
    arguments += ["--kernel", str(tmp_path / "kernel.tsv")]
    arguments += ["--resolution", str(tmp_path / "run.json")]
    arguments += ["--kernel-diagnostics", str(tmp_path / "diagnostics.tsv")]
    assert main(arguments) == 0

    kernel, _ = _read_kernel(tmp_path / "kernel.tsv")
    [row] = read_tsv(tmp_path / "diagnostics.tsv")
    assert minimum_resolvable_separation(kernel, 0.5) == 4
    assert row["minimum_resolvable_separation"] == "4"
