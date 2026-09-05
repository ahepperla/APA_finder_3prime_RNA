"""Command-line entry points used by people, tests, and Nextflow processes."""

from __future__ import annotations

import argparse
import json
import platform
import tempfile
from collections.abc import Iterator, Mapping
from dataclasses import asdict
from heapq import merge
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pysam

from . import __version__
from .alignments import inspect_alignment, prepare_alignment, validate_contigs
from .annotation import annotate_candidates, load_known_pacs, load_motif_catalog
from .calibration import (
    CalibrationMetrics,
    calculate_metrics,
    classify_metrics,
    empirical_kernel,
    kernel_correlations,
    observation_offsets,
    pooled_kernel,
)
from .clustering import (
    candidates_as_rows,
    cluster_exact_boundaries,
    discover_proximal_pacs,
)
from .errors import PacusageError
from .evidence import (
    extract_evidence,
    infer_strandedness,
    write_bedgraphs,
    write_evidence,
)
from .models import EvidenceObservation, PacCandidate
from .parameters import load_parameters, write_resolved_parameters
from .profiles import ProtocolProfile, get_profile, resolve_profile_defaults
from .quantification import build_count_outputs, quantify_exact, quantify_proximal
from .reference import (
    annotation_contigs,
    iter_annotation_features,
    parse_annotation,
    prepare_reference,
    transcript_ends,
)
from .report import build_report
from .samples import control_mapping_rows, read_and_validate_samples, write_normalized_samples
from .statistics import add_bh_fdr, cmh_kmer_test, filter_testable_features, motif_usage_scores
from .tableio import iter_tsv, read_tsv, sha256_file, write_json, write_tsv


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        arguments.function(arguments)
    except PacusageError as error:
        parser.exit(2, f"PACusage error: {error}\n")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pacusage",
        description="PAC discovery and differential usage workflow utilities.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    validate = commands.add_parser("validate", help="validate parameters and sample sheet")
    validate.add_argument("--params", required=True)
    validate.add_argument("--samples")
    validate.add_argument("--output-dir", default=".")
    validate.add_argument("--skip-alignment-open", action="store_true")
    validate.set_defaults(function=command_validate)

    reference = commands.add_parser("prepare-reference", help="copy and index a FASTA")
    reference.add_argument("--fasta", required=True)
    reference.add_argument("--fai")
    reference.add_argument("--output", required=True)
    reference.add_argument("--metadata", required=True)
    reference.set_defaults(function=command_prepare_reference)

    alignment = commands.add_parser("prepare-alignment", help="sort and index one BAM/CRAM")
    alignment.add_argument("--sample-id", required=True)
    alignment.add_argument("--alignment", required=True)
    alignment.add_argument("--reference", required=True)
    alignment.add_argument("--output", required=True)
    alignment.add_argument("--metadata", required=True)
    alignment.add_argument("--threads", type=int, default=1)
    alignment.set_defaults(function=command_prepare_alignment)

    strand = commands.add_parser("infer-strandedness", help="infer layout and strandedness")
    strand.add_argument("--sample-id", required=True)
    strand.add_argument("--alignment", required=True)
    strand.add_argument("--reference", required=True)
    strand.add_argument("--annotation", required=True)
    strand.add_argument("--profile", required=True)
    strand.add_argument("--layout", default="auto")
    strand.add_argument("--strandedness", default="auto")
    strand.add_argument("--evidence-source", default="auto")
    strand.add_argument("--endpoint-model", default="auto")
    strand.add_argument("--minimum-informative", type=int, default=10000)
    strand.add_argument("--maximum-sampled", type=int, default=200000)
    strand.add_argument("--decision-fraction", type=float, default=0.8)
    strand.add_argument("--random-seed", type=int, default=1729)
    strand.add_argument("--output", required=True)
    strand.add_argument("--qc", required=True)
    strand.set_defaults(function=command_infer_strandedness)

    calibrate = commands.add_parser("calibrate", help="select one compatible evidence model")
    calibrate.add_argument("--samples", required=True)
    calibrate.add_argument("--alignments", nargs="+", required=True)
    calibrate.add_argument("--resolutions", nargs="+", required=True)
    calibrate.add_argument("--reference", required=True)
    calibrate.add_argument("--annotation", required=True)
    calibrate.add_argument("--params", required=True)
    calibrate.add_argument("--output", required=True)
    calibrate.add_argument("--kernel", required=True)
    calibrate.add_argument("--resolution", required=True)
    calibrate.add_argument("--threads", type=int, default=1)
    calibrate.set_defaults(function=command_calibrate)

    calibration_reference = commands.add_parser(
        "calibration-reference",
        help="build compact annotated transcript ends for calibration",
    )
    calibration_reference.add_argument("--annotation", required=True)
    calibration_reference.add_argument("--output", required=True)
    calibration_reference.set_defaults(function=command_calibration_reference)

    calibrate_sample = commands.add_parser(
        "calibrate-sample",
        help="calculate calibration summaries for one sample",
    )
    calibrate_sample.add_argument("--sample-id", required=True)
    calibrate_sample.add_argument("--alignment", required=True)
    calibrate_sample.add_argument("--resolution", required=True)
    calibrate_sample.add_argument("--reference", required=True)
    calibrate_sample.add_argument("--transcript-ends", required=True)
    calibrate_sample.add_argument("--params", required=True)
    calibrate_sample.add_argument("--output", required=True)
    calibrate_sample.add_argument("--threads", type=int, default=1)
    calibrate_sample.set_defaults(function=command_calibrate_sample)

    aggregate_calibration = commands.add_parser(
        "aggregate-calibration",
        help="combine per-sample calibration summaries",
    )
    aggregate_calibration.add_argument("--samples", required=True)
    aggregate_calibration.add_argument("--calibrations", nargs="+", required=True)
    aggregate_calibration.add_argument("--params", required=True)
    aggregate_calibration.add_argument("--output", required=True)
    aggregate_calibration.add_argument("--kernel", required=True)
    aggregate_calibration.add_argument("--resolution", required=True)
    aggregate_calibration.set_defaults(function=command_aggregate_calibration)

    evidence = commands.add_parser("extract-evidence", help="extract one sample's evidence")
    evidence.add_argument("--sample-id", required=True)
    evidence.add_argument("--alignment", required=True)
    evidence.add_argument("--reference", required=True)
    evidence.add_argument("--resolution", required=True)
    evidence.add_argument("--params", required=True)
    evidence.add_argument("--tsv", required=True)
    evidence.add_argument("--parquet", required=True)
    evidence.add_argument("--plus-track", required=True)
    evidence.add_argument("--minus-track", required=True)
    evidence.add_argument("--qc", required=True)
    evidence.add_argument("--threads", type=int, default=1)
    evidence.set_defaults(function=command_extract_evidence)

    cluster = commands.add_parser("cluster", help="build the condition-blind PAC atlas")
    cluster.add_argument("--evidence", nargs="+", required=True)
    cluster.add_argument("--samples", required=True)
    cluster.add_argument("--resolution", required=True)
    cluster.add_argument("--kernel", required=True)
    cluster.add_argument("--params", required=True)
    cluster.add_argument("--known-pacs")
    cluster.add_argument("--accepted", required=True)
    cluster.add_argument("--rejected", required=True)
    cluster.add_argument("--qc", required=True)
    cluster.set_defaults(function=command_cluster)

    annotate = commands.add_parser("annotate", help="annotate accepted PAC candidates")
    annotate.add_argument("--candidates", required=True)
    annotate.add_argument("--reference", required=True)
    annotate.add_argument("--annotation", required=True)
    annotate.add_argument("--resolution", required=True)
    annotate.add_argument("--params", required=True)
    annotate.add_argument("--known-pacs")
    annotate.add_argument("--metadata", required=True)
    annotate.add_argument("--bed", required=True)
    annotate.add_argument("--motifs")
    annotate.add_argument("--checksum", required=True)
    annotate.set_defaults(function=command_annotate)

    quantify = commands.add_parser("quantify", help="quantify a frozen atlas in one sample")
    quantify.add_argument("--sample-id", required=True)
    quantify.add_argument("--evidence", required=True)
    quantify.add_argument("--atlas", required=True)
    quantify.add_argument("--resolution", required=True)
    quantify.add_argument("--kernel", required=True)
    quantify.add_argument("--params", required=True)
    quantify.add_argument("--counts", required=True)
    quantify.add_argument("--qc", required=True)
    quantify.set_defaults(function=command_quantify)

    merge = commands.add_parser("merge-counts", help="merge sample PAC count tables")
    merge.add_argument("--counts", nargs="+", required=True)
    merge.add_argument("--atlas", required=True)
    merge.add_argument("--wide", required=True)
    merge.add_argument("--long", required=True)
    merge.add_argument("--gene-totals", required=True)
    merge.add_argument("--pau", required=True)
    merge.add_argument("--params", required=True)
    merge.add_argument("--testable", required=True)
    merge.add_argument("--filtering", required=True)
    merge.set_defaults(function=command_merge_counts)

    motifs = commands.add_parser("motif-scores", help="build equal-gene motif usage scores")
    motifs.add_argument("--pau", required=True)
    motifs.add_argument("--atlas", required=True)
    motifs.add_argument("--output", required=True)
    motifs.add_argument("--include-known-rescue", action="store_true")
    motifs.set_defaults(function=command_motif_scores)

    kmers = commands.add_parser(
        "kmer-enrichment", help="test upstream k-mers among gained/increased PACs"
    )
    kmers.add_argument("--statistics", required=True)
    kmers.add_argument("--atlas", required=True)
    kmers.add_argument("--k", type=int, default=6)
    kmers.add_argument("--output-dir", required=True)
    kmers.set_defaults(function=command_kmer_enrichment)

    report = commands.add_parser("report", help="create the portable HTML report")
    report.add_argument("--results", required=True)
    report.add_argument("--output", required=True)
    report.set_defaults(function=command_report)
    return parser


def command_validate(args: argparse.Namespace) -> None:
    params = load_parameters(args.params)
    if args.samples:
        params["input"] = str(Path(args.samples).resolve())
    else:
        params["input"] = str(Path(params["input"]).expanduser().resolve())
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    for key in ("fasta", "gtf"):
        path = Path(params[key]).expanduser().resolve()
        if not path.is_file():
            raise PacusageError(f"Required {key} file does not exist: {path}")
        params[key] = str(path)
    for key in ("known_pacs", "chromosome_aliases"):
        if params.get(key):
            path = Path(params[key]).expanduser().resolve()
            if not path.is_file():
                raise PacusageError(f"Optional {key} file does not exist: {path}")
            params[key] = str(path)
    samples, checks = read_and_validate_samples(
        params["input"], params, check_files=not args.skip_alignment_open
    )
    if not args.skip_alignment_open:
        with tempfile.TemporaryDirectory(prefix="pacusage-validation-") as temporary:
            validation_fasta = Path(temporary) / "genome.fa"
            prepare_reference(params["fasta"], validation_fasta)
            with pysam.FastaFile(str(validation_fasta)) as fasta_handle:
                fasta_contigs = dict(
                    zip(fasta_handle.references, fasta_handle.lengths, strict=True)
                )
            reference_annotation_contigs = annotation_contigs(params["gtf"])
            aliases = _read_aliases(params.get("chromosome_aliases"))
            for sample in samples:
                details = inspect_alignment(sample.alignment, validation_fasta)
                if details["layout"] == "mixed":
                    raise PacusageError(
                        f"Sample {sample.sample_id} contains a meaningful mixture of paired "
                        "and unpaired records. Split it before analysis."
                    )
                checks.append(
                    {
                        "sample_id": sample.sample_id,
                        "check": "alignment_open",
                        "status": "PASS",
                        "detail": f"{details['layout']}; sort={details['sort_order']}",
                        "exploratory_insufficient_replicates": "",
                    }
                )
                validate_contigs(
                    details["contigs"],
                    fasta_contigs,
                    reference_annotation_contigs,
                    aliases,
                )
    write_normalized_samples(samples, output / "normalized_samples.tsv")
    write_resolved_parameters(params, output / "resolved_params.yaml")
    write_tsv(checks, output / "input_validation.tsv")
    write_tsv(control_mapping_rows(samples), output / "control_mapping.tsv")
    checksums = [
        {
            "role": "sample_sheet",
            "path": str(Path(params["input"]).resolve()),
            "sha256": sha256_file(params["input"]),
        },
        {"role": "fasta", "path": params["fasta"], "sha256": sha256_file(params["fasta"])},
        {"role": "annotation", "path": params["gtf"], "sha256": sha256_file(params["gtf"])},
    ]
    checksums.extend(
        {"role": "alignment", "path": sample.alignment, "sha256": sha256_file(sample.alignment)}
        for sample in samples
    )
    write_tsv(checksums, output / "input_checksums.tsv")
    versions = [
        {"software": "pacusage", "version": __version__},
        {"software": "python", "version": platform.python_version()},
        {"software": "pysam", "version": pysam.__version__},
    ]
    write_tsv(versions, output / "software_versions.tsv")
    write_json(
        {
            "pacusage_version": __version__,
            "samples": len(samples),
            "conditions": sorted({sample.condition for sample in samples}),
            "parameters_sha256": sha256_file(output / "resolved_params.yaml"),
        },
        output / "run_manifest.json",
    )


def command_prepare_reference(args: argparse.Namespace) -> None:
    metadata = prepare_reference(args.fasta, args.output, args.fai)
    write_tsv([metadata], args.metadata)


def command_prepare_alignment(args: argparse.Namespace) -> None:
    metadata = prepare_alignment(args.alignment, args.output, args.reference, threads=args.threads)
    metadata["sample_id"] = args.sample_id
    write_json(metadata, Path(args.metadata).with_suffix(".json"))
    write_tsv([metadata], args.metadata)


def command_infer_strandedness(args: argparse.Namespace) -> None:
    inspection = inspect_alignment(args.alignment, args.reference)
    if args.layout == "auto":
        layout = inspection["layout"]
        if layout not in {"SE", "PE"}:
            raise PacusageError(
                f"Sample {args.sample_id}: layout could not be inferred ({layout}). "
                "Set layout explicitly after inspecting the file."
            )
    else:
        layout = args.layout
        if inspection["layout"] in {"SE", "PE"} and inspection["layout"] != layout:
            raise PacusageError(
                f"Sample {args.sample_id}: requested layout {layout} conflicts with "
                f"observed {inspection['layout']} records."
            )
    diagnostics = infer_strandedness(
        args.alignment,
        args.reference,
        iter_annotation_features(args.annotation, {"exon"}),
        args.minimum_informative,
        args.maximum_sampled,
        args.decision_fraction,
        args.random_seed,
    )
    defaults = resolve_profile_defaults(
        args.profile, layout, args.strandedness, args.evidence_source, args.endpoint_model
    )
    requested_strand = defaults["strandedness"]
    inferred = diagnostics["inferred_strandedness"]
    if requested_strand == "auto":
        if inferred == "ambiguous":
            raise PacusageError(
                f"Sample {args.sample_id}: strandedness is ambiguous ({diagnostics['reason']}). "
                "Set strandedness explicitly if external protocol evidence supports it."
            )
        resolved_strand = inferred
    else:
        resolved_strand = requested_strand
        if inferred != "ambiguous" and inferred != resolved_strand:
            raise PacusageError(
                f"Sample {args.sample_id}: profile/requested strandedness {resolved_strand} "
                f"conflicts with BAM evidence ({inferred})."
            )
    resolution = {
        "sample_id": args.sample_id,
        **defaults,
        "layout": layout,
        "strandedness": resolved_strand,
        "alignment_basename": Path(args.alignment).name,
    }
    write_json(resolution, args.output)
    write_tsv([{"sample_id": args.sample_id, **diagnostics, **resolution}], args.qc)


def command_calibrate(args: argparse.Namespace) -> None:
    params = load_parameters(args.params)
    sample_rows = {row["sample_id"]: row for row in read_tsv(args.samples)}
    resolutions = {
        value["sample_id"]: value for value in (_read_json(path) for path in args.resolutions)
    }
    alignments = {_sample_id_from_alignment(path): path for path in args.alignments}
    if set(sample_rows) != set(resolutions) or set(sample_rows) != set(alignments):
        raise PacusageError(
            "Calibration inputs do not contain the same sample IDs in the normalized sheet, "
            "alignment files, and resolution files."
        )
    profiles = {value["library_profile"] for value in resolutions.values()}
    if len(profiles) != 1:
        raise PacusageError(
            "Samples use incompatible library profiles. Analyze each protocol in a separate run."
        )
    profile = get_profile(next(iter(profiles)))
    ends = transcript_ends(
        iter_annotation_features(args.annotation, {"exon", "transcript", "mrna"})
    )
    common_sources = set(profile.candidate_sources(next(iter(resolutions.values()))["layout"]))
    for resolution in resolutions.values():
        common_sources.intersection_update(profile.candidate_sources(resolution["layout"]))
        if resolution["evidence_source"] != "auto":
            common_sources.intersection_update({resolution["evidence_source"]})
    if not common_sources:
        raise PacusageError("No evidence source is compatible with every sample.")

    source_metrics: dict[str, list[CalibrationMetrics]] = {}
    source_offsets: dict[str, list[Mapping[int, int]]] = {}
    for source in sorted(common_sources):
        metrics = []
        offset_lists = []
        for sample_id in sorted(sample_rows):
            resolution = resolutions[sample_id]
            observations, _ = extract_evidence(
                sample_id,
                alignments[sample_id],
                args.reference,
                resolution["layout"],
                resolution["strandedness"],
                source,
                min_mapq=int(params["min_mapq"]),
                require_unique=bool(params["require_unique"]),
                require_proper_pair=bool(params["require_proper_pair"]),
                exclude_duplicates=bool(params["exclude_duplicates"]),
                excluded_contigs=params["excluded_contigs"],
                contig_aliases=_read_aliases(params.get("chromosome_aliases")),
                threads=args.threads,
            )
            offsets, genes, clips = observation_offsets(
                observations,
                ends,
                int(params["calibration_max_distance"]),
                int(params["pac_min_sample_count"]),
            )
            offset_lists.append(offsets)
            metrics.append(
                calculate_metrics(
                    sample_id,
                    source,
                    offsets,
                    genes,
                    clips,
                    float(params["calibration_quantile_low"]),
                    float(params["calibration_quantile_high"]),
                )
            )
        kernels = [
            empirical_kernel(
                offsets,
                -int(params["calibration_max_distance"]),
                int(params["calibration_max_distance"]),
            )
            for offsets in offset_lists
        ]
        correlations = kernel_correlations(kernels)
        source_metrics[source] = [
            type(metric)(**{**metric.__dict__, "reproducibility": correlation})
            for metric, correlation in zip(metrics, correlations, strict=True)
        ]
        source_offsets[source] = offset_lists

    _write_calibration_outputs(
        source_metrics,
        source_offsets,
        resolutions,
        profile,
        params,
        args.output,
        args.kernel,
        args.resolution,
    )


def command_calibration_reference(args: argparse.Namespace) -> None:
    ends = transcript_ends(
        iter_annotation_features(args.annotation, {"exon", "transcript", "mrna"})
    )
    write_tsv(
        (
            {
                "contig": contig,
                "strand": strand,
                "coordinate": coordinate,
                "gene_id": gene_id,
            }
            for (contig, strand), values in sorted(ends.items())
            for coordinate, gene_id in values
        ),
        args.output,
        ["contig", "strand", "coordinate", "gene_id"],
    )


def command_calibrate_sample(args: argparse.Namespace) -> None:
    params = load_parameters(args.params)
    resolution = _read_json(args.resolution)
    if resolution.get("sample_id") != args.sample_id:
        raise PacusageError(
            f"Calibration resolution belongs to {resolution.get('sample_id')!r}, "
            f"not {args.sample_id!r}."
        )
    profile = get_profile(resolution["library_profile"])
    sources = set(profile.candidate_sources(resolution["layout"]))
    if resolution["evidence_source"] != "auto":
        sources.intersection_update({resolution["evidence_source"]})
    if not sources:
        raise PacusageError(f"Sample {args.sample_id} has no compatible calibration source.")

    ends = _read_transcript_ends(args.transcript_ends)
    source_results: dict[str, dict[str, object]] = {}
    for source in sorted(sources):
        observations, _ = extract_evidence(
            args.sample_id,
            args.alignment,
            args.reference,
            resolution["layout"],
            resolution["strandedness"],
            source,
            min_mapq=int(params["min_mapq"]),
            require_unique=bool(params["require_unique"]),
            require_proper_pair=bool(params["require_proper_pair"]),
            exclude_duplicates=bool(params["exclude_duplicates"]),
            excluded_contigs=params["excluded_contigs"],
            contig_aliases=_read_aliases(params.get("chromosome_aliases")),
            threads=args.threads,
        )
        offsets, genes, clips = observation_offsets(
            observations,
            ends,
            int(params["calibration_max_distance"]),
            int(params["pac_min_sample_count"]),
        )
        metrics = calculate_metrics(
            args.sample_id,
            source,
            offsets,
            genes,
            clips,
            float(params["calibration_quantile_low"]),
            float(params["calibration_quantile_high"]),
        )
        source_results[source] = {
            "metrics": asdict(metrics),
            "offset_counts": {str(offset): count for offset, count in sorted(offsets.items())},
        }
    write_json(
        {
            "sample_id": args.sample_id,
            "resolution": resolution,
            "sources": source_results,
        },
        args.output,
    )


def command_aggregate_calibration(args: argparse.Namespace) -> None:
    params = load_parameters(args.params)
    sample_rows = {row["sample_id"]: row for row in read_tsv(args.samples)}
    payloads: dict[str, dict[str, Any]] = {}
    for path in args.calibrations:
        payload = _read_json(path)
        sample_id = str(payload.get("sample_id", ""))
        if not sample_id:
            raise PacusageError(f"Per-sample calibration summary {path} lacks a sample_id.")
        if sample_id in payloads:
            raise PacusageError(f"Duplicate per-sample calibration summary for {sample_id}.")
        payloads[sample_id] = payload
    if set(sample_rows) != set(payloads):
        raise PacusageError(
            "Per-sample calibration summaries do not match the normalized sample sheet."
        )
    resolutions = {
        sample_id: payloads[sample_id]["resolution"] for sample_id in sorted(payloads)
    }
    profiles = {value["library_profile"] for value in resolutions.values()}
    if len(profiles) != 1:
        raise PacusageError(
            "Samples use incompatible library profiles. Analyze each protocol in a separate run."
        )
    profile = get_profile(next(iter(profiles)))
    common_sources = set.intersection(
        *(set(payload["sources"]) for payload in payloads.values())
    )
    if not common_sources:
        raise PacusageError("No evidence source is compatible with every sample.")

    source_metrics: dict[str, list[CalibrationMetrics]] = {}
    source_offsets: dict[str, list[dict[int, int]]] = {}
    for source in sorted(common_sources):
        metrics = []
        offsets = []
        for sample_id in sorted(payloads):
            result = payloads[sample_id]["sources"][source]
            metrics.append(CalibrationMetrics(**result["metrics"]))
            offsets.append(
                {int(offset): int(count) for offset, count in result["offset_counts"].items()}
            )
        kernels = [
            empirical_kernel(
                values,
                -int(params["calibration_max_distance"]),
                int(params["calibration_max_distance"]),
            )
            for values in offsets
        ]
        correlations = kernel_correlations(kernels)
        source_metrics[source] = [
            type(metric)(**{**metric.__dict__, "reproducibility": correlation})
            for metric, correlation in zip(metrics, correlations, strict=True)
        ]
        source_offsets[source] = offsets

    _write_calibration_outputs(
        source_metrics,
        source_offsets,
        resolutions,
        profile,
        params,
        args.output,
        args.kernel,
        args.resolution,
    )


def _write_calibration_outputs(
    source_metrics: dict[str, list[CalibrationMetrics]],
    source_offsets: dict[str, list[Mapping[int, int]]],
    resolutions: dict[str, dict[str, Any]],
    profile: ProtocolProfile,
    params: dict[str, Any],
    output_path: str | Path,
    kernel_path: str | Path,
    resolution_path: str | Path,
) -> None:
    selected_source = _select_source(source_metrics, resolutions, profile, params)
    selected_metrics = source_metrics[selected_source]
    requested_models = {value["endpoint_model"] for value in resolutions.values()}
    if len(requested_models) != 1:
        raise PacusageError("All samples must request the same endpoint_model.")
    requested_model = next(iter(requested_models))
    classified = [classify_metrics(metric, params) for metric in selected_metrics]
    if requested_model in {"exact_boundary", "proximal_tag"}:
        endpoint_model = requested_model
        if any(
            metric.calibration_genes < int(params["calibration_min_genes"])
            or metric.reproducibility < float(params["calibration_min_kernel_correlation"])
            for metric in classified
        ):
            raise PacusageError(
                f"Explicit endpoint_model={requested_model} failed minimum calibration "
                "gene or reproducibility requirements."
            )
    else:
        classes = {metric.classification for metric in classified}
        if classes == {"exact"}:
            endpoint_model = "exact_boundary"
        elif classes == {"proximal"}:
            endpoint_model = "proximal_tag"
        else:
            detail = ", ".join(
                f"{metric.sample_id}={metric.classification}" for metric in classified
            )
            raise PacusageError(
                f"Samples do not resolve to one compatible endpoint model: {detail}."
            )
    kernels = [
        empirical_kernel(
            offsets,
            -int(params["calibration_max_distance"]),
            int(params["calibration_max_distance"]),
        )
        for offsets in source_offsets[selected_source]
    ]
    pooled = pooled_kernel(kernels)
    kernel_minimum = -int(params["calibration_max_distance"])
    write_tsv([asdict(metric) for metric in classified], output_path)
    write_tsv(
        (
            {"offset": index + kernel_minimum, "weight": value}
            for index, value in enumerate(pooled)
            if value > 0
        ),
        kernel_path,
        ["offset", "weight"],
    )
    write_json(
        {
            "library_profile": profile.name,
            "evidence_source": selected_source,
            "endpoint_model": endpoint_model,
            "kernel_minimum_offset": kernel_minimum,
            "kernel_maximum_offset": int(params["calibration_max_distance"]),
            "calibration_version": "PACusage-0.1",
            "samples": {
                sample_id: {
                    "layout": value["layout"],
                    "strandedness": value["strandedness"],
                }
                for sample_id, value in resolutions.items()
            },
        },
        resolution_path,
    )


def _select_source(
    source_metrics: dict[str, list[CalibrationMetrics]],
    resolutions: dict[str, dict[str, Any]],
    profile: ProtocolProfile,
    params: dict[str, Any],
) -> str:
    explicit = {value["evidence_source"] for value in resolutions.values()}
    explicit.discard("auto")
    if explicit:
        if len(explicit) != 1:
            raise PacusageError("Samples request conflicting evidence_source values.")
        return next(iter(explicit))
    if profile.name == "plasmidsaurus_3prime":
        return "read_3p"
    if profile.name == "exact_boundary":
        poly_a = source_metrics.get("polyA_junction", [])
        if poly_a and all(
            metric.calibration_genes >= int(params["calibration_min_genes"]) for metric in poly_a
        ):
            return "polyA_junction"
        return "fragment_3p" if "fragment_3p" in source_metrics else "read_3p"
    ranked = sorted(
        (
            (float(np.median([metric.reproducibility for metric in metrics])), source)
            for source, metrics in source_metrics.items()
        ),
        reverse=True,
    )
    threshold = float(params["calibration_min_kernel_correlation"])
    margin = float(params["calibration_min_model_margin"])
    if not ranked or ranked[0][0] < threshold:
        raise PacusageError("No evidence source passed the calibration reproducibility threshold.")
    if len(ranked) > 1 and ranked[0][0] - ranked[1][0] < margin:
        raise PacusageError(
            "Evidence-source calibration is ambiguous: the best source does not exceed "
            f"the runner-up by calibration_min_model_margin={margin}."
        )
    return ranked[0][1]


def command_extract_evidence(args: argparse.Namespace) -> None:
    params = load_parameters(args.params)
    resolution = _read_json(args.resolution)
    sample_resolution = resolution.get("samples", {}).get(args.sample_id, resolution)
    if "layout" not in sample_resolution or "strandedness" not in sample_resolution:
        raise PacusageError(
            f"Run resolution does not contain layout/strandedness for sample {args.sample_id}."
        )
    observations, qc = extract_evidence(
        args.sample_id,
        args.alignment,
        args.reference,
        sample_resolution["layout"],
        sample_resolution["strandedness"],
        resolution["evidence_source"],
        min_mapq=int(params["min_mapq"]),
        require_unique=bool(params["require_unique"]),
        require_proper_pair=bool(params["require_proper_pair"]),
        exclude_duplicates=bool(params["exclude_duplicates"]),
        excluded_contigs=params["excluded_contigs"],
        contig_aliases=_read_aliases(params.get("chromosome_aliases")),
        threads=args.threads,
    )
    write_evidence(observations, args.tsv, args.parquet)
    write_bedgraphs(observations, args.plus_track, args.minus_track)
    write_tsv([qc], args.qc)


def command_cluster(args: argparse.Namespace) -> None:
    params = load_parameters(args.params)
    resolution = _read_json(args.resolution)
    sample_conditions = {
        row["sample_id"]: row["condition"] for row in iter_tsv(args.samples)
    }
    observations = _iter_observations(args.evidence, merge_sorted=True)
    known = load_known_pacs(args.known_pacs)
    if resolution["endpoint_model"] == "exact_boundary":
        accepted, rejected = cluster_exact_boundaries(
            observations,
            int(params["pac_seed_radius"]),
            int(params["pac_cluster_radius"]),
            int(params["pac_min_total_count"]),
            int(params["pac_min_sample_count"]),
            float(params["pac_min_supporting_samples"]),
            known,
            int(params["known_pac_rescue_total"]),
            int(params["known_pac_match_radius"]),
            observations_sorted=True,
            sample_conditions=sample_conditions,
        )
        minimum_resolution = int(params["pac_cluster_radius"])
    else:
        kernel, minimum = _read_kernel(args.kernel)
        accepted, rejected, minimum_resolution = discover_proximal_pacs(
            observations,
            kernel,
            minimum,
            float(params["proximal_kernel_overlap_threshold"]),
            int(params["pac_min_total_count"]),
            int(params["pac_min_sample_count"]),
            float(params["pac_min_supporting_samples"]),
            int(params["proximal_bin_size"]),
            float(params["proximal_assignment_likelihood_ratio"]),
            observations_sorted=True,
            sample_conditions=sample_conditions,
        )
    write_tsv(candidates_as_rows(accepted), args.accepted, compresslevel=1)
    write_tsv(candidates_as_rows(rejected), args.rejected, compresslevel=1)
    write_tsv(
        [
            {
                "endpoint_model": resolution["endpoint_model"],
                "accepted_pacs": len(accepted),
                "rejected_candidates": len(rejected),
                "minimum_resolvable_separation": minimum_resolution,
                "support_requirement": _support_requirement_description(
                    float(params["pac_min_supporting_samples"])
                ),
            }
        ],
        args.qc,
    )


def _support_requirement_description(threshold: float) -> str:
    if threshold < 1:
        return (
            f"{threshold:g} of samples within one condition, rounded up"
        )
    return f"{int(threshold)} samples within one condition"


def command_annotate(args: argparse.Namespace) -> None:
    params = load_parameters(args.params)
    resolution = _read_json(args.resolution)
    candidates = [
        PacCandidate(
            contig=row["contig"],
            strand=row["strand"],
            coordinate=int(row["coordinate"]),
            total_count=int(row["total_count"]),
            supporting_samples=int(row["supporting_samples"]),
            capped_support=int(row["capped_support"]),
            member_coordinates=tuple(
                int(value) for value in row["member_coordinates"].split(",") if value
            ),
            status=row["status"],
            rejection_reason=row.get("rejection_reason", ""),
            fraction_within_2nt=float(row.get("fraction_within_2nt", 0)),
            width_90=int(row.get("width_90", 0)),
            poly_a_clip_fraction=float(row.get("poly_a_clip_fraction", 0)),
            local_enrichment=float(row.get("local_enrichment", 0)),
            coordinate_interval_low=_optional_int(row.get("coordinate_interval_low")),
            coordinate_interval_high=_optional_int(row.get("coordinate_interval_high")),
            coordinate_bootstrap_successes=int(row.get("coordinate_bootstrap_successes", 0)),
            region_start=_optional_int(row.get("region_start")),
            region_end=_optional_int(row.get("region_end")),
            resolution_nt=int(row.get("resolution_nt", 0)),
            total_supporting_samples=int(row.get("total_supporting_samples", 0)),
            supporting_condition=row.get("supporting_condition", ""),
        )
        for row in read_tsv(args.candidates)
    ]
    known = load_known_pacs(args.known_pacs)
    rows = annotate_candidates(
        candidates,
        parse_annotation(args.annotation, {"gene", "exon"}),
        args.reference,
        str(params["assembly"]),
        resolution["endpoint_model"],
        resolution["evidence_source"],
        int(params["max_downstream_distance"]),
        int(params["pas_scan_upstream_far"]),
        int(params["pas_scan_upstream_near"]),
        int(params["pas_core_upstream_far"]),
        int(params["pas_core_upstream_near"]),
        int(params["internal_priming_window"]),
        int(params["internal_priming_max_a_run"]),
        float(params["internal_priming_max_a_fraction"]),
        known,
        int(params["known_pac_match_radius"]),
        load_motif_catalog(params.get("pas_motif_catalog")),
    )
    write_tsv(rows, args.metadata)
    _write_bed(rows, args.bed)
    if args.motifs:
        motif_columns = [
            "pac_id",
            "gene_id",
            "strand",
            "upstream_sequence",
            "all_pas_motifs",
            "primary_pas_motif",
            "primary_pas_motif_rna",
            "primary_motif_class",
            "primary_motif_position",
            "primary_motif_in_core",
            "known_rescue_only",
        ]
        write_tsv(rows, args.motifs, motif_columns)
    checksum = sha256_file(args.metadata)
    Path(args.checksum).write_text(f"{checksum}  {Path(args.metadata).name}\n")


def command_quantify(args: argparse.Namespace) -> None:
    params = load_parameters(args.params)
    resolution = _read_json(args.resolution)
    observations = _iter_observations([args.evidence])
    atlas = read_tsv(args.atlas)
    if resolution["endpoint_model"] == "exact_boundary":
        rows, qc = quantify_exact(observations, atlas, int(params["pac_cluster_radius"]))
    else:
        kernel, minimum = _read_kernel(args.kernel)
        rows, qc = quantify_proximal(
            observations,
            atlas,
            kernel,
            minimum,
            float(params["proximal_assignment_likelihood_ratio"]),
        )
    for row in rows:
        row["sample_id"] = args.sample_id
    qc["sample_id"] = args.sample_id
    write_tsv(rows, args.counts, ["sample_id", "pac_id", "count"])
    write_tsv([qc], args.qc)


def command_merge_counts(args: argparse.Namespace) -> None:
    params = load_parameters(args.params)
    atlas = read_tsv(args.atlas)
    per_sample: dict[str, pd.Series] = {}
    for path in args.counts:
        frame = pd.read_csv(path, sep="\t", usecols=["sample_id", "pac_id", "count"])
        if frame.empty:
            continue
        sample_ids = frame["sample_id"].astype(str).unique()
        if len(sample_ids) != 1:
            raise PacusageError(f"Count table {path} contains multiple sample IDs.")
        sample_id = sample_ids[0]
        per_sample[sample_id] = frame.set_index(frame["pac_id"].astype(str))["count"]
    wide, long, totals, pau = build_count_outputs(per_sample, atlas)
    wide.to_csv(args.wide, sep="\t", index=False, compression="infer")
    long.to_parquet(args.long, index=False)
    totals.to_csv(args.gene_totals, sep="\t", index=False, compression="infer")
    pau.to_csv(args.pau, sep="\t", index=False, compression="infer")
    _, filtering = filter_testable_features(
        wide,
        int(params["min_gene_total"]),
        int(params["min_site_count"]),
        float(params["min_site_usage"]),
        int(params["min_test_supporting_samples"]),
    )
    model_input = wide[
        (wide["gene_id"].astype(str) != "")
        & ~wide["gene_id"].astype(str).str.contains(",", regex=False)
    ]
    model_input.to_csv(args.testable, sep="\t", index=False, compression="infer")
    filtering.to_csv(args.filtering, sep="\t", index=False)


def command_motif_scores(args: argparse.Namespace) -> None:
    pau = pd.read_csv(args.pau, sep="\t")
    atlas = pd.read_csv(args.atlas, sep="\t")
    scores = motif_usage_scores(pau, atlas, excluded_rescue_sites=not args.include_known_rescue)
    scores.to_csv(args.output, sep="\t", index=False)


def command_kmer_enrichment(args: argparse.Namespace) -> None:
    statistics_directory = Path(args.statistics)
    output_directory = Path(args.output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)
    atlas = pd.read_csv(args.atlas, sep="\t").fillna("")
    atlas = atlas[
        ~atlas["known_rescue_only"].map(_as_bool)
        & (atlas["candidate_status"].astype(str) != "motif_assisted_rescue")
        & ~atlas["ambiguous_gene_assignment"].map(_as_bool)
        & (atlas["gene_id"].astype(str) != "")
    ].copy()
    completed = []
    for event_path in sorted(statistics_directory.glob("*.events.tsv.gz")):
        events = pd.read_csv(event_path, sep="\t").fillna("")
        event_id_column = "pac_id" if "pac_id" in events else "feature_id"
        selected = set(
            events.loc[
                events["event_type"].isin(["gained", "increased_usage"]),
                event_id_column,
            ].astype(str)
        )
        rows = _kmer_rows(atlas, selected, args.k)
        if rows:
            adjusted = add_bh_fdr(row["p_value"] for row in rows)
            for row, fdr in zip(rows, adjusted, strict=True):
                row["fdr"] = fdr
        comparison = event_path.name.removesuffix(".events.tsv.gz")
        destination = output_directory / f"{comparison}.kmer_enrichment.tsv.gz"
        write_tsv(
            rows,
            destination,
            [
                "kmer",
                "common_odds_ratio",
                "ci_low",
                "ci_high",
                "p_value",
                "fdr",
                "event_pacs",
                "background_pacs",
                "informative_genes",
                "analysis_label",
            ],
        )
        completed.append({"comparison": comparison, "tested_kmers": len(rows)})
    write_tsv(completed, output_directory / "kmer_enrichment_status.tsv")


def command_report(args: argparse.Namespace) -> None:
    build_report(args.results, args.output)


def _read_json(path: str | Path) -> dict[str, Any]:
    with Path(path).open() as handle:
        return json.load(handle)


def _read_aliases(path: str | Path | None) -> dict[str, str]:
    if not path:
        return {}
    aliases: dict[str, str] = {}
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            values = line.rstrip("\n").split("\t")
            if len(values) < 2:
                raise PacusageError(
                    f"Chromosome alias file {path}, line {line_number}, "
                    "must contain at least two tab-separated columns."
                )
            if line_number == 1 and values[0].lower() in {"alias", "source", "from"}:
                continue
            aliases[values[0]] = values[1]
    return aliases


def _optional_int(value: object) -> int | None:
    return None if value in (None, "", "None", "nan") else int(value)


def _as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def _kmer_rows(
    atlas: pd.DataFrame, selected_pac_ids: set[str], kmer_length: int
) -> list[dict[str, object]]:
    event_presence: dict[str, dict[str, list[bool]]] = {}
    background_presence: dict[str, dict[str, list[bool]]] = {}
    all_kmers: set[str] = set()
    for gene_id, gene_rows in atlas.groupby("gene_id"):
        event_rows = gene_rows[gene_rows["pac_id"].astype(str).isin(selected_pac_ids)]
        background_rows = gene_rows[~gene_rows["pac_id"].astype(str).isin(selected_pac_ids)]
        if event_rows.empty or background_rows.empty:
            continue
        sequences = gene_rows["upstream_sequence"].astype(str).str.upper()
        gene_kmers = {
            sequence[index : index + kmer_length]
            for sequence in sequences
            for index in range(max(0, len(sequence) - kmer_length + 1))
            if set(sequence[index : index + kmer_length]) <= set("ACGT")
        }
        for kmer in gene_kmers:
            event_presence.setdefault(kmer, {})[str(gene_id)] = [
                kmer in sequence
                for sequence in event_rows["upstream_sequence"].astype(str).str.upper()
            ]
            background_presence.setdefault(kmer, {})[str(gene_id)] = [
                kmer in sequence
                for sequence in background_rows["upstream_sequence"].astype(str).str.upper()
            ]
        all_kmers.update(gene_kmers)
    rows = []
    for kmer in sorted(all_kmers):
        result = cmh_kmer_test(event_presence[kmer], background_presence[kmer])
        if int(result["informative_genes"]) == 0:
            continue
        rows.append(
            {
                "kmer": kmer,
                **result,
                "event_pacs": sum(len(values) for values in event_presence[kmer].values()),
                "background_pacs": sum(
                    len(values) for values in background_presence[kmer].values()
                ),
                "analysis_label": "exploratory",
            }
        )
    return rows


def _sample_id_from_alignment(path: str) -> str:
    name = Path(path).name
    if name.endswith(".bam") or name.endswith(".cram"):
        return name.rsplit(".", 1)[0]
    return name


def _iter_observations(
    paths: list[str],
    merge_sorted: bool = False,
) -> Iterator[EvidenceObservation]:
    streams = [_iter_observation_file(path, require_sorted=merge_sorted) for path in paths]
    if merge_sorted:
        yield from merge(*streams, key=_observation_sort_key)
        return
    for stream in streams:
        yield from stream


def _iter_observation_file(
    path: str | Path,
    require_sorted: bool,
) -> Iterator[EvidenceObservation]:
    path = Path(path)
    rows: Iterator[Mapping[str, object]]
    if path.suffix.lower() == ".parquet":
        rows = _iter_parquet_rows(path)
    else:
        rows = iter_tsv(path)
    previous: tuple[str, str, int, str] | None = None
    for row in rows:
        observation = EvidenceObservation(
            sample_id=str(row["sample_id"]),
            contig=str(row["contig"]),
            strand=str(row["strand"]),
            coordinate=int(row["coordinate"]),
            count=int(row["count"]),
            poly_a_clip_count=int(row.get("poly_a_clip_count", 0)),
            evidence_source=str(row["evidence_source"]),
        )
        key = _observation_sort_key(observation)
        if require_sorted and previous is not None and key < previous:
            raise PacusageError(f"Evidence table is not coordinate-sorted: {path}")
        previous = key
        yield observation


def _iter_parquet_rows(path: Path) -> Iterator[dict[str, object]]:
    columns = [
        "sample_id",
        "contig",
        "strand",
        "coordinate",
        "count",
        "poly_a_clip_count",
        "evidence_source",
    ]
    parquet = pq.ParquetFile(path)
    for batch in parquet.iter_batches(batch_size=100_000, columns=columns):
        values = batch.to_pydict()
        for row in zip(*(values[column] for column in columns), strict=True):
            yield dict(zip(columns, row, strict=True))


def _observation_sort_key(
    observation: EvidenceObservation,
) -> tuple[str, str, int, str]:
    return (
        observation.contig,
        observation.strand,
        observation.coordinate,
        observation.sample_id,
    )


def _read_kernel(path: str | Path) -> tuple[np.ndarray, int]:
    rows = read_tsv(path)
    if not rows:
        return np.asarray([1.0]), 0
    minimum = min(int(row["offset"]) for row in rows)
    maximum = max(int(row["offset"]) for row in rows)
    kernel = np.zeros(maximum - minimum + 1)
    for row in rows:
        kernel[int(row["offset"]) - minimum] = float(row["weight"])
    return kernel, minimum


def _read_transcript_ends(
    path: str | Path,
) -> dict[tuple[str, str], list[tuple[int, str]]]:
    ends: dict[tuple[str, str], list[tuple[int, str]]] = {}
    for row in iter_tsv(path):
        ends.setdefault((row["contig"], row["strand"]), []).append(
            (int(row["coordinate"]), row["gene_id"])
        )
    for key in ends:
        ends[key].sort()
    return ends


def _write_bed(rows: list[dict[str, object]], path: str | Path) -> None:
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as handle:
        for row in rows:
            coordinate = int(row["coordinate"])
            if row["endpoint_model"] == "proximal_tag":
                start = int(row["region_start"])
                end = int(row["region_end"])
            else:
                start = coordinate
                end = coordinate + 1
            handle.write(
                "\t".join(
                    [
                        str(row["contig"]),
                        str(start),
                        str(end),
                        str(row["pac_id"]),
                        str(row["total_count"]),
                        str(row["strand"]),
                    ]
                )
                + "\n"
            )


if __name__ == "__main__":
    raise SystemExit(main())
