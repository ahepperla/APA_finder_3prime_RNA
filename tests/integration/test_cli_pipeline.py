import json
from pathlib import Path

import pysam
import pytest

from pacusage.annotation import annotate_candidates
from pacusage.cli import main
from pacusage.clustering import candidates_as_rows, cluster_exact_boundaries
from pacusage.evidence import write_evidence, write_splice_continuations
from pacusage.models import EvidenceObservation, PacCandidate, SpliceContinuation
from pacusage.quantification import build_count_outputs, quantify_exact
from pacusage.reference import parse_annotation, prepare_reference
from pacusage.tableio import read_tsv, write_tsv


def test_small_exact_boundary_pipeline(tmp_path: Path) -> None:
    fasta = tmp_path / "source.fa"
    fasta.write_text(">chr1\n" + "C" * 500 + "\n")
    prepared = tmp_path / "prepared.fa"
    prepare_reference(fasta, prepared)
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\ttest\tgene\t101\t400\t.\t+\t.\tgene_id "g1"; gene_name "G1";\n'
        'chr1\ttest\texon\t101\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    observations = [
        EvidenceObservation("c1", "chr1", "+", 300, 4),
        EvidenceObservation("c2", "chr1", "+", 300, 5),
        EvidenceObservation("t1", "chr1", "+", 300, 2),
        EvidenceObservation("t1", "chr1", "+", 350, 6),
        EvidenceObservation("t2", "chr1", "+", 300, 2),
        EvidenceObservation("t2", "chr1", "+", 350, 7),
    ]
    conditions = {"c1": "control", "c2": "control", "t1": "treatment", "t2": "treatment"}
    candidates, _ = cluster_exact_boundaries(
        observations, 2, 12, 5, 2, 2, sample_conditions=conditions
    )
    atlas = annotate_candidates(
        candidates,
        parse_annotation(gtf),
        prepared,
        "test",
        "exact_boundary",
        "read_3p",
        5000,
        50,
        5,
        35,
        10,
        20,
        6,
        0.6,
    )
    per_sample = {}
    for sample_id in ["c1", "c2", "t1", "t2"]:
        rows, qc = quantify_exact(
            [row for row in observations if row.sample_id == sample_id], atlas, 12
        )
        assert qc["accepted_fragments"] == qc["assigned_fragments"] + qc["unassigned_fragments"]
        per_sample[sample_id] = rows
    wide, _, totals, pau = build_count_outputs(per_sample, atlas)
    assert len(wide) == 2
    assert totals["gene_total"].min() > 0
    assert all(
        abs(value - 1) < 1e-9 for value in pau.groupby(["gene_id", "sample_id"])["pau"].sum()
    )


def test_proximal_cluster_cli_streams_parquet_evidence(tmp_path: Path, resolved_params) -> None:
    evidence_paths = []
    splice_continuation_paths = []
    for sample_id in ("a", "b"):
        observations = [
            EvidenceObservation(sample_id, "chr1", "+", 90, 3),
            EvidenceObservation(sample_id, "chr1", "-", 210, 4),
        ]
        tsv = tmp_path / f"{sample_id}.tsv.gz"
        parquet = tmp_path / f"{sample_id}.parquet"
        write_evidence(observations, tsv, parquet)
        evidence_paths.append(str(parquet))
        splice_continuations = tmp_path / f"{sample_id}.splice_continuations.tsv.gz"
        write_splice_continuations([], splice_continuations)
        splice_continuation_paths.append(str(splice_continuations))

    params = tmp_path / "resolved_params.yaml"
    resolved_params(
        params,
        pac_min_total_count=1,
        pac_min_sample_count=1,
        pac_min_supporting_samples=0.75,
        proximal_bin_size=1,
    )
    resolution = tmp_path / "resolution.json"
    resolution.write_text('{"endpoint_model":"proximal_tag","evidence_source":"read_3p"}\n')
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_id\tcondition\n"
        "a\ttreatment\n"
        "b\ttreatment\n"
    )
    kernel = tmp_path / "kernel.tsv"
    kernel.write_text("offset\tweight\n10\t1.0\n")
    annotation = single_exon_annotation(tmp_path)
    accepted = tmp_path / "accepted.tsv"
    rejected = tmp_path / "rejected.tsv"
    qc = tmp_path / "qc.tsv"

    assert (
        main(
            [
                "cluster",
                "--evidence",
                *evidence_paths,
                "--splice-continuations",
                *splice_continuation_paths,
                "--annotation",
                str(annotation),
                "--samples",
                str(samples),
                "--resolution",
                str(resolution),
                "--kernel",
                str(kernel),
                "--params",
                str(params),
                "--accepted",
                str(accepted),
                "--rejected",
                str(rejected),
                "--qc",
                str(qc),
            ]
        )
        == 0
    )
    rows = read_tsv(accepted)
    assert [(row["strand"], int(row["coordinate"])) for row in rows] == [
        ("+", 100),
        ("-", 200),
    ]
    assert all(int(row["resolution_nt"]) == 1 for row in rows)
    qc_rows = read_tsv(qc)
    assert qc_rows[0]["support_requirement"] == (
        "0.75 of samples within one condition, rounded up"
    )
    assert qc_rows[0]["internal_exon_end_filter_enabled"] == "True"
    assert qc_rows[0]["internal_exon_end_rejected"] == "0"
    assert [row["supporting_conditions"] for row in rows] == ["treatment", "treatment"]


def test_proximal_cluster_cli_rejects_constitutive_readthrough(
    tmp_path: Path, resolved_params
) -> None:
    evidence_paths = []
    continuation_paths = []
    sample_ids = ("control_1", "control_2", "treatment_1", "treatment_2")
    for sample_id in sample_ids:
        observations = [EvidenceObservation(sample_id, "chr1", "+", 90, 3)]
        tsv = tmp_path / f"{sample_id}.tsv.gz"
        parquet = tmp_path / f"{sample_id}.parquet"
        write_evidence(observations, tsv, parquet)
        evidence_paths.append(str(parquet))
        continuation_path = tmp_path / f"{sample_id}.splice_continuations.tsv.gz"
        write_splice_continuations(
            [
                SpliceContinuation(
                    sample_id,
                    "chr1",
                    "+",
                    50,
                    150,
                    200,
                    300,
                    2,
                )
            ],
            continuation_path,
        )
        continuation_paths.append(str(continuation_path))

    params = tmp_path / "resolved_params.yaml"
    resolved_params(
        params,
        pac_min_total_count=1,
        pac_min_sample_count=1,
        pac_min_supporting_samples=1,
        proximal_bin_size=1,
        constitutive_readthrough_min_junction_count=2,
    )
    resolution = tmp_path / "resolution.json"
    resolution.write_text('{"endpoint_model":"proximal_tag","evidence_source":"read_3p"}\n')
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_id\tcondition\n"
        "control_1\tcontrol\n"
        "control_2\tcontrol\n"
        "treatment_1\ttreatment\n"
        "treatment_2\ttreatment\n"
    )
    kernel = tmp_path / "kernel.tsv"
    kernel.write_text("offset\tweight\n10\t1.0\n")
    annotation = single_exon_annotation(tmp_path)
    accepted = tmp_path / "accepted.tsv"
    rejected = tmp_path / "rejected.tsv"
    qc = tmp_path / "qc.tsv"

    assert (
        main(
            [
                "cluster",
                "--evidence",
                *evidence_paths,
                "--splice-continuations",
                *continuation_paths,
                "--annotation",
                str(annotation),
                "--samples",
                str(samples),
                "--resolution",
                str(resolution),
                "--kernel",
                str(kernel),
                "--params",
                str(params),
                "--accepted",
                str(accepted),
                "--rejected",
                str(rejected),
                "--qc",
                str(qc),
            ]
        )
        == 0
    )
    assert read_tsv(accepted) == []
    assert [row["rejection_reason"] for row in read_tsv(rejected)] == [
        "constitutive_readthrough"
    ]
    assert read_tsv(rejected)[0]["supporting_conditions"] == "control;treatment"
    assert read_tsv(qc)[0]["constitutive_readthrough_rejected"] == "1"


def single_exon_annotation(directory: Path) -> Path:
    """An annotation whose only gene has no internal exon."""
    path = directory / "single_exon.gtf"
    path.write_text('chr1\ttest\texon\t1\t1000\t.\t+\t.\tgene_id "g"; transcript_id "t1";\n')
    return path


TWO_EXON_GENE = (
    'chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id "g"; transcript_id "t1";\n'
    'chr1\ttest\texon\t201\t300\t.\t+\t.\tgene_id "g"; transcript_id "t1";\n'
)


def exon_end_cluster(
    tmp_path: Path, resolved_params, annotation: str | None, **params: object
) -> tuple[list[str], dict[str, str], list[dict[str, str]], list[dict[str, str]]]:
    """Cluster reads piled at exon 1's donor (interbase 100) and near the end
    of exon 2 of TWO_EXON_GENE. The kernel peaks 10 nt downstream, so the
    donor's pile makes a candidate at 110. Returns the arguments, the QC row,
    and the accepted and rejected rows."""
    evidence_paths = []
    continuation_paths = []
    for sample_id in ("a", "b"):
        observations = [
            EvidenceObservation(sample_id, "chr1", "+", 100, 3),
            EvidenceObservation(sample_id, "chr1", "+", 290, 3),
        ]
        parquet = tmp_path / f"{sample_id}.parquet"
        write_evidence(observations, tmp_path / f"{sample_id}.tsv.gz", parquet)
        evidence_paths.append(str(parquet))
        continuation_path = tmp_path / f"{sample_id}.splice_continuations.tsv.gz"
        write_splice_continuations([], continuation_path)
        continuation_paths.append(str(continuation_path))
    parameter_file = tmp_path / "resolved_params.yaml"
    resolved_params(
        parameter_file,
        pac_min_total_count=1,
        pac_min_sample_count=1,
        pac_min_supporting_samples=1,
        proximal_bin_size=1,
        **params,
    )
    resolution = tmp_path / "resolution.json"
    resolution.write_text('{"endpoint_model":"proximal_tag","evidence_source":"read_3p"}\n')
    samples = tmp_path / "samples.tsv"
    samples.write_text("sample_id\tcondition\na\ttreatment\nb\ttreatment\n")
    kernel = tmp_path / "kernel.tsv"
    kernel.write_text("offset\tweight\n10\t1.0\n")
    accepted = tmp_path / "accepted.tsv"
    rejected = tmp_path / "rejected.tsv"
    qc = tmp_path / "qc.tsv"
    arguments = ["cluster", "--evidence", *evidence_paths]
    arguments += ["--splice-continuations", *continuation_paths]
    if annotation is not None:
        path = tmp_path / "genes.gtf"
        path.write_text(annotation)
        arguments += ["--annotation", str(path)]
    arguments += ["--samples", str(samples), "--resolution", str(resolution)]
    arguments += ["--kernel", str(kernel), "--params", str(parameter_file)]
    arguments += ["--accepted", str(accepted), "--rejected", str(rejected), "--qc", str(qc)]
    assert main(arguments) == 0
    return arguments, read_tsv(qc)[0], read_tsv(accepted), read_tsv(rejected)


def test_proximal_cluster_cli_rejects_internal_exon_ends(
    tmp_path: Path, resolved_params, capsys: pytest.CaptureFixture[str]
) -> None:
    arguments, summary, accepted, rejected = exon_end_cluster(
        tmp_path, resolved_params, TWO_EXON_GENE
    )
    assert [int(row["coordinate"]) for row in accepted] == [300]
    assert [(int(row["coordinate"]), row["rejection_reason"]) for row in rejected] == [
        (110, "internal_exon_end")
    ]
    assert summary["internal_exon_end_filter_enabled"] == "True"
    assert summary["accepted_before_internal_exon_end_filter"] == "2"
    assert summary["internal_exon_end_rejected"] == "1"
    assert summary["accepted_before_constitutive_readthrough_filter"] == "1"
    # Without the annotation the filter cannot run.
    position = arguments.index("--annotation")
    del arguments[position : position + 2]
    with pytest.raises(SystemExit):
        main(arguments)
    assert "--annotation is required" in capsys.readouterr().err


def test_proximal_cluster_cli_keeps_a_donor_two_bins_from_a_transcript_end(
    tmp_path: Path, resolved_params
) -> None:
    # Another gene ends 2 nt past the donor, within two 1-nt bins, so a
    # transcript may end there and the donor is left out.
    nearby = TWO_EXON_GENE + (
        'chr1\ttest\texon\t51\t102\t.\t+\t.\tgene_id "h"; transcript_id "h1";\n'
    )
    _, summary, accepted, rejected = exon_end_cluster(tmp_path, resolved_params, nearby)
    assert [int(row["coordinate"]) for row in accepted] == [110, 300]
    assert rejected == []
    assert summary["internal_exon_end_rejected"] == "0"


def test_proximal_cluster_cli_runs_without_an_annotation_when_the_filter_is_off(
    tmp_path: Path, resolved_params
) -> None:
    _, summary, accepted, rejected = exon_end_cluster(
        tmp_path, resolved_params, None, internal_exon_end_filter=False
    )
    assert [int(row["coordinate"]) for row in accepted] == [110, 300]
    assert rejected == []
    assert summary["internal_exon_end_filter_enabled"] == "False"
    assert summary["internal_exon_end_rejected"] == "0"


def test_merge_comparison_family_statistics(tmp_path: Path) -> None:
    family_a = tmp_path / "family-1"
    family_b = tmp_path / "family-2"
    family_a.mkdir()
    family_b.mkdir()
    write_tsv(
        [{"gene_id": "g1", "precision": 10, "family": "control_a"}],
        family_a / "gene_precision.tsv.gz",
    )
    write_tsv(
        [{"gene_id": "g2", "precision": 20, "family": "control_b"}],
        family_b / "gene_precision.tsv.gz",
    )
    write_tsv(
        [{"gene_id": "g1", "pac_id": "p1", "delta_pau": 0.25}],
        family_a / "fitted_pau.tsv.gz",
    )
    write_tsv(
        [{"gene_id": "g2", "pac_id": "p2", "delta_pau": -0.30}],
        family_b / "fitted_pau.tsv.gz",
    )
    write_tsv(
        [{"gene_id": "g1", "pvalue": 0.01}],
        family_a / "treatment_a_vs_control_a.genes.tsv.gz",
    )
    write_tsv(
        [{"gene_id": "g2", "pvalue": 0.02}],
        family_b / "treatment_b_vs_control_b.genes.tsv.gz",
    )

    output = tmp_path / "merged"
    assert (
        main(
            [
                "merge-statistics",
                "--inputs",
                str(family_a),
                str(family_b),
                "--output-dir",
                str(output),
            ]
        )
        == 0
    )
    assert read_tsv(output / "gene_precision.tsv.gz") == [
        {"gene_id": "g1", "precision": "10", "family": "control_a"},
        {"gene_id": "g2", "precision": "20", "family": "control_b"},
    ]
    assert [row["pac_id"] for row in read_tsv(output / "fitted_pau.tsv.gz")] == [
        "p1",
        "p2",
    ]
    assert (output / "treatment_a_vs_control_a.genes.tsv.gz").is_file()
    assert (output / "treatment_b_vs_control_b.genes.tsv.gz").is_file()


def test_kmer_enrichment_tolerates_unavailable_model_statistics(
    tmp_path: Path, resolved_params
) -> None:
    statistics = tmp_path / "statistics"
    statistics.mkdir()
    write_tsv(
        [
            {
                "pac_id": "pac-event",
                "event_type": "gained",
                "pvalue_pac": None,
                "gene_fdr": None,
                "delta_pau_ci_low": None,
                "delta_pau_ci_high": None,
                "bootstrap_successes": None,
            },
            {
                "pac_id": "pac-background",
                "event_type": None,
                "pvalue_pac": None,
                "gene_fdr": None,
                "delta_pau_ci_low": None,
                "delta_pau_ci_high": None,
                "bootstrap_successes": None,
            },
        ],
        statistics / "treatment_vs_control.events.tsv.gz",
    )
    write_tsv(
        [
            {"pac_id": "pac-event"},
            {"pac_id": "pac-background"},
        ],
        statistics / "treatment_vs_control.pacs.tsv.gz",
    )
    atlas = tmp_path / "atlas.tsv.gz"
    write_tsv(
        [
            {
                "pac_id": "pac-event",
                "gene_id": "gene-1",
                "upstream_sequence": "AAAT",
                "known_rescue_only": False,
                "candidate_status": "",
                "ambiguous_gene_assignment": False,
            },
            {
                "pac_id": "pac-background",
                "gene_id": "gene-1",
                "upstream_sequence": "GGGT",
                "known_rescue_only": False,
                "candidate_status": "",
                "ambiguous_gene_assignment": False,
            },
        ],
        atlas,
    )

    params = tmp_path / "resolved_params.yaml"
    resolved_params(params, motif_kmer_length=2)
    output = tmp_path / "kmer"
    arguments = ["kmer-enrichment", "--statistics", str(statistics), "--atlas", str(atlas)]
    assert main([*arguments, "--params", str(params), "--output-dir", str(output)]) == 0
    assert read_tsv(output / "kmer_enrichment_status.tsv") == [
        {"comparison": "treatment_vs_control", "tested_kmers": "4"}
    ]
    assert (output / "treatment_vs_control.kmer_enrichment.tsv.gz").is_file()

    resolved_params(params, run_kmer_enrichment=False)
    skipped = tmp_path / "skipped"
    assert main([*arguments, "--params", str(params), "--output-dir", str(skipped)]) == 0
    assert read_tsv(skipped / "kmer_enrichment_status.tsv") == [{"status": "not_requested"}]
    assert not list(skipped.glob("*.kmer_enrichment.tsv.gz"))


def test_kmer_enrichment_restricts_background_to_tested_pacs(
    tmp_path: Path, resolved_params
) -> None:
    """Background PACs not in the .pacs file should be excluded from k-mer analysis."""
    statistics = tmp_path / "statistics"
    statistics.mkdir()
    write_tsv(
        [
            {
                "pac_id": "pac-event",
                "event_type": "gained",
                "pvalue_pac": None,
                "gene_fdr": None,
                "delta_pau_ci_low": None,
                "delta_pau_ci_high": None,
                "bootstrap_successes": None,
            },
            {
                "pac_id": "pac-background",
                "event_type": None,
                "pvalue_pac": None,
                "gene_fdr": None,
                "delta_pau_ci_low": None,
                "delta_pau_ci_high": None,
                "bootstrap_successes": None,
            },
        ],
        statistics / "treatment_vs_control.events.tsv.gz",
    )
    write_tsv(
        [
            {"pac_id": "pac-event"},
            {"pac_id": "pac-background"},
        ],
        statistics / "treatment_vs_control.pacs.tsv.gz",
    )

    atlas = tmp_path / "atlas.tsv.gz"
    write_tsv(
        [
            {
                "pac_id": "pac-event",
                "gene_id": "gene-1",
                "upstream_sequence": "AAAT",
                "known_rescue_only": False,
                "candidate_status": "",
                "ambiguous_gene_assignment": False,
            },
            {
                "pac_id": "pac-background",
                "gene_id": "gene-1",
                "upstream_sequence": "GGGT",
                "known_rescue_only": False,
                "candidate_status": "",
                "ambiguous_gene_assignment": False,
            },
            {
                "pac_id": "pac-untested",
                "gene_id": "gene-1",
                "upstream_sequence": "TTTT",
                "known_rescue_only": False,
                "candidate_status": "",
                "ambiguous_gene_assignment": False,
            },
        ],
        atlas,
    )

    params = tmp_path / "resolved_params.yaml"
    resolved_params(params, motif_kmer_length=2)

    def enrichment(atlas_path: Path, name: str) -> tuple[list[dict], list[dict]]:
        output = tmp_path / name
        arguments = ["kmer-enrichment", "--statistics", str(statistics), "--atlas", str(atlas_path)]
        assert main([*arguments, "--params", str(params), "--output-dir", str(output)]) == 0
        return (
            read_tsv(output / "kmer_enrichment_status.tsv"),
            read_tsv(output / "treatment_vs_control.kmer_enrichment.tsv.gz"),
        )

    status, rows = enrichment(atlas, "with-untested")
    assert status == [{"comparison": "treatment_vs_control", "tested_kmers": "4"}]
    # TT occurs only in the untested PAC.
    assert "TT" not in {row["kmer"] for row in rows}
    # The untested PAC changes nothing: the result equals a run without it.
    tested_only = tmp_path / "atlas-tested.tsv.gz"
    write_tsv([row for row in read_tsv(atlas) if row["pac_id"] != "pac-untested"], tested_only)
    assert enrichment(tested_only, "tested-only") == (status, rows)


def test_kmer_enrichment_raises_on_missing_pacs_file(
    tmp_path: Path, resolved_params, capsys
) -> None:
    """Missing .pacs file should raise PacusageError."""
    statistics = tmp_path / "statistics"
    statistics.mkdir()
    write_tsv(
        [
            {
                "pac_id": "pac-event",
                "event_type": "gained",
                "pvalue_pac": None,
                "gene_fdr": None,
                "delta_pau_ci_low": None,
                "delta_pau_ci_high": None,
                "bootstrap_successes": None,
            },
            {
                "pac_id": "pac-background",
                "event_type": None,
                "pvalue_pac": None,
                "gene_fdr": None,
                "delta_pau_ci_low": None,
                "delta_pau_ci_high": None,
                "bootstrap_successes": None,
            },
        ],
        statistics / "treatment_vs_control.events.tsv.gz",
    )

    atlas = tmp_path / "atlas.tsv.gz"
    write_tsv(
        [
            {
                "pac_id": "pac-event",
                "gene_id": "gene-1",
                "upstream_sequence": "AAAT",
                "known_rescue_only": False,
                "candidate_status": "",
                "ambiguous_gene_assignment": False,
            },
            {
                "pac_id": "pac-background",
                "gene_id": "gene-1",
                "upstream_sequence": "GGGT",
                "known_rescue_only": False,
                "candidate_status": "",
                "ambiguous_gene_assignment": False,
            },
        ],
        atlas,
    )

    params = tmp_path / "resolved_params.yaml"
    resolved_params(params, motif_kmer_length=2)
    output = tmp_path / "kmer"
    arguments = ["kmer-enrichment", "--statistics", str(statistics), "--atlas", str(atlas)]

    with pytest.raises(SystemExit) as exc_info:
        main([*arguments, "--params", str(params), "--output-dir", str(output)])
    assert exc_info.value.code == 2
    assert "treatment_vs_control.pacs.tsv.gz" in capsys.readouterr().err


def test_merge_statistics_concatenates_precision_shards_with_model_status(tmp_path: Path) -> None:
    family_a = tmp_path / "family-A"
    family_b = tmp_path / "family-B"
    family_a.mkdir()
    family_b.mkdir()
    write_tsv(
        [
            {
                "gene_id": "g1",
                "precision": "10.5",
                "family": "A",
                "model_status": "drimseq",
            },
            {
                "gene_id": "g2",
                "precision": "20",
                "family": "A",
                "model_status": "drimseq_add_uniform",
            },
        ],
        family_a / "gene_precision.tsv.gz",
    )
    write_tsv(
        [
            {
                "gene_id": "g3",
                "precision": "30",
                "family": "B",
                "model_status": "fit_unavailable",
            },
        ],
        family_b / "gene_precision.tsv.gz",
    )

    output = tmp_path / "merged"
    assert (
        main(
            [
                "merge-statistics",
                "--inputs",
                str(family_a),
                str(family_b),
                "--output-dir",
                str(output),
            ]
        )
        == 0
    )
    merged = read_tsv(output / "gene_precision.tsv.gz")
    assert len(merged) == 3
    assert merged[0]["gene_id"] == "g1"
    assert merged[0]["precision"] == "10.5"
    assert merged[0]["family"] == "A"
    assert merged[0]["model_status"] == "drimseq"
    assert merged[1]["gene_id"] == "g2"
    assert merged[1]["precision"] == "20"
    assert merged[1]["family"] == "A"
    assert merged[1]["model_status"] == "drimseq_add_uniform"
    assert merged[2]["gene_id"] == "g3"
    assert merged[2]["precision"] == "30"
    assert merged[2]["family"] == "B"
    assert merged[2]["model_status"] == "fit_unavailable"


def test_merge_statistics_rejects_mismatched_precision_headers(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    family_a = tmp_path / "family-A"
    family_b = tmp_path / "family-B"
    family_a.mkdir()
    family_b.mkdir()
    write_tsv(
        [
            {
                "gene_id": "g1",
                "precision": "10.5",
                "family": "A",
                "model_status": "drimseq",
            },
            {
                "gene_id": "g2",
                "precision": "20",
                "family": "A",
                "model_status": "drimseq_add_uniform",
            },
        ],
        family_a / "gene_precision.tsv.gz",
    )
    # Family B lacks the model_status column
    write_tsv(
        [
            {"gene_id": "g3", "precision": "30", "family": "B"},
        ],
        family_b / "gene_precision.tsv.gz",
    )

    output = tmp_path / "merged"
    with pytest.raises(SystemExit) as exc_info:
        main(
            [
                "merge-statistics",
                "--inputs",
                str(family_a),
                str(family_b),
                "--output-dir",
                str(output),
            ]
        )
    assert exc_info.value.code == 2
    assert "TSV headers differ between families" in capsys.readouterr().err


def test_merge_statistics_orders_shards_by_family_not_task_number(tmp_path: Path) -> None:
    # Task numbers follow task-creation order, which can change between runs.
    shards = {"family-1": "Vehicle", "family-2": "DMSO", "family-3": "TreatmentA"}
    for directory, family in shards.items():
        (tmp_path / directory).mkdir()
        precision_row = {
            "gene_id": f"{family}_gene",
            "precision": "1",
            "family": family,
            "model_status": "drimseq",
        }
        fitted_row = {
            "gene_id": f"{family}_gene",
            "condition": f"{family}_treated",
            "control_condition": family,
        }
        write_tsv([precision_row], tmp_path / directory / "gene_precision.tsv.gz")
        write_tsv([fitted_row], tmp_path / directory / "fitted_pau.tsv.gz")

    output = tmp_path / "merged"
    arguments = ["merge-statistics", "--inputs"]
    arguments += [str(tmp_path / directory) for directory in shards]
    assert main([*arguments, "--output-dir", str(output)]) == 0
    precision = read_tsv(output / "gene_precision.tsv.gz")
    assert [row["family"] for row in precision] == ["DMSO", "TreatmentA", "Vehicle"]
    fitted = read_tsv(output / "fitted_pau.tsv.gz")
    assert [row["control_condition"] for row in fitted] == ["DMSO", "TreatmentA", "Vehicle"]


def _annotate_arguments(
    tmp_path: Path, candidates: list[PacCandidate], model: str, resolved_params
) -> list[str]:
    """Inputs for `pacusage annotate` around a one-gene reference."""
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "C" * 500 + "\n")
    pysam.faidx(str(fasta))
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\ttest\tgene\t101\t400\t.\t+\t.\tgene_id "g1"; gene_name "G1";\n'
        'chr1\ttest\texon\t101\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    params = tmp_path / "resolved_params.yaml"
    resolved_params(params, fasta=str(fasta), gtf=str(gtf))
    resolution = tmp_path / "resolution.json"
    resolution.write_text(json.dumps({"endpoint_model": model, "evidence_source": "read_3p"}))
    # CLUSTER_PACS writes an empty table, with no header, when nothing passes.
    accepted = tmp_path / "accepted_candidates.tsv.gz"
    write_tsv(candidates_as_rows(candidates), accepted)
    return [
        "annotate",
        "--candidates",
        str(accepted),
        "--reference",
        str(fasta),
        "--annotation",
        str(gtf),
        "--resolution",
        str(resolution),
        "--params",
        str(params),
        "--metadata",
        str(tmp_path / "pacs.v1.metadata.tsv.gz"),
        "--bed",
        str(tmp_path / "pacs.v1.bed.gz"),
        "--motifs",
        str(tmp_path / "pac_motifs.tsv.gz"),
        "--checksum",
        str(tmp_path / "pacs.v1.sha256"),
    ]


def test_annotate_rejects_an_empty_candidate_table(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], resolved_params
) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(_annotate_arguments(tmp_path, [], "proximal_tag", resolved_params))
    assert exc_info.value.code == 2
    error = capsys.readouterr().err
    assert "No PAC candidates passed discovery" in error
    assert "qc/pac_discovery.tsv" in error
    assert "atlas/rejected_candidates.tsv.gz" in error
    assert not (tmp_path / "pacs.v1.metadata.tsv.gz").exists()


def test_annotate_writes_an_atlas_for_one_candidate(tmp_path: Path, resolved_params) -> None:
    candidate = PacCandidate(
        contig="chr1",
        strand="+",
        coordinate=300,
        total_count=9,
        supporting_samples=2,
        capped_support=2,
        member_coordinates=(300,),
        status="primary",
    )
    assert main(_annotate_arguments(tmp_path, [candidate], "exact_boundary", resolved_params)) == 0
    metadata = read_tsv(tmp_path / "pacs.v1.metadata.tsv.gz")
    assert [row["pac_id"] for row in metadata] == ["PACv1.test.chr1.+.300"]


def test_motif_scores_cli_writes_motif_and_class_levels(tmp_path: Path, resolved_params) -> None:
    # Create a small atlas and PAU table
    atlas = tmp_path / "atlas.tsv.gz"
    write_tsv(
        [
            {
                "pac_id": "p1",
                "gene_id": "g1",
                "primary_motif_class": "other_variant",
                "known_rescue_only": False,
                "primary_pas_motif_rna": "UAUAAA",
            },
            {
                "pac_id": "p2",
                "gene_id": "g1",
                "primary_motif_class": "other_variant",
                "known_rescue_only": False,
                "primary_pas_motif_rna": "AGUAAA",
            },
        ],
        atlas,
    )
    pau = tmp_path / "pau.tsv"
    write_tsv(
        [
            {
                "gene_id": "g1",
                "pac_id": "p1",
                "sample_id": "s1",
                "count": 3,
                "gene_total": 4,
                "pau": 0.75,
            },
            {
                "gene_id": "g1",
                "pac_id": "p2",
                "sample_id": "s1",
                "count": 1,
                "gene_total": 4,
                "pau": 0.25,
            },
        ],
        pau,
    )
    params = tmp_path / "resolved_params.yaml"
    resolved_params(params, min_gene_total=1)

    motif_output = tmp_path / "motif_scores.tsv"
    class_output = tmp_path / "class_scores.tsv"

    # Run with --class-output
    assert (
        main(
            [
                "motif-scores",
                "--pau",
                str(pau),
                "--atlas",
                str(atlas),
                "--output",
                str(motif_output),
                "--class-output",
                str(class_output),
                "--params",
                str(params),
            ]
        )
        == 0
    )

    # Check motif-level output has hexamer columns
    motif_scores = read_tsv(motif_output)
    assert len(motif_scores) == 2  # Two hexamers
    assert all("primary_pas_motif_rna" in row for row in motif_scores)
    assert set(row["primary_pas_motif_rna"] for row in motif_scores) == {
        "UAUAAA",
        "AGUAAA",
    }
    expected_motif_keys = {
        "sample_id",
        "primary_motif_class",
        "motif_usage",
        "transformed_motif_usage",
        "informative_genes",
        "primary_pas_motif_rna",
    }
    assert all(set(row.keys()) >= expected_motif_keys for row in motif_scores)

    # Check class-level output has exactly the five columns
    class_scores = read_tsv(class_output)
    assert len(class_scores) == 1  # One class pooling the hexamers
    assert class_scores[0]["primary_motif_class"] == "other_variant"
    expected_class_keys = {
        "sample_id",
        "primary_motif_class",
        "motif_usage",
        "transformed_motif_usage",
        "informative_genes",
    }
    assert set(class_scores[0].keys()) == expected_class_keys
    assert class_scores[0]["informative_genes"] == "1"
