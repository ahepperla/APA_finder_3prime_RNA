from pathlib import Path

from pacusage.annotation import annotate_candidates
from pacusage.cli import main
from pacusage.clustering import cluster_exact_boundaries
from pacusage.evidence import write_evidence, write_splice_continuations
from pacusage.models import EvidenceObservation, SpliceContinuation
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
    candidates, _ = cluster_exact_boundaries(observations, 2, 12, 5, 2, 2)
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


def test_proximal_cluster_cli_streams_parquet_evidence(tmp_path: Path) -> None:
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

    params = tmp_path / "params.yaml"
    params.write_text(
        "\n".join(
            [
                "input: samples.tsv",
                "assembly: test",
                "fasta: genome.fa",
                "gtf: genes.gtf",
                "pac_min_total_count: 1",
                "pac_min_sample_count: 1",
                "pac_min_supporting_samples: 0.75",
                "proximal_bin_size: 1",
            ]
        )
        + "\n"
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


def test_proximal_cluster_cli_rejects_constitutive_readthrough(tmp_path: Path) -> None:
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

    params = tmp_path / "params.yaml"
    params.write_text(
        "\n".join(
            [
                "input: samples.tsv",
                "assembly: test",
                "fasta: genome.fa",
                "gtf: genes.gtf",
                "pac_min_total_count: 1",
                "pac_min_sample_count: 1",
                "pac_min_supporting_samples: 1",
                "proximal_bin_size: 1",
                "constitutive_readthrough_min_junction_count: 2",
                "constitutive_readthrough_min_replicate_support: all",
            ]
        )
        + "\n"
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
    assert read_tsv(qc)[0]["constitutive_readthrough_rejected"] == "1"


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
        [{"gene_id": "g1", "feature_id": "p1", "delta_pau": 0.25}],
        family_a / "fitted_pau.tsv.gz",
    )
    write_tsv(
        [{"gene_id": "g2", "feature_id": "p2", "delta_pau": -0.30}],
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
    assert [row["feature_id"] for row in read_tsv(output / "fitted_pau.tsv.gz")] == [
        "p1",
        "p2",
    ]
    assert (output / "treatment_a_vs_control_a.genes.tsv.gz").is_file()
    assert (output / "treatment_b_vs_control_b.genes.tsv.gz").is_file()
