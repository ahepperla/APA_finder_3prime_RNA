from pathlib import Path

from pacusage.annotation import annotate_candidates
from pacusage.cli import main
from pacusage.clustering import cluster_exact_boundaries
from pacusage.evidence import write_evidence
from pacusage.models import EvidenceObservation
from pacusage.quantification import build_count_outputs, quantify_exact
from pacusage.reference import parse_annotation, prepare_reference
from pacusage.tableio import read_tsv


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
    for sample_id in ("a", "b"):
        observations = [
            EvidenceObservation(sample_id, "chr1", "+", 90, 3),
            EvidenceObservation(sample_id, "chr1", "-", 210, 4),
        ]
        tsv = tmp_path / f"{sample_id}.tsv.gz"
        parquet = tmp_path / f"{sample_id}.parquet"
        write_evidence(observations, tsv, parquet)
        evidence_paths.append(str(parquet))

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
