from pathlib import Path

from pacusage.annotation import annotate_candidates
from pacusage.clustering import cluster_exact_boundaries
from pacusage.models import EvidenceObservation
from pacusage.quantification import build_count_outputs, quantify_exact
from pacusage.reference import parse_annotation, prepare_reference


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
