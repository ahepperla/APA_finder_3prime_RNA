from pathlib import Path

import pysam

from pacusage.annotation import annotate_candidates, find_motifs
from pacusage.models import PacCandidate
from pacusage.reference import parse_annotation


def test_motif_position_and_minus_strand_annotation(tmp_path: Path) -> None:
    sequence = list("C" * 300)
    sequence[70:76] = list("AATAAA")
    sequence[224:230] = list("TTTATT")
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "".join(sequence) + "\n")
    pysam.faidx(str(fasta))
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\ttest\tgene\t51\t101\t.\t+\t.\tgene_id "plus"; gene_name "Plus";\n'
        'chr1\ttest\texon\t51\t101\t.\t+\t.\tgene_id "plus"; transcript_id "p1";\n'
        'chr1\ttest\tgene\t201\t251\t.\t-\t.\tgene_id "minus"; gene_name "Minus";\n'
        'chr1\ttest\texon\t201\t251\t.\t-\t.\tgene_id "minus"; transcript_id "m1";\n'
    )
    candidates = [
        PacCandidate("chr1", "+", 101, 10, 2, 6, (101,)),
        PacCandidate("chr1", "-", 200, 10, 2, 6, (200,)),
    ]
    rows = annotate_candidates(
        candidates,
        parse_annotation(gtf),
        fasta,
        "test",
        "exact_boundary",
        "read_3p",
        100,
        50,
        5,
        35,
        10,
        20,
        6,
        0.6,
    )
    assert rows[0]["gene_id"] == "plus"
    assert rows[0]["primary_pas_motif"] == "AATAAA"
    assert rows[0]["primary_pas_motif_rna"] == "AAUAAA"
    assert rows[1]["primary_pas_motif_rna"] == "AAUAAA"
    assert rows[1]["gene_id"] == "minus"


def test_all_motif_matches_are_retained() -> None:
    matches = find_motifs("AATAAACCCAATAAA", 20, 18, 5)
    assert len([match for match in matches if match[0] == "AATAAA"]) == 2


def test_gff3_parent_relationships_are_resolved(tmp_path: Path) -> None:
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=g1;Name=Gene1\n"
        "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1\t100\t.\t+\t.\tParent=t1\n"
    )
    features = parse_annotation(gff)
    exon = next(feature for feature in features if feature.feature_type == "exon")
    assert exon.gene_id == "g1"
    assert exon.transcript_id == "t1"
