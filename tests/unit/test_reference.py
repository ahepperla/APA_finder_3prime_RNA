"""GFF3 parentage in the annotation parser, in the Ensembl and NCBI styles."""

from __future__ import annotations

from pathlib import Path

import pytest

from pacusage.annotation import FeatureIndex
from pacusage.errors import PacusageError
from pacusage.evidence import ExonBinIndex
from pacusage.reference import gene_name_map, parse_annotation, transcript_ends


def gff3(source: str, rows: list[tuple[str, int, int, str, str]]) -> str:
    lines = ["##gff-version 3"]
    for feature_type, start, end, strand, attributes in rows:
        fields = ["chr1", source, feature_type, str(start), str(end), ".", strand, ".", attributes]
        lines.append("\t".join(fields))
    return "\n".join(lines) + "\n"


# Ensembl GFF3 prefixes gene and transcript IDs and gives exons only a Parent.
# Gene A (+) has two transcripts ending at 600 and 800 and three exons; gene B
# (-) is an lncRNA, which Ensembl writes as ncRNA_gene with an lnc_RNA child.
ENSEMBL_GFF3 = gff3(
    "ensembl",
    [
        ("gene", 101, 800, "+", "ID=gene:ENSG01;Name=GENEA;gene_id=ENSG01"),
        ("mRNA", 101, 600, "+",
         "ID=transcript:ENST01;Parent=gene:ENSG01;Name=GENEA-201;transcript_id=ENST01"),
        ("exon", 101, 200, "+", "Parent=transcript:ENST01;Name=ENSE01;exon_id=ENSE01"),
        ("exon", 301, 600, "+", "Parent=transcript:ENST01;Name=ENSE02;exon_id=ENSE02"),
        ("mRNA", 101, 800, "+",
         "ID=transcript:ENST02;Parent=gene:ENSG01;Name=GENEA-202;transcript_id=ENST02"),
        ("exon", 101, 200, "+", "Parent=transcript:ENST02;Name=ENSE01;exon_id=ENSE01"),
        ("exon", 301, 800, "+", "Parent=transcript:ENST02;Name=ENSE03;exon_id=ENSE03"),
        ("ncRNA_gene", 2001, 3000, "-", "ID=gene:ENSG02;Name=GENEB;biotype=lncRNA;gene_id=ENSG02"),
        ("lnc_RNA", 2101, 2900, "-",
         "ID=transcript:ENST03;Parent=gene:ENSG02;Name=GENEB-201;transcript_id=ENST03"),
        ("exon", 2101, 2300, "-", "Parent=transcript:ENST03;Name=ENSE04;exon_id=ENSE04"),
        ("exon", 2501, 2900, "-", "Parent=transcript:ENST03;Name=ENSE05;exon_id=ENSE05"),
    ],
)

# NCBI GFF3 carries the gene symbol as gene= on most records. The second exon
# has neither gene= nor transcript_id=, so only its Parent chain names them.
NCBI_GFF3 = gff3(
    "RefSeq",
    [
        ("gene", 1001, 2000, "+", "ID=gene-OR4F5;Name=OR4F5;gbkey=Gene;gene=OR4F5"),
        ("mRNA", 1001, 2000, "+",
         "ID=rna-NM_01.2;Parent=gene-OR4F5;gene=OR4F5;product=olfactory receptor;"
         "transcript_id=NM_01.2"),
        ("exon", 1001, 1100, "+",
         "ID=exon-NM_01.2-1;Parent=rna-NM_01.2;gene=OR4F5;transcript_id=NM_01.2"),
        ("exon", 1501, 2000, "+", "ID=exon-NM_01.2-2;Parent=rna-NM_01.2;gbkey=mRNA"),
    ],
)


def write(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "genes.gff3"
    path.write_text(text)
    return path


def test_ensembl_exons_resolve_their_transcript_and_gene(tmp_path: Path) -> None:
    features = parse_annotation(write(tmp_path, ENSEMBL_GFF3))
    exons = [
        (feature.gene_id, feature.transcript_id, feature.start, feature.end)
        for feature in features
        if feature.feature_type == "exon"
    ]
    assert exons == [
        ("ENSG01", "ENST01", 100, 200),
        ("ENSG01", "ENST01", 300, 600),
        ("ENSG01", "ENST02", 100, 200),
        ("ENSG01", "ENST02", 300, 800),
        ("ENSG02", "ENST03", 2100, 2300),
        ("ENSG02", "ENST03", 2500, 2900),
    ]
    transcripts = {
        (feature.feature_type, feature.transcript_id): feature.gene_id
        for feature in features
        if feature.feature_type in {"mrna", "lnc_rna"}
    }
    assert transcripts == {
        ("mrna", "ENST01"): "ENSG01",
        ("mrna", "ENST02"): "ENSG01",
        ("lnc_rna", "ENST03"): "ENSG02",
    }
    # Ensembl's ncRNA_gene is a gene record, and only gene records' Name
    # values are gene names: exon and transcript Name values are not.
    genes = [
        (feature.gene_id, feature.gene_name)
        for feature in features
        if feature.feature_type == "gene"
    ]
    assert genes == [("ENSG01", "GENEA"), ("ENSG02", "GENEB")]
    assert all(feature.gene_name == "" for feature in features if feature.feature_type != "gene")
    assert gene_name_map(features) == {"ENSG01": "GENEA", "ENSG02": "GENEB"}


def test_ensembl_transcripts_give_calibration_ends(tmp_path: Path) -> None:
    ends = transcript_ends(parse_annotation(write(tmp_path, ENSEMBL_GFF3)))
    # Gene A's transcripts end at 600 and 800, so it has no single end; they
    # would merge into one end at 800 if the exons lost their transcripts.
    assert dict(ends) == {("chr1", "-"): [(2100, "ENSG02")]}


def test_ensembl_genes_are_assigned_and_indexed(tmp_path: Path) -> None:
    features = parse_annotation(write(tmp_path, ENSEMBL_GFF3), {"gene", "exon"})
    index = FeatureIndex(features)

    def assigned(strand: str, coordinate: int) -> list[tuple[str, str]]:
        return [
            (label, feature.gene_id)
            for label, feature in index.assign("chr1", strand, coordinate, 5000)
        ]

    assert assigned("+", 800) == [("last_exon", "ENSG01")]
    assert assigned("+", 250) == [("intron", "ENSG01")]
    assert assigned("-", 2100) == [("last_exon", "ENSG02")]
    # The lncRNA's intron is inside its ncRNA_gene record.
    assert assigned("-", 2400) == [("intron", "ENSG02")]
    exons = ExonBinIndex(features).query("chr1", 150, 160)
    assert [feature.gene_id for feature in exons] == ["ENSG01"]


def test_ncbi_exons_resolve_through_gene_symbols_or_parents(tmp_path: Path) -> None:
    features = parse_annotation(write(tmp_path, NCBI_GFF3))
    records = [
        (feature.feature_type, feature.gene_id, feature.transcript_id) for feature in features
    ]
    assert records == [
        ("gene", "OR4F5", ""),
        ("mrna", "OR4F5", "NM_01.2"),
        ("exon", "OR4F5", "NM_01.2"),
        ("exon", "OR4F5", "NM_01.2"),
    ]
    assert gene_name_map(features) == {"OR4F5": "OR4F5"}


def test_a_gff3_exon_whose_parent_is_missing_is_rejected(tmp_path: Path) -> None:
    path = write(
        tmp_path,
        "##gff-version 3\n"
        "chr1\tensembl\texon\t101\t200\t.\t+\t.\tParent=transcript:ENST09;exon_id=ENSE09\n",
    )
    with pytest.raises(PacusageError, match="feature lacks a gene identifier"):
        parse_annotation(path)


def test_a_gff3_parent_cycle_is_rejected(tmp_path: Path) -> None:
    path = write(
        tmp_path,
        "##gff-version 3\n"
        "chr1\tsource\tmRNA\t101\t200\t.\t+\t.\tID=t1;Parent=t2\n"
        "chr1\tsource\tmRNA\t101\t200\t.\t+\t.\tID=t2;Parent=t1\n"
        "chr1\tsource\texon\t101\t200\t.\t+\t.\tParent=t1\n",
    )
    with pytest.raises(PacusageError, match="feature lacks a gene identifier"):
        parse_annotation(path)
