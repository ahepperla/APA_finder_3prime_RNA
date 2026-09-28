from collections import defaultdict
from pathlib import Path

import pysam
import pytest

from pacusage.annotation import (
    FeatureIndex,
    KnownSiteIndex,
    annotate_candidates,
    find_motifs,
    load_known_pacs,
)
from pacusage.models import PacCandidate
from pacusage.reference import GenomicFeature, annotation_contigs, parse_annotation


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
    assert rows[0]["region_start"] == 101
    assert rows[0]["region_end"] == 102
    assert rows[0]["resolution_nt"] == 0
    assert rows[0]["total_supporting_samples"] == 0
    assert rows[0]["supporting_condition"] == ""
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


def test_annotation_contigs_scans_without_materializing_features(tmp_path: Path) -> None:
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "# comment\n"
        'chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "g1";\n'
        'chr2\ttest\texon\t1\t100\t.\t-\t.\tgene_id "g2"; transcript_id "t2";\n'
        'chr3\ttest\tregion\t1\t100\t.\t.\t.\tID=region3\n'
    )
    assert annotation_contigs(gtf) == {"chr1", "chr2"}


def _reference_assign(
    contig: str,
    strand: str,
    coordinate: int,
    maximum_downstream: int,
    features: list[GenomicFeature],
) -> list[tuple[str, GenomicFeature]]:
    """Reference implementation of FeatureIndex.assign for equivalence testing."""
    # Extract genes and exons
    genes: dict[tuple[str, str, str], GenomicFeature] = {}
    exons: list[GenomicFeature] = []
    terminal_exons: set[tuple[str, str, int, int, str]] = set()

    for feature in features:
        if feature.feature_type == "gene":
            genes[(feature.contig, feature.strand, feature.gene_id)] = feature
        elif feature.feature_type == "exon":
            exons.append(feature)

    # Build genes from exons if needed
    if not genes:
        gene_parts: dict[tuple[str, str, str], list[GenomicFeature]] = defaultdict(list)
        for exon in exons:
            gene_parts[(exon.contig, exon.strand, exon.gene_id)].append(exon)
        for key, values in gene_parts.items():
            first = values[0]
            genes[key] = GenomicFeature(
                contig=first.contig,
                start=min(item.start for item in values),
                end=max(item.end for item in values),
                strand=first.strand,
                feature_type="gene",
                gene_id=first.gene_id,
                gene_name=first.gene_name,
            )

    # Find terminal exons
    transcript_exons: dict[tuple[str, str, str], list[GenomicFeature]] = defaultdict(list)
    for exon in exons:
        key = (exon.contig, exon.strand, exon.transcript_id or exon.gene_id)
        transcript_exons[key].append(exon)

    for exons_list in transcript_exons.values():
        terminal = (
            max(exons_list, key=lambda item: item.end)
            if exons_list[0].strand == "+"
            else min(exons_list, key=lambda item: item.start)
        )
        terminal_exons.add(
            (terminal.contig, terminal.strand, terminal.start, terminal.end, terminal.gene_id)
        )

    # Implement the matching logic (inclusive rule)
    exon_matches = [
        feature
        for feature in exons
        if feature.contig == contig
        and feature.strand == strand
        and feature.start <= coordinate <= feature.end
    ]
    terminal = [
        feature
        for feature in exon_matches
        if (feature.contig, feature.strand, feature.start, feature.end, feature.gene_id)
        in terminal_exons
    ]
    if terminal:
        return [("terminal_exon", item) for item in terminal]
    if exon_matches:
        return [("other_exon", item) for item in exon_matches]

    intronic = [
        feature
        for key, feature in genes.items()
        if key[0] == contig
        and key[1] == strand
        and feature.start <= coordinate <= feature.end
    ]
    if intronic:
        return [("intronic", item) for item in intronic]

    # Find downstream genes
    same_strand_genes = [
        feature for key, feature in genes.items() if key[0] == contig and key[1] == strand
    ]
    downstream = []
    for gene in same_strand_genes:
        distance = coordinate - gene.end if strand == "+" else gene.start - coordinate
        if 0 < distance <= maximum_downstream:
            intervening = any(
                other.gene_id != gene.gene_id
                and (
                    gene.end < other.start <= coordinate
                    if strand == "+"
                    else coordinate <= other.end < gene.start
                )
                for other in same_strand_genes
            )
            if not intervening:
                downstream.append(("downstream", gene))

    if downstream:
        nearest = min(
            abs(coordinate - (item[1].end if strand == "+" else item[1].start))
            for item in downstream
        )
        return [
            item
            for item in downstream
            if abs(coordinate - (item[1].end if strand == "+" else item[1].start)) == nearest
        ]
    return []


def _gtf(rows: list[tuple[str, int, int, str, str, str]]) -> str:
    """GTF text from (feature, start, end, strand, gene_id, transcript_id), 1-based."""
    lines = []
    for feature, start, end, strand, gene_id, transcript_id in rows:
        attributes = f'gene_id "{gene_id}";'
        if transcript_id:
            attributes += f' transcript_id "{transcript_id}";'
        lines.append(f"chr1\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n")
    return "".join(lines)


WITH_GENE_RECORDS = [
    # Overlapping plus-strand genes, one with two transcripts sharing exons,
    # spanning several 50-nt bins.
    ("gene", 10, 100, "+", "g1", ""),
    ("exon", 10, 40, "+", "g1", "g1_t1"),
    ("exon", 60, 100, "+", "g1", "g1_t1"),
    ("exon", 10, 50, "+", "g1", "g1_t2"),
    ("exon", 70, 100, "+", "g1", "g1_t2"),
    ("gene", 80, 150, "+", "g2", ""),
    ("exon", 80, 150, "+", "g2", "g2_t1"),
    # Two genes ending at the same coordinate, then a gene between an upstream
    # gene and the positions after it.
    ("gene", 170, 200, "+", "g3", ""),
    ("exon", 170, 200, "+", "g3", "g3_t1"),
    ("gene", 185, 200, "+", "g4", ""),
    ("exon", 185, 200, "+", "g4", "g4_t1"),
    ("gene", 230, 240, "+", "g5", ""),
    ("exon", 230, 240, "+", "g5", "g5_t1"),
    # Exons of a gene without a gene record: not a gene here.
    ("exon", 300, 320, "+", "g6", "g6_t1"),
    # Minus strand, including shared starts and a gene nested downstream.
    ("gene", 110, 200, "-", "g7", ""),
    ("exon", 110, 150, "-", "g7", "g7_t1"),
    ("exon", 160, 200, "-", "g7", "g7_t1"),
    ("gene", 130, 250, "-", "g8", ""),
    ("exon", 130, 250, "-", "g8", "g8_t1"),
    ("gene", 280, 300, "-", "g9", ""),
    ("exon", 280, 300, "-", "g9", "g9_t1"),
    ("gene", 280, 290, "-", "g10", ""),
    ("exon", 280, 290, "-", "g10", "g10_t1"),
]

WITHOUT_GENE_RECORDS = [row for row in WITH_GENE_RECORDS if row[0] == "exon"]


@pytest.mark.parametrize("rows", [WITH_GENE_RECORDS, WITHOUT_GENE_RECORDS])
def test_feature_index_matches_brute_force_assignment(tmp_path: Path, rows) -> None:
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(_gtf(rows))
    features = parse_annotation(gtf)
    index = FeatureIndex(features, bin_size=50)

    def key(assignments):
        return sorted(
            (label, feature.gene_id, feature.start, feature.end) for label, feature in assignments
        )

    classes = set()
    for coordinate in range(0, 400):
        for strand in ("+", "-"):
            for maximum_downstream in (5, 30, 100):
                expected = _reference_assign(
                    "chr1", strand, coordinate, maximum_downstream, features
                )
                observed = index.assign("chr1", strand, coordinate, maximum_downstream)
                assert key(observed) == key(expected), (coordinate, strand, maximum_downstream)
                classes.update(label for label, _ in expected)
    assert classes == {"terminal_exon", "other_exon", "intronic", "downstream"}


def test_known_site_index_matches() -> None:
    """Test KnownSiteIndex.matches with exact ±radius boundaries on both strands."""
    sites = [
        ("chr1", "+", 100),
        ("chr1", "+", 110),
        ("chr1", "-", 200),
        ("chr1", "-", 210),
        ("chr2", "+", 50),
    ]
    index = KnownSiteIndex(sites)

    # Test exact matches at boundaries
    assert index.matches("chr1", "+", 100, 0) == [100]
    assert index.matches("chr1", "+", 100, 5) == [100]
    assert index.matches("chr1", "+", 100, 10) == [100, 110]
    assert index.matches("chr1", "+", 105, 5) == [100, 110]
    assert index.matches("chr1", "+", 111, 0) == []
    assert index.matches("chr1", "+", 111, 1) == [110]

    # Test minus strand
    assert index.matches("chr1", "-", 200, 10) == [200, 210]
    assert index.matches("chr1", "-", 205, 5) == [200, 210]

    # Test other contigs
    assert index.matches("chr2", "+", 50, 0) == [50]
    assert index.matches("chr3", "+", 100, 10) == []

    # Test deduplication
    dup_index = KnownSiteIndex([("chr1", "+", 100), ("chr1", "+", 100)])
    assert dup_index.matches("chr1", "+", 100, 0) == [100]


def test_known_site_index_bool() -> None:
    """Test KnownSiteIndex.__bool__."""
    assert not KnownSiteIndex([])
    assert KnownSiteIndex([("chr1", "+", 100)])


def test_load_known_pacs_from_bed(tmp_path: Path) -> None:
    """Test load_known_pacs on a small BED6 file.

    Validates that + strand uses field 2 (end) and - strand uses field 1 (start).
    """
    bed = tmp_path / "known.bed"
    bed.write_text(
        "chr1\t95\t101\tpac1\t100\t+\n"
        "chr1\t199\t205\tpac2\t100\t-\n"
        "chr2\t45\t51\tpac3\t100\t+\n"
    )
    index = load_known_pacs(bed)

    # Plus strand: uses end (field 2)
    assert index.matches("chr1", "+", 101, 0) == [101]
    # Minus strand: uses start (field 1)
    assert index.matches("chr1", "-", 199, 0) == [199]
    # Other contig
    assert index.matches("chr2", "+", 51, 0) == [51]

    # Empty path
    empty_index = load_known_pacs(None)
    assert not empty_index
