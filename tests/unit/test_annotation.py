from collections import defaultdict
from pathlib import Path

import pysam
import pytest

from pacusage.annotation import (
    FeatureIndex,
    KnownSiteIndex,
    annotate_candidates,
    find_motifs,
    internal_exon_donors,
    load_known_pacs,
)
from pacusage.models import PacCandidate
from pacusage.reference import (
    GenomicFeature,
    annotation_contigs,
    gene_name_map,
    parse_annotation,
)


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
    assert rows[0]["gene_name"] == "Plus"
    assert rows[1]["gene_name"] == "Minus"
    assert rows[0]["primary_pas_motif"] == "AATAAA"
    assert rows[0]["primary_pas_motif_rna"] == "AAUAAA"
    # An exact PAC's browser interval is its boundary base, in BED and 1-based.
    assert (rows[0]["chrom"], rows[0]["start"], rows[0]["end"]) == ("chr1", 101, 102)
    assert rows[0]["locus"] == "chr1:102-102"
    assert (rows[1]["start"], rows[1]["end"], rows[1]["locus"]) == (200, 201, "chr1:201-201")
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


def test_gene_name_map_prefers_gene_records_and_falls_back_to_ids() -> None:
    features = [
        GenomicFeature("chr1", 0, 100, "+", "exon", "g1", "t1", "ExonLineName"),
        GenomicFeature("chr1", 0, 100, "+", "gene", "g1", "", "Gene1"),
        GenomicFeature("chr1", 200, 300, "+", "exon", "g2", "t2", "Gene2"),
        GenomicFeature("chr1", 400, 500, "+", "gene", "g3", "", ""),
        GenomicFeature("chr1", 400, 500, "+", "exon", "g3", "t3", ""),
    ]
    assert gene_name_map(features) == {"g1": "Gene1", "g2": "Gene2", "g3": "g3"}


def test_gff3_exon_and_transcript_names_are_not_gene_names(tmp_path: Path) -> None:
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=g1;Name=Gene1\n"
        "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1;Name=Gene1-201\n"
        "chr1\ttest\texon\t1\t100\t.\t+\t.\tParent=t1;Name=exon-1\n"
    )
    features = parse_annotation(gff)
    names = {feature.feature_type: feature.gene_name for feature in features}
    assert names == {"gene": "Gene1", "mrna": "", "exon": ""}
    assert gene_name_map(features) == {"g1": "Gene1"}


def test_atlas_gene_names_line_up_with_gene_ids(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "C" * 400 + "\n")
    pysam.faidx(str(fasta))
    gtf = tmp_path / "genes.gtf"
    # Three genes end at the PAC. Their names sort in a different order than
    # their IDs, and g_c has no name.
    gtf.write_text(
        'chr1\ttest\tgene\t51\t201\t.\t+\t.\tgene_id "g_a"; gene_name "Zeta";\n'
        'chr1\ttest\texon\t51\t201\t.\t+\t.\tgene_id "g_a"; transcript_id "a1";\n'
        'chr1\ttest\tgene\t101\t201\t.\t+\t.\tgene_id "g_b"; gene_name "Alpha";\n'
        'chr1\ttest\texon\t101\t201\t.\t+\t.\tgene_id "g_b"; transcript_id "b1";\n'
        'chr1\ttest\tgene\t101\t201\t.\t+\t.\tgene_id "g_c";\n'
        'chr1\ttest\texon\t101\t201\t.\t+\t.\tgene_id "g_c"; transcript_id "c1";\n'
    )
    rows = annotate_candidates(
        [PacCandidate("chr1", "+", 201, 10, 2, 6, (201,))],
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
    assert rows[0]["gene_id"] == "g_a,g_b,g_c"
    assert rows[0]["gene_name"] == "Zeta,Alpha,g_c"
    assert rows[0]["ambiguous_gene_assignment"] is True
    # A PAC shared by several genes belongs to no one gene's last exon.
    assert rows[0]["last_exon"] == ""


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

    # A transcript's final exon is terminal unless it overlaps an internal
    # exon of another transcript of the gene, or it is a single-exon model
    # apart from every exon of the gene's multi-exon transcripts.
    def overlap(first: GenomicFeature, second: GenomicFeature) -> bool:
        return first.start < second.end and second.start < first.end

    by_gene: dict[tuple[str, str, str], list[list[GenomicFeature]]] = defaultdict(list)
    for exons_list in transcript_exons.values():
        first = exons_list[0]
        by_gene[(first.contig, first.strand, first.gene_id)].append(exons_list)
    for transcripts in by_gene.values():
        finals = [
            max(exons_list, key=lambda item: item.end)
            if exons_list[0].strand == "+"
            else min(exons_list, key=lambda item: item.start)
            for exons_list in transcripts
        ]
        spliced = [exon for exons_list in transcripts if len(exons_list) > 1 for exon in exons_list]
        for index, final in enumerate(finals):
            other_internal = [
                exon
                for other, exons_list in enumerate(transcripts)
                if other != index
                for exon in exons_list
                if exon is not finals[other]
            ]
            if any(overlap(final, exon) for exon in other_internal):
                continue
            if len(transcripts[index]) == 1 and spliced and not any(
                overlap(final, exon) for exon in spliced
            ):
                continue
            terminal_exons.add((final.contig, final.strand, final.start, final.end, final.gene_id))

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
    # A retained-intron model whose final exon covers an internal exon, and a
    # single-exon fragment inside an intron: neither ends the gene.
    ("gene", 330, 390, "+", "g11", ""),
    ("exon", 330, 340, "+", "g11", "g11_t1"),
    ("exon", 350, 360, "+", "g11", "g11_t1"),
    ("exon", 370, 390, "+", "g11", "g11_t1"),
    ("exon", 330, 340, "+", "g11", "g11_retained"),
    ("exon", 350, 390, "+", "g11", "g11_retained"),
    ("exon", 363, 367, "+", "g11", "g11_fragment"),
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


def _index(rows: list[tuple[str, int, int, str, str, str]], tmp_path: Path) -> FeatureIndex:
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(_gtf(rows))
    return FeatureIndex(parse_annotation(gtf))


def test_retained_intron_and_fragment_ends_are_not_terminal_exons(tmp_path: Path) -> None:
    index = _index(
        [
            ("gene", 101, 1000, "+", "g1", ""),
            ("exon", 101, 200, "+", "g1", "main"),
            ("exon", 401, 500, "+", "g1", "main"),
            ("exon", 801, 1000, "+", "g1", "main"),
            # Its final exon keeps the intron after exon 2.
            ("exon", 101, 200, "+", "g1", "retained"),
            ("exon", 401, 1000, "+", "g1", "retained"),
            # 3'-incomplete: it stops inside exon 2.
            ("exon", 101, 200, "+", "g1", "incomplete"),
            ("exon", 401, 450, "+", "g1", "incomplete"),
            # A single-exon model inside intron 1.
            ("exon", 251, 300, "+", "g1", "fragment"),
        ],
        tmp_path,
    )

    def classify(coordinate: int) -> tuple[str, list[tuple[int, int]], str]:
        assignments = index.assign("chr1", "+", coordinate, 100)
        spans = sorted((feature.start, feature.end) for _, feature in assignments)
        return assignments[0][0], spans, index.last_exon(assignments)

    # The retained intron is an exon of one model, not a transcript end.
    assert classify(700) == ("other_exon", [(400, 1000)], "")
    assert classify(450) == ("other_exon", [(400, 450), (400, 500), (400, 1000)], "")
    assert classify(280) == ("other_exon", [(250, 300)], "")
    assert classify(1000) == ("terminal_exon", [(800, 1000)], "chr1:801-1000")
    assert classify(350) == ("intronic", [(100, 1000)], "")


def test_last_exons_merge_overlapping_ends_and_separate_alternative_ones(
    tmp_path: Path,
) -> None:
    index = _index(
        [
            ("gene", 101, 900, "+", "g1", ""),
            ("exon", 101, 200, "+", "g1", "short"),
            ("exon", 401, 500, "+", "g1", "short"),
            # A longer 3' UTR on the same last exon.
            ("exon", 101, 200, "+", "g1", "long"),
            ("exon", 401, 650, "+", "g1", "long"),
            # A last exon that starts where the longer one ends.
            ("exon", 101, 200, "+", "g1", "touching"),
            ("exon", 651, 700, "+", "g1", "touching"),
            # An alternative last exon.
            ("exon", 101, 200, "+", "g1", "alternative"),
            ("exon", 801, 900, "+", "g1", "alternative"),
        ],
        tmp_path,
    )

    def last_exon(coordinate: int) -> str:
        return index.last_exon(index.assign("chr1", "+", coordinate, 100))

    assert last_exon(500) == "chr1:401-700"
    assert last_exon(650) == "chr1:401-700"
    assert last_exon(700) == "chr1:401-700"
    assert last_exon(900) == "chr1:801-900"
    # A PAC past the gene's end belongs to its 3'-most last exon.
    assert index.assign("chr1", "+", 950, 100)[0][0] == "downstream"
    assert last_exon(950) == "chr1:801-900"


def test_minus_strand_downstream_pacs_take_the_three_prime_last_exon(tmp_path: Path) -> None:
    index = _index(
        [
            ("gene", 2001, 2600, "-", "g2", ""),
            ("exon", 2501, 2600, "-", "g2", "proximal"),
            ("exon", 2201, 2300, "-", "g2", "proximal"),
            ("exon", 2501, 2600, "-", "g2", "distal"),
            ("exon", 2001, 2100, "-", "g2", "distal"),
        ],
        tmp_path,
    )
    downstream = index.assign("chr1", "-", 1950, 100)
    assert downstream[0][0] == "downstream"
    assert index.last_exon(downstream) == "chr1:2001-2100"
    assert index.last_exon(index.assign("chr1", "-", 2200, 100)) == "chr1:2201-2300"
    assert index.last_exon(index.assign("chr1", "-", 2550, 100)) == ""


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


def exon(
    start: int, end: int, transcript_id: str, gene_id: str = "g", strand: str = "+"
) -> GenomicFeature:
    return GenomicFeature("chr1", start, end, strand, "exon", gene_id, transcript_id)


def test_internal_exon_donors_are_the_three_prime_ends_of_spliced_exons() -> None:
    exons = [
        # t1 has three exons; t2 skips the middle one and ends further on.
        exon(100, 200, "t1"), exon(300, 400, "t1"), exon(500, 600, "t1"),
        exon(100, 200, "t2"), exon(500, 650, "t2"),
        # Single-exon genes have no donors.
        exon(2000, 2100, "s1", "single"),
        # On the minus strand, the 3' end of an exon spliced onward is its start.
        exon(1000, 1100, "m1", "minus", "-"), exon(1300, 1400, "m1", "minus", "-"),
    ]
    assert internal_exon_donors(exons, 25) == {("chr1", "+"): [200, 400], ("chr1", "-"): [1300]}


def test_internal_exon_donors_keep_clear_of_transcript_ends() -> None:
    exons = [exon(100, 200, "t1"), exon(300, 400, "t1"), exon(500, 600, "t1")]
    # Another gene on the strand ends 10 nt past the donor at 400, so a
    # transcript may end there.
    ending = exons + [exon(300, 410, "o1", "other")]
    assert internal_exon_donors(ending, 25) == {("chr1", "+"): [200]}
    assert internal_exon_donors(ending, 5) == {("chr1", "+"): [200, 400]}
    # A gene on the other strand does not protect it.
    opposite = exons + [exon(300, 410, "o1", "other", "-")]
    assert internal_exon_donors(opposite, 25) == {("chr1", "+"): [200, 400]}


def test_genes_without_transcript_ids_have_no_donors() -> None:
    # Without transcript IDs, alternative last exons would look spliced onward.
    exons = [exon(100, 200, ""), exon(300, 400, ""), exon(500, 600, "")]
    assert internal_exon_donors(exons, 25) == {}
