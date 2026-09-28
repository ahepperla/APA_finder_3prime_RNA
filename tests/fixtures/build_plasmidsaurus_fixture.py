#!/usr/bin/env python3
"""Generate the Plasmidsaurus-like proximal-tag fixture.

Plasmidsaurus 3' tag reads end upstream of the polyadenylation site, not on
it. Here every read's aligned 3' end sits 20-280 nt upstream of its PAC, at a
deterministic quantile of a triangular distribution, so each sample yields the
same calibration kernel.

Calibration compares read ends with annotated transcript ends, and only genes
with a single annotated end take part. So the 24 single-PAC genes calibrate
the kernel. The 18 multi-PAC genes annotate one transcript per PAC, keeping
their alternative ends out of calibration. Their PACs are 400 nt apart, farther
than the kernel reaches, and sit on multiples of the 25-nt discovery bin, so
each read belongs to exactly one PAC.

The fixture has its own FASTA, annotation, and sample sheet. The exact-boundary
fixture next to it is unchanged.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pysam

ROOT = Path(__file__).resolve().parent / "plasmidsaurus"
CONTIG = "chr1"
READ_LENGTH = 50
PAC_SPACING = 400
GENE_SPACING = 3000
FIRST_ORIGIN = 2000
GENE_LENGTH = 1500
# Triangular read-end offsets upstream of the PAC: minimum, mode, maximum.
OFFSET_LOW, OFFSET_MODE, OFFSET_HIGH = 20, 150, 280

# (sample_id, condition, control), in sample-sheet order.
ROWS = [
    ("DMSO_1", "DMSO", ""),
    ("DMSO_2", "DMSO", ""),
    ("TRA_1", "TreatmentA", "DMSO"),
    ("TRA_2", "TreatmentA", "DMSO"),
]
SAMPLES = [sample_id for sample_id, _, _ in ROWS]
CONDITION = {sample_id: condition for sample_id, condition, _ in ROWS}

# The test profile requires 20 calibration genes; these 24 are the only
# genes with one annotated end, so keep at least 20 of them.
SINGLE_PAC_GENES = 24
TWO_PAC_GENES = 14
THREE_PAC_GENES = 4
# Two-PAC genes whose usage moves from the distal to the proximal PAC in
# TreatmentA. Proportions are listed distal PAC first.
SHIFTED_GENES = (0, 5, 10)
SHIFT_PROPORTIONS = {"DMSO": (0.75, 0.25), "TreatmentA": (0.25, 0.75)}
# Two-PAC genes used exactly 3:2 in every sample, at different depths, so they
# must show no change. The other null genes are random draws, which can differ
# between conditions by chance.
EXACT_NULL_GENES = (7, 12)
EXACT_NULL_DEPTHS = (250, 300, 200, 350)
NULL_TWO_PAC_BASES = [(0.6, 0.4), (0.7, 0.3)]
NULL_THREE_PAC_BASE = (0.5, 0.3, 0.2)


def gene_layout() -> list[dict[str, object]]:
    """Genes in genomic order: id, strand, origin, and PACs from distal to proximal."""
    kinds = (
        [("single", index) for index in range(SINGLE_PAC_GENES)]
        + [("two", index) for index in range(TWO_PAC_GENES)]
        + [("three", index) for index in range(THREE_PAC_GENES)]
    )
    genes = []
    for position, (kind, index) in enumerate(kinds):
        origin = FIRST_ORIGIN + GENE_SPACING * position
        strand = "+" if position % 2 == 0 else "-"
        pac_count = {"single": 1, "two": 2, "three": 3}[kind]
        # Plus-strand genes end at origin + GENE_LENGTH and minus-strand genes
        # at origin; each further PAC lies PAC_SPACING upstream.
        distal = origin + GENE_LENGTH if strand == "+" else origin
        step = -PAC_SPACING if strand == "+" else PAC_SPACING
        genes.append(
            {
                "gene_id": f"{kind}{index:02d}",
                "kind": kind,
                "index": index,
                "strand": strand,
                "origin": origin,
                "pacs": [distal + step * rank for rank in range(pac_count)],
            }
        )
    return genes


def design(gene: dict[str, object]) -> str:
    """The gene's role: calibration, shifted, exact_null, or null."""
    if gene["kind"] == "single":
        return "calibration"
    if gene["kind"] == "two" and gene["index"] in SHIFTED_GENES:
        return "shifted"
    if gene["kind"] == "two" and gene["index"] in EXACT_NULL_GENES:
        return "exact_null"
    return "null"


def site_counts(random: np.random.RandomState) -> dict[tuple[str, int], dict[str, int]]:
    """Reads per sample for each (strand, PAC coordinate)."""
    sites: dict[tuple[str, int], dict[str, int]] = {}
    for gene in gene_layout():
        pacs = gene["pacs"]
        per_sample: dict[str, tuple[int, ...]] = {}
        for sample_number, sample_id in enumerate(SAMPLES):
            if gene["kind"] == "single":
                per_sample[sample_id] = (int(120 + 20 * (gene["index"] % 4) + 5 * sample_number),)
                continue
            if design(gene) == "exact_null":
                depth = EXACT_NULL_DEPTHS[sample_number]
                per_sample[sample_id] = (depth * 3 // 5, depth * 2 // 5)
                continue
            if design(gene) == "shifted":
                proportions = np.asarray(SHIFT_PROPORTIONS[CONDITION[sample_id]], dtype=float)
                depth = 280 + 20 * (sample_number % 2)
            else:
                base = (
                    NULL_TWO_PAC_BASES[gene["index"] % len(NULL_TWO_PAC_BASES)]
                    if gene["kind"] == "two"
                    else NULL_THREE_PAC_BASE
                )
                proportions = random.dirichlet(np.asarray(base, dtype=float) * 100)
                depth = int(random.negative_binomial(20, 20 / (20 + 250)))
            per_sample[sample_id] = tuple(
                int(value) for value in random.multinomial(depth, proportions)
            )
        for rank, coordinate in enumerate(pacs):
            sites[(gene["strand"], coordinate)] = {
                sample_id: per_sample[sample_id][rank] for sample_id in SAMPLES
            }
    return sites


def check_support(sites: dict[tuple[str, int], dict[str, int]]) -> None:
    """Every PAC must pass default discovery: 10 reads in total and 2 or more
    reads in both replicates of at least one condition."""
    for key, counts in sites.items():
        total = sum(counts.values())
        supported = any(
            all(counts[sample_id] >= 2 for sample_id in SAMPLES if CONDITION[sample_id] == group)
            for group in set(CONDITION.values())
        )
        if total < 10 or not supported:
            raise SystemExit(f"Plasmidsaurus fixture PAC {key} lacks discovery support: {counts}")


def read_offsets(count: int) -> list[int]:
    """Deterministic quantiles of the triangular offset distribution."""
    span = OFFSET_HIGH - OFFSET_LOW
    left = OFFSET_MODE - OFFSET_LOW
    offsets = []
    for rank in range(count):
        quantile = (rank + 0.5) / count
        if quantile < left / span:
            value = OFFSET_LOW + math.sqrt(quantile * span * left)
        else:
            value = OFFSET_HIGH - math.sqrt((1 - quantile) * span * (OFFSET_HIGH - OFFSET_MODE))
        offsets.append(int(round(value)))
    return offsets


def contig_length() -> int:
    return FIRST_ORIGIN + GENE_SPACING * len(gene_layout()) + FIRST_ORIGIN


def write_reference() -> None:
    sequence = list("C" * contig_length())
    gtf = []
    for gene in gene_layout():
        origin, strand, gene_id = gene["origin"], gene["strand"], gene["gene_id"]
        for coordinate in gene["pacs"]:
            # Canonical AATAAA 20 nt upstream of each PAC (TTTATT on minus).
            if strand == "+":
                sequence[coordinate - 20 : coordinate - 14] = list("AATAAA")
            else:
                sequence[coordinate + 14 : coordinate + 20] = list("TTTATT")
        start, end = origin + 1, origin + GENE_LENGTH
        gtf.append(
            f"{CONTIG}\tfixture\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
            f'gene_id "{gene_id}"; gene_name "{gene_id.upper()}";'
        )
        # One transcript per PAC: each ends at its PAC, so a multi-PAC gene has
        # several annotated ends and stays out of calibration.
        for rank, coordinate in enumerate(gene["pacs"]):
            if strand == "+":
                exon_start, exon_end = start, coordinate
            else:
                exon_start, exon_end = coordinate + 1, end
            gtf.append(
                f"{CONTIG}\tfixture\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t"
                f'gene_id "{gene_id}"; transcript_id "{gene_id}_tx{rank + 1}";'
            )
    (ROOT / "genome.fa").write_text(f">{CONTIG}\n" + "".join(sequence) + "\n")
    (ROOT / "genes.gtf").write_text("\n".join(gtf) + "\n")


def aligned_read(name: str, strand: str, end_coordinate: int) -> pysam.AlignedSegment:
    """A forward-stranded read whose aligned 3' end is ``end_coordinate``."""
    read = pysam.AlignedSegment()
    read.query_name = name
    read.flag = 16 if strand == "-" else 0
    read.reference_id = 0
    # Interbase: the plus-strand 3' end is the reference end, the minus-strand
    # 3' end the reference start.
    read.reference_start = end_coordinate if strand == "-" else end_coordinate - READ_LENGTH
    read.mapping_quality = 60
    read.cigar = [(0, READ_LENGTH)]
    read.query_sequence = ("CGT" * (READ_LENGTH // 3 + 1))[:READ_LENGTH]
    read.query_qualities = pysam.qualitystring_to_array("I" * READ_LENGTH)
    read.set_tag("NH", 1)
    return read


def write_alignment(sample_id: str, sites: dict[tuple[str, int], dict[str, int]]) -> str:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": contig_length()}],
    }
    reads = []
    for (strand, coordinate), counts in sites.items():
        for serial, offset in enumerate(read_offsets(counts[sample_id])):
            end = coordinate - offset if strand == "+" else coordinate + offset
            name = f"{sample_id}_{strand}{coordinate}_{serial:04d}"
            reads.append(aligned_read(name, strand, end))
    reads.sort(key=lambda read: (read.reference_start, read.query_name))
    path = ROOT / f"{sample_id}.bam"
    with pysam.AlignmentFile(path, "wb", header=header) as handle:
        for read in reads:
            handle.write(read)
    pysam.index(str(path))
    return path.name


def write_expected(sites: dict[tuple[str, int], dict[str, int]]) -> None:
    """Designed PACs: gene, design role, rank from distal, and reads per sample."""
    site_details = {
        (gene["strand"], coordinate): (gene["gene_id"], design(gene), rank)
        for gene in gene_layout()
        for rank, coordinate in enumerate(gene["pacs"])
    }
    header = ["gene_id", "design", "strand", "coordinate", "rank_from_distal", *SAMPLES]
    lines = ["\t".join(header)]
    for strand, coordinate in sorted(sites, key=lambda item: (item[1], item[0])):
        gene_id, role, rank = site_details[(strand, coordinate)]
        values = [str(sites[(strand, coordinate)][sample_id]) for sample_id in SAMPLES]
        lines.append("\t".join([gene_id, role, strand, str(coordinate), str(rank), *values]))
    (ROOT / "expected_pacs.tsv").write_text("\n".join(lines) + "\n")


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    for pattern in ("*.bam", "*.bai", "genome.fa.fai"):
        for path in ROOT.glob(pattern):
            path.unlink()
    # RandomState keeps a frozen stream across NumPy versions. It is separate
    # from the exact-boundary fixture's stream, which it must not disturb.
    random = np.random.RandomState(20260928)
    sites = site_counts(random)
    check_support(sites)
    write_reference()
    paths = {sample_id: write_alignment(sample_id, sites) for sample_id in SAMPLES}
    write_expected(sites)
    with (ROOT / "samples.tsv").open("w") as handle:
        handle.write("sample_id\talignment\tcondition\tcontrol\treplicate\tlibrary_profile\n")
        for index, (sample_id, condition, control) in enumerate(ROWS, start=1):
            handle.write(
                f"{sample_id}\t{paths[sample_id]}\t{condition}\t{control}\t{index}\t"
                "plasmidsaurus_3prime\n"
            )
    print(f"Wrote the Plasmidsaurus-like fixture to {ROOT}")


if __name__ == "__main__":
    sys.exit(main())
