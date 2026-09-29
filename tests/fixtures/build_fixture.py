#!/usr/bin/env python3
"""Generate the BAM/CRAM integration fixture used by the test profile.

chr1 carries the original hand-made genes: gene_plus with PACs at 300, 350,
and 400 (350 is gained in TreatmentA) and single-PAC gene_minus at 800. chr2
adds 32 background genes so that DRIMSeq's cross-gene precision moderation
behaves as it does on real data. Designed chr2 genes have known effects, and
the remaining null genes are seeded Dirichlet-multinomial draws. chr3 holds
two multi-exon genes for the APA patterns: ipa01 gains an intronic PAC in
TreatmentA, and ale01 switches between two alternative last exons.
"""

from __future__ import annotations

import sys
from pathlib import Path

import build_plasmidsaurus_fixture
import numpy as np
import pysam

ROOT = Path(__file__).resolve().parent
CONTIGS = {"chr1": 2000, "chr2": 34000, "chr3": 8000}
READ_LENGTH = 30

# (sample_id, condition, control), in sample-sheet order.
ROWS = [
    ("DMSO_1", "DMSO", ""),
    ("DMSO_2", "DMSO", ""),
    ("TRA_1", "TreatmentA", "DMSO"),
    ("TRA_2", "TreatmentA", "DMSO"),
    ("VEH_1", "Vehicle", ""),
    ("VEH_2", "Vehicle", ""),
    ("TRB_1", "TreatmentB", "Vehicle"),
    ("TRB_2", "TreatmentB", "Vehicle"),
    ("RES_1", "Rescue", "TreatmentA"),
    ("RES_2", "Rescue", "TreatmentA"),
]
SAMPLES = [sample_id for sample_id, _, _ in ROWS]
CONDITION = {sample_id: condition for sample_id, condition, _ in ROWS}

# chr1 counts per sample for (+300, +350, +400, -800).
CHR1_COUNTS = {
    "DMSO_1": (40, 0, 80, 80),
    "DMSO_2": (30, 0, 90, 70),
    "TRA_1": (40, 70, 40, 80),
    "TRA_2": (30, 80, 50, 90),
    "VEH_1": (40, 0, 70, 80),
    "VEH_2": (40, 0, 80, 80),
    "TRB_1": (80, 0, 30, 80),
    "TRB_2": (80, 0, 40, 80),
    "RES_1": (20, 80, 70, 80),
    "RES_2": (20, 90, 80, 80),
}

BACKGROUND_GENES = 32
NULL_BASES = [(0.5, 0.3, 0.2), (0.4, 0.35, 0.25), (0.6, 0.25, 0.15)]
TWO_PAC_BASES = [(0.6, 0.4), (0.7, 0.3)]
NULL_DEPTHS = [60, 150, 300]


def gene_origin(index: int) -> int:
    return 1000 * index


def gene_strand(index: int) -> str:
    return "+" if index % 2 else "-"


def pac_coordinates(index: int, pac_count: int) -> list[int]:
    """Interbase PAC coordinates in genomic order.

    Plus-strand genes end at origin + 400 and minus-strand genes at origin,
    so the distal PAC sits on the annotated transcript end.
    """
    origin = gene_origin(index)
    if gene_strand(index) == "+":
        coordinates = [origin + 300, origin + 350, origin + 400]
        return coordinates[-pac_count:]
    return [origin, origin + 50, origin + 100][:pac_count]


def designed_proportions(index: int) -> dict[str, tuple[float, ...]] | None:
    """Proportions by condition, in genomic PAC order, for designed genes."""
    everywhere = dict.fromkeys(CONDITION.values())
    if index == 1:  # middle PAC lost in TreatmentA, so gained in Rescue vs TreatmentA
        values = {condition: (0.35, 0.30, 0.35) for condition in everywhere}
        values["TreatmentA"] = (0.5, 0.0, 0.5)
        return values
    if index == 2:  # middle PAC gained only in TreatmentB
        values = {condition: (0.6, 0.0, 0.4) for condition in everywhere}
        values["TreatmentB"] = (0.4, 0.3, 0.3)
        return values
    if 3 <= index <= 6:  # usage shift in TreatmentA only
        values = {condition: (0.5, 0.3, 0.2) for condition in everywhere}
        values["TreatmentA"] = (0.25, 0.3, 0.45)
        return values
    if index == 9:  # gained in TreatmentA at an internal-priming-like site
        values = {condition: (0.6, 0.0, 0.4) for condition in everywhere}
        values["TreatmentA"] = (0.4, 0.3, 0.3)
        return values
    return None


def exact_counts(index: int) -> dict[str, tuple[int, ...]] | None:
    """Literal counts for the exact null genes."""
    if index == 7:  # TreatmentB is exactly three times Vehicle
        vehicle = {"VEH_1": (60, 30, 30), "VEH_2": (50, 35, 35)}
        counts = {}
        for sample_id in SAMPLES:
            replicate = 1 if sample_id.endswith("_1") else 2
            base = vehicle[f"VEH_{replicate}"]
            factor = 3 if sample_id.startswith("TRB") else 1
            counts[sample_id] = tuple(factor * value for value in base)
        return counts
    if index == 8:  # exactly 3:2 in every sample
        return {sample_id: (120, 80) for sample_id in SAMPLES}
    return None


def pac_count_for(index: int) -> int:
    if index == 8 or (index > 9 and index % 5 == 0):
        return 2
    return 3


def background_counts(random: np.random.RandomState) -> dict[tuple[str, str, int], dict[str, int]]:
    """Counts keyed by (contig, strand, coordinate) for the chr2 genes."""
    sites: dict[tuple[str, str, int], dict[str, int]] = {}
    for index in range(1, BACKGROUND_GENES + 1):
        strand = gene_strand(index)
        pac_count = pac_count_for(index)
        coordinates = pac_coordinates(index, pac_count)
        literal = exact_counts(index)
        designed = designed_proportions(index)
        per_sample: dict[str, tuple[int, ...]] = {}
        for sample_number, sample_id in enumerate(SAMPLES):
            if literal is not None:
                per_sample[sample_id] = literal[sample_id]
            elif designed is not None:
                depth = 280 + 20 * (sample_number % 3)
                proportions = np.asarray(designed[CONDITION[sample_id]], dtype=float)
                draw = random.multinomial(depth, proportions)
                per_sample[sample_id] = tuple(int(value) for value in draw)
            else:
                bases = NULL_BASES if pac_count == 3 else TWO_PAC_BASES
                base = np.asarray(bases[index % len(bases)], dtype=float)
                mean = NULL_DEPTHS[index % len(NULL_DEPTHS)]
                total = int(random.negative_binomial(20, 20 / (20 + mean)))
                proportions = random.dirichlet(base * 100)
                draw = random.multinomial(total, proportions)
                per_sample[sample_id] = tuple(int(value) for value in draw)
        for position, coordinate in enumerate(coordinates):
            sites[("chr2", strand, coordinate)] = {
                sample_id: per_sample[sample_id][position] for sample_id in SAMPLES
            }
    return sites


# chr3 genes, plus strand: (gene_id, transcripts as 1-based exons, interbase
# PAC coordinates, proportions in controls and in TreatmentA). ipa01's
# intronic PAC lies 1,250 nt from its annotated end, outside calibration's
# 1,000-nt window, and ale01 has two ends, so neither changes the kernel.
CHR3_GENES = [
    (
        "ipa01",
        {"ipa01_tx": [(1001, 1100), (2201, 2400)]},
        (1150, 2400),
        (0.1, 0.9),
        (0.6, 0.4),
    ),
    (
        "ale01",
        {
            "ale01_short": [(4001, 4100), (4801, 5000)],
            "ale01_long": [(4001, 4100), (5401, 5600)],
        },
        (5000, 5600),
        (0.8, 0.2),
        (0.2, 0.8),
    ),
]


def chr3_counts(random: np.random.RandomState) -> dict[tuple[str, str, int], dict[str, int]]:
    """Counts for the chr3 genes, drawn after chr2's so its draws stay the same."""
    sites: dict[tuple[str, str, int], dict[str, int]] = {}
    for _, _, coordinates, baseline, treated in CHR3_GENES:
        per_sample = {}
        for sample_number, sample_id in enumerate(SAMPLES):
            depth = 280 + 20 * (sample_number % 3)
            proportions = treated if CONDITION[sample_id] == "TreatmentA" else baseline
            draw = random.multinomial(depth, np.asarray(proportions, dtype=float))
            per_sample[sample_id] = tuple(int(value) for value in draw)
        for position, coordinate in enumerate(coordinates):
            sites[("chr3", "+", coordinate)] = {
                sample_id: per_sample[sample_id][position] for sample_id in SAMPLES
            }
    return sites


def chr1_counts() -> dict[tuple[str, str, int], dict[str, int]]:
    keys = [("chr1", "+", 300), ("chr1", "+", 350), ("chr1", "+", 400), ("chr1", "-", 800)]
    return {
        key: {sample_id: CHR1_COUNTS[sample_id][position] for sample_id in SAMPLES}
        for position, key in enumerate(keys)
    }


def check_support(sites: dict[tuple[str, str, int], dict[str, int]]) -> None:
    """Every PAC must pass default discovery: 10 reads in total and 2 or more
    reads in both replicates of at least one condition."""
    for key, counts in sites.items():
        total = sum(counts.values())
        supported = any(
            all(counts[sample_id] >= 2 for sample_id in SAMPLES if CONDITION[sample_id] == group)
            for group in set(CONDITION.values())
        )
        if total < 10 or not supported:
            raise SystemExit(f"Fixture PAC {key} lacks discovery support: {counts}")


def write_reference() -> None:
    chr1 = list("C" * CONTIGS["chr1"])
    chr1[360:366] = list("AATAAA")
    chr1[830:836] = list("TTTATT")
    chr2 = list("C" * CONTIGS["chr2"])
    gtf = [
        'chr1\tfixture\tgene\t101\t400\t.\t+\t.\tgene_id "gene_plus"; gene_name "GenePlus";',
        'chr1\tfixture\texon\t101\t400\t.\t+\t.\tgene_id "gene_plus"; transcript_id "plus_tx";',
        'chr1\tfixture\tgene\t801\t1100\t.\t-\t.\tgene_id "gene_minus"; gene_name "GeneMinus";',
        'chr1\tfixture\texon\t801\t1100\t.\t-\t.\tgene_id "gene_minus"; transcript_id "minus_tx";',
    ]
    for index in range(1, BACKGROUND_GENES + 1):
        origin = gene_origin(index)
        strand = gene_strand(index)
        coordinates = pac_coordinates(index, pac_count_for(index))
        # Canonical AATAAA 20 nt upstream of the distal PAC and, for every third
        # gene, ATTAAA upstream of the proximal PAC (reverse complements on -).
        if strand == "+":
            distal, proximal = coordinates[-1], coordinates[0]
            chr2[distal - 20 : distal - 14] = list("AATAAA")
            if index % 3 == 0 and len(coordinates) == 3:
                chr2[proximal - 20 : proximal - 14] = list("ATTAAA")
        else:
            distal, proximal = coordinates[0], coordinates[-1]
            chr2[distal + 14 : distal + 20] = list("TTTATT")
            if index % 3 == 0 and len(coordinates) == 3:
                chr2[proximal + 14 : proximal + 20] = list("TTTAAT")
        if index == 9:
            # An A-rich stretch just downstream of the middle PAC mimics
            # internal priming.
            chr2[origin + 350 : origin + 360] = list("A" * 10)
        gene_id = f"bg{index:02d}"
        fields = f"chr2\tfixture\t{{}}\t{origin + 1}\t{origin + 400}\t.\t{strand}\t."
        gtf.append(fields.format("gene") + f'\tgene_id "{gene_id}"; gene_name "{gene_id.upper()}";')
        gtf.append(fields.format("exon") + f'\tgene_id "{gene_id}"; transcript_id "{gene_id}_tx";')
    chr3 = list("C" * CONTIGS["chr3"])
    for gene_id, transcripts, coordinates, _, _ in CHR3_GENES:
        for coordinate in coordinates:
            chr3[coordinate - 20 : coordinate - 14] = list("AATAAA")
        exons = [exon for parts in transcripts.values() for exon in parts]
        start = min(first for first, _ in exons)
        end = max(last for _, last in exons)
        gtf.append(
            f"chr3\tfixture\tgene\t{start}\t{end}\t.\t+\t.\t"
            f'gene_id "{gene_id}"; gene_name "{gene_id.upper()}";'
        )
        for transcript_id, parts in transcripts.items():
            for first, last in parts:
                gtf.append(
                    f"chr3\tfixture\texon\t{first}\t{last}\t.\t+\t.\t"
                    f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
                )
    fasta = ROOT / "genome.fa"
    fasta.write_text(
        ">chr1\n" + "".join(chr1) + "\n>chr2\n" + "".join(chr2)
        + "\n>chr3\n" + "".join(chr3) + "\n"
    )
    pysam.faidx(str(fasta))
    (ROOT / "genes.gtf").write_text("\n".join(gtf) + "\n")


def aligned_read(
    name: str, contig: str, coordinate: int, strand: str, serial: int
) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment()
    read.query_name = f"{name}_{serial:05d}"
    read.flag = 16 if strand == "-" else 0
    read.reference_id = list(CONTIGS).index(contig)
    # The aligned 3' end sits on the PAC coordinate: the reference end on the
    # plus strand, the reference start on the minus strand.
    read.reference_start = coordinate if strand == "-" else coordinate - READ_LENGTH
    read.mapping_quality = 60
    read.cigar = [(0, READ_LENGTH)]
    read.query_sequence = "CGT" * (READ_LENGTH // 3)
    read.query_qualities = pysam.qualitystring_to_array("I" * READ_LENGTH)
    read.set_tag("NH", 1)
    return read


def write_alignment(
    sample_id: str,
    sites: dict[tuple[str, str, int], dict[str, int]],
    cram: bool = False,
    unsorted: bool = False,
) -> str:
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted" if unsorted else "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in CONTIGS.items()],
    }
    reads = []
    serial = 0
    for (contig, strand, coordinate), counts in sites.items():
        for _ in range(counts[sample_id]):
            serial += 1
            reads.append(aligned_read(sample_id, contig, coordinate, strand, serial))
    reads.sort(
        key=lambda read: (read.reference_id, read.reference_start, read.query_name),
        reverse=unsorted,
    )
    suffix = ".cram" if cram else ".bam"
    path = ROOT / f"{sample_id}{suffix}"
    mode = "wc" if cram else "wb"
    kwargs = {"reference_filename": str(ROOT / "genome.fa")} if cram else {}
    with pysam.AlignmentFile(path, mode, header=header, **kwargs) as handle:
        for read in reads:
            handle.write(read)
    if not unsorted:
        pysam.index(str(path))
    return path.name


def write_expected_counts(sites: dict[tuple[str, str, int], dict[str, int]]) -> None:
    gene_by_site: dict[tuple[str, str, int], str] = {
        ("chr1", "+", 300): "gene_plus",
        ("chr1", "+", 350): "gene_plus",
        ("chr1", "+", 400): "gene_plus",
        ("chr1", "-", 800): "gene_minus",
    }
    for index in range(1, BACKGROUND_GENES + 1):
        for coordinate in pac_coordinates(index, pac_count_for(index)):
            gene_by_site[("chr2", gene_strand(index), coordinate)] = f"bg{index:02d}"
    for gene_id, _, coordinates, _, _ in CHR3_GENES:
        for coordinate in coordinates:
            gene_by_site[("chr3", "+", coordinate)] = gene_id
    lines = ["\t".join(["gene_id", "pac_id", *SAMPLES])]
    for key in sorted(sites, key=lambda item: (item[0], item[2], item[1])):
        contig, strand, coordinate = key
        pac_id = f"PACv1.synthetic.{contig}.{strand}.{coordinate}"
        values = [str(sites[key][sample_id]) for sample_id in SAMPLES]
        lines.append("\t".join([gene_by_site[key], pac_id, *values]))
    (ROOT / "expected_pac_counts.tsv").write_text("\n".join(lines) + "\n")


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    for pattern in ("*.bam", "*.bai", "*.cram", "*.crai", "genome.fa.fai"):
        for path in ROOT.glob(pattern):
            path.unlink()
    # RandomState keeps a frozen stream across NumPy versions.
    random = np.random.RandomState(20260927)
    sites = {**chr1_counts(), **background_counts(random), **chr3_counts(random)}
    check_support(sites)
    write_reference()
    layout = {"DMSO_1": (False, True), "VEH_2": (True, False)}
    paths = {
        sample_id: write_alignment(sample_id, sites, *layout.get(sample_id, (False, False)))
        for sample_id in SAMPLES
    }
    write_expected_counts(sites)
    with (ROOT / "samples.tsv").open("w") as handle:
        handle.write("sample_id\talignment\tcondition\tcontrol\treplicate\n")
        for index, (sample_id, condition, control) in enumerate(ROWS, start=1):
            handle.write(f"{sample_id}\t{paths[sample_id]}\t{condition}\t{control}\t{index}\n")
    with (ROOT / "samples_incompatible.tsv").open("w") as handle:
        handle.write("sample_id\talignment\tcondition\tcontrol\tlibrary_profile\n")
        for index, (sample_id, condition, control) in enumerate(ROWS[:4]):
            profile = "exact_boundary" if index < 2 else "plasmidsaurus_3prime"
            handle.write(f"{sample_id}\t{paths[sample_id]}\t{condition}\t{control}\t{profile}\n")
    (ROOT / "genome.fa.fai").unlink(missing_ok=True)
    print(f"Wrote PACusage integration fixture to {ROOT}")
    # The Plasmidsaurus-like run needs proximal-tag reads, which these
    # exact-boundary alignments cannot provide; it has its own fixture.
    build_plasmidsaurus_fixture.main()


if __name__ == "__main__":
    sys.exit(main())
