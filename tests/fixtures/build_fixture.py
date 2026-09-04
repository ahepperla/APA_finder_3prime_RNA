#!/usr/bin/env python3
"""Generate the tiny BAM/CRAM integration fixture used by the test profile."""

from __future__ import annotations

import sys
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parent


def write_reference() -> None:
    sequence = list("C" * 2000)
    sequence[360:366] = list("AATAAA")
    sequence[830:836] = list("TTTATT")
    fasta = ROOT / "genome.fa"
    fasta.write_text(">chr1\n" + "".join(sequence) + "\n")
    pysam.faidx(str(fasta))
    (ROOT / "genes.gtf").write_text(
        'chr1\tfixture\tgene\t101\t400\t.\t+\t.\tgene_id "gene_plus"; gene_name "GenePlus";\n'
        'chr1\tfixture\texon\t101\t400\t.\t+\t.\tgene_id "gene_plus"; transcript_id "plus_tx";\n'
        'chr1\tfixture\tgene\t801\t1100\t.\t-\t.\tgene_id "gene_minus"; gene_name "GeneMinus";\n'
        'chr1\tfixture\texon\t801\t1100\t.\t-\t.\tgene_id "gene_minus"; transcript_id "minus_tx";\n'
    )


def aligned_read(name: str, coordinate: int, strand: str, serial: int) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment()
    read.query_name = f"{name}_{serial:04d}"
    read.flag = 16 if strand == "-" else 0
    read.reference_id = 0
    read.reference_start = coordinate if strand == "-" else coordinate - 30
    read.mapping_quality = 60
    read.cigar = [(0, 30)]
    read.query_sequence = "CGT" * 10
    read.query_qualities = pysam.qualitystring_to_array("I" * 30)
    read.set_tag("NH", 1)
    return read


def write_alignment(
    sample_id: str,
    sites: dict[tuple[int, str], int],
    cram: bool = False,
    unsorted: bool = False,
) -> str:
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted" if unsorted else "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 2000}],
    }
    reads = []
    serial = 0
    for (coordinate, strand), count in sites.items():
        for _ in range(count):
            serial += 1
            reads.append(aligned_read(sample_id, coordinate, strand, serial))
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


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    for pattern in ("*.bam", "*.bai", "*.cram", "*.crai", "genome.fa.fai"):
        for path in ROOT.glob(pattern):
            path.unlink()
    write_reference()
    definitions = {
        "DMSO_1": ({(400, "+"): 8, (300, "+"): 4, (800, "-"): 8}, False, True),
        "DMSO_2": ({(400, "+"): 9, (300, "+"): 3, (800, "-"): 7}, False, False),
        "TRA_1": (
            {(400, "+"): 4, (300, "+"): 4, (350, "+"): 7, (800, "-"): 8},
            False,
            False,
        ),
        "TRA_2": (
            {(400, "+"): 5, (300, "+"): 3, (350, "+"): 8, (800, "-"): 9},
            False,
            False,
        ),
        "VEH_1": ({(400, "+"): 7, (300, "+"): 4, (800, "-"): 8}, False, False),
        "VEH_2": ({(400, "+"): 8, (300, "+"): 4, (800, "-"): 8}, True, False),
        "TRB_1": ({(400, "+"): 3, (300, "+"): 8, (800, "-"): 8}, False, False),
        "TRB_2": ({(400, "+"): 4, (300, "+"): 8, (800, "-"): 8}, False, False),
        "RES_1": (
            {(400, "+"): 7, (300, "+"): 2, (350, "+"): 8, (800, "-"): 8},
            False,
            False,
        ),
        "RES_2": (
            {(400, "+"): 8, (300, "+"): 2, (350, "+"): 9, (800, "-"): 8},
            False,
            False,
        ),
    }
    paths = {
        sample_id: write_alignment(sample_id, sites, cram, unsorted)
        for sample_id, (sites, cram, unsorted) in definitions.items()
    }
    rows = [
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
    with (ROOT / "samples.tsv").open("w") as handle:
        handle.write("sample_id\talignment\tcondition\tcontrol\treplicate\n")
        for index, (sample_id, condition, control) in enumerate(rows, start=1):
            handle.write(f"{sample_id}\t{paths[sample_id]}\t{condition}\t{control}\t{index}\n")
    with (ROOT / "samples_plasmidsaurus.tsv").open("w") as handle:
        handle.write("sample_id\talignment\tcondition\tcontrol\tlibrary_profile\n")
        for sample_id, condition, control in rows[:4]:
            handle.write(
                f"{sample_id}\t{paths[sample_id]}\t{condition}\t{control}\tplasmidsaurus_3prime\n"
            )
    with (ROOT / "samples_incompatible.tsv").open("w") as handle:
        handle.write("sample_id\talignment\tcondition\tcontrol\tlibrary_profile\n")
        for index, (sample_id, condition, control) in enumerate(rows[:4]):
            profile = "exact_boundary" if index < 2 else "plasmidsaurus_3prime"
            handle.write(f"{sample_id}\t{paths[sample_id]}\t{condition}\t{control}\t{profile}\n")
    (ROOT / "genome.fa.fai").unlink(missing_ok=True)
    print(f"Wrote PACusage integration fixture to {ROOT}")


if __name__ == "__main__":
    sys.exit(main())
