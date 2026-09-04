"""Gene, transcript, sequence, PAS motif, and internal-priming annotation."""

from __future__ import annotations

import gzip
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

import pysam

from .errors import PacusageError
from .evidence import reverse_complement
from .models import PacCandidate
from .reference import GenomicFeature

DEFAULT_MOTIFS = (
    ("AATAAA", "canonical", 1),
    ("ATTAAA", "common_variant", 2),
    ("TATAAA", "other_variant", 3),
    ("AGTAAA", "other_variant", 4),
    ("AAGAAA", "other_variant", 5),
    ("AATACA", "other_variant", 6),
    ("CATAAA", "other_variant", 7),
    ("GATAAA", "other_variant", 8),
    ("AATATA", "other_variant", 9),
    ("AATAGA", "other_variant", 10),
)


def load_known_pacs(path: str | Path | None) -> set[tuple[str, str, int]]:
    if not path:
        return set()
    opener = gzip.open if str(path).endswith(".gz") else open
    sites: set[tuple[str, str, int]] = set()
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            strand = fields[5] if len(fields) > 5 and fields[5] in {"+", "-"} else "+"
            coordinate = int(fields[2]) if strand == "+" else int(fields[1])
            sites.add((fields[0], strand, coordinate))
    return sites


def annotate_candidates(
    candidates: Iterable[PacCandidate],
    features: list[GenomicFeature],
    fasta_path: str | Path,
    assembly: str,
    endpoint_model: str,
    evidence_source: str,
    maximum_downstream_distance: int,
    scan_upstream_far: int,
    scan_upstream_near: int,
    core_upstream_far: int,
    core_upstream_near: int,
    internal_priming_window: int,
    internal_priming_max_a_run: int,
    internal_priming_max_a_fraction: float,
    known_sites: set[tuple[str, str, int]] | None = None,
    known_match_radius: int = 12,
    motif_catalog: tuple[tuple[str, str, int], ...] = DEFAULT_MOTIFS,
) -> list[dict[str, object]]:
    known_sites = known_sites or set()
    feature_index = FeatureIndex(features)
    fasta = pysam.FastaFile(str(fasta_path))
    rows: list[dict[str, object]] = []
    try:
        ordered_candidates = sorted(
            candidates, key=lambda item: (item.contig, item.coordinate, item.strand)
        )
        for candidate in ordered_candidates:
            assignments = feature_index.assign(
                candidate.contig,
                candidate.strand,
                candidate.coordinate,
                maximum_downstream_distance,
            )
            assignment_class = assignments[0][0] if assignments else "intergenic"
            gene_ids = sorted({item[1].gene_id for item in assignments})
            gene_names = sorted({item[1].gene_name for item in assignments if item[1].gene_name})
            ambiguous = len(gene_ids) > 1
            upstream = oriented_interval_sequence(
                fasta,
                candidate.contig,
                candidate.strand,
                candidate.coordinate,
                scan_upstream_far,
                scan_upstream_near,
                upstream=True,
            )
            downstream = oriented_interval_sequence(
                fasta,
                candidate.contig,
                candidate.strand,
                candidate.coordinate,
                0,
                internal_priming_window,
                upstream=False,
            )
            motif_matches = find_motifs(
                upstream,
                scan_upstream_far,
                core_upstream_far,
                core_upstream_near,
                motif_catalog,
            )
            primary = choose_primary_motif(motif_matches)
            longest_a = longest_run(downstream, "A")
            a_fraction = downstream.count("A") / len(downstream) if downstream else 0.0
            internal_priming = (
                longest_a > internal_priming_max_a_run
                or a_fraction > internal_priming_max_a_fraction
            )
            known_matches = sorted(
                coordinate
                for contig, strand, coordinate in known_sites
                if contig == candidate.contig
                and strand == candidate.strand
                and abs(coordinate - candidate.coordinate) <= known_match_radius
            )
            pac_id = (
                f"PACv1.{assembly}.{candidate.contig}.{candidate.strand}.{candidate.coordinate}"
            )
            rows.append(
                {
                    "pac_id": pac_id,
                    "contig": candidate.contig,
                    "coordinate": candidate.coordinate,
                    "strand": candidate.strand,
                    "endpoint_model": endpoint_model,
                    "coordinate_precision": (
                        "exact" if endpoint_model == "exact_boundary" else "estimated"
                    ),
                    "resolved_evidence_source": evidence_source,
                    "primary_evidence_type": (
                        "exact_boundary" if endpoint_model == "exact_boundary" else "endpoint_peak"
                    ),
                    "resolution_group": pac_id,
                    "merged_candidate_coordinates": ",".join(
                        map(str, candidate.member_coordinates)
                    ),
                    "total_count": candidate.total_count,
                    "supporting_samples": candidate.supporting_samples,
                    "fraction_within_2nt": round(candidate.fraction_within_2nt, 6),
                    "width_90": candidate.width_90,
                    "local_strand_enrichment": round(candidate.local_enrichment, 6),
                    "poly_a_clip_fraction": round(candidate.poly_a_clip_fraction, 6),
                    "coordinate_interval_low": candidate.coordinate_interval_low,
                    "coordinate_interval_high": candidate.coordinate_interval_high,
                    "coordinate_bootstrap_successes": (candidate.coordinate_bootstrap_successes),
                    "calibration_profile": "PACusage-0.1",
                    "candidate_status": candidate.status,
                    "known_pac": bool(known_matches),
                    "known_pac_coordinates": ",".join(map(str, known_matches)),
                    "known_rescue_only": candidate.status == "known_rescue_only",
                    "gene_id": ",".join(gene_ids),
                    "gene_name": ",".join(gene_names),
                    "assignment_class": assignment_class,
                    "ambiguous_gene_assignment": ambiguous,
                    "upstream_sequence": upstream,
                    "downstream_sequence": downstream,
                    "all_pas_motifs": ";".join(
                        f"{motif}:{position}:{motif_class}:{int(in_core)}"
                        for motif, position, motif_class, _, in_core in motif_matches
                    ),
                    "primary_pas_motif": primary[0],
                    "primary_pas_motif_rna": (
                        primary[0].replace("T", "U") if primary[0] else "none"
                    ),
                    "primary_motif_class": primary[2],
                    "primary_motif_position": primary[1],
                    "primary_motif_in_core": primary[4],
                    "downstream_a_fraction": round(a_fraction, 6),
                    "downstream_longest_a_run": longest_a,
                    "internal_priming_flag": internal_priming,
                    "confidence": confidence_tier(
                        candidate, ambiguous, internal_priming, bool(known_matches)
                    ),
                }
            )
    finally:
        fasta.close()
    add_proximal_distal_ranks(rows)
    return rows


class FeatureIndex:
    def __init__(self, features: Iterable[GenomicFeature], bin_size: int = 10000):
        self.bin_size = bin_size
        self.genes: dict[tuple[str, str, str], GenomicFeature] = {}
        self.exons: dict[tuple[str, str], list[GenomicFeature]] = defaultdict(list)
        self.transcript_exons: dict[tuple[str, str, str], list[GenomicFeature]] = defaultdict(list)
        for feature in features:
            if feature.feature_type == "gene":
                self.genes[(feature.contig, feature.strand, feature.gene_id)] = feature
            if feature.feature_type == "exon":
                self.exons[(feature.contig, feature.strand)].append(feature)
                key = (feature.contig, feature.strand, feature.transcript_id or feature.gene_id)
                self.transcript_exons[key].append(feature)
        if not self.genes:
            gene_parts: dict[tuple[str, str, str], list[GenomicFeature]] = defaultdict(list)
            for exon_list in self.exons.values():
                for exon in exon_list:
                    gene_parts[(exon.contig, exon.strand, exon.gene_id)].append(exon)
            for key, values in gene_parts.items():
                first = values[0]
                self.genes[key] = GenomicFeature(
                    contig=first.contig,
                    start=min(item.start for item in values),
                    end=max(item.end for item in values),
                    strand=first.strand,
                    feature_type="gene",
                    gene_id=first.gene_id,
                    gene_name=first.gene_name,
                )
        self.terminal_exons: set[tuple[str, str, int, int, str]] = set()
        for exons in self.transcript_exons.values():
            terminal = (
                max(exons, key=lambda item: item.end)
                if exons[0].strand == "+"
                else min(exons, key=lambda item: item.start)
            )
            self.terminal_exons.add(
                (terminal.contig, terminal.strand, terminal.start, terminal.end, terminal.gene_id)
            )

    def assign(
        self, contig: str, strand: str, coordinate: int, maximum_downstream: int
    ) -> list[tuple[str, GenomicFeature]]:
        exons = [
            feature
            for feature in self.exons.get((contig, strand), [])
            if feature.start <= coordinate <= feature.end
        ]
        terminal = [
            feature
            for feature in exons
            if (feature.contig, feature.strand, feature.start, feature.end, feature.gene_id)
            in self.terminal_exons
        ]
        if terminal:
            return [("terminal_exon", item) for item in terminal]
        if exons:
            return [("other_exon", item) for item in exons]
        intronic = [
            feature
            for key, feature in self.genes.items()
            if key[0] == contig and key[1] == strand and feature.start <= coordinate <= feature.end
        ]
        if intronic:
            return [("intronic", item) for item in intronic]
        downstream = []
        same_strand_genes = [
            feature for key, feature in self.genes.items() if key[0] == contig and key[1] == strand
        ]
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


def oriented_interval_sequence(
    fasta: pysam.FastaFile,
    contig: str,
    strand: str,
    coordinate: int,
    far: int,
    near: int,
    upstream: bool,
) -> str:
    contig_length = fasta.get_reference_length(contig)
    if upstream:
        if strand == "+":
            start, end = coordinate - far, coordinate - near
        else:
            start, end = coordinate + near, coordinate + far
    else:
        if strand == "+":
            start, end = coordinate + far, coordinate + near
        else:
            start, end = coordinate - near, coordinate - far
    start, end = max(0, min(start, end)), min(contig_length, max(start, end))
    sequence = fasta.fetch(contig, start, end).upper()
    return sequence if strand == "+" else reverse_complement(sequence)


def find_motifs(
    upstream_sequence: str,
    upstream_far: int,
    core_far: int,
    core_near: int,
    motif_catalog: tuple[tuple[str, str, int], ...] = DEFAULT_MOTIFS,
) -> list[tuple[str, int, str, int, bool]]:
    matches: list[tuple[str, int, str, int, bool]] = []
    for motif, motif_class, priority in motif_catalog:
        start = 0
        while True:
            index = upstream_sequence.find(motif, start)
            if index < 0:
                break
            distance = upstream_far - index
            in_core = core_near <= distance <= core_far
            matches.append((motif, distance, motif_class, priority, in_core))
            start = index + 1
    return sorted(matches, key=lambda item: (item[3], not item[4], abs(item[1] - 20), item[1]))


def load_motif_catalog(
    path: str | Path | None,
) -> tuple[tuple[str, str, int], ...]:
    if not path:
        return DEFAULT_MOTIFS
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line_number == 1 and fields[0].lower() == "motif":
                continue
            if len(fields) < 2:
                raise PacusageError(
                    f"Motif catalog {path}, line {line_number}: "
                    "expected motif, class, and optional priority."
                )
            motif = fields[0].upper()
            if len(motif) != 6 or set(motif).difference("ACGT"):
                raise PacusageError(
                    f"Motif catalog {path}, line {line_number}: {motif!r} "
                    "must be a six-base DNA sequence."
                )
            priority = int(fields[2]) if len(fields) > 2 else len(rows) + 1
            rows.append((motif, fields[1], priority))
    return tuple(sorted(rows, key=lambda row: (row[2], row[0])))


def choose_primary_motif(
    matches: list[tuple[str, int, str, int, bool]],
) -> tuple[str, int | str, str, int, bool]:
    return matches[0] if matches else ("", "", "no_recognized_motif", 999, False)


def longest_run(sequence: str, base: str) -> int:
    longest = current = 0
    for value in sequence:
        current = current + 1 if value == base else 0
        longest = max(longest, current)
    return longest


def confidence_tier(
    candidate: PacCandidate,
    ambiguous_assignment: bool,
    internal_priming: bool,
    known_match: bool,
) -> str:
    if ambiguous_assignment or internal_priming:
        return "low"
    if candidate.status == "known_rescue_only":
        return "moderate"
    if known_match and candidate.supporting_samples >= 2:
        return "high"
    return "moderate"


def add_proximal_distal_ranks(rows: list[dict[str, object]]) -> None:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        if row["gene_id"] and not row["ambiguous_gene_assignment"]:
            grouped[(str(row["gene_id"]), str(row["strand"]))].append(row)
    for (_, strand), values in grouped.items():
        ordered = sorted(
            values,
            key=lambda row: int(row["coordinate"]),
            reverse=strand == "-",
        )
        for rank, row in enumerate(ordered, start=1):
            row["proximal_distal_rank"] = rank
            row["proximal_distal_label"] = (
                "proximal" if rank == 1 else "distal" if rank == len(ordered) else "middle"
            )
