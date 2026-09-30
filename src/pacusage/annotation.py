"""Gene, transcript, sequence, PAS motif, and internal-priming annotation."""

from __future__ import annotations

import gzip
from bisect import bisect_left, bisect_right
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

import pysam

from .errors import PacusageError
from .evidence import reverse_complement
from .models import IDENTITY_COLUMNS, PacCandidate
from .reference import GenomicFeature, gene_name_map

ATLAS_COLUMNS = [
    *IDENTITY_COLUMNS,
    "coordinate",
    "coordinate_precision",
    "resolution_nt",
    "merged_candidate_coordinates",
    "gene_region",
    "last_exon_locus",
    "ambiguous_gene_assignment",
    "proximal_distal_rank",
    "proximal_distal_label",
    "total_count",
    "supporting_samples",
    "total_supporting_samples",
    "best_supporting_condition",
    "candidate_status",
    "confidence",
    "known_pac",
    "known_pac_coordinates",
    "known_rescue_only",
    "internal_priming_flag",
    "downstream_a_fraction",
    "downstream_longest_a_run",
    "primary_pas_motif",
    "primary_pas_motif_rna",
    "primary_motif_class",
    "primary_motif_position",
    "primary_motif_in_core",
    "all_pas_motifs",
    "upstream_sequence",
    "downstream_sequence",
    "fraction_within_2nt",
    "width_90",
    "local_strand_enrichment",
    "poly_a_clip_fraction",
    "endpoint_model",
    "resolved_evidence_source",
    "primary_evidence_type",
    "calibration_profile",
]

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


class KnownSiteIndex:
    """Efficient lookup of known PAC sites within a specified radius."""

    def __init__(self, sites: Iterable[tuple[str, str, int]] = ()) -> None:
        """Build index from (contig, strand, coordinate) tuples.

        Sites are deduplicated and sorted by coordinate for each contig/strand pair.
        """
        self._sites: dict[tuple[str, str], list[int]] = defaultdict(list)
        for contig, strand, coordinate in sites:
            self._sites[(contig, strand)].append(coordinate)
        for key in self._sites:
            self._sites[key] = sorted(set(self._sites[key]))

    def matches(
        self, contig: str, strand: str, coordinate: int, radius: int
    ) -> list[int]:
        """Return sorted known coordinates within [coordinate - radius, coordinate + radius]."""
        sites = self._sites.get((contig, strand), [])
        if not sites:
            return []
        lower_idx = bisect_left(sites, coordinate - radius)
        upper_idx = bisect_right(sites, coordinate + radius)
        return sites[lower_idx:upper_idx]

    def __bool__(self) -> bool:
        """Return True if any site exists."""
        return bool(self._sites)


def load_known_pacs(path: str | Path | None) -> KnownSiteIndex:
    """Load known PAC sites from a BED6 file.

    Returns an empty KnownSiteIndex if path is falsy.
    Uses field 2 (end) for + strand, field 1 (start) for - strand.
    """
    sites: list[tuple[str, str, int]] = []
    if path:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                strand = fields[5] if len(fields) > 5 and fields[5] in {"+", "-"} else "+"
                coordinate = int(fields[2]) if strand == "+" else int(fields[1])
                sites.append((fields[0], strand, coordinate))
    return KnownSiteIndex(sites)


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
    known_sites: KnownSiteIndex | None = None,
    known_match_radius: int = 12,
    motif_catalog: tuple[tuple[str, str, int], ...] = DEFAULT_MOTIFS,
) -> list[dict[str, object]]:
    known_sites = known_sites or KnownSiteIndex()
    feature_index = FeatureIndex(features)
    names = gene_name_map(features)
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
            gene_region = assignments[0][0] if assignments else "intergenic"
            gene_ids = sorted({item[1].gene_id for item in assignments})
            # Names line up with the IDs, one per gene.
            gene_names = [names.get(gene_id, gene_id) for gene_id in gene_ids]
            ambiguous = len(gene_ids) > 1
            last_exon = "" if ambiguous else feature_index.last_exon(assignments)
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
            known_matches = known_sites.matches(
                candidate.contig,
                candidate.strand,
                candidate.coordinate,
                known_match_radius,
            )
            pac_id = (
                f"PACv1.{assembly}.{candidate.contig}.{candidate.strand}.{candidate.coordinate}"
            )
            start, end = browser_interval(candidate, endpoint_model)
            rows.append(
                {
                    "pac_id": pac_id,
                    "chrom": candidate.contig,
                    "start": start,
                    "end": end,
                    "locus": f"{candidate.contig}:{start + 1}-{end}",
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
                    "merged_candidate_coordinates": ",".join(
                        map(str, candidate.member_coordinates)
                    ),
                    "total_count": candidate.total_count,
                    "supporting_samples": candidate.supporting_samples,
                    "total_supporting_samples": candidate.total_supporting_samples,
                    "best_supporting_condition": candidate.best_supporting_condition,
                    "fraction_within_2nt": round(candidate.fraction_within_2nt, 6),
                    "width_90": candidate.width_90,
                    "local_strand_enrichment": round(candidate.local_enrichment, 6),
                    "poly_a_clip_fraction": round(candidate.poly_a_clip_fraction, 6),
                    "resolution_nt": candidate.resolution_nt,
                    "calibration_profile": "PACusage-0.1",
                    "candidate_status": candidate.status,
                    "known_pac": bool(known_matches),
                    "known_pac_coordinates": ",".join(map(str, known_matches)),
                    "known_rescue_only": candidate.status == "known_rescue_only",
                    "gene_id": ",".join(gene_ids),
                    "gene_name": ",".join(gene_names),
                    "gene_region": gene_region,
                    "last_exon_locus": last_exon,
                    "ambiguous_gene_assignment": ambiguous,
                    # Ranked below for PACs of a single gene; blank otherwise.
                    "proximal_distal_rank": "",
                    "proximal_distal_label": "",
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


def browser_interval(candidate: PacCandidate, endpoint_model: str) -> tuple[int, int]:
    """The PAC's BED interval: the boundary base for an exact PAC, the
    resolution region for a proximal-tag PAC."""
    if endpoint_model == "proximal_tag":
        return int(candidate.region_start), int(candidate.region_end)
    return candidate.coordinate, candidate.coordinate + 1


class FeatureIndex:
    """Gene assignment for PAC coordinates (design section 6).

    Containment is inclusive at both ends: in interbase coordinates, a 3-prime
    boundary equals its exon's ``end``. Exons and genes are binned so a query
    reads one bin; downstream genes are found by bisection.
    """

    def __init__(self, features: Iterable[GenomicFeature], bin_size: int = 10000):
        self.bin_size = bin_size
        genes: dict[tuple[str, str, str], GenomicFeature] = {}
        exons: list[GenomicFeature] = []
        for feature in features:
            if feature.feature_type == "gene":
                genes[(feature.contig, feature.strand, feature.gene_id)] = feature
            if feature.feature_type == "exon":
                exons.append(feature)
        if not genes:
            # Without gene records, each gene spans its exons.
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
        self.terminal_exons: set[tuple[str, str, int, int, str]] = set()
        # Each terminal exon's last-exon cluster, and each gene's 3'-most
        # cluster, as 1-based loci.
        self.last_exon_of: dict[tuple[str, str, int, int, str], str] = {}
        self.three_prime_last_exon: dict[tuple[str, str, str], str] = {}
        for gene_key, transcripts in transcripts_by_gene(exons).items():
            terminals = gene_terminal_exons(transcripts)
            clusters = last_exon_clusters(terminals)
            for exon in terminals:
                key = (exon.contig, exon.strand, exon.start, exon.end, exon.gene_id)
                self.terminal_exons.add(key)
                self.last_exon_of[key] = next(
                    locus for start, end, locus in clusters if start <= exon.start <= end
                )
            if clusters:
                three_prime = clusters[-1] if gene_key[1] == "+" else clusters[0]
                self.three_prime_last_exon[gene_key] = three_prime[2]
        self.exon_bins = self._bins(exons)
        self.gene_bins = self._bins(genes.values())
        by_strand: dict[tuple[str, str], list[GenomicFeature]] = defaultdict(list)
        for gene in genes.values():
            by_strand[(gene.contig, gene.strand)].append(gene)
        self.genes_by_end = {
            key: sorted(values, key=lambda gene: (gene.end, gene.gene_id))
            for key, values in by_strand.items()
        }
        self.gene_ends = {
            key: [gene.end for gene in values] for key, values in self.genes_by_end.items()
        }
        self.genes_by_start = {
            key: sorted(values, key=lambda gene: (gene.start, gene.gene_id))
            for key, values in by_strand.items()
        }
        self.gene_starts = {
            key: [gene.start for gene in values] for key, values in self.genes_by_start.items()
        }

    def _bins(
        self, features: Iterable[GenomicFeature]
    ) -> dict[tuple[str, str, int], list[GenomicFeature]]:
        bins: dict[tuple[str, str, int], list[GenomicFeature]] = defaultdict(list)
        for feature in features:
            for number in range(feature.start // self.bin_size, feature.end // self.bin_size + 1):
                bins[(feature.contig, feature.strand, number)].append(feature)
        return bins

    def _containing(
        self,
        bins: dict[tuple[str, str, int], list[GenomicFeature]],
        contig: str,
        strand: str,
        coordinate: int,
    ) -> list[GenomicFeature]:
        candidates = bins.get((contig, strand, coordinate // self.bin_size), [])
        return [item for item in candidates if item.start <= coordinate <= item.end]

    def last_exon(self, assignments: list[tuple[str, GenomicFeature]]) -> str:
        """The last exon of one gene that holds a PAC: its terminal exon's
        cluster, or for a PAC downstream of the gene, the gene's 3'-most
        cluster. Other PACs have none."""
        if not assignments:
            return ""
        label = assignments[0][0]
        if label == "last_exon":
            loci = {
                self.last_exon_of[
                    (exon.contig, exon.strand, exon.start, exon.end, exon.gene_id)
                ]
                for _, exon in assignments
            }
        elif label == "downstream_of_gene":
            loci = {
                self.three_prime_last_exon.get((gene.contig, gene.strand, gene.gene_id), "")
                for _, gene in assignments
            }
        else:
            return ""
        return loci.pop() if len(loci) == 1 else ""

    def assign(
        self, contig: str, strand: str, coordinate: int, maximum_downstream: int
    ) -> list[tuple[str, GenomicFeature]]:
        exons = self._containing(self.exon_bins, contig, strand, coordinate)
        terminal = [
            exon
            for exon in exons
            if (exon.contig, exon.strand, exon.start, exon.end, exon.gene_id)
            in self.terminal_exons
        ]
        if terminal:
            return [("last_exon", exon) for exon in terminal]
        if exons:
            return [("internal_exon", exon) for exon in exons]
        intronic = self._containing(self.gene_bins, contig, strand, coordinate)
        if intronic:
            return [("intron", gene) for gene in intronic]
        # Downstream: the nearest same-strand gene(s) ending upstream of the
        # PAC in transcript orientation, within maximum_downstream. No other
        # same-strand gene can lie between them: it would either contain the
        # PAC or end nearer to it.
        key = (contig, strand)
        if strand == "+":
            ends = self.gene_ends.get(key, [])
            upper = bisect_left(ends, coordinate)
            if upper and coordinate - ends[upper - 1] <= maximum_downstream:
                lower = bisect_left(ends, ends[upper - 1])
                genes = self.genes_by_end[key][lower:upper]
                return [("downstream_of_gene", gene) for gene in genes]
        else:
            starts = self.gene_starts.get(key, [])
            lower = bisect_right(starts, coordinate)
            if lower < len(starts) and starts[lower] - coordinate <= maximum_downstream:
                upper = bisect_right(starts, starts[lower])
                genes = self.genes_by_start[key][lower:upper]
                return [("downstream_of_gene", gene) for gene in genes]
        return []


def transcripts_by_gene(
    exons: Iterable[GenomicFeature],
) -> dict[tuple[str, str, str], list[list[GenomicFeature]]]:
    """Each gene's transcripts as lists of exons, keyed by (contig, strand,
    gene_id), in annotation order. Exons without a transcript ID form one
    transcript per gene."""
    transcript_exons: dict[tuple[str, str, str], list[GenomicFeature]] = defaultdict(list)
    for exon in exons:
        transcript_exons[(exon.contig, exon.strand, exon.transcript_id or exon.gene_id)].append(
            exon
        )
    grouped: dict[tuple[str, str, str], list[list[GenomicFeature]]] = defaultdict(list)
    for transcript in transcript_exons.values():
        first = transcript[0]
        grouped[(first.contig, first.strand, first.gene_id)].append(transcript)
    return grouped


def three_prime_end(exon: GenomicFeature) -> int:
    """An exon's 3' end in transcript orientation, as an interbase coordinate:
    its end on the plus strand, its start on the minus strand."""
    return exon.end if exon.strand == "+" else exon.start


def internal_exon_donors(
    exons: Iterable[GenomicFeature],
    tolerance: int,
) -> dict[tuple[str, str], list[int]]:
    """Each strand's internal exon 3' ends, sorted: the 3' end of every exon
    that a transcript splices onward from.

    A donor within ``tolerance`` of a terminal exon's 3' end, in any gene on
    the strand, is left out, since a transcript can end there. Genes without
    transcript IDs have no donors: their exons form one pseudo-transcript, so
    alternative last exons would look spliced onward.
    """
    donors: dict[tuple[str, str], set[int]] = defaultdict(set)
    ends: dict[tuple[str, str], list[int]] = defaultdict(list)
    for (contig, strand, _), transcripts in transcripts_by_gene(exons).items():
        ends[(contig, strand)].extend(
            three_prime_end(exon) for exon in gene_terminal_exons(transcripts)
        )
        if any(not exon.transcript_id for transcript in transcripts for exon in transcript):
            continue
        for transcript in transcripts:
            final = final_exon(transcript)
            donors[(contig, strand)].update(
                three_prime_end(exon) for exon in transcript if exon is not final
            )
    result: dict[tuple[str, str], list[int]] = {}
    for key, values in donors.items():
        terminal = sorted(ends[key])
        kept = [
            donor
            for donor in sorted(values)
            if bisect_right(terminal, donor + tolerance) == bisect_left(terminal, donor - tolerance)
        ]
        if kept:
            result[key] = kept
    return result


def final_exon(transcript: list[GenomicFeature]) -> GenomicFeature:
    """A transcript's 3'-most exon."""
    if transcript[0].strand == "+":
        return max(transcript, key=lambda item: item.end)
    return min(transcript, key=lambda item: item.start)


def exons_overlap(first: GenomicFeature, second: GenomicFeature) -> bool:
    return first.start < second.end and second.start < first.end


def gene_terminal_exons(transcripts: list[list[GenomicFeature]]) -> list[GenomicFeature]:
    """The final exons of a gene's transcripts that are true transcript ends.

    A final exon is left out when it overlaps an internal exon of another
    transcript of the gene, as the final exons of retained-intron and
    3'-incomplete models do, or when it is a single-exon model apart from
    every exon of the gene's multi-exon transcripts, such as a fragment inside
    an intron. Without transcript IDs a gene is one transcript, so only its
    3'-most exon is terminal.
    """
    finals = [final_exon(transcript) for transcript in transcripts]
    internal = [
        [exon for exon in transcript if exon is not finals[index]]
        for index, transcript in enumerate(transcripts)
    ]
    spliced = [exon for transcript in transcripts if len(transcript) > 1 for exon in transcript]
    terminals = []
    for index, final in enumerate(finals):
        inside_internal = any(
            exons_overlap(final, exon)
            for other, exons in enumerate(internal)
            if other != index
            for exon in exons
        )
        stray = (
            len(transcripts[index]) == 1
            and bool(spliced)
            and not any(exons_overlap(final, exon) for exon in spliced)
        )
        if not inside_internal and not stray:
            terminals.append(final)
    return terminals


def last_exon_clusters(terminals: list[GenomicFeature]) -> list[tuple[int, int, str]]:
    """A gene's terminal exons merged into last exons, in genomic order, as
    (start, end, 1-based locus). Exons sharing a coordinate merge, because
    containment is inclusive."""
    clusters: list[list[int]] = []
    for exon in sorted(terminals, key=lambda item: (item.start, item.end)):
        if clusters and exon.start <= clusters[-1][1]:
            clusters[-1][1] = max(clusters[-1][1], exon.end)
        else:
            clusters.append([exon.start, exon.end])
    contig = terminals[0].contig if terminals else ""
    return [(start, end, f"{contig}:{start + 1}-{end}") for start, end in clusters]


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
