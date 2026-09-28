"""Reference preparation and annotation parsing."""

from __future__ import annotations

import gzip
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path

import pysam

from .errors import PacusageError
from .tableio import sha256_file

# Ensembl GFF3 writes gene and transcript IDs with these prefixes, and writes
# gene records under the gene's biotype family.
GENE_ID_PREFIX = "gene:"
TRANSCRIPT_ID_PREFIX = "transcript:"
GENE_FEATURE_TYPES = {"gene", "ncrna_gene", "pseudogene"}
# Parts of a transcript, which never parent another feature.
TRANSCRIPT_PART_TYPES = {
    "exon",
    "cds",
    "utr",
    "five_prime_utr",
    "three_prime_utr",
    "start_codon",
    "stop_codon",
}


@dataclass(frozen=True)
class GenomicFeature:
    contig: str
    start: int
    end: int
    strand: str
    feature_type: str
    gene_id: str
    transcript_id: str = ""
    gene_name: str = ""


def prepare_reference(
    fasta: str | Path,
    output_fasta: str | Path,
) -> dict[str, str]:
    """Link the FASTA and index it in the task directory.

    The source is never copied or written to: ``faidx`` runs on the link, so the
    generated index lands next to ``output_fasta``, not next to the source.
    """
    source = Path(fasta).resolve()
    destination = Path(output_fasta)
    if not source.is_file():
        raise PacusageError(f"Genome FASTA does not exist: {source}")

    # Detect compression by magic bytes.
    with open(source, "rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        raise PacusageError(
            f"FASTA {source} is compressed. Decompress it (gunzip or bgzip -d) "
            "and set fasta to the uncompressed file."
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination_fai = Path(f"{destination}.fai")
    for path in (destination, destination_fai):
        if path.is_symlink():
            path.unlink()
        elif path.exists():
            if path.samefile(source) or Path(f"{source}.fai").resolve() == path.resolve():
                raise PacusageError(f"Reference preparation would replace its source: {path}")
            path.unlink()
    destination.symlink_to(source)

    try:
        pysam.faidx(str(destination))
    except Exception as error:
        raise PacusageError(
            f"FASTA {source} could not be indexed. Ensure it is uncompressed "
            "and contains valid FASTA records."
        ) from error

    if not destination_fai.is_file():
        raise PacusageError(f"Reference preparation did not create {destination_fai}.")
    return {
        "source_fasta": str(source),
        "prepared_fasta": str(destination),
        "prepared_fai": str(destination_fai),
        "action": "generated_index",
        "fasta_sha256": sha256_file(destination),
        "fai_sha256": sha256_file(destination_fai),
    }


def annotation_contigs(path: str | Path) -> set[str]:
    return {
        fields[0]
        for _, fields in _iter_annotation_rows(path)
        if fields[6] in {"+", "-"}
    }


def parse_annotation(
    path: str | Path, feature_types: set[str] | None = None
) -> list[GenomicFeature]:
    return list(iter_annotation_features(path, feature_types))


def iter_annotation_features(
    path: str | Path, feature_types: set[str] | None = None
) -> Iterator[GenomicFeature]:
    parents = _gff3_parents(path)
    selected_types = {_feature_type(value) for value in feature_types} if feature_types else None
    for line_number, fields in _iter_annotation_rows(path):
        contig, _, raw_type, start, end, _, strand, _, raw_attributes = fields
        feature_type = _feature_type(raw_type)
        if strand not in {"+", "-"} or (
            selected_types is not None and feature_type not in selected_types
        ):
            continue
        attributes = parse_attributes(raw_attributes)
        record_id = attributes.get("ID", "")
        parent = attributes.get("Parent", "")
        transcript_id = attributes.get("transcript_id", "")
        if not transcript_id and feature_type in {"transcript", "mrna"}:
            transcript_id = record_id.removeprefix(TRANSCRIPT_ID_PREFIX)
        if not transcript_id and feature_type == "exon" and parent:
            transcript_id = _parent_transcript(parent, parents)
        gene_id = attributes.get("gene_id") or attributes.get("gene", "")
        if not gene_id and feature_type == "gene":
            gene_id = record_id.removeprefix(GENE_ID_PREFIX)
        if not gene_id and parent:
            gene_id = _parent_gene(parent, parents)
        # A GFF3 Name is the gene's name only on a gene record; on an exon or
        # transcript it names that feature. NCBI annotations use gene.
        gene_name = attributes.get("gene_name") or attributes.get("gene", "")
        if not gene_name and feature_type == "gene":
            gene_name = attributes.get("Name", "")
        if not gene_id and feature_type in {"gene", "exon", "transcript", "mrna"}:
            raise PacusageError(
                f"Annotation {path}, line {line_number}: feature lacks a gene identifier."
            )
        yield GenomicFeature(
            contig=contig,
            start=int(start) - 1,
            end=int(end),
            strand=strand,
            feature_type=feature_type,
            gene_id=gene_id,
            transcript_id=transcript_id,
            gene_name=gene_name,
        )


@dataclass(frozen=True)
class _ParentRecord:
    """A GFF3 feature that other features can name as their Parent."""

    feature_type: str
    parent: str
    gene_id: str
    transcript_id: str


def _gff3_parents(path: str | Path) -> dict[str, _ParentRecord]:
    """Every GFF3 feature that can be a Parent, by ID. GTF lines carry no ID."""
    parents: dict[str, _ParentRecord] = {}
    for _, fields in _iter_annotation_rows(path):
        feature_type = _feature_type(fields[2])
        if feature_type in TRANSCRIPT_PART_TYPES:
            continue
        attributes = parse_attributes(fields[8])
        record_id = attributes.get("ID", "")
        if record_id:
            parents[record_id] = _ParentRecord(
                feature_type=feature_type,
                parent=attributes.get("Parent", ""),
                gene_id=attributes.get("gene_id") or attributes.get("gene", ""),
                transcript_id=attributes.get("transcript_id", ""),
            )
    return parents


def _parent_gene(record_id: str, parents: Mapping[str, _ParentRecord]) -> str:
    """The gene a GFF3 Parent chain leads to, or "" when it leads nowhere."""
    seen: set[str] = set()
    while record_id in parents and record_id not in seen:
        seen.add(record_id)
        record = parents[record_id]
        if record.gene_id:
            return record.gene_id
        if record.feature_type == "gene":
            return record_id.removeprefix(GENE_ID_PREFIX)
        record_id = record.parent
    return ""


def _parent_transcript(record_id: str, parents: Mapping[str, _ParentRecord]) -> str:
    record = parents.get(record_id)
    if record is not None and record.transcript_id:
        return record.transcript_id
    return record_id.removeprefix(TRANSCRIPT_ID_PREFIX)


def _feature_type(raw: str) -> str:
    feature_type = raw.lower()
    return "gene" if feature_type in GENE_FEATURE_TYPES else feature_type


def _iter_annotation_rows(path: str | Path) -> Iterator[tuple[int, list[str]]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise PacusageError(
                    f"Annotation {path}, line {line_number}: expected 9 tab-separated fields."
                )
            yield line_number, fields


def parse_attributes(raw: str) -> dict[str, str]:
    values: dict[str, str] = {}
    if "=" in raw and ' "' not in raw:
        for item in raw.strip().strip(";").split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                values[key.strip()] = value.strip().split(",")[0]
        return values
    for item in raw.strip().strip(";").split(";"):
        parts = item.strip().split(None, 1)
        if len(parts) == 2:
            values[parts[0]] = parts[1].strip().strip('"')
    return values


def gene_name_map(features: Iterable[GenomicFeature]) -> dict[str, str]:
    """One name per gene: its gene record's name, else the first name on its
    other records, else the gene ID itself."""
    from_gene_records: dict[str, str] = {}
    from_other_records: dict[str, str] = {}
    gene_ids: dict[str, None] = {}
    for feature in features:
        if not feature.gene_id:
            continue
        gene_ids.setdefault(feature.gene_id)
        names = from_gene_records if feature.feature_type == "gene" else from_other_records
        if feature.gene_name:
            names.setdefault(feature.gene_id, feature.gene_name)
    return {
        gene_id: from_gene_records.get(gene_id) or from_other_records.get(gene_id) or gene_id
        for gene_id in gene_ids
    }


def transcript_ends(
    features: Iterable[GenomicFeature],
) -> dict[tuple[str, str], list[tuple[int, str]]]:
    transcript_ranges: dict[tuple[str, str, str], list[int]] = {}
    transcript_gene: dict[tuple[str, str, str], str] = {}
    for feature in features:
        if feature.feature_type not in {"exon", "transcript", "mrna"}:
            continue
        transcript = feature.transcript_id or feature.gene_id
        key = (feature.contig, feature.strand, transcript)
        if key not in transcript_ranges:
            transcript_ranges[key] = [feature.start, feature.end]
        else:
            transcript_ranges[key][0] = min(transcript_ranges[key][0], feature.start)
            transcript_ranges[key][1] = max(transcript_ranges[key][1], feature.end)
        transcript_gene[key] = feature.gene_id
    gene_ends: dict[tuple[str, str, str], set[int]] = defaultdict(set)
    for key, (start, end) in transcript_ranges.items():
        contig, strand, _ = key
        coordinate = end if strand == "+" else start
        gene_ends[(contig, strand, transcript_gene[key])].add(coordinate)
    result: dict[tuple[str, str], list[tuple[int, str]]] = defaultdict(list)
    for (contig, strand, gene_id), coordinates in gene_ends.items():
        if len(coordinates) == 1:
            result[(contig, strand)].append((next(iter(coordinates)), gene_id))
    for key in result:
        result[key] = sorted(set(result[key]))
    return result

