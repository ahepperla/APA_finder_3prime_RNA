"""Reference preparation and annotation parsing."""

from __future__ import annotations

import gzip
import shutil
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import pysam

from .errors import PacusageError
from .tableio import sha256_file


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
    supplied_fai: str | Path | None = None,
) -> dict[str, str]:
    source = Path(fasta).resolve()
    destination = Path(output_fasta)
    if not source.is_file():
        raise PacusageError(f"Genome FASTA does not exist: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    destination_fai = Path(f"{destination}.fai")

    action = "generated_index"
    if supplied_fai and Path(supplied_fai).is_file():
        if validate_fai(source, supplied_fai):
            shutil.copyfile(supplied_fai, destination_fai)
            action = "reused_index"
        else:
            raise PacusageError(
                f"Supplied FASTA index {supplied_fai} does not match {source}. "
                "Remove it or provide the matching index."
            )
    else:
        sibling = Path(f"{source}.fai")
        if sibling.is_file() and validate_fai(source, sibling):
            shutil.copyfile(sibling, destination_fai)
            action = "reused_index"
        else:
            try:
                pysam.faidx(str(destination))
            except Exception as error:
                raise PacusageError(
                    f"FASTA {source} could not be indexed. Ensure it is uncompressed "
                    "or BGZF-compressed and contains valid FASTA records."
                ) from error
    if not destination_fai.is_file():
        raise PacusageError(f"Reference preparation did not create {destination_fai}.")
    return {
        "source_fasta": str(source),
        "prepared_fasta": str(destination),
        "prepared_fai": str(destination_fai),
        "action": action,
        "fasta_sha256": sha256_file(destination),
        "fai_sha256": sha256_file(destination_fai),
    }


def validate_fai(fasta: str | Path, fai: str | Path) -> bool:
    try:
        expected = _scan_fasta_lengths(fasta)
        observed: dict[str, int] = {}
        with Path(fai).open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    return False
                observed[fields[0]] = int(fields[1])
        return observed == expected
    except (OSError, ValueError):
        return False


def _scan_fasta_lengths(path: str | Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name: str | None = None
    length = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                if not name or name in lengths:
                    raise ValueError("Invalid or duplicate FASTA sequence name")
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        lengths[name] = length
    if not lengths:
        raise ValueError("No FASTA records")
    return lengths


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
    transcript_to_gene: dict[str, str] = {}
    for _, fields in _iter_annotation_rows(path):
        feature_type = fields[2].lower()
        if feature_type not in {"transcript", "mrna"}:
            continue
        raw_attributes = fields[8]
        attributes = parse_attributes(raw_attributes)
        transcript_id = attributes.get("transcript_id") or attributes.get("ID", "")
        parent_gene = attributes.get("gene_id") or attributes.get("Parent", "")
        if transcript_id and parent_gene:
            transcript_to_gene[transcript_id] = parent_gene

    selected_types = {value.lower() for value in feature_types} if feature_types else None
    for line_number, fields in _iter_annotation_rows(path):
        contig, _, feature_type, start, end, _, strand, _, raw_attributes = fields
        feature_type = feature_type.lower()
        if strand not in {"+", "-"} or (
            selected_types is not None and feature_type not in selected_types
        ):
            continue
        attributes = parse_attributes(raw_attributes)
        transcript_id = attributes.get("transcript_id", "")
        if not transcript_id and feature_type in {"transcript", "mrna"}:
            transcript_id = attributes.get("ID", "")
        if not transcript_id and feature_type == "exon":
            transcript_id = attributes.get("Parent", "").split(",")[0]
        gene_id = attributes.get("gene_id") or attributes.get("gene", "")
        if not gene_id and feature_type == "gene":
            gene_id = attributes.get("ID", "")
        if not gene_id and transcript_id:
            gene_id = transcript_to_gene.get(transcript_id, "")
        gene_name = attributes.get("gene_name") or attributes.get("Name", "")
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


def transcript_ends(features: list[GenomicFeature]) -> dict[tuple[str, str], list[tuple[int, str]]]:
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


def iter_exons(features: list[GenomicFeature]) -> Iterator[GenomicFeature]:
    return (feature for feature in features if feature.feature_type == "exon")
