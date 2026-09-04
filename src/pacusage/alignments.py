"""Alignment inspection, sorting, indexing, and source preservation."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import pysam

from .errors import PacusageError
from .tableio import sha256_file


def inspect_alignment(path: str | Path, reference: str | Path | None = None) -> dict[str, Any]:
    path = Path(path)
    mode = "rc" if path.suffix.lower() == ".cram" else "rb"
    try:
        with pysam.AlignmentFile(
            str(path), mode, reference_filename=str(reference) if reference else None
        ) as handle:
            header = handle.header.to_dict()
            sort_order = header.get("HD", {}).get("SO", "unknown")
            contigs = dict(zip(handle.references, handle.lengths, strict=True))
            paired = 0
            unpaired = 0
            for index, record in enumerate(handle.fetch(until_eof=True)):
                if record.is_paired:
                    paired += 1
                else:
                    unpaired += 1
                if index >= 9999:
                    break
    except Exception as error:
        raise PacusageError(
            f"Cannot read alignment {path}. Confirm it is an intact BAM/CRAM and "
            "that CRAM uses the supplied reference."
        ) from error
    total = paired + unpaired
    if total == 0:
        layout = "unknown"
    elif paired / total >= 0.95:
        layout = "PE"
    elif unpaired / total >= 0.95:
        layout = "SE"
    else:
        layout = "mixed"
    return {
        "sort_order": sort_order,
        "contigs": contigs,
        "layout": layout,
        "sampled_records": total,
    }


def prepare_alignment(
    source: str | Path,
    output: str | Path,
    reference: str | Path,
    threads: int = 1,
) -> dict[str, Any]:
    source = Path(source).resolve()
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    before = sha256_file(source)
    details = inspect_alignment(source, reference)
    is_cram = source.suffix.lower() == ".cram"
    if destination.exists() and destination.samefile(source):
        # Nextflow may stage the input under the requested output name as a
        # symlink. Remove only that work-directory link before creating output.
        destination.unlink()
    if details["sort_order"] == "coordinate":
        shutil.copyfile(source, destination)
        action = "reused_alignment"
    else:
        args = ["-@", str(max(1, threads)), "-o", str(destination)]
        if is_cram:
            args.extend(["-O", "CRAM", "--reference", str(reference)])
        args.append(str(source))
        try:
            pysam.sort(*args, catch_stdout=False)
        except Exception as error:
            raise PacusageError(f"Failed to coordinate-sort alignment {source}.") from error
        action = "sorted_and_indexed"

    index_path = Path(f"{destination}.crai" if is_cram else f"{destination}.bai")
    sibling_index = _find_index(source)
    reusable_index = (
        action == "reused_alignment"
        and sibling_index
        and _index_is_usable(source, sibling_index, reference)
    )
    if reusable_index:
        shutil.copyfile(sibling_index, index_path)
        action = "reused_alignment_and_index"
    else:
        try:
            pysam.index("-@", str(max(1, threads)), str(destination), str(index_path))
        except Exception as error:
            raise PacusageError(f"Failed to index prepared alignment {destination}.") from error
        if action == "reused_alignment":
            action = "indexed"

    after = sha256_file(source)
    if before != after:
        raise PacusageError(f"Source alignment changed during preparation: {source}")
    prepared = inspect_alignment(destination, reference)
    if prepared["sort_order"] != "coordinate":
        raise PacusageError(f"Prepared alignment is not coordinate sorted: {destination}")
    return {
        "source_alignment": str(source),
        "prepared_alignment": str(destination),
        "prepared_index": str(index_path),
        "format": "CRAM" if is_cram else "BAM",
        "action": action,
        "layout_detected": prepared["layout"],
        "source_sha256": before,
        "prepared_sha256": sha256_file(destination),
    }


def _find_index(path: Path) -> Path | None:
    candidates = (
        [Path(f"{path}.crai"), path.with_suffix(".crai")]
        if path.suffix.lower() == ".cram"
        else [Path(f"{path}.bai"), path.with_suffix(".bai")]
    )
    return next((candidate for candidate in candidates if candidate.is_file()), None)


def _index_is_usable(path: Path, index: Path, reference: str | Path) -> bool:
    mode = "rc" if path.suffix.lower() == ".cram" else "rb"
    try:
        with pysam.AlignmentFile(
            str(path),
            mode,
            index_filename=str(index),
            reference_filename=str(reference),
        ) as handle:
            next(iter(handle.fetch()), None)
        return True
    except Exception:
        return False


def validate_contigs(
    alignment_contigs: dict[str, int],
    fasta_contigs: dict[str, int],
    annotation_contigs: set[str],
    aliases: dict[str, str] | None = None,
) -> None:
    aliases = aliases or {}
    normalized = {aliases.get(name, name): length for name, length in alignment_contigs.items()}
    absent = sorted(set(normalized).difference(fasta_contigs))
    length_mismatch = sorted(
        name
        for name, length in normalized.items()
        if name in fasta_contigs and fasta_contigs[name] != length
    )
    annotation_overlap = set(normalized).intersection(annotation_contigs)
    if absent or length_mismatch or not annotation_overlap:
        details = []
        if absent:
            details.append(f"alignment contigs absent from FASTA: {', '.join(absent[:10])}")
        if length_mismatch:
            details.append(f"length mismatches: {', '.join(length_mismatch[:10])}")
        if not annotation_overlap:
            details.append("no alignment contigs overlap annotation contigs")
        raise PacusageError("Reference/alignment contig incompatibility: " + "; ".join(details))
