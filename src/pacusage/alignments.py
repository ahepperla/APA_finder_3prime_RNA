"""Alignment inspection, sorting, indexing, and source preservation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pysam

from .errors import PacusageError
from .tableio import read_tsv, sha256_file


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


def source_fingerprint(path: str | Path) -> dict[str, int]:
    """Size and modification time, which also key Nextflow's task cache."""
    status = Path(path).stat()
    return {"source_size": status.st_size, "source_mtime_ns": status.st_mtime_ns}


def prepare_alignment(
    source: str | Path,
    output: str | Path,
    reference: str | Path,
    threads: int = 1,
) -> dict[str, Any]:
    """Link a coordinate-sorted source, or sort it, and index the result.

    A sorted source is not copied: the output is an absolute symlink to it, so
    the source must stay in place and unmodified for the rest of the analysis.
    """
    source = Path(source).resolve()
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fingerprint = source_fingerprint(source)
    source_sha256 = sha256_file(source)
    details = inspect_alignment(source, reference)
    is_cram = source.suffix.lower() == ".cram"
    # Nextflow stages the input under the requested output name. Remove that
    # work-directory link so nothing is ever written through it.
    _remove_staged_link(destination, source)
    if details["sort_order"] == "coordinate":
        destination.symlink_to(source)
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
    _remove_staged_link(index_path, None)
    sibling_index = _find_index(source)
    reusable_index = (
        action == "reused_alignment"
        and sibling_index
        and _index_is_usable(source, sibling_index, reference)
    )
    if reusable_index:
        index_path.symlink_to(sibling_index.resolve())
        action = "reused_alignment_and_index"
    else:
        try:
            pysam.index("-@", str(max(1, threads)), str(destination), str(index_path))
        except Exception as error:
            raise PacusageError(f"Failed to index prepared alignment {destination}.") from error
        if action == "reused_alignment":
            action = "indexed"

    if source_fingerprint(source) != fingerprint:
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
        "source_sha256": source_sha256,
        "prepared_sha256": (
            source_sha256 if destination.is_symlink() else sha256_file(destination)
        ),
        **fingerprint,
    }


def check_prepared_source(alignment: str | Path, metadata: str | Path) -> None:
    """Fail if a linked source alignment moved or changed after preparation."""
    rows = read_tsv(metadata)
    if len(rows) != 1:
        raise PacusageError(f"Alignment preparation record {metadata} must contain one row.")
    row = rows[0]
    source = row["source_alignment"]
    if row.get("action") == "sorted_and_indexed":
        # A sorted copy in the work directory, which the pipeline owns. Every
        # other action links the source, whose path may resolve differently
        # inside a container, so it is checked by its recorded path.
        return
    if not Path(source).is_file():
        raise PacusageError(
            f"Source alignment {source} is not readable here. Prepared alignments link to "
            "their sources, so sources must stay in place, and container runs need their "
            "directories in bind_paths."
        )
    try:
        recorded = {key: int(row[key]) for key in ("source_size", "source_mtime_ns")}
    except (KeyError, ValueError) as error:
        raise PacusageError(
            f"Alignment preparation record {metadata} lacks the source fingerprint."
        ) from error
    if source_fingerprint(source) != recorded:
        raise PacusageError(
            f"Source alignment {source} changed after preparation. Restore it, or rerun "
            "with -resume so PREPARE_ALIGNMENT records the current file."
        )


def _remove_staged_link(path: Path, source: Path | None) -> None:
    if path.is_symlink():
        path.unlink()
    elif path.exists():
        if source is not None and path.samefile(source):
            path.unlink()
        else:
            raise PacusageError(f"Refusing to replace an existing file: {path}")


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
