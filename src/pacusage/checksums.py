"""Complete the input checksum manifest with per-sample alignment hashes."""

from __future__ import annotations

from .errors import PacusageError


def record_input_checksums(
    base_rows: list[dict[str, str]],
    sample_rows: list[dict[str, str]],
    metadata_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Append one alignment row per sample, in sample-sheet order.

    ``base_rows`` are VALIDATE's sample-sheet, FASTA, and GTF rows. The hashes
    come from PREPARE_ALIGNMENT, which reads each alignment once.
    """
    if not base_rows:
        raise PacusageError("Base input checksums are empty.")
    if any(row.get("role") == "alignment" for row in base_rows):
        raise PacusageError(
            "Base input checksums already contain alignment rows. VALIDATE_INPUTS was "
            "probably reused from an older PACusage version; rerun without -resume."
        )

    metadata_by_sample: dict[str, dict[str, str]] = {}
    for metadata in metadata_rows:
        sample_id = metadata.get("sample_id", "").strip()
        if not sample_id:
            raise PacusageError("Metadata row has an empty sample_id.")
        for field in ("source_alignment", "source_sha256"):
            if not metadata.get(field, "").strip():
                raise PacusageError(f"Metadata for sample {sample_id} has an empty {field}.")
        if sample_id in metadata_by_sample:
            raise PacusageError(f"Duplicate metadata entries for sample {sample_id}.")
        metadata_by_sample[sample_id] = metadata

    sample_ids = {sample.get("sample_id", "") for sample in sample_rows}
    for sample_id in metadata_by_sample:
        if sample_id not in sample_ids:
            raise PacusageError(
                f"Metadata row names sample {sample_id}, which is not in the sample sheet."
            )

    output = list(base_rows)
    for sample in sample_rows:
        sample_id = sample.get("sample_id", "")
        alignment = sample.get("alignment", "")
        metadata = metadata_by_sample.get(sample_id)
        if metadata is None:
            raise PacusageError(f"No metadata row found for sample {sample_id}.")
        if metadata["source_alignment"] != alignment:
            raise PacusageError(
                f"Sample {sample_id}: source_alignment mismatch. "
                f"Metadata has {metadata['source_alignment']!r}, "
                f"but sample has {alignment!r}."
            )
        output.append(
            {"role": "alignment", "path": alignment, "sha256": metadata["source_sha256"]}
        )
    return output
