"""Preparation links sources instead of copying them, and never writes beside them."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest
from alignment_builders import aligned_segment, sorted_copy, write_alignment, write_reference

from pacusage import alignments
from pacusage.alignments import check_prepared_source, prepare_alignment
from pacusage.cli import main
from pacusage.errors import PacusageError
from pacusage.reference import prepare_reference
from pacusage.tableio import read_tsv, sha256_file, write_tsv

CONTIGS = [("chr1", 2000)]


def directory_state(directory: Path) -> dict[str, tuple[str, int]]:
    """Every file's hash and modification time, to prove nothing changed."""
    return {
        path.name: (sha256_file(path), path.stat().st_mtime_ns)
        for path in sorted(directory.iterdir())
    }


def reads() -> list:
    return [aligned_segment(f"r{index}", 0, 0, 100 + index, "30M") for index in (5, 1, 3)]


@pytest.fixture()
def sources(tmp_path: Path) -> dict[str, Path]:
    source = tmp_path / "source"
    source.mkdir()
    reference = write_reference(source / "genome.fa", CONTIGS)
    unsorted = write_alignment(source / "unsorted.bam", CONTIGS, reads())
    ordered = sorted_copy(unsorted, source / "sorted.bam")
    os.remove(f"{reference}.fai")
    return {"directory": source, "reference": reference, "unsorted": unsorted, "sorted": ordered}


def test_sorted_alignment_is_linked_and_hashed_once(tmp_path, sources, monkeypatch) -> None:
    import pysam

    pysam.index(str(sources["sorted"]))
    before = directory_state(sources["directory"])
    hashed = []
    real_sha256 = alignments.sha256_file
    monkeypatch.setattr(
        alignments, "sha256_file", lambda path: hashed.append(Path(path)) or real_sha256(path)
    )
    output = tmp_path / "work" / "S.bam"
    reference = tmp_path / "work" / "genome.fa"
    prepare_reference(sources["reference"], reference)
    metadata = prepare_alignment(sources["sorted"], output, reference)

    assert output.is_symlink()
    assert os.readlink(output) == str(sources["sorted"].resolve())
    index = Path(f"{output}.bai")
    assert index.is_symlink()
    assert os.readlink(index) == str(Path(f"{sources['sorted']}.bai").resolve())
    assert metadata["action"] == "reused_alignment_and_index"
    assert hashed == [sources["sorted"].resolve()]
    assert metadata["prepared_sha256"] == metadata["source_sha256"]
    assert metadata["source_sha256"] == sha256_file(sources["sorted"])
    status = sources["sorted"].stat()
    assert (metadata["source_size"], metadata["source_mtime_ns"]) == (
        status.st_size,
        status.st_mtime_ns,
    )
    assert directory_state(sources["directory"]) == before


def test_sorted_alignment_without_index_is_indexed_in_the_task_directory(
    tmp_path, sources
) -> None:
    before = directory_state(sources["directory"])
    output = tmp_path / "work" / "S.bam"
    metadata = prepare_alignment(sources["sorted"], output, sources["reference"])
    assert output.is_symlink()
    index = Path(f"{output}.bai")
    assert index.is_file() and not index.is_symlink()
    assert metadata["action"] == "indexed"
    assert directory_state(sources["directory"]) == before


def test_unsorted_alignment_is_sorted_into_a_real_file(tmp_path, sources) -> None:
    before = directory_state(sources["directory"])
    output = tmp_path / "work" / "S.bam"
    metadata = prepare_alignment(sources["unsorted"], output, sources["reference"])
    assert output.is_file() and not output.is_symlink()
    assert metadata["action"] == "sorted_and_indexed"
    assert metadata["prepared_sha256"] == sha256_file(output)
    assert metadata["prepared_sha256"] != metadata["source_sha256"]
    assert directory_state(sources["directory"]) == before


def test_staged_input_link_is_replaced_not_written_through(tmp_path, sources) -> None:
    # Nextflow stages the input under the output name as a symlink.
    work = tmp_path / "work"
    work.mkdir()
    staged = work / "S.bam"
    staged.symlink_to(sources["unsorted"])
    before = directory_state(sources["directory"])
    prepare_alignment(staged, staged, sources["reference"])
    assert staged.is_file() and not staged.is_symlink()
    assert directory_state(sources["directory"]) == before


def test_reference_is_linked_and_indexed_beside_the_link(tmp_path, sources) -> None:
    before = directory_state(sources["directory"])
    output = tmp_path / "work" / "genome.fa"
    metadata = prepare_reference(sources["reference"], output)
    assert output.is_symlink()
    assert os.readlink(output) == str(sources["reference"].resolve())
    assert Path(f"{output}.fai").is_file()
    assert metadata["action"] == "generated_index"
    assert metadata["fasta_sha256"] == sha256_file(sources["reference"])
    assert directory_state(sources["directory"]) == before


def test_reference_reuses_a_valid_source_index(tmp_path, sources) -> None:
    import pysam

    pysam.faidx(str(sources["reference"]))
    before = directory_state(sources["directory"])
    output = tmp_path / "work" / "genome.fa"
    metadata = prepare_reference(sources["reference"], output)
    assert metadata["action"] == "reused_index"
    assert Path(f"{output}.fai").read_bytes() == Path(f"{sources['reference']}.fai").read_bytes()
    assert directory_state(sources["directory"]) == before


def test_reference_refuses_to_replace_its_source(sources) -> None:
    with pytest.raises(PacusageError, match="would replace its source"):
        prepare_reference(sources["reference"], sources["reference"])


def _record(tmp_path: Path, metadata: dict) -> Path:
    path = tmp_path / "S.alignment_preparation.tsv"
    write_tsv([metadata], path)
    return path


def test_source_check_accepts_an_unchanged_link(tmp_path, sources) -> None:
    output = tmp_path / "work" / "S.bam"
    metadata = prepare_alignment(sources["sorted"], output, sources["reference"])
    check_prepared_source(output, _record(tmp_path, metadata))


def test_source_check_rejects_a_changed_source(tmp_path, sources) -> None:
    output = tmp_path / "work" / "S.bam"
    metadata = prepare_alignment(sources["sorted"], output, sources["reference"])
    record = _record(tmp_path, metadata)
    status = sources["sorted"].stat()
    os.utime(sources["sorted"], ns=(status.st_atime_ns, status.st_mtime_ns + 1_000_000_000))
    with pytest.raises(PacusageError, match="changed after preparation"):
        check_prepared_source(output, record)


def test_source_check_explains_a_missing_source(tmp_path, sources) -> None:
    output = tmp_path / "work" / "S.bam"
    metadata = prepare_alignment(sources["sorted"], output, sources["reference"])
    record = _record(tmp_path, metadata)
    sources["sorted"].rename(sources["directory"] / "moved.bam")
    with pytest.raises(PacusageError, match="bind_paths"):
        check_prepared_source(output, record)


def test_source_check_ignores_a_sorted_copy(tmp_path, sources) -> None:
    output = tmp_path / "work" / "S.bam"
    metadata = prepare_alignment(sources["unsorted"], output, sources["reference"])
    record = _record(tmp_path, metadata)
    # The source may change: the work directory holds its own sorted copy.
    sources["unsorted"].write_bytes(b"changed")
    check_prepared_source(output, record)


def test_validate_hashes_no_alignment_and_writes_nothing_beside_the_fasta(
    tmp_path, sources, monkeypatch
) -> None:
    import pysam

    cram = sources["directory"] / "S.cram"
    pysam.view(
        "-C", "-T", str(sources["reference"]), "-o", str(cram), str(sources["sorted"]),
        catch_stdout=False,
    )
    os.remove(f"{sources['reference']}.fai")
    for name in ("treated_1.bam", "treated_2.bam"):
        (sources["directory"] / name).write_bytes(sources["sorted"].read_bytes())
    sheet = sources["directory"] / "samples.tsv"
    sheet.write_text(
        "sample_id\talignment\tcondition\tcontrol\n"
        "S1\tsorted.bam\tA\t\n"
        "S2\tS.cram\tA\t\n"
        "S3\ttreated_1.bam\tB\tA\n"
        "S4\ttreated_2.bam\tB\tA\n"
    )
    annotation = sources["directory"] / "genes.gtf"
    annotation.write_text(
        'chr1\ttest\texon\t101\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    params = tmp_path / "params.json"
    params.write_text(
        json.dumps(
            {
                "input": str(sheet),
                "assembly": "test",
                "fasta": str(sources["reference"]),
                "gtf": str(annotation),
            }
        )
    )
    before = directory_state(sources["directory"])
    hashed = []
    real_sha256 = sha256_file
    monkeypatch.setattr(
        "pacusage.cli.sha256_file", lambda path: hashed.append(Path(path).name) or real_sha256(path)
    )
    output = tmp_path / "validated"
    assert main(["validate", "--params", str(params), "--output-dir", str(output)]) == 0

    assert directory_state(sources["directory"]) == before
    assert not any(name.endswith((".bam", ".cram")) for name in hashed)
    rows = read_tsv(output / "input_checksums.tsv")
    assert [row["role"] for row in rows] == ["sample_sheet", "fasta", "annotation"]
    assert rows[1]["sha256"] == sha256_file(sources["reference"])


def test_source_check_uses_the_recorded_action_not_the_resolved_path(tmp_path, sources) -> None:
    # Inside a container a linked source can resolve to a different path; the
    # record's action, not path equality, says whether a source is linked.
    output = tmp_path / "work" / "S.bam"
    metadata = prepare_alignment(sources["sorted"], output, sources["reference"])
    record = _record(tmp_path, {**metadata, "source_alignment": str(tmp_path / "gone.bam")})
    with pytest.raises(PacusageError, match="bind_paths"):
        check_prepared_source(output, record)
