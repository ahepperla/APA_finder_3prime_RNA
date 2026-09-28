"""Tests for recording per-sample alignment hashes."""

import pytest

from pacusage.checksums import record_input_checksums
from pacusage.cli import main
from pacusage.errors import PacusageError


def test_record_input_checksums_writes_validate_format(tmp_path) -> None:
    """Test that record_input_checksums produces the expected output format."""
    # Write base checksums TSV (sample_sheet, fasta, annotation)
    base_file = tmp_path / "base_checksums.tsv"
    base_file.write_text(
        "role\tpath\tsha256\n"
        "sample_sheet\t/path/to/samples.tsv\taaa111aaa111aaa111aaa111aaa111aaa111aaa111aaa111aaa111aaa111aaa1\n"
        "fasta\t/path/to/genome.fa\tbbb222bbb222bbb222bbb222bbb222bbb222bbb222bbb222bbb222bbb222bbb2\n"
        "annotation\t/path/to/genes.gtf\tccc333ccc333ccc333ccc333ccc333ccc333ccc333ccc333ccc333ccc333ccc3\n"
    )

    # Write normalized samples TSV with two samples
    samples_file = tmp_path / "normalized_samples.tsv"
    samples_file.write_text(
        "sample_id\talignment\tcondition\tcontrol\tcontrol_condition\treplicate\tbatch\tdonor\tlayout\tstrandedness\tlibrary_profile\tevidence_source\n"
        "S1\t/path/to/S1.bam\ttreatment\tcontrol\tcontrol\t1\t1\tdonor1\tSE\tunstranded\tprofile1\tauto\n"
        "S2\t/path/to/S2.bam\tcontrol\t\t\t1\t1\tdonor2\tSE\tunstranded\tprofile1\tauto\n"
    )

    # Write metadata for S1
    metadata_s1 = tmp_path / "S1_metadata.tsv"
    metadata_s1.write_text(
        "source_alignment\tprepared_alignment\tprepared_index\tformat\taction\tlayout_detected\tsource_sha256\tprepared_sha256\tsample_id\n"
        "/path/to/S1.bam\t/cache/S1.sorted.bam\t/cache/S1.sorted.bam.bai\tBAM\tsorted_and_indexed\tSE\tddd444ddd444ddd444ddd444ddd444ddd444ddd444ddd444ddd444ddd444ddd4\teee555eee555eee555eee555eee555eee555eee555eee555eee555eee555eee5\tS1\n"
    )

    # Write metadata for S2
    metadata_s2 = tmp_path / "S2_metadata.tsv"
    metadata_s2.write_text(
        "source_alignment\tprepared_alignment\tprepared_index\tformat\taction\tlayout_detected\tsource_sha256\tprepared_sha256\tsample_id\n"
        "/path/to/S2.bam\t/cache/S2.sorted.bam\t/cache/S2.sorted.bam.bai\tBAM\treused_alignment\tSE\tfff666fff666fff666fff666fff666fff666fff666fff666fff666fff666fff6\tggg777ggg777ggg777ggg777ggg777ggg777ggg777ggg777ggg777ggg777ggg7\tS2\n"
    )

    output_file = tmp_path / "output_checksums.tsv"

    # Run the command with metadata files in reverse order (S2 first)
    ret = main(
        [
            "record-input-checksums",
            "--base", str(base_file),
            "--samples", str(samples_file),
            "--alignment-metadata", str(metadata_s2), str(metadata_s1),
            "--output", str(output_file),
        ]
    )

    assert ret == 0
    output_text = output_file.read_text()

    # Expected output: header + base rows + alignment rows
    # (in sample_rows order, not metadata order)
    expected = (
        "role\tpath\tsha256\n"
        "sample_sheet\t/path/to/samples.tsv\taaa111aaa111aaa111aaa111aaa111aaa111aaa111aaa111aaa111aaa111aaa1\n"
        "fasta\t/path/to/genome.fa\tbbb222bbb222bbb222bbb222bbb222bbb222bbb222bbb222bbb222bbb222bbb2\n"
        "annotation\t/path/to/genes.gtf\tccc333ccc333ccc333ccc333ccc333ccc333ccc333ccc333ccc333ccc333ccc3\n"
        "alignment\t/path/to/S1.bam\tddd444ddd444ddd444ddd444ddd444ddd444ddd444ddd444ddd444ddd444ddd4\n"
        "alignment\t/path/to/S2.bam\tfff666fff666fff666fff666fff666fff666fff666fff666fff666fff666fff6\n"
    )

    assert output_text == expected


def test_record_input_checksums_empty_base_rows() -> None:
    """Test that empty base_rows raises PacusageError."""
    with pytest.raises(PacusageError, match="Base input checksums are empty"):
        record_input_checksums([], [], [])


def test_record_input_checksums_mismatched_alignment_path() -> None:
    """Test that mismatched source_alignment raises PacusageError naming the sample."""
    base = [
        {"role": "sample_sheet", "path": "/path", "sha256": "abc123"},
        {"role": "fasta", "path": "/path", "sha256": "def456"},
        {"role": "annotation", "path": "/path", "sha256": "ghi789"},
    ]
    samples = [
        {"sample_id": "S1", "alignment": "/path/to/S1.bam"},
    ]
    metadata = [
        {
            "sample_id": "S1",
            "source_alignment": "/path/to/different.bam",
            "source_sha256": "jkl000",
        },
    ]
    with pytest.raises(
        PacusageError,
        match="Sample S1: source_alignment mismatch. Metadata has '/path/to/different.bam', "
        "but sample has '/path/to/S1.bam'",
    ):
        record_input_checksums(base, samples, metadata)


def test_record_input_checksums_sample_without_metadata() -> None:
    """Test that a sample without metadata raises PacusageError."""
    base = [
        {"role": "sample_sheet", "path": "/path", "sha256": "abc123"},
        {"role": "fasta", "path": "/path", "sha256": "def456"},
        {"role": "annotation", "path": "/path", "sha256": "ghi789"},
    ]
    samples = [
        {"sample_id": "S1", "alignment": "/path/to/S1.bam"},
    ]
    metadata = []
    with pytest.raises(PacusageError, match="No metadata row found for sample S1"):
        record_input_checksums(base, samples, metadata)


def test_record_input_checksums_metadata_for_unknown_sample() -> None:
    """Test that metadata for an unknown sample raises PacusageError."""
    base = [
        {"role": "sample_sheet", "path": "/path", "sha256": "abc123"},
        {"role": "fasta", "path": "/path", "sha256": "def456"},
        {"role": "annotation", "path": "/path", "sha256": "ghi789"},
    ]
    samples = [
        {"sample_id": "S1", "alignment": "/path/to/S1.bam"},
    ]
    metadata = [
        {
            "sample_id": "S2",
            "source_alignment": "/path/to/S2.bam",
            "source_sha256": "jkl000",
        },
    ]
    with pytest.raises(PacusageError, match="Metadata row names sample S2"):
        record_input_checksums(base, samples, metadata)


def test_record_input_checksums_duplicate_sample_in_metadata() -> None:
    """Test that duplicate sample_id in metadata raises PacusageError."""
    base = [
        {"role": "sample_sheet", "path": "/path", "sha256": "abc123"},
        {"role": "fasta", "path": "/path", "sha256": "def456"},
        {"role": "annotation", "path": "/path", "sha256": "ghi789"},
    ]
    samples = [
        {"sample_id": "S1", "alignment": "/path/to/S1.bam"},
    ]
    metadata = [
        {
            "sample_id": "S1",
            "source_alignment": "/path/to/S1.bam",
            "source_sha256": "jkl000",
        },
        {
            "sample_id": "S1",
            "source_alignment": "/path/to/S1.bam",
            "source_sha256": "mno111",
        },
    ]
    with pytest.raises(PacusageError, match="Duplicate metadata entries for sample S1"):
        record_input_checksums(base, samples, metadata)


def test_record_input_checksums_metadata_file_with_multiple_rows(tmp_path) -> None:
    """Test that metadata file with multiple rows raises PacusageError."""
    # Write base checksums TSV
    base_file = tmp_path / "base_checksums.tsv"
    base_file.write_text(
        "role\tpath\tsha256\n"
        "sample_sheet\t/path/to/samples.tsv\taaa111aaa111aaa111aaa111aaa111aaa111aaa111aaa111aaa111aaa111aaa1\n"
        "fasta\t/path/to/genome.fa\tbbb222bbb222bbb222bbb222bbb222bbb222bbb222bbb222bbb222bbb222bbb2\n"
        "annotation\t/path/to/genes.gtf\tccc333ccc333ccc333ccc333ccc333ccc333ccc333ccc333ccc333ccc333ccc3\n"
    )

    # Write samples TSV
    samples_file = tmp_path / "samples.tsv"
    samples_file.write_text(
        "sample_id\talignment\tcondition\tcontrol\tcontrol_condition\treplicate\tbatch\tdonor\tlayout\tstrandedness\tlibrary_profile\tevidence_source\n"
        "S1\t/path/to/S1.bam\ttreatment\tcontrol\tcontrol\t1\t1\tdonor1\tSE\tunstranded\tprofile1\tauto\n"
    )

    # Write metadata with two data rows (bad)
    metadata_file = tmp_path / "metadata.tsv"
    metadata_file.write_text(
        "source_alignment\tprepared_alignment\tprepared_index\tformat\taction\tlayout_detected\tsource_sha256\tprepared_sha256\tsample_id\n"
        "/path/to/S1.bam\t/cache/S1.sorted.bam\t/cache/S1.sorted.bam.bai\tBAM\tsorted\tSE\taaa\tbbb\tS1\n"
        "/path/to/S2.bam\t/cache/S2.sorted.bam\t/cache/S2.sorted.bam.bai\tBAM\tsorted\tSE\tccc\tddd\tS2\n"
    )

    output_file = tmp_path / "output.tsv"

    with pytest.raises(SystemExit) as exc_info:
        main(
            [
                "record-input-checksums",
                "--base", str(base_file),
                "--samples", str(samples_file),
                "--alignment-metadata", str(metadata_file),
                "--output", str(output_file),
            ]
        )

    assert exc_info.value.code == 2


def test_record_input_checksums_metadata_empty_sample_id(tmp_path) -> None:
    """Test that metadata with empty sample_id raises PacusageError."""
    base = [
        {"role": "sample_sheet", "path": "/path", "sha256": "abc123"},
        {"role": "fasta", "path": "/path", "sha256": "def456"},
        {"role": "annotation", "path": "/path", "sha256": "ghi789"},
    ]
    samples = [
        {"sample_id": "S1", "alignment": "/path/to/S1.bam"},
    ]
    metadata = [
        {
            "sample_id": "",
            "source_alignment": "/path/to/S1.bam",
            "source_sha256": "jkl000",
        },
    ]
    with pytest.raises(PacusageError, match="empty sample_id"):
        record_input_checksums(base, samples, metadata)


def test_record_input_checksums_metadata_empty_source_alignment() -> None:
    """Test that metadata with empty source_alignment raises PacusageError."""
    base = [
        {"role": "sample_sheet", "path": "/path", "sha256": "abc123"},
        {"role": "fasta", "path": "/path", "sha256": "def456"},
        {"role": "annotation", "path": "/path", "sha256": "ghi789"},
    ]
    samples = [
        {"sample_id": "S1", "alignment": "/path/to/S1.bam"},
    ]
    metadata = [
        {
            "sample_id": "S1",
            "source_alignment": "",
            "source_sha256": "jkl000",
        },
    ]
    with pytest.raises(PacusageError, match="empty source_alignment"):
        record_input_checksums(base, samples, metadata)


def test_record_input_checksums_metadata_empty_source_sha256() -> None:
    """Test that metadata with empty source_sha256 raises PacusageError."""
    base = [
        {"role": "sample_sheet", "path": "/path", "sha256": "abc123"},
        {"role": "fasta", "path": "/path", "sha256": "def456"},
        {"role": "annotation", "path": "/path", "sha256": "ghi789"},
    ]
    samples = [
        {"sample_id": "S1", "alignment": "/path/to/S1.bam"},
    ]
    metadata = [
        {
            "sample_id": "S1",
            "source_alignment": "/path/to/S1.bam",
            "source_sha256": "",
        },
    ]
    with pytest.raises(PacusageError, match="empty source_sha256"):
        record_input_checksums(base, samples, metadata)
