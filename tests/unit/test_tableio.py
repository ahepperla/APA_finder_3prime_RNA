import gzip
import time

import pandas as pd

from pacusage.cli import _write_bed
from pacusage.evidence import write_bedgraphs
from pacusage.models import EvidenceObservation
from pacusage.tableio import gzip_compression, iter_tsv, write_tsv


def test_tsv_round_trip_accepts_a_generator(tmp_path) -> None:
    path = tmp_path / "rows.tsv.gz"
    write_tsv(
        ({"name": f"row-{index}", "value": index} for index in range(3)),
        path,
        ["name", "value"],
    )
    assert list(iter_tsv(path)) == [
        {"name": "row-0", "value": "0"},
        {"name": "row-1", "value": "1"},
        {"name": "row-2", "value": "2"},
    ]


def test_gzip_tables_have_zero_timestamp(tmp_path) -> None:
    path = tmp_path / "rows.tsv.gz"
    write_tsv([{"col": "value"}], path, ["col"])
    data = path.read_bytes()
    assert data[:3] == b"\x1f\x8b\x08"
    assert data[4:8] == b"\x00\x00\x00\x00"


def test_gzip_tables_are_byte_identical_across_times(tmp_path, monkeypatch) -> None:
    path_a = tmp_path / "a" / "rows.tsv.gz"
    path_a.parent.mkdir(parents=True)
    path_b = tmp_path / "b" / "rows.tsv.gz"
    path_b.parent.mkdir(parents=True)

    monkeypatch.setattr(time, "time", lambda: 1_000_000_000.0)
    write_tsv([{"col": "val1"}, {"col": "val2"}], path_a, ["col"])

    monkeypatch.setattr(time, "time", lambda: 1_500_000_000.0)
    write_tsv([{"col": "val1"}, {"col": "val2"}], path_b, ["col"])

    data_a = path_a.read_bytes()
    data_b = path_b.read_bytes()
    assert data_a == data_b
    assert gzip.decompress(data_a) == b"col\nval1\nval2\n"


def test_bedgraphs_have_zero_timestamp(tmp_path) -> None:
    plus_path = tmp_path / "plus.bedGraph.gz"
    minus_path = tmp_path / "minus.bedGraph.gz"

    observations = [
        EvidenceObservation(
            sample_id="sample1",
            contig="chr1",
            strand="+",
            coordinate=100,
            count=5,
        ),
        EvidenceObservation(
            sample_id="sample1",
            contig="chr1",
            strand="-",
            coordinate=200,
            count=3,
        ),
        EvidenceObservation(
            sample_id="sample1",
            contig="chr1",
            strand="+",
            coordinate=150,
            count=7,
        ),
    ]

    write_bedgraphs(observations, plus_path, minus_path)

    plus_bytes = plus_path.read_bytes()
    minus_bytes = minus_path.read_bytes()

    assert plus_bytes[:3] == b"\x1f\x8b\x08"
    assert plus_bytes[4:8] == b"\x00\x00\x00\x00"
    assert minus_bytes[:3] == b"\x1f\x8b\x08"
    assert minus_bytes[4:8] == b"\x00\x00\x00\x00"

    plus_decompressed = gzip.decompress(plus_bytes).decode()
    minus_decompressed = gzip.decompress(minus_bytes).decode()

    assert plus_decompressed == "chr1\t100\t101\t5\nchr1\t150\t151\t7\n"
    assert minus_decompressed == "chr1\t200\t201\t3\n"


def test_bed_has_zero_timestamp(tmp_path) -> None:
    path = tmp_path / "test.bed.gz"

    rows = [
        {
            "contig": "chr1",
            "coordinate": 100,
            "endpoint_model": "exact_boundary",
            "region_start": 0,
            "region_end": 0,
            "pac_id": "pac1",
            "total_count": 5,
            "strand": "+",
        },
        {
            "contig": "chr2",
            "coordinate": 200,
            "endpoint_model": "proximal_tag",
            "region_start": 195,
            "region_end": 205,
            "pac_id": "pac2",
            "total_count": 10,
            "strand": "-",
        },
    ]

    _write_bed(rows, path)

    data = path.read_bytes()
    assert data[:3] == b"\x1f\x8b\x08"
    assert data[4:8] == b"\x00\x00\x00\x00"

    decompressed = gzip.decompress(data).decode()
    assert decompressed == "chr1\t100\t101\tpac1\t5\t+\nchr2\t195\t205\tpac2\t10\t-\n"


def test_pandas_gzip_compression_has_zero_timestamp(tmp_path) -> None:
    gz_path = tmp_path / "x.tsv.gz"
    tsv_path = tmp_path / "x.tsv"

    df = pd.DataFrame({"a": [1, 2]})

    df.to_csv(gz_path, sep="\t", index=False, compression=gzip_compression(gz_path))
    gz_data = gz_path.read_bytes()
    assert gz_data[:3] == b"\x1f\x8b\x08"
    assert gz_data[4:8] == b"\x00\x00\x00\x00"
    assert gzip.decompress(gz_data) == b"a\n1\n2\n"

    df.to_csv(tsv_path, sep="\t", index=False, compression=gzip_compression(tsv_path))
    tsv_data = tsv_path.read_bytes()
    assert tsv_data[:4] != b"\x1f\x8b\x08\x08"
    assert tsv_data == b"a\n1\n2\n"
