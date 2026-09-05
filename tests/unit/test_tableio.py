from pacusage.tableio import iter_tsv, write_tsv


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
