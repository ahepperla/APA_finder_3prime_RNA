"""Small deterministic table and checksum helpers."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
from collections.abc import Iterable
from itertools import chain
from pathlib import Path
from typing import Any


def open_text(
    path: str | Path,
    mode: str = "rt",
    compresslevel: int | None = None,
):
    path = Path(path)
    if path.suffix == ".gz":
        # Zero timestamp makes identical content produce identical bytes,
        # and the atlas checksum seeds the statistics.
        if any(op in mode for op in ("w", "a", "x")):
            binary = gzip.GzipFile(
                filename=path,
                mode=mode.replace("t", ""),
                compresslevel=9 if compresslevel is None else compresslevel,
                mtime=0,
            )
            if "b" in mode:
                return binary
            return io.TextIOWrapper(binary, newline="")
        # Reading: unchanged behavior.
        options = {"newline": ""}
        if compresslevel is not None:
            options["compresslevel"] = compresslevel
        return gzip.open(path, mode, **options)
    return path.open(mode, newline="")


def gzip_compression(path: str | Path) -> str | dict[str, object]:
    """Return compression argument for pandas to_csv with zero gzip timestamp."""
    if str(path).endswith(".gz"):
        return {"method": "gzip", "mtime": 0}
    return "infer"


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    return list(iter_tsv(path))


def iter_tsv(path: str | Path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(
    rows: Iterable[dict[str, Any]],
    path: str | Path,
    fieldnames: list[str] | None = None,
    compresslevel: int | None = None,
) -> None:
    iterator = iter(rows)
    if fieldnames is None:
        first = next(iterator, None)
        fieldnames = list(first) if first is not None else []
        if first is not None:
            iterator = chain([first], iterator)
    with open_text(path, "wt", compresslevel=compresslevel) as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        if fieldnames:
            writer.writeheader()
            writer.writerows(iterator)


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(value: Any, path: str | Path) -> None:
    with Path(path).open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
