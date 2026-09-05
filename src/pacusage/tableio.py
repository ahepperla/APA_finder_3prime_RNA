"""Small deterministic table and checksum helpers."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections.abc import Iterable
from itertools import chain
from pathlib import Path
from typing import Any


def open_text(path: str | Path, mode: str = "rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    return list(iter_tsv(path))


def iter_tsv(path: str | Path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(
    rows: Iterable[dict[str, Any]], path: str | Path, fieldnames: list[str] | None = None
) -> None:
    iterator = iter(rows)
    if fieldnames is None:
        first = next(iterator, None)
        fieldnames = list(first) if first is not None else []
        if first is not None:
            iterator = chain([first], iterator)
    with open_text(path, "wt") as handle:
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
