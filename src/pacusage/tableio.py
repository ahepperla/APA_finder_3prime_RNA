"""Small deterministic table and checksum helpers."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any


def open_text(path: str | Path, mode: str = "rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(
    rows: Iterable[dict[str, Any]], path: str | Path, fieldnames: list[str] | None = None
) -> None:
    materialized = list(rows)
    if fieldnames is None:
        fieldnames = list(materialized[0]) if materialized else []
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
            writer.writerows(materialized)


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
