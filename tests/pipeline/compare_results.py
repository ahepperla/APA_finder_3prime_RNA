#!/usr/bin/env python3
"""Check that two PACusage result directories hold byte-identical files.

tests/pipeline/run_nextflow.sh uses this to show that a fresh run reproduces an earlier
one. Only the output directory may differ: it appears in
manifest/resolved_params.yaml and in that file's checksum in
manifest/run_manifest.json. Nextflow's own reports in pipeline_info/ record
the run's times, so they are not compared.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def comparable(path: Path, name: str) -> bytes:
    data = path.read_bytes()
    if name == "manifest/resolved_params.yaml":
        lines = data.decode().splitlines(keepends=True)
        return "".join(line for line in lines if not line.startswith("outdir:")).encode()
    if name == "manifest/run_manifest.json":
        manifest = json.loads(data)
        manifest.pop("parameters_sha256", None)
        return json.dumps(manifest, sort_keys=True).encode()
    return data


def main(left: str, right: str) -> int:
    roots = (Path(left), Path(right))
    names = [
        {
            path.relative_to(root).as_posix()
            for path in root.rglob("*")
            if path.is_file() and path.relative_to(root).parts[0] != "pipeline_info"
        }
        for root in roots
    ]
    problems = [f"only in {roots[0]}: {name}" for name in sorted(names[0] - names[1])]
    problems += [f"only in {roots[1]}: {name}" for name in sorted(names[1] - names[0])]
    shared = sorted(names[0] & names[1])
    problems += [
        f"differs: {name}"
        for name in shared
        if comparable(roots[0] / name, name) != comparable(roots[1] / name, name)
    ]
    for problem in problems:
        print(problem, file=sys.stderr)
    print(f"Compared {len(shared)} files: {len(problems)} differences.")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main(*sys.argv[1:3]))
