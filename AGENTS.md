# PACusage Repository Guide

## Codebase Discovery

Use Codebase Memory MCP for cross-file structure, callers, and impact analysis
before changing shared pipeline or Python modules. Prefer `search_graph`,
`trace_path`, and `get_code_snippet` for structural discovery. Use `rg` for
literal strings, configuration values, and non-code files. Check index
coverage before making broad or negative structural claims.

## Architecture And Conventions

- PACusage is a Nextflow DSL2 pipeline with a Python library and CLI in
  `src/pacusage`. The installed `pacusage` command resolves to
  `pacusage.cli:main`.
- `main.nf` is the workflow entry point. `workflows/pacusage.nf` orchestrates
  `VALIDATE_INPUTS`, `PREPARATION`, `DISCOVERY`, `QUANTIFICATION`,
  `STATISTICS`, and `BUILD_REPORT`, in that order.
- Preserve the analytical separation between protocol calibration and atlas
  construction. Discovery uses a condition-blind frozen atlas; differential
  testing compares each treatment only with its declared direct control.
- Keep workflow modules small and use their declared channels rather than
  bypassing stages with ad hoc files. Explicit Nextflow CLI parameters override
  values supplied by `analysis.yaml`.
- Do not modify source FASTA or alignment inputs. The pipeline records
  checksums, resolved parameters, software versions, preparation actions, and
  the frozen-atlas checksum for reproducibility.

## Input And Data Invariants

- `analysis.yaml` requires `input`, `assembly`, `fasta`, and `gtf`; `outdir`
  should be supplied for normal analysis runs.
- Sample sheets are TSV files with required `sample_id`, `alignment`,
  `condition`, and `control` columns. Alignment paths are relative to the
  sample sheet. A root control has a genuinely empty `control` field, never
  the literal value `NA`.
- Each condition has one consistent direct control. Control relationships must
  be acyclic, and modeled conditions normally need at least two biological
  replicates.
- PAC IDs and internal coordinates are zero-based interbase coordinates.
  Exact-boundary BED records are `[coordinate, coordinate + 1)`; proximal-tag
  BED records are `[region_start, region_end)`.
- Keep exact-boundary and proximal-tag discovery semantics distinct. Do not
  represent a proximal-tag estimate as nucleotide-resolution cleavage evidence.
- Counts are raw assigned fragment counts. PAU is each PAC count divided by
  total assigned PAC counts for its gene and sample; do not add pseudocounts.

## Outputs And Generated State

- Normal results contain `manifest/`, `qc/`, `evidence/`, `atlas/`, `counts/`,
  `statistics/`, `motifs/`, `tracks/`, and `report/index.html`.
- `work/`, `.nextflow/`, `results-test/`, caches, generated reports, and
  Nextflow logs are generated state. Do not treat them as maintained source.
- Retain Nextflow work directories during an analysis so `-resume` remains
  reliable. Prepared references, alignments, and per-sample count tables stay
  in the cache unless their corresponding `save_*` option is enabled.

## Validation

Use the smallest relevant check first, then run broader validation for
pipeline, schema, data-contract, or statistics changes:

```bash
pytest
ruff check src tests
```

```bash
tests/run_nextflow.sh
```

The integration script runs the R bootstrap test, rebuilds fixtures, executes
the test Nextflow profile, checks required result artifacts, and validates
their contents. When running the fixture manually with Conda:

```bash
python tests/fixtures/build_fixture.py
nextflow run . -profile test,conda -resume
```
