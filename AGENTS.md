# PACusage Repository Guide

## Codebase Discovery

Use Codebase Memory MCP for cross-file structure, callers, and impact analysis
before changing shared pipeline or Python modules. Prefer `search_graph`,
`trace_path`, and `get_code_snippet` for structural discovery. Use `rg` for
literal strings, configuration values, and non-code files. Check index
coverage before making broad or negative structural claims.

## Architecture And Conventions

- PACusage is a Nextflow DSL2 pipeline with a Python package and CLI in
  `src/pacusage`, the statistics in `scripts/fit_usage_model.R`, and the
  figures in `scripts/plot_usage_figures.R`. Nextflow
  tasks run `bin/pacusage`, which uses the installed package (as in the
  Apptainer image) or else this checkout's `src/`.
- `main.nf` is the workflow entry point. `workflows/pacusage.nf` orchestrates
  `VALIDATE_INPUTS`, `PREPARATION`, `DISCOVERY`, `QUANTIFICATION`,
  `STATISTICS`, `PLOT_FIGURES`, and `BUILD_REPORT`, in that order.
- Preserve the analytical separation between protocol calibration and atlas
  construction. The frozen atlas is built without the treatment–control
  contrasts:
  - candidates come from every sample's read ends pooled;
  - a PAC needs replicate support within the condition that supports it;
  - the readthrough filter also works per condition.

  Differential testing compares each treatment only with its declared direct
  control.
- Keep workflow modules small and use their declared channels rather than
  bypassing stages with ad hoc files. Explicit Nextflow CLI parameters override
  values supplied by `analysis.yaml`.
- Parameters are defined once, in `nextflow_schema.json`: types, ranges,
  defaults, and descriptions.
  - `nextflow.config` repeats the defaults, and `tests/unit/test_config.py`
    keeps the two in step.
  - VALIDATE_INPUTS resolves and validates the parameters once, against the
    staged schema.
  - Later steps read `resolved_params.yaml` as it is. The exception is
    INFER_STRANDEDNESS, which takes its few parameters, including
    `excluded_contigs` and the aliases, on the command line so that other
    parameter changes don't rerun it. Nextflow applies the
    resource, scheduling, and publication parameters itself.
  - RECORD_SOFTWARE_VERSIONS loads the R packages of both R scripts right
    after validation. Keep the R scripts out of VALIDATE_INPUTS, or every edit
    to one reruns the whole pipeline on `-resume`.
- Do not modify source FASTA or alignment inputs.
  - Prepared alignments, indexes, and the prepared FASTA are symlinks to the
    sources. Never write through them.
  - Keep any generated index in the work directory. htslib writes a missing
    `.fai` beside whatever FASTA path it is given, so pass it a link in the
    task directory, never the source path.
- The pipeline records checksums, resolved parameters, software versions,
  preparation actions, and the frozen-atlas checksum for reproducibility.
  Published gzip files and PDFs carry no timestamps, so reruns are
  byte-identical.
  The atlas checksum seeds the statistics, so any change to the atlas file's
  contents or columns changes the seeded results.
- Each alignment is read in full twice: once by INFER_STRANDEDNESS, and once
  by SCAN_ALIGNMENT (`src/pacusage/scan.py`) for every candidate evidence
  source together. EXTRACT_3PRIME_EVIDENCE writes the chosen source's
  evidence from the scan.
- Statistical filtering belongs to `fit_usage_model.R`: it filters each
  comparison family's own samples and writes the reasons, per PAC, to
  `FAMILY.statistical_filtering.tsv.gz`.
- Calls belong to `fit_usage_model.R` too. PLOT_FIGURES draws the calls in
  the `.pacs` and `.genes` tables and never makes its own.

## Input And Data Invariants

- `analysis.yaml` requires `input`, `assembly`, `fasta`, and `gtf`; `outdir`
  should be supplied for normal analysis runs. The FASTA must be uncompressed.
- Sample sheets are TSV files with required `sample_id`, `alignment`,
  `condition`, and `control` columns. Alignment paths are relative to the
  sample sheet. A root control has a genuinely empty `control` field, never
  the literal value `NA`.
- Each condition has one consistent direct control. Control relationships must
  be acyclic, and modeled conditions normally need at least two biological
  replicates.
- PAC IDs and internal coordinates are zero-based interbase coordinates. A
  PAC's `start`/`end` (and its atlas BED record) are `[coordinate,
  coordinate + 1)` for an exact-boundary PAC and its resolution region for a
  proximal-tag PAC; `locus` is the same interval, 1-based. Gene assignment treats an
  exon or gene as containing coordinates `start` through `end` inclusive,
  because a 3′ boundary equals its exon's `end`.
- Keep exact-boundary and proximal-tag discovery semantics distinct. Do not
  represent a proximal-tag estimate as nucleotide-resolution cleavage evidence.
- Counts are raw assigned fragment counts. PAU is each PAC count divided by
  total assigned PAC counts for its gene and sample; do not add pseudocounts.

## Outputs And Generated State

- Normal results contain `manifest/`, `qc/`, `evidence/`, `atlas/`, `counts/`,
  `statistics/`, `motifs/`, `tracks/`, `figures/`, and `report/index.html`.
- Every PAC-level table starts with the identity block `pac_id, gene_id,
  gene_name, chrom, start, end, strand, locus`. It is `IDENTITY_COLUMNS` in
  `src/pacusage/models.py` and in the R script, and those names are reserved as
  sample IDs. Gene-level tables start with `gene_id, gene_name`.
  - A gene's name comes from `reference.gene_name_map`, falling back to its
    ID.
  - `feature_id` is DRIMSeq's internal name for a PAC and is never published.
- `excluded_contigs` defaults to the mitochondrial names `chrM`, `MT`, and
  `chrMT`, and every step, strandedness inference included, drops them.
- `pipeline_info/` holds Nextflow's execution report, timeline, trace, and
  DAG. `nextflow.config` sets them after the profiles, so a profile's `outdir`
  applies. Each run replaces them, and they record the run's times, so the
  rerun comparison skips them. Only `.nextflow.log`, `.nextflow/`, and the
  default `work/` go to the launch directory.
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
Rscript tests/r/test_usage_model.R scripts/fit_usage_model.R
Rscript tests/r/test_usage_model_simulation.R scripts/fit_usage_model.R
Rscript tests/r/test_usage_figures.R scripts/plot_usage_figures.R scripts/fit_usage_model.R
tests/pipeline/run_nextflow.sh
```

The integration script:
- runs the R statistics and figure tests and rebuilds the fixtures;
- executes the test Nextflow profile, checks the required result artifacts,
  and validates their contents (`tests/pipeline/verify_integration.py`);
- checks that the fixture files are unchanged;
- checks that a fresh rerun reproduces every published file outside
  `pipeline_info/`;
- checks that a mixed-protocol sample sheet fails before discovery;
- runs the Plasmidsaurus-like proximal-tag fixture in
  `tests/fixtures/plasmidsaurus/`, which has its own reference and annotation;
- checks that exact-boundary reads run under the Plasmidsaurus profile warn
  at calibration and stop before quantification.

Keep that fixture's calibration genes single-ended. Its multi-PAC genes
annotate one transcript per PAC, so alternative ends stay out of the
calibration kernel.

When running the fixture manually with Conda:

```bash
python tests/fixtures/build_fixture.py
nextflow run . -profile test,conda -resume
```
