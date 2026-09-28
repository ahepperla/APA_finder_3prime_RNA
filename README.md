# PACusage

PACusage is a Nextflow DSL2 pipeline for discovering polyadenylation-site
clusters (PACs) from deduplicated BAM or CRAM files and testing differential
PAC usage. It keeps protocol calibration separate from atlas construction,
uses a condition-blind frozen atlas, quantifies raw fragment counts, and tests
each treatment against its declared control.

The implementation follows the design in
[`pacusage_hpc_implementation_plan.md`](pacusage_hpc_implementation_plan.md).
Version 0.1 supports exact-boundary and calibrated proximal-tag evidence,
single- and paired-end alignments, GTF/GFF3 annotation, known-PAC rescue,
DRIMSeq/stageR testing, motif summaries, browser tracks, and a portable HTML
report.

## Requirements

- Nextflow 24.04 or newer
- Java 17 or newer
- Conda/Mamba for the `conda` profile, or Apptainer for the `apptainer` profile
- A POSIX-like execution environment

Read trimming, alignment, and deduplication happen before PACusage.

## Quick Start

Create a tab-separated sample sheet. Control rows have a genuinely empty final
field, not the word `NA`.

```text
sample_id	alignment	condition	control
DMSO_1	/data/a.bam	DMSO	
DMSO_2	/data/b.bam	DMSO	
TRA_1	/data/c.bam	TreatmentA	DMSO
TRA_2	/data/d.bam	TreatmentA	DMSO
```

Create `analysis.yaml`:

```yaml
input: samples.tsv
outdir: results
assembly: GRCh38
fasta: /reference/GRCh38.fa
gtf: /reference/gencode.gtf
```

Run locally:

```bash
nextflow run /path/to/pacusage \
  -profile local,conda \
  -params-file analysis.yaml \
  -resume
```

Run on Slurm:

```bash
nextflow run /path/to/pacusage \
  -profile slurm,apptainer \
  -params-file analysis.yaml \
  -c institution.config \
  -resume
```

The Slurm profile separates CPU-parallel BAM work from serial, memory-heavy
steps. A job that runs out of memory or time is retried up to twice, with more
memory and more time on each attempt. Completed tasks remain reusable with
`-resume`; the statistics steps rerun when `scripts/fit_usage_model.R` changes.
Each direct comparison family first runs its whole-family DRIMSeq fit, then
scatters selected bootstrap genes into independent batches before a
deterministic family-level merge.

`FIT_USAGE_MODEL` receives eight CPUs by default. Increase only that process
with `--statistics_cpus 16`; those cores parallelize DRIMSeq's family-wide
fits and tests.

Genes with a PAC that has no reads in some condition cannot be fitted as they
are. For those genes, the family is refitted `dm_zero_sensitivity_repeats`
times (5 by default). Each refit replaces their zero counts with small seeded
values between 0 and 0.1, following DRIMSeq's `addUniform` rule. The reported
result is the median of the refits. PACs whose effect changes direction or
spreads by more than `dm_zero_max_delta_pau_spread`, or that fewer than 80%
of the refits can fit, are flagged `zero_boundary_unstable` and left untested. These perturbed values are used
only for fitting and never appear in output tables. A condition with no reads
at all for a gene leaves that gene untested in the comparisons that use the
condition (`model_status` is `group_without_counts`).

Bootstrap intervals hold each gene's precision at the family-fit estimate, so
every replicate needs one refit of the proportions, or two when a simulated
zero group needs the same perturbation. One draw gives every
comparison in the family its interval. In simulations with four replicates per
group, nominal 95% intervals covered the true change in PAU 80-93% of the
time (about 88% on average), because they do not include uncertainty in the
precision.
`bootstrap_status` records why an interval is missing.

Bootstrap work uses separate four-core jobs by default; tune their allocation
with `--statistics_bootstrap_cpus 4` and genes per job with
`--statistics_bootstrap_batch_size 500`. At most eight batch jobs are submitted
at once by default; tune that cap with `--statistics_bootstrap_max_forks 8`.
Each statistics attempt requests 12 hours of walltime and each bootstrap
attempt 6 hours, multiplied by the attempt number on retry. To skip intervals
in an exploratory run, set `--dm_bootstrap_replicates 0`. Setting
`--dm_bootstrap_include_candidates false` leaves p-values and event calls
unchanged, but omits bootstrap intervals for non-significant candidate genes.

stageR's stage-wise procedure confirms both PACs of a two-PAC gene together.
When such a gene passes the screen, both PACs report a `pac_fdr` of 0.

Build the default Apptainer image once on a networked system:

```bash
apptainer build --fakeroot containers/pacusage.sif containers/Apptainer.def
```

The build host does not need Conda, Mamba, or micromamba. The Miniforge base
image supplies Conda inside the build. Apptainer, network access, and either
fakeroot support or privileged build access are required. Omit `--fakeroot`
when using a privileged build service.

The SIF and all analysis inputs can then be moved to an offline cluster. Set
`container` in YAML when the image is stored elsewhere.

Some HPC installations do not automatically expose shared filesystems inside
Apptainer. Add their host roots to `analysis.yaml`:

```yaml
container: /project/software/PACusage/containers/pacusage.sif
bind_paths:
  - /vast
```

PACusage passes each entry as an Apptainer bind mount with the same host and
container path. Nextflow work directories are mounted separately.

Prepared alignments, their indexes, and the prepared FASTA are symbolic links
to your source files, not copies. Only an unsorted alignment is sorted into a
new file, and a missing index is generated in the work directory. Therefore:
- keep the source files in place and unmodified until the analysis is
  accepted;
- with Apptainer, make sure `bind_paths` covers their directories.

Each step that reads a linked alignment first checks that the file's size and
modification time still match what `PREPARE_ALIGNMENT` recorded. If they
don't, the step fails and asks for a rerun with `-resume`. The
`save_prepared_alignments` and `save_prepared_reference` options still
publish real copies.

Explicit command-line values override YAML values:

```bash
nextflow run /path/to/pacusage \
  -profile local,conda \
  -params-file analysis.yaml \
  --min_mapq 30 \
  -resume
```

Use `nextflow run /path/to/pacusage --help` for the required inputs and
[`nextflow_schema.json`](nextflow_schema.json) for all parameters and defaults.

## Library Profiles

`generic_3prime` is the default. It tests compatible evidence sources and
requires a reproducible, clearly preferred calibration model.

`plasmidsaurus_3prime` defaults to single-end, forward-stranded `read_3p`
evidence and a `proximal_tag` endpoint model. The BAM evidence must remain
compatible with those defaults.

`exact_boundary` is for assays preserving the transcript/poly(A) junction. It
prefers `polyA_junction` evidence when enough genes support it and otherwise
uses the profile-compatible aligned edge.

Calibration is parallelized by sample. Nextflow first builds one compact table
of annotated transcript ends, then submits one `SCAN_ALIGNMENT` task per
BAM/CRAM, and finally combines the small per-sample summaries into the
run-level calibration files. On Slurm, a 48-sample run can therefore schedule
up to 48 independent jobs instead of scanning all alignments in one
large-memory job. Each scan reads its alignment once, name-collating a
paired-end file once, for every candidate evidence source together:
- calibration summaries;
- aggregated end observations;
- splice continuations;
- filtering counts.

After calibration picks the run's evidence source, `EXTRACT_3PRIME_EVIDENCE`
writes that source's evidence from the scan without reading the alignment
again.

Calibration also checks the pooled offset kernel and writes the result to
`qc/calibration_kernel_diagnostics.tsv`. Read ends scattered around one site
give a kernel with one mode. A proximal-tag run warns when:
- the kernel, smoothed over its own minimum resolvable separation, has more
  than one separated mode; or
- the samples' median central interval is more than 6 times that separation.

A kernel of PAC spacings looks like this, for example when exact-boundary
reads run under a proximal-tag profile, or when unannotated alternative
polyadenylation sites lie near the annotated ends that calibration uses. The
warning appears in the Nextflow log and at the top of the report, and the run
goes on. Check that the library profile suits the data before trusting its
proximal-tag PACs. Exact-boundary runs record the same values with status
`not_applicable`, because exact discovery does not assign reads through the
kernel.

When multiple library chemistries resolve differently, run them separately.
A batch term cannot recover information lost through incompatible endpoint
definitions.

## Sample Sheet

Required columns are `sample_id`, `alignment`, `condition`, and `control`.
Paths are resolved relative to the sample sheet.

Optional columns include `replicate`, `batch`, `donor`, `layout`,
`strandedness`, `library_profile`, and `evidence_source`. Additional columns
can be included in the model by naming them in `model_covariates`.

Every condition must have one consistent direct control value. A blank control
marks a root control condition. Nested comparisons are supported: a condition
may be tested against its own control and also serve as the control for another
condition, such as `WT -> disease_vehicle -> disease_drug`. Control
relationships must be acyclic. By default every modeled condition needs two
biological replicates.

## Coordinates and Interpretation

PACusage uses zero-based interbase coordinates internally and in PAC IDs.
Exact-boundary BED entries represent `[coordinate, coordinate + 1)`;
proximal-tag BED entries span `[region_start, region_end)`.

An **end observation** is a transcript-oriented aligned boundary. It is called
an exact PAC only after the selected protocol and calibration support
nucleotide resolution. Proximal-tag libraries produce estimated PAC
coordinates and resolution groups.

The two discovery modes remain separate:

- `exact_boundary` clusters observed cleavage boundaries directly and reports
  nucleotide-scale representatives.
- `proximal_tag` bins endpoints (25 nt by default), convolves each
  chromosome/strand signal with the calibrated offset kernel, finds regional
  score peaks, and merges peaks closer than the calibrated minimum resolvable
  separation. The atlas reports `region_start`, `region_end`, and
  `resolution_nt`; its representative coordinate is not a claimed
  single-nucleotide cleavage site.

Discovery and quantification stream the per-sample Parquet evidence files.
This avoids expanding every endpoint over every kernel offset and keeps
large, many-sample runs bounded by one chromosome/strand group at a time.

Candidate support is replicate-coherent: `pac_min_supporting_samples` must be
met by samples from at least one condition. Evidence from unrelated conditions
cannot be combined merely to pass the discovery threshold. The atlas remains
condition-blind otherwise: direction and treatment effect are not used during
discovery, and the final atlas is the union of candidates supported by any
condition. Whole-number values such as `2` require that many samples. Values
strictly between zero and one are fractions of the samples in a condition,
rounded up; for example, `0.5` requires three of five replicates.

For `proximal_tag` libraries, PACusage additionally guards against ordinary
aligned endpoints in splice-continued exons. It derives immediate
exon-to-next-exon edges from full-read CIGAR `N` operations, then applies the
post-discovery `constitutive_readthrough` filter. A candidate is rejected only
when its representative coordinate is in an upstream CIGAR exon block and
direct downstream continuation meets the threshold independently in every
condition. This does not use PAS annotation, endpoint prominence, or a
treatment effect to create candidates.

The default is intentionally strict:
`constitutive_readthrough_min_junction_count: 2` requires two direct
continuation reads in a qualifying sample, and
`constitutive_readthrough_min_replicate_support: all` requires every replicate
within every condition. This avoids overloading `1`: as with
`pac_min_supporting_samples`, a numeric `1` means one replicate. The setting
can instead be a fraction in `(0, 1)`, rounded up within each condition, or an
absolute count. A condition that
does not meet this continuation requirement protects its candidate, allowing a
condition-specific terminal exon to remain eligible for downstream testing.
The per-sample continuation tables are written to `evidence/` for audit.

PAU is the raw PAC count divided by all assigned PAC counts for that gene in
one sample. No pseudocount is added to count or observed-PAU matrices.

## Outputs

The main output directories are:

- `manifest/`: resolved parameters, normalized samples, checksums, versions
- `qc/`: validation, preparation, calibration, filtering, and conservation
- `evidence/`: aggregate sample-level end observations
- `atlas/`: versioned PAC BED and metadata
- `counts/`: raw PAC counts, gene totals, and observed PAU
- `statistics/`: family omnibus and treatment-versus-control tests and events
- `motifs/`: motif preference results
- `tracks/`: separate non-negative plus/minus bedGraph files
- `report/index.html`: self-contained report

Motif outputs retain the broad `primary_motif_class` and the exact matched PAS
in both genomic DNA (`primary_pas_motif`) and RNA (`primary_pas_motif_rna`)
notation.

Prepared references, alignments, and per-sample count tables stay in the
Nextflow cache unless their `save_*` options are enabled.

## Development

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[dev]'
pytest
```

The statistics tests need the pipeline's R packages (DRIMSeq, stageR, limma,
and yaml). Point `R_LIBS` at a local library if they are not installed
globally:

```bash
Rscript tests/test_usage_model.R scripts/fit_usage_model.R
Rscript tests/test_usage_model_simulation.R scripts/fit_usage_model.R
```

The integration fixture can be generated with:

```bash
python tests/fixtures/build_fixture.py
nextflow run . -profile test,conda -resume
```

`tests/run_nextflow.sh` runs the R tests and the fixture pipeline, then
checks the results. It also checks that:
- the fixture files are unchanged afterwards;
- a fresh run in a new work directory reproduces every published file;
- a sample sheet mixing two protocols fails at calibration, before discovery;
- the Plasmidsaurus-like fixture runs end to end and passes
  `tests/verify_plasmidsaurus.py`;
- exact-boundary reads run under the Plasmidsaurus profile warn at
  calibration, then stop at ANNOTATE_PACS.

`build_fixture.py` also writes that Plasmidsaurus-like fixture, to
`tests/fixtures/plasmidsaurus/`. It has its own reference and annotation, and
its reads end 20-280 nt upstream of their PACs, so it exercises proximal-tag
calibration, discovery, and quantification. To run it by hand, pass absolute
paths:

```bash
nextflow run . -profile test,conda \
  --input "$PWD/tests/fixtures/plasmidsaurus/samples.tsv" \
  --fasta "$PWD/tests/fixtures/plasmidsaurus/genome.fa" \
  --gtf "$PWD/tests/fixtures/plasmidsaurus/genes.gtf" \
  --endpoint_model auto --outdir results-plasmidsaurus
```

If proximal-tag discovery accepts no PAC at all, the run stops at
ANNOTATE_PACS. The error names `qc/pac_discovery.tsv` and
`atlas/rejected_candidates.tsv.gz`, which explain why. A calibration kernel
warning earlier in the log often points to the cause.

The statistical process requires DRIMSeq, stageR, and limma. The supplied
Conda files install them from Bioconda and pin DRIMSeq to 1.38.0, the version
the statistics were validated with. Runs are network-independent after the
environment or container has been prepared.

## Reproducibility

PACusage never modifies source FASTA or alignment files, and never writes
beside them. It records all of the following:
- input checksums, with each alignment hashed once by its own
  `PREPARE_ALIGNMENT` task;
- resolved parameters;
- software versions;
- preparation actions;
- a checksum of the frozen atlas.

Gzip files are written without timestamps. So a rerun with the same inputs,
parameters, and software reproduces every published file byte for byte,
including the atlas checksum that seeds the statistics.

Nextflow work directories provide resumability; retain them until the
analysis is accepted.
