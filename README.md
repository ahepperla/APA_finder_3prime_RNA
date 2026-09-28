# PACusage

PACusage is a Nextflow pipeline for 3′-end RNA-seq. It finds polyadenylation
site clusters (PACs) in BAM or CRAM alignments, counts the reads at each PAC,
and tests every treatment against its own control for changes in PAC usage:
the share of a gene's reads at each of its PACs.

It keeps a few things strictly separate:

- **Protocol calibration and discovery.** The library's 3′-end behavior is
  calibrated against annotated transcript ends before any PAC is called.
- **Discovery and testing.** One PAC atlas is built from every sample's read
  ends, without the treatment–control contrasts, then frozen, and only then
  counted and tested.
- **Comparisons.** Each treatment is compared only with the control it names
  in the sample sheet.

Counts are raw read counts; PAU (PAC usage) is a PAC's count divided by its
gene's total, with no pseudocounts. Differential usage is modeled with DRIMSeq
and adjusted with stageR.

The design, including every statistical choice, is in
[docs/design.md](docs/design.md). Decisions made while building it are in
[docs/decisions.md](docs/decisions.md).

## Contents

- [Requirements](#requirements)
- [Set up on the cluster](#set-up-on-the-cluster)
- [Inputs](#inputs)
- [Run](#run)
- [Outputs](#outputs)
- [How it works](#how-it-works)
- [Statistics](#statistics)
- [Reproducibility and -resume](#reproducibility-and--resume)
- [Troubleshooting](#troubleshooting)
- [Development](#development)

## Requirements

- Nextflow 24.04 or newer, with Java 17 or newer. CI runs the tests with
  Nextflow 25.04.7.
- Apptainer, for Slurm runs. Conda works too, mainly for local runs.
- Deduplicated BAM or CRAM files. Trimming, alignment, and duplicate marking
  happen before PACusage.
- An uncompressed genome FASTA, and a GTF or GFF3 annotation for the same
  assembly.

## Set up on the cluster

Clone the pipeline and build its image once, on a machine with network access:

```bash
git clone https://github.com/ahepperla/APA_finder_3prime_RNA.git pacusage
cd pacusage
apptainer build --fakeroot containers/pacusage.sif containers/Apptainer.def
```

`--fakeroot` needs user-namespace support; omit it on a build service that
runs as root. The build needs only Apptainer and network access: the image
installs its own Conda environment, from `envs/pacusage.yml`.

The image contains the PACusage Python package. **Rebuild it after every
`git pull`**, so the Python code in the image matches the R and Nextflow code
in the checkout. After that, runs need no network access; the image and inputs
can be copied to an offline cluster, with `container:` pointing at the image.

## Inputs

### Sample sheet

A tab-separated file with one row per sample:

```text
sample_id	alignment	condition	control
DMSO_1	bams/dmso_1.bam	DMSO
DMSO_2	bams/dmso_2.bam	DMSO
TRA_1	bams/treatment_a_1.bam	TreatmentA	DMSO
TRA_2	bams/treatment_a_2.bam	TreatmentA	DMSO
```

| Column | Meaning |
|---|---|
| `sample_id` | Unique; letters, digits, `.`, `_`, and `-` only. |
| `alignment` | BAM or CRAM, relative to the sample sheet. Sorted or not; PACusage sorts and indexes as needed. |
| `condition` | The sample's condition. |
| `control` | The condition it is compared with. Leave it **empty** (not `NA`) for a control condition. |

The rules:
- Every condition names one control, consistently across its samples.
- Every control is used by at least one treatment.
- Controls can chain, as in `WT → disease_vehicle → disease_drug`: a
  condition can be tested against its control and serve as another
  condition's control. There can be no cycles.
- Each condition needs at least two biological replicates. For an explicitly
  exploratory run, set `insufficient_replicates_policy: warn`.

Optional columns are `replicate`, `batch`, `donor`, `layout`, `strandedness`,
`library_profile`, and `evidence_source`; the last four override the matching
parameter for that sample. Any other column can be used as a model covariate
by naming it in `model_covariates`. Covariates are categorical and must be
filled in for every sample.

### analysis.yaml

Parameters go in a YAML file passed with `-params-file`. Four are required:

```yaml
input: samples.tsv
assembly: GRCh38            # used in PAC IDs
fasta: /refs/GRCh38.primary_assembly.genome.fa
gtf: /refs/gencode.v46.annotation.gtf
outdir: results
```

A typical HPC file adds the Slurm settings and the directories Apptainer must
see:

```yaml
library_profile: plasmidsaurus_3prime
slurm_account: my_lab
slurm_partition: general
scratch: /work/users/me       # work directory goes to <scratch>/pacusage-work
bind_paths:                   # every directory holding inputs, as a YAML list
  - /proj/my_lab
  - /work/users/me
```

`examples/analysis.yaml` has both. Relative paths in the YAML are resolved
against the directory you launch from. Values given on the command line (for
example `--min_mapq 30`) override the YAML. `nextflow run . --help` lists
every parameter with its default, and
[nextflow_schema.json](nextflow_schema.json) holds their full definitions.

List parameters (`bind_paths`, `excluded_contigs`, `model_covariates`) must be
YAML lists; a command-line value is not split.

### Mitochondrial reads are excluded by default

`excluded_contigs` defaults to `[chrM, MT, chrMT]`, the usual names of the
mitochondrial genome. Reads on these contigs are dropped by every step,
including strandedness inference, so mitochondrial genes get no PACs, counts,
or tests. There are two reasons:
- Calibration weights read ends by read, and mitochondrial transcripts are
  among the most abundant in most cells. Their packed, single-end genes would
  dominate the offset kernel.
- Their 3′ ends come from tRNA processing, not from the cleavage and
  polyadenylation that PACusage models.

A contig is excluded when its name in the alignment, or its
`chromosome_aliases` target, is on the list. Two things to know:
- **Other names aren't caught.** A mitochondrial contig named something else,
  such as `M` or `NC_012920.1`, isn't excluded. Either alias it to `chrM`, or
  list it.
- **A supplied list replaces the default.** It doesn't add to it, so
  `excluded_contigs: [chrY]` keeps mitochondrial reads, and
  `excluded_contigs: []` keeps every contig. To add a contig, repeat the
  mitochondrial names: `[chrM, MT, chrMT, chrY]`.

### Optional reference files

- **`known_pacs`**: a BED6 file of known PACs, such as a PolyASite atlas. Each
  record's PAC is its 3′ edge in transcript orientation: `end` on the plus
  strand, `start` on the minus strand. Known PACs mark atlas matches, and a
  known PAC with fewer reads than a novel one needs can be kept as
  `known_rescue_only`.
- **`chromosome_aliases`**: a two-column TSV mapping alignment contig names
  to FASTA and annotation names, for example `1` to `chr1`.
- **`pas_motif_catalog`**: a TSV of `motif`, `class`, and optional
  `priority`, replacing the built-in poly(A)-signal catalog.

### What PACusage never does to your files

Sources are never modified, and nothing is written beside them. The prepared
FASTA and every sorted alignment are **symbolic links** to your files, and any
index is built in the work directory. Only an unsorted alignment is copied,
into a sorted file in the work directory. Keep the sources in place and
unchanged until the analysis is finished; each step that reads an alignment
checks that its size and modification time still match. With Apptainer,
`bind_paths` must cover every directory the links point into.

## Run

On Slurm, with Apptainer:

```bash
nextflow run /path/to/pacusage \
  -profile slurm,apptainer \
  -params-file analysis.yaml \
  -resume
```

Run Nextflow itself from a login or interactive node, for example inside
`tmux` or a long-running batch job; it submits every step as its own Slurm
job. Site-wide settings, such as a default cluster queue or module loads, can
go in a Nextflow config file passed with `-c institution.config`.

Locally, with Conda:

```bash
nextflow run /path/to/pacusage -profile local,conda -params-file analysis.yaml -resume
```

The profiles are:
- `slurm`: submits to Slurm with `slurm_account`, `slurm_partition`, and
  `slurm_qos`;
- `local`: runs on this machine;
- `apptainer`: uses `containers/pacusage.sif`, or the image named by
  `container`;
- `conda`: builds `envs/pacusage.yml`;
- `test`: runs the small built-in fixture.

### Resources

Each step has a resource label in `conf/base.config`. A step killed for memory
or time (exit status 137, 140, or 143) is retried up to twice, with its memory
and time limits multiplied by the attempt number. The statistics have their
own settings:

| Parameter | Default | Controls |
|---|---|---|
| `statistics_cpus` | 8 | CPUs for each comparison family's model fit |
| `statistics_bootstrap_cpus` | 4 | CPUs for each bootstrap batch job |
| `statistics_bootstrap_batch_size` | 500 | genes per bootstrap batch job |
| `statistics_bootstrap_max_forks` | 8 | bootstrap batch jobs running at once |
| `dm_bootstrap_replicates` | 200 | bootstrap replicates per gene; `0` skips intervals |

## Outputs

```text
results/
  manifest/    resolved_params.yaml, normalized_samples.tsv, input_checksums.tsv,
               software_versions.tsv, run_manifest.json, library_resolution.json,
               calibration_kernel.tsv
  qc/          input validation, control mapping, reference and alignment
               preparation, strandedness, library calibration, kernel
               diagnostics, fragment filtering, PAC discovery, quantification
  evidence/    per-sample read-end tables (TSV and Parquet) and splice
               continuations
  tracks/      per-sample plus- and minus-strand bedGraphs of read ends
  atlas/       pacs.v1.bed.gz, pacs.v1.metadata.tsv.gz, pacs.v1.sha256,
               rejected_candidates.tsv.gz
  counts/      pac_counts.tsv.gz (samples as columns), pac_counts.long.parquet,
               gene_totals.tsv.gz, observed_pau.tsv.gz
  statistics/  per comparison: CONDITION_vs_CONTROL.genes / .pacs / .events;
               per family: FAMILY.gene_omnibus and FAMILY.statistical_filtering;
               fitted_pau.tsv.gz, gene_precision.tsv.gz
  motifs/      pac_motifs.tsv.gz, motif_scores.tsv, per-comparison motif
               preference and exploratory k-mer enrichment tables
  report/      index.html, a self-contained report
  pipeline_info/
               Nextflow's execution report, timeline, trace, and DAG; each
               run, resumed or not, replaces them
```

`prepared_reference/`, `prepared_alignments/`, and `counts/per_sample/` also
appear when `save_prepared_reference`, `save_prepared_alignments`, or
`save_intermediates` is set. The `motifs/` folder also holds class-level
versions of the scores and preference tables (`motif_class_scores*.tsv`,
`*.preference_class*.tsv.gz`).

A PAC ID such as `PACv1.GRCh38.chr1.+.1234567` holds the atlas version,
assembly, contig, strand, and the PAC's representative coordinate. Coordinates
are zero-based and interbase, so a plus-strand PAC at 1234567 is the boundary
after the 1,234,567th base.

### Columns

Every table with one row per PAC starts with the same identity columns:

| Column | Meaning |
|---|---|
| `pac_id` | The PAC's ID |
| `gene_id`, `gene_name` | Its gene, and the gene's name from the annotation (comma-separated, in the same order, for a PAC assigned to several genes) |
| `chrom`, `start`, `end` | Its interval in the BED convention (0-based start, exclusive end), the same interval as `atlas/pacs.v1.bed.gz`: `[coordinate, coordinate + 1)` for an exact-boundary PAC, and the resolution region for a proximal-tag PAC |
| `strand` | `+` or `-` |
| `locus` | The same interval, 1-based, such as `chr1:1234568-1234568`, to paste into IGV or the UCSC browser |

Gene-level tables start with `gene_id` and `gene_name`.
- **Where names come from:** a gene's name is its `gene_name` attribute, or
  `Name` on a GFF3 gene record, or NCBI's `gene` attribute.
- **Unnamed genes:** a gene without a name is named by its `gene_id`.
- **Counts:** `pac_counts.tsv.gz` has the identity columns, then one column per
  sample, with rows in genomic order.

After the identity columns, the results tables run from the answer to the
evidence. In `.pacs` and `.events` the groups are:
1. **the comparison and the call:** condition, control, `event_type`;
2. **the effect:** fitted PAU per group, the change in PAU, and its interval;
3. **significance:** `pac_fdr`, `gene_fdr`, the p-values, and the test
   statistics;
4. **the PAC's annotation;**
5. **the data behind the call:** supporting samples, gene totals, raw counts,
   and observed PAU;
6. **model diagnostics:** precision, stabilization, and the bootstrap.

The `.pacs` tables are the main results: one row per tested PAC. The
`.events` tables keep only the PACs with an event.

| Event | When |
|---|---|
| `gained`, `lost` | A PAC is detected (reads in at least two samples) in only one group, and that group's usage passes the fitted-PAU thresholds. The group without it has at least `min_gene_total` reads at the gene. The gene and the PAC pass their FDRs. The PAC is neither ambiguous nor flagged for internal priming. |
| `gained_candidate`, `lost_candidate` | The detection rules hold, but significance, stability, coverage, or confidence does not, or the comparison is exploratory. These are not confirmed calls. |
| `increased_usage`, `decreased_usage` | The PAC is detected in both groups, passes both FDRs, and its usage changes by at least `min_abs_delta_pau`. |
| `dominant_switch` | The most-used PAC differs between the groups, in a gene that passes `gene_fdr`. |
| `complexity_gain`, `complexity_loss` | The number of PACs with fitted PAU of at least `event_min_treatment_pau` differs, in a gene that passes `gene_fdr`. |

## How it works

1. **Validate.** The parameters, sample sheet, references, and alignment
   contigs are checked. The R statistics packages are loaded alongside, so a
   missing one stops the run within minutes.
2. **Prepare.** The FASTA and each alignment are linked and indexed, or
   sorted if needed. Each alignment is hashed once, for
   `manifest/input_checksums.tsv`.
3. **Strandedness.** It is inferred per sample from reads over unambiguous
   exons, unless the sample sheet or the profile sets it.
4. **Scan and calibrate.** Each alignment is read once to collect read ends
   for every candidate evidence source. Calibration then compares them with
   the annotated ends of genes whose transcripts share one end. It picks one
   evidence source and endpoint model for the whole run, and builds the
   offset kernel. Samples that calibrate differently stop the run: analyze
   different protocols separately.
5. **Discover.** Candidates are found in every sample's read ends pooled,
   without the treatment–control contrasts. A PAC must be supported by
   replicates within the condition that supports it:
   - **Exact-boundary** libraries cluster read ends directly, at nucleotide
     resolution.
   - **Proximal-tag** libraries find peaks of the kernel-matched signal, and
     report each PAC with its resolution interval rather than as a cleavage
     site. Candidates inside exon blocks that are spliced onward in every
     condition are rejected as readthrough.
6. **Annotate.** Genes, poly(A) signals, internal-priming flags, and known-PAC
   matches are added, and the atlas is frozen and checksummed.
7. **Count.** Each sample's read ends are assigned to the frozen atlas, giving
   raw counts, gene totals, and PAU.
8. **Test.** Each comparison family (a control and the treatments that name
   it) is fitted and tested. Motif-class preference and an exploratory k-mer
   enrichment follow.
9. **Report.** Everything is summarized in `report/index.html`.

### Library profiles

| `library_profile` | Use for | Behavior |
|---|---|---|
| `generic_3prime` (default) | Any 3′-end protocol | Calibrates every compatible evidence source and requires one to be clearly best. |
| `plasmidsaurus_3prime` | Plasmidsaurus 3′ tag-seq | Single-end, forward-stranded `read_3p` ends, modeled as proximal tags. |
| `exact_boundary` | Assays that keep the transcript/poly(A) junction | Uses poly(A)-junction ends when enough genes have them, otherwise the aligned 3′ edge. |

`endpoint_model` and `evidence_source` can override what calibration picks.

For proximal-tag runs, calibration also checks the pooled kernel and writes
`qc/calibration_kernel_diagnostics.tsv`. Read ends scattered around one site
give a kernel with a single mode. A kernel with several separated modes, or
one far wider than its resolution, warns in the Nextflow log and at the top
of the report. The run continues, but the warning usually means the profile
does not suit the data, or that unannotated alternative ends lie near the
annotated ones calibration uses.

## Statistics

- **Filtering.** Each comparison family filters its own samples' counts.
  - A gene needs `min_gene_total` reads and at least two PACs that pass.
  - A PAC needs `min_site_count` reads, `min_site_usage` of its gene's reads,
    and reads in `min_test_supporting_samples` samples.
  - `statistics/FAMILY.statistical_filtering.tsv.gz` lists every PAC as
    tested, or with the reasons it was not.
- **Model.** DRIMSeq fits a Dirichlet-multinomial model per gene with the
  design `~ model_covariates + condition`, with the family's control as the
  reference.
  - The family-wide test is descriptive.
  - Each treatment's own gene-level test is the screen. stageR confirms PACs
    within screened genes at `site_fdr`.
  - A reported significant PAC passes both `gene_fdr` and `site_fdr`. stageR
    confirms both PACs of a two-PAC gene together, so both report a
    `pac_fdr` of 0 when the gene passes.
- **Genes with a group without reads.** A gene with a PAC that has no reads
  in some group can't be fitted as it is.
  - It is refitted `dm_zero_sensitivity_repeats` times, with its zeros
    replaced by small seeded values (DRIMSeq's `addUniform` rule), and the
    median result is reported. These values never appear in any count or
    PAU table.
  - A PAC whose effect changes direction, spreads by more than
    `dm_zero_max_delta_pau_spread`, or fits in fewer than 80% of the refits
    is flagged `zero_boundary_unstable` and left untested.
  - These refits are seeded from `random_seed` and the atlas. With another
    seed, a stabilized gene's precision and p-values can move by orders of
    magnitude while its change in PAU barely moves, so read them alongside
    `model_status` (`drimseq_add_uniform`).
  - A condition with no reads at all for a gene leaves that gene untested in
    its comparisons (`model_status` `group_without_counts`).
- **Intervals.** A parametric bootstrap gives a 95% interval for the change
  in PAU in screened genes and detection candidates.
  - The intervals hold each gene's precision at its fitted value. In
    simulations with four replicates per group, they covered the true change
    in about 88% of cases (80-93%), so read them as approximate.
  - `bootstrap_status` says why an interval is missing.
- **Motif preference.** A limma model on the same design tests whether a
  treatment shifts usage toward PACs with a given poly(A) signal. It uses
  genes with at least `min_gene_total` reads in every sample, and each gene
  counts equally.
  - One table tests each primary hexamer (`.preference`).
  - A second table tests each motif class (`.preference_class`): canonical,
    the common `ATTAAA` variant, other variants, and no motif.
- **k-mer enrichment** is exploratory. It is a Cochran-Mantel-Haenszel test
  per k-mer, stratified by gene. It compares gained and increased PACs with
  the other PACs of their genes that were tested in the same comparison.

## Reproducibility and -resume

Every run records:
- the resolved parameters;
- the input checksums;
- the software versions, including the R packages;
- a checksum of the frozen atlas.

Gzip files carry no timestamps and every random step is seeded, so a fresh
rerun with the same inputs, parameters, and software reproduces every
published file byte for byte. The exception is Nextflow's own reports in
`pipeline_info/`, which record the run's times.

Keep the Nextflow work directory until the analysis is final: `-resume` reuses
finished steps from it.

- **Which steps rerun.** Every step from the alignment scan onward reads the
  resolved parameters, so changing any parameter reruns them all. That
  includes `outdir`, the `save_*` flags, and the Slurm and CPU settings.
  Alignment preparation and strandedness inference are reused unless their
  own inputs change. Set those parameters before the first run.
- **After updating PACusage.** Nextflow does not track the Python package, so
  start a fresh run: use a new work directory, or leave out `-resume`. The R
  statistics script is tracked, and changes to it rerun the statistics.

## Troubleshooting

- **Mitochondrial genes are missing from the results**: that's the
  `excluded_contigs` default. See
  [Mitochondrial reads are excluded by default](#mitochondrial-reads-are-excluded-by-default).
  A mitochondrial contig under another name, such as `M`, is not excluded until
  you alias or list it.
- **`Invalid parameters: ...`**: the message names each bad parameter, what
  it must be, and what it means. Check `analysis.yaml` against `--help`.
- **`strandedness is ambiguous`**: too few reads overlap unambiguous exons, or
  the orientation is mixed. Set `strandedness` for the samples if the
  protocol is known.
- **`Samples use incompatible library profiles`**, **`do not resolve to one
  compatible endpoint model`**, or **`No evidence source passed the
  calibration reproducibility threshold`**: the samples don't calibrate
  alike. Run each protocol separately, or set `library_profile`,
  `endpoint_model`, or `evidence_source` when the protocol is known.
- **`No PAC candidates passed discovery`**:
  - `qc/pac_discovery.tsv` and `atlas/rejected_candidates.tsv.gz` say why
    each candidate was rejected;
  - a calibration kernel warning earlier in the log often names the cause.
- **`Source alignment ... is not readable here`**: the linked source is not
  visible inside the container. Add its directory to `bind_paths`.
- **`Source alignment ... changed after preparation`**: a source file was
  modified or replaced during the run. Restore it, or rerun with `-resume` so
  that preparation records the new file.
- **`FASTA ... is compressed`**: decompress it (`gunzip` or `bgzip -d`) and
  point `fasta` at the uncompressed file.
- **A step fails or runs slowly**: `pipeline_info/execution_report.html` shows
  each task's run time and memory. In `pipeline_info/execution_trace.txt`, a
  task's `native_id` is its Slurm job ID, and its `hash`, such as `3f/a1b2c3`,
  begins its path under the work directory, where `.command.log` and
  `.command.err` hold its output.

## Development

```bash
python -m venv .venv
.venv/bin/pip install -e '.[dev]'
.venv/bin/ruff check src tests
.venv/bin/python -m pytest
```

The R tests need DRIMSeq, stageR, limma, BiocParallel, and yaml. Set
`R_LIBS` to a library that has them:

```bash
Rscript tests/r/test_usage_model.R scripts/fit_usage_model.R
Rscript tests/r/test_usage_model_simulation.R scripts/fit_usage_model.R
```

`tests/pipeline/run_nextflow.sh` runs everything end to end, with `.venv/bin`
first on `PATH`. It builds the fixtures (`tests/fixtures/build_fixture.py`)
and runs the R tests, then runs the `test` profile and checks:
- the results, and that a fresh rerun reproduces every file outside
  `pipeline_info/`;
- that a mixed-protocol sample sheet stops before discovery;
- the Plasmidsaurus-like fixture;
- that exact-boundary reads run as Plasmidsaurus tags warn at calibration.

| Path | Contents |
|---|---|
| `main.nf`, `workflows/`, `subworkflows/local/`, `modules/local/` | The Nextflow pipeline, one process per module file |
| `nextflow.config`, `conf/`, `nextflow_schema.json` | Parameters, profiles, resources, and the parameter schema |
| `src/pacusage/` | The Python package behind every `pacusage` step |
| `scripts/fit_usage_model.R` | The statistics: DRIMSeq, stageR, bootstrap, events, motif preference |
| `bin/pacusage` | The command Nextflow tasks run |
| `containers/`, `envs/` | The Apptainer recipe and the Conda environment |
| `tests/unit/`, `tests/integration/` | pytest tests |
| `tests/r/` | R tests for the statistics |
| `tests/pipeline/` | The end-to-end script and its result checks |
| `tests/fixtures/` | Synthetic references and alignments, and their builders |
| `docs/` | The design and the decisions log |
