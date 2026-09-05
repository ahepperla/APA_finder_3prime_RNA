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
steps. Low, medium, and high-memory jobs automatically retry up to twice after
an OOM-style exit, increasing their requested memory on each attempt. Completed
tasks remain reusable with `-resume`.

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
of annotated transcript ends, submits one calibration task per BAM/CRAM, and
then combines the small per-sample summaries into the run-level calibration
files. On Slurm, a 48-sample run can therefore schedule up to 48 independent
calibration jobs instead of scanning all alignments in one large-memory job.

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
condition.

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

The integration fixture can be generated with:

```bash
python tests/fixtures/build_fixture.py
nextflow run . -profile test,conda -resume
```

The statistical process requires DRIMSeq, stageR, and limma. The supplied
Conda files install them from Bioconda. Runs are network-independent after the
environment or container has been prepared.

## Reproducibility

PACusage never modifies source FASTA or alignment files. It records input
checksums, resolved parameters, software versions, preparation actions, and a
checksum of the frozen atlas. Nextflow work directories provide resumability;
retain them until the analysis is accepted.
