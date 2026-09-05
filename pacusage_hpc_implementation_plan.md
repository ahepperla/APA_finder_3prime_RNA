# PACusage: Nextflow Implementation Plan

## Goal

Build a readable, reproducible command-line pipeline that:

1. Starts from deduplicated BAM or CRAM files.
2. Prepares any missing FASTA or alignment indexes.
3. Sorts alignments when necessary.
4. Extracts protocol-appropriate 3-prime evidence from each accepted read or
   fragment.
5. Discovers and annotates candidate polyadenylation site clusters at the
   resolution supported by the library protocol.
6. Counts PAC usage within each gene.
7. Tests every non-control condition against its declared control condition.
8. Runs locally or on a Slurm HPC cluster through Nextflow.

Read trimming, alignment, and deduplication are outside the pipeline.

Use these terms consistently:

- **End observation:** A transcript-oriented read or fragment coordinate used
  as 3-prime evidence. It is not automatically a cleavage coordinate.
- **Candidate PAC:** A polyadenylation site inferred from clustered end
  observations or a calibrated 3-prime-proximal signal.
- **Exact PAC:** A nucleotide-resolution site supported by a protocol that
  preserves the transcript/poly(A) junction.
- **PAS motif:** An upstream sequence motif such as AATAAA.
- **PAU:** A PAC's fraction of all PAC counts assigned to its gene in one
  sample.

## User Experience

The standard user action must be a single command.

Local:

```text
nextflow run /path/to/pacusage \
  -profile local,conda \
  -params-file analysis.yaml \
  -resume
```

Slurm:

```text
nextflow run /path/to/pacusage \
  -profile slurm,apptainer \
  -params-file analysis.yaml \
  -resume
```

Do not require separate initialization, validation, testing, or reporting
commands. Validation runs first and reporting runs last.

Follow the nf-core parameter model:

- checked-in configuration defines documented defaults;
- `-params-file analysis.yaml` overrides those defaults;
- explicit command-line values such as `--min_mapq 30` override the YAML;
- an optional `-c institution.config` supplies site-specific Nextflow settings;
- the resolved parameters are printed and saved;
- `--help` is generated from a parameter schema.

The run must work without network access after environments or containers have
been prepared.

## Inputs

### Sample Sheet

Require a tab-separated file with four columns:

```text
sample_id  alignment   condition   control
DMSO_1     /data/a.bam DMSO
DMSO_2     /data/b.bam DMSO
TRA_1      /data/c.bam TreatmentA DMSO
TRA_2      /data/d.bam TreatmentA DMSO
TRB_1      /data/e.bam TreatmentB DMSO
TRB_2      /data/f.bam TreatmentB DMSO
```

The `control` cells for DMSO must be truly empty TSV fields.

Required meanings:

- `sample_id`: unique filesystem-safe sample name.
- `alignment`: BAM or CRAM path, resolved relative to the sample sheet.
- `condition`: experimental condition.
- `control`: control condition name, or blank when the row belongs to a control
  condition.

Validation rules:

1. Every sample ID is unique.
2. Every condition has one consistent `control` value across all its samples.
3. A blank `control` marks that row's condition as a control condition.
4. Every nonblank `control` exactly matches a condition in the sample sheet.
5. Every referenced condition is a control condition with blank values on all
   its rows.
6. Every control condition is referenced by at least one non-control condition.
7. A condition cannot reference itself or another treatment condition.
8. Every modeled control and treatment condition has at least two biological
   replicates by default.

Create a normalized `control_condition` field:

- control rows use their own `condition`;
- non-control rows use their declared `control`.

Conditions sharing a `control_condition` form a comparison family. For example,
`DMSO`, `TreatmentA`, and `TreatmentB` form the `DMSO` family. Multiple
comparison families are allowed.

The control mapping is condition-level. It does not imply paired samples.

### Optional Sample Columns

Support:

- `replicate`: descriptive replicate label;
- `batch`: technical batch;
- `donor`: biological individual or paired source;
- `layout`: `SE`, `PE`, or `auto`;
- `strandedness`: `forward`, `reverse`, or `auto`;
- `library_profile`: a supported protocol profile or `generic_3prime`;
- `evidence_source`: `auto`, `read_3p`, `read_5p`, `fragment_3p`, or
  `polyA_junction`;
- additional model covariates.

Behavior:

- missing `layout` inherits the pipeline default, which is `auto`;
- missing `strandedness` inherits the pipeline default, which is `auto`;
- missing `library_profile` inherits the pipeline-level profile, which defaults
  to `generic_3prime`;
- missing `evidence_source` inherits `auto` and is resolved by the selected
  library profile and calibration;
- `batch`, `donor`, and custom covariates enter the model only when named in
  `model_covariates`;
- every covariate used in a model must be complete for the samples in that
  comparison family;
- `replicate` does not create statistical pairing by itself.

### Reference Files

Require:

- genome assembly name;
- genome FASTA;
- GTF or GFF3 annotation.

Optionally accept:

- a known-PAC BED file;
- chromosome-name aliases.

The FASTA does not need to be indexed. A `PREPARE_REFERENCE` process must:

1. Validate a supplied `.fai` when available.
2. Generate an index with `samtools faidx` when needed.
3. Write generated files in the Nextflow work directory.
4. Never modify the source FASTA or its directory.
5. Emit the prepared FASTA and index together.
6. Fail clearly if the FASTA is not indexable.

CRAM inputs must be compatible with the supplied FASTA.

## Parameters

Use flat, descriptive names in `nextflow_schema.json`. A typical user YAML:

```yaml
input: samples.tsv
outdir: results

assembly: GRCh38
fasta: /reference/GRCh38.fa
gtf: /reference/gencode.gtf
known_pacs: null

layout: auto
strandedness: auto
library_profile: generic_3prime
endpoint_model: auto
evidence_source: auto
min_mapq: 20
require_unique: true
require_proper_pair: true
exclude_duplicates: true
excluded_contigs: []

min_replicates_per_condition: 2
insufficient_replicates_policy: error

strand_min_informative_fragments: 10000
strand_max_sampled_fragments: 200000
strand_decision_fraction: 0.80

calibration_min_genes: 100
calibration_max_distance: 1000
calibration_quantile_low: 0.05
calibration_quantile_high: 0.95
calibration_min_kernel_correlation: 0.80
calibration_min_model_margin: 0.10
exact_max_median_abs_offset: 2
exact_max_central_width: 12
exact_min_boundary_fraction: 0.50
proximal_min_median_upstream_offset: 3

pac_seed_radius: 2
pac_cluster_radius: 12
pac_min_total_count: 10
pac_min_sample_count: 2
pac_min_supporting_samples: 2
known_pac_rescue_total: 5
known_pac_match_radius: 12
proximal_kernel_overlap_threshold: 0.50
proximal_assignment_likelihood_ratio: 3.0
proximal_bin_size: 25

max_downstream_distance: 5000
pas_scan_upstream_far: 50
pas_scan_upstream_near: 5
pas_core_upstream_far: 35
pas_core_upstream_near: 10
pas_motif_catalog: null
internal_priming_window: 20
internal_priming_max_a_run: 6
internal_priming_max_a_fraction: 0.60

min_gene_total: 20
min_site_count: 5
min_site_usage: 0.01
min_test_supporting_samples: 2

model_covariates: []
gene_fdr: 0.05
site_fdr: 0.05
min_abs_delta_pau: 0.10
dm_bootstrap_replicates: 200
dm_bootstrap_min_success_fraction: 0.80
dm_zero_sensitivity_repeats: 5
dm_zero_max_delta_pau_spread: 0.02
random_seed: 1729
event_min_treatment_pau: 0.05
event_max_control_pau: 0.01
event_min_supporting_samples: 2
motif_preference_min_genes: 50
motif_kmer_length: 6
run_kmer_enrichment: true

save_prepared_reference: false
save_prepared_alignments: false
save_intermediates: false
```

Defaults must permit a standard run with a shorter YAML containing only
`input`, `outdir`, `assembly`, `fasta`, and `gtf`.

`insufficient_replicates_policy` accepts `error` or `warn`. The default
`error` stops before modeling. `warn` permits an explicitly exploratory run,
marks every affected comparison in all statistical outputs, and suppresses
claims of confirmed gained or lost PACs.

For Plasmidsaurus data, the user explicitly overrides the generic default:

```yaml
library_profile: plasmidsaurus_3prime
```

## Workflow

### 1. Validate Inputs

Validate:

- parameters against `nextflow_schema.json`;
- sample-sheet structure and control relationships;
- file existence and readability;
- BAM/CRAM format;
- FASTA, annotation, and alignment contig compatibility;
- optional covariate completeness;
- output paths.

Write a normalized sample sheet and a resolved parameter file before expensive
work begins.

### 2. Prepare Reference and Alignments

For each BAM or CRAM:

1. Check file integrity.
2. Detect whether it is coordinate sorted.
3. Reuse a valid index when possible.
4. Index a sorted but unindexed file.
5. Sort and index an unsorted file with `samtools sort` and `samtools index`.
6. Preserve BAM as BAM and CRAM as CRAM.
7. Never modify the source file.
8. Record whether the file was reused, indexed, or sorted and indexed.

Prepared files remain in the Nextflow cache unless publication is requested.

### 3. Resolve Protocol, Layout, Strandedness, and Calibration

Separate protocol interpretation from downstream PAC analysis. Implement a
small strategy interface so additional 3-prime protocols can be added without
changing annotation, quantification, or statistics.

Support these initial profiles:

#### `plasmidsaurus_3prime`

Defaults:

- single-end layout;
- forward strandedness;
- `read_3p` evidence source;
- reads sequence toward the transcript 3-prime end;
- `proximal_tag` endpoint model;
- expected signal concentrated within roughly the terminal 400 nt.

Treat these as defaults that are verified from the BAM, not unquestionable
truths. Allow sample-level layout or strandedness overrides.

#### `exact_boundary`

Use when the assay preserves the transcript/poly(A) junction or another exact
3-prime boundary. Resolve one source for the whole sample:
`polyA_junction` when it is observed in at least `calibration_min_genes`, with
at least `pac_min_sample_count` observations in each contributing gene;
otherwise use the profile-defined read or fragment edge. Do not mix source
definitions within a sample. The transcript-oriented boundary may be clustered
directly as nucleotide-level PAC evidence.

#### `generic_3prime`

Infer layout and strandedness, evaluate the compatible evidence sources, then
run calibration against high-confidence annotated transcript ends. For SE data,
evaluate `read_3p`, `read_5p`, and `polyA_junction` when available. For PE data,
also evaluate `fragment_3p`. Select a source and endpoint model only when its
reproducibility score passes `calibration_min_kernel_correlation` and exceeds
the next-best model by `calibration_min_model_margin`. Fail with a clear
message when the evidence is ambiguous. Never silently choose the first
candidate source.

Allow `endpoint_model` to be explicitly set to:

```text
auto
exact_boundary
proximal_tag
```

The explicit setting overrides the profile default but must still produce
calibration diagnostics.

#### Evidence Sources

Protocol profiles resolve one evidence source before PAC discovery:

- `read_3p`: aligned transcript-oriented 3-prime edge of the profile-selected
  read;
- `read_5p`: aligned transcript-oriented 5-prime edge of the profile-selected
  read;
- `fragment_3p`: outer transcript-oriented 3-prime edge of a reconstructed
  proper pair;
- `polyA_junction`: genomic boundary adjacent to a compatible terminal
  poly(A)-like soft clip.

The profile strategy must define which read is selected when a paired protocol
uses a read-level source. The default behavior remains the SE read 3-prime edge
or PE fragment 3-prime edge. Adding a protocol should require a new profile
definition, not changes to quantification or statistics.

#### Calibration

Use unambiguous, well-covered, single-terminal-site genes to compare observed
3-prime boundaries with annotated transcript ends. Define a positive offset as
an observation upstream of the annotated transcript end. Estimate:

- strand-oriented endpoint offset distribution;
- median and modal offset;
- interval containing the configured central fraction of observations;
- fraction ending at or immediately adjacent to the annotated end;
- fraction carrying compatible terminal poly(A)-like soft clipping;
- effective spatial resolution.

For each candidate source, smooth and normalize a per-sample offset histogram
as described under proximal discovery. Define its reproducibility score as the
median Pearson correlation between each sample kernel and the equal-weight
leave-one-sample-out pooled kernel. Choose the source and endpoint model once
for the run using the median sample score; every individual sample must still
pass the compatibility threshold.

Classify calibration deterministically:

- `exact`: reproducibility passes the threshold, absolute median offset is at
  most `exact_max_median_abs_offset`, the configured central interval width is
  at most `exact_max_central_width`, and the adjacent-boundary fraction is at
  least `exact_min_boundary_fraction`;
- `proximal`: reproducibility passes the threshold, the median upstream offset
  is at least `proximal_min_median_upstream_offset`, and the configured central
  interval lies within `calibration_max_distance`;
- `ambiguous`: too few calibration genes, inadequate reproducibility,
  conflicting source scores, or neither quantitative definition passes.

Store the inferred endpoint model, calibration interval, and resolution for
every sample. Do not silently call nucleotide-resolution PACs from an
ambiguous or proximal profile.

Version 1 constructs one frozen atlas, so all samples in a run must resolve to
the same evidence source and endpoint model. Pool calibration kernels with
equal weight per sample and require every sample to correlate with the pooled
kernel by at least `calibration_min_kernel_correlation`. Fail the run when
protocol, evidence source, or calibration resolution is incompatible across
samples, and instruct the user to analyze incompatible protocols in separate
runs. Do not attempt to statistically correct incompatible library chemistries
with a batch term.

#### Layout

When `layout=auto`, inspect alignment flags to classify the sample as SE or PE.
Fail if the sample contains a meaningful incompatible mixture.

#### Strandedness Definitions

For SE:

- `forward`: alignment strand equals transcript strand;
- `reverse`: alignment strand is opposite transcript strand.

For PE:

- `forward`: read 1 strand equals transcript strand;
- `reverse`: read 1 strand is opposite transcript strand.

#### Automatic Strandedness Inference

Construct diagnostic intervals from annotated exonic sequence that:

- belongs unambiguously to one gene;
- has no opposite-strand overlap;
- prioritizes terminal exons.

Reservoir-sample eligible alignments so coordinate order does not bias the
result. Report:

```text
informative_fragments
forward_count
reverse_count
forward_fraction
reverse_fraction
inferred_strandedness
```

Default decision:

- forward fraction at least 0.80: `forward`;
- reverse fraction at least 0.80: `reverse`;
- otherwise: fail as ambiguous unless the user supplied a value.

Warn prominently when samples disagree with one another or with their selected
library profile.

### 4. Extract 3-Prime Evidence

Apply these default exclusions:

- unmapped;
- secondary;
- supplementary;
- QC failed;
- duplicate marked;
- below the MAPQ threshold;
- excluded contig;
- non-unique when a reliable `NH` or equivalent tag is available.

For paired-end input, count one fragment per pair and require by default:

- both mates mapped;
- both mates passing filters;
- both mates on the same contig;
- primary alignments;
- proper-pair flag;
- consistent query name and mate designation.

Use a name-collated stream, such as `samtools collate`, to reconstruct pairs.
Do not repeatedly search a coordinate-sorted BAM for mates. Do not impose a
maximum genomic insert size because spliced fragments may span large introns.

Use zero-based interbase coordinates internally:

- plus-strand transcript: maximum outer aligned reference end;
- minus-strand transcript: minimum outer aligned reference start.

Soft-clipped bases do not extend the genomic boundary. Retain aggregate
soft-clip information for poly(A)-tail evidence.

With the default evidence source, each accepted SE read or PE pair contributes
one transcript-oriented end observation:

- SE: the aligned 3-prime edge of the read;
- PE: the outer transcript-oriented 3-prime edge of the reconstructed fragment.

When a profile resolves another evidence source, extract the corresponding
read edge or soft-clip junction using the same strand and coordinate
conventions. Record the resolved source with every aggregate observation.

Aggregate identical chromosome, strand, and boundary combinations per sample.
Store the observation as exact-boundary evidence only for an
`exact_boundary` profile. For `proximal_tag`, retain it as an observed endpoint
whose offset distribution is interpreted by the calibrated profile.

### 5. Build the PAC Atlas

Build one condition-blind atlas from all samples. Do not use CPM, condition
labels, or PAS motifs during primary discovery.

Use separate discovery strategies:

- `exact_boundary`: apply the support-ranked clustering below directly to
  observed endpoints;
- `proximal_tag`: apply the deterministic matched-kernel procedure below and
  report an estimated PAC coordinate with an assay-resolution interval.

#### Deterministic Proximal-Tag Discovery

For the compatible samples in the run:

1. Represent every calibration observation as a signed transcript-oriented
   offset from its high-confidence annotated transcript end.
2. Build a one-nucleotide empirical offset histogram for each sample over the
   calibrated range. Smooth it once with normalized triangular weights
   `[1, 2, 3, 2, 1]`, normalize it to sum to one, and average sample kernels
   with equal sample weight.
3. Bin observed endpoints and the pooled calibration kernel at
   `proximal_bin_size` resolution. For every implied regional PAC coordinate
   \(x\), calculate the matched score
   \(L_s(x)=\sum_e \min(y_{es},3)K_s(e-x)\), where \(e\) is an observed
   endpoint, \(y_{es}\) is its sample count, and \(K_s\) is the sample kernel.
4. Identify local maxima of the summed matched score. Rank maxima by supporting
   sample count, summed capped support, matched score, raw support, and genomic
   coordinate, in that order.
5. Apply the same total-count, per-sample-count, and supporting-sample
   requirements used for exact-boundary candidates.

Search only coordinates implied by observed endpoints shifted over the
nonzero calibrated kernel support. Implement scoring as streaming,
FFT-backed convolution by chromosome and strand so the algorithm does not
expand every endpoint across every kernel offset or retain all samples'
evidence in memory.

Define kernel overlap at separation \(\delta\) as:

\[
O(\delta)=\sum_d \min(K(d),K(d-\delta))
\]

The minimum resolvable separation is the smallest positive \(\delta\) for which
`O(delta) <= proximal_kernel_overlap_threshold`. Merge maxima closer than this
distance into one PAC resolution group. Store the strongest maximum as the
representative coordinate, retain the merged maxima, and report
`region_start`, `region_end`, and `resolution_nt`. Never report unresolved
maxima as independently quantified PACs or interpret the representative as a
nucleotide-resolution cleavage site.

For `proximal_tag`, known transcript ends may be used for calibration and
annotation, but PAS motif sequence must not be used to select or reposition
primary candidates. This prevents circular motif-preference results. An
optional motif-assisted rescue pass may be produced separately, and rescued
sites must be excluded from motif-preference testing.

#### Exact-Boundary Clustering

For every observed coordinate, calculate support within
`coordinate +/- pac_seed_radius`:

- number of supporting samples;
- capped per-sample support, using `min(local_count, 3)` by default;
- raw local fragment count.

Rank candidate seeds by:

1. number of supporting samples;
2. capped support;
3. raw support;
4. coordinate for deterministic tie-breaking.

For each chromosome and strand:

1. Select the strongest unassigned seed.
2. Assign unassigned endpoints within `seed +/- pac_cluster_radius`.
3. Store the seed as the representative coordinate.
4. Repeat until all candidates are assigned.

Default unannotated PAC requirements:

- at least 10 total fragments;
- at least two supporting samples;
- at least two fragments in each supporting sample.

A known PAC may be retained with:

- at least five total fragments;
- at least two supporting samples.

Mark a known PAC that passes only the relaxed known-site thresholds as
`known_rescue_only`. Exclude these sites, as well as motif-assisted rescue
sites, from primary motif-preference and k-mer analyses. Produce a sensitivity
table that repeats motif summaries with known-rescue-only sites included.

Keep rejected candidates and their reasons. For exact-boundary data, the
default 12 nt radius defines the initial assay resolution and should be tested
at 8, 12, and 18 nt. For proximal-tag data, derive the assignment window and
minimum resolvable PAC separation from calibration rather than pretending the
same 12 nt resolution applies.

### 6. Identify and Annotate PACs

Assign IDs using genomic identity:

```text
PACv1.GRCh38.chr1.+.123456789
```

The ID contains atlas version, assembly, contig, strand, and representative
interbase coordinate. Gene names and proximal/distal ranks are annotations, not
primary identifiers.

Also store:

- `endpoint_model`;
- `coordinate_precision`: `exact` or `estimated`;
- assay-resolution region and width;
- PAC resolution-group identifier and merged candidate coordinates;
- calibration profile and version;
- resolved evidence source;
- primary evidence type, such as exact boundary, endpoint peak, coverage edge,
  or known-site match.

Freeze and checksum the atlas before quantification and testing. Rebuilding
with new samples or parameters creates a new atlas version.

For each PAC, record:

- total count and sample support;
- fraction of ends within 2 nt of the representative;
- width containing 90 percent of ends;
- local strand-specific enrichment;
- compatible terminal poly(A)-like soft clipping;
- strand-oriented upstream and downstream sequence;
- every recognized upstream PAS motif, its position, and motif class;
- one primary PAS motif selected by catalog priority and expected position;
- downstream A fraction and longest A run;
- same-strand known-PAC matches;
- gene and transcript assignments;
- terminal-exon, other-exon, intronic, downstream, or intergenic class;
- proximal-to-distal rank;
- confidence tier and explicit reasons.

Internal-priming evidence is a flag by default, not an automatic exclusion.
Strict filtering may be exposed as an optional parameter.

Assign genes in this order:

1. Same-strand terminal exon.
2. Same-strand nonterminal exon.
3. Same-strand intron.
4. Downstream of a transcript end within the configured distance, without
   crossing an intervening same-strand gene.
5. Intergenic.

Retain ambiguous assignments in the atlas but exclude them from primary
within-gene testing.

### 7. Quantify the Frozen Atlas

Recount every sample after the atlas is frozen.

For exact-boundary profiles, define nonoverlapping PAC territories by:

- capping assignment distance at `pac_cluster_radius`;
- splitting nearby territories at the midpoint between representatives;
- assigning each end observation to at most one PAC.

For proximal-tag profiles, construct assignment regions from the calibrated
strand-oriented endpoint distribution. Unresolved nearby candidates have
already been merged into one PAC resolution group. When regions for distinct,
resolvable PACs overlap, assign an observation to the maximum-likelihood PAC
only when its likelihood is at least
`proximal_assignment_likelihood_ratio` times the next-best likelihood;
otherwise mark it ambiguous. Keep final statistical counts as integers and
report ambiguous counts by gene and sample.

For gene \(g\), PAC \(p\), and sample \(s\):

\[
y_{gps} = \text{raw PAC count}
\]

\[
n_{gs} = \sum_p y_{gps}
\]

\[
PAU_{gps} = \frac{y_{gps}}{n_{gs}}
\]

PAU is missing when `n_gs=0`. Do not add pseudocounts to count or observed-PAU
matrices.

Verify:

- accepted fragments equal assigned plus unassigned fragments;
- PAC counts sum to gene totals;
- PAU sums to 1 whenever the gene total is positive.

### 8. Filter Testable Genes

Keep discovery filtering separate from statistical filtering.

Suggested defaults:

```text
at least 2 PACs per gene
gene total >= 20
PAC count >= 5
PAC usage >= 1 percent
support in >= 2 samples
```

Apply filters across all samples in a comparison family without requiring
support in a particular condition. Record why every excluded gene or PAC was
not tested.

### 9. Model Differential PAC Usage

Use raw PAC counts for inference and PAU for interpretation.

For gene \(g\), sample \(s\), and \(K_g\) PACs:

\[
\mathbf p_{gs} \sim
Dirichlet(\boldsymbol\alpha_{gs})
\]

\[
\mathbf y_{gs} \mid \mathbf p_{gs}, n_{gs}
\sim Multinomial(n_{gs}, \mathbf p_{gs})
\]

Parameterize:

\[
\boldsymbol\alpha_{gs} = \kappa_g\boldsymbol\pi_{gs}
\]

where:

- \(\boldsymbol\pi_{gs}\) is expected PAC usage;
- \(\kappa_g\) is gene-level precision shared across conditions in version 1;
- lower precision represents more variation among biological replicates.

Control counts estimate baseline probabilities and precision. They are not used
directly as alpha values.

Within each comparison family, build:

```text
full model:    ~ model_covariates + condition
reduced model: ~ model_covariates
```

When `model_covariates` is empty, these become:

```text
full model:    ~ condition
reduced model: ~ 1
```

Set the family's declared control condition as the condition reference level.
If `donor` or `batch` should adjust the model, the user adds it to
`model_covariates`.

Validate that every design is full rank. Fail with a readable explanation of
missing or confounded covariates.

Use DRIMSeq for:

- a family-wide gene-level omnibus likelihood-ratio test;
- a contrast-specific gene-level test for every treatment versus its declared
  control;
- condition coefficients for every treatment versus its declared control;
- fitted PAC proportions;
- gene-level precision.

Do not generate treatment-versus-treatment comparisons in version 1.

The family-wide omnibus result is descriptive and must not gate a particular
treatment-control event. For each declared comparison, use its
contrast-specific gene p-value as the screening p-value and its PAC-level
p-values as confirmation p-values. Apply `stageR` two-stage adjustment within
that comparison:

1. Screen genes for the selected treatment-control contrast.
2. Confirm contributing PAC effects within screened genes.
3. Report BH-adjusted gene-level screening FDR using `gene_fdr`.
4. Run `stageR` at `site_fdr` and report its stage-wise PAC adjusted p-values.
   A reported significant PAC must pass both the gene and PAC thresholds.
5. State explicitly that the correction family is one declared
   treatment-control comparison.

For every treatment-control result, report:

- fitted control PAU;
- fitted treatment PAU;
- delta PAU, defined as treatment minus control;
- a 95 percent parametric-bootstrap confidence interval when the gene is
  bootstrap-eligible;
- raw replicate counts and observed PAU;
- gene-level precision;
- reconstructed alpha values
  `alpha = fitted_probability * precision`;
- omnibus and PAC-level p-values and FDR;
- whether the effect exceeds `min_abs_delta_pau`.

For example:

```text
DMSO:      original PAC 100%, new PAC 0%
Treatment: original PAC  50%, new PAC 50%
```

Report:

```text
original PAC delta PAU = -0.50
new PAC delta PAU      = +0.50
```

Never alter the stored raw counts or observed PAU. Fit the unmodified count
matrix first. When an otherwise testable PAC is all zero in one group and
DRIMSeq returns a boundary or non-finite fit, use DRIMSeq's documented
`add_uniform` behavior only as a model-fitting fallback:

1. Derive a deterministic seed from `random_seed`, atlas checksum, comparison,
   gene ID, and repeat number.
2. Keep PAC ordering deterministic by genomic coordinate.
3. Record that stabilization was used; never write the perturbed values to a
   count or PAU output.
4. Repeat the stabilized fit `dm_zero_sensitivity_repeats` times.
5. Flag the PAC as `zero_boundary_unstable` when effect direction changes,
   fewer than 80 percent of fits converge, or fitted delta PAU spans more than
   `dm_zero_max_delta_pau_spread`.
6. Do not assign a PAC-level p-value when stable finite fits cannot be obtained.
   Keep the detection evidence and label the inferential result unavailable
   rather than inventing a value.

For PACs in genes that pass the contrast-specific screen or satisfy the
pre-statistical gained/lost event criteria, estimate delta-PAU uncertainty with
a parametric bootstrap:

1. Simulate sample count vectors from the fitted Dirichlet-multinomial model at
   each sample's observed gene total and design row.
2. Refit the identical model and recompute fitted treatment-control delta PAU.
3. Use `dm_bootstrap_replicates` successful replicates and deterministic
   per-gene seeds.
4. Report percentile 95 percent intervals and bootstrap success counts.
5. Leave the interval unavailable and flag the reason when the success
   fraction is below `dm_bootstrap_min_success_fraction`.

#### Researcher-Facing PAC Events

Create one gene-level event summary for every treatment-control comparison.
Classify PACs using both statistical evidence and detection evidence:

- **gained:** not detectably used in adequately covered controls, reproducibly
  detected in treatment, positive delta PAU, and significant PAC-level FDR;
- **lost:** the reverse of gained;
- **increased usage:** detected in both groups with significant positive delta
  PAU;
- **decreased usage:** detected in both groups with significant negative delta
  PAU;
- **dominant switch:** the highest-usage PAC differs between treatment and
  control;
- **complexity gain or loss:** the number of reproducibly detected PACs changes.

When detection and effect-size requirements for `gained` or `lost` pass but a
stable PAC-level p-value cannot be obtained, report `gained_candidate` or
`lost_candidate`. These labels are detection-supported but not statistically
confirmed and must never be merged with significant gained/lost counts.

Default gained-PAC requirements:

- the parent gene passes the relevant contrast-specific gene test;
- control samples meet the gene-coverage threshold;
- fitted control PAU is at most `event_max_control_pau`;
- fitted treatment PAU is at least `event_min_treatment_pau`;
- treatment support meets `event_min_supporting_samples`;
- delta PAU is at least `min_abs_delta_pau`;
- PAC-level FDR passes `site_fdr`;
- the PAC is not low confidence, `zero_boundary_unstable`, or flagged as likely
  internal priming.

Use the label `gained` or `not detected in control`, not `de novo`, unless
independent evidence establishes biological absence in the control.

The event table should contain:

```text
gene_id
condition
control_condition
event_type
pac_id
site_class
known_pac
known_rescue_only
control_supporting_samples
treatment_supporting_samples
control_gene_total
treatment_gene_total
fitted_control_pau
fitted_treatment_pau
delta_pau
gene_fdr
pac_fdr
confidence
internal_priming_flag
model_status
zero_boundary_unstable
delta_pau_ci_low
delta_pau_ci_high
bootstrap_successes
dominant_pac_control
dominant_pac_treatment
```

#### Motif Preference

Classify each PAC into one primary motif class:

- canonical `AATAAA`;
- common `ATTAAA` variant;
- other recognized PAS variant;
- no recognized motif.

Use DNA alphabet in output and provide the strand-oriented sequence so it reads
5-prime to 3-prime in transcript orientation. Scan the configurable upstream
window and separately flag whether a match lies in the expected core window.
Retain all matches even though one primary class is selected for summaries.

For every treatment-control comparison:

1. Select before looking at condition effects a fixed set of adequately covered
   genes with at least two primary-discovery PACs.
2. Remove motif-assisted and `known_rescue_only` PACs, then renormalize usage
   across the remaining PACs within each sample and gene. Exclude a gene when
   fewer than two primary PACs remain or their retained total is zero.
3. For each sample and gene, sum the renormalized PAU for PACs belonging to each
   motif class.
4. Give every gene equal weight when averaging motif-class usage within a
   sample.
5. Transform each sample-level score with
   `asin(sqrt(score))`, which accepts exact zero and one without pseudocounts.
6. Fit a `limma` linear model for each motif class using the same condition and
   covariate design as the main model.
7. Test the treatment-control coefficient and adjust across motif classes
   within the comparison using BH.

This tests whether treatment shifts transcript usage toward PACs carrying a
particular motif class without allowing highly expressed genes to dominate.
Report untransformed group means and delta motif usage for interpretation,
along with transformed-scale coefficients and p-values. Report the number of
informative genes and do not test when it is below
`motif_preference_min_genes`. Repeat the descriptive summaries, but not the
primary significance claim, with `known_rescue_only` sites included as a
sensitivity analysis.

As a secondary exploratory analysis, test presence or absence of every upstream
k-mer among gained or increased-usage PACs against testable background PACs
from the same genes. Use a Cochran-Mantel-Haenszel test stratified by gene,
retaining only genes containing both event and background PACs and within-gene
variation for that k-mer. Report the common odds ratio, confidence interval,
raw PAC counts, informative-gene count, and BH-adjusted p-value. Label this
analysis exploratory because sequence context, PAC position, and event
selection can remain confounded.

### 10. Report Results

Produce a static HTML report with:

- input and reference preparation;
- sorting and indexing actions;
- condition-to-control mapping;
- fragment exclusions;
- proper-pair rates;
- strandedness evidence;
- protocol, evidence-source, and calibration compatibility;
- end-observation and PAC counts;
- proximal-kernel resolution and merged resolution groups;
- PAC width, peakiness, annotation, motif, and internal-priming summaries;
- motif-class usage by sample and condition;
- treatment-versus-control motif-preference tests;
- enriched upstream k-mers for gained or increased PACs;
- known-PAC overlap;
- testable-gene counts;
- PAU correlation and PCA for adequately covered genes;
- model convergence and p-value distributions;
- zero-boundary stabilization and bootstrap success summaries;
- searchable gained, lost, switched, and redistributed PAC events;
- top genes for every condition-versus-control test;
- per-gene plots of genomic PAC positions, raw counts, observed replicate PAU,
  and fitted PAU.

The report must remain usable after being copied off the cluster.

## Outputs

```text
results/
  manifest/
    resolved_params.yaml
    normalized_samples.tsv
    software_versions.tsv
    input_checksums.tsv
    run_manifest.json
  qc/
    input_validation.tsv
    reference_preparation.tsv
    alignment_preparation.tsv
    control_mapping.tsv
    library_calibration.tsv
    strandedness.tsv
    fragment_filtering.tsv
    pac_discovery.tsv
    quantification.tsv
  prepared_reference/          # optional publication
  prepared_alignments/         # optional publication
  evidence/
    SAMPLE.3prime_evidence.parquet
    SAMPLE.3prime_evidence.tsv.gz
  atlas/
    pacs.v1.bed.gz
    pacs.v1.metadata.tsv.gz
    rejected_candidates.tsv.gz
  counts/
    pac_counts.tsv.gz
    pac_counts.long.parquet
    gene_totals.tsv.gz
    observed_pau.tsv.gz
  statistics/
    FAMILY.gene_omnibus.tsv.gz
    CONDITION_vs_CONTROL.genes.tsv.gz
    CONDITION_vs_CONTROL.pacs.tsv.gz
    CONDITION_vs_CONTROL.events.tsv.gz
    fitted_pau.tsv.gz
    gene_precision.tsv.gz
  motifs/
    pac_motifs.tsv.gz
    CONDITION_vs_CONTROL.preference.tsv.gz
    CONDITION_vs_CONTROL.preference_known_rescue_sensitivity.tsv.gz
    CONDITION_vs_CONTROL.kmer_enrichment.tsv.gz
  tracks/
    SAMPLE.plus.3prime_evidence.bedGraph.gz
    SAMPLE.minus.3prime_evidence.bedGraph.gz
  report/
    index.html
```

BedGraph values are non-negative counts. Strand is encoded by the separate
plus/minus filenames rather than by signed values.

## Implementation Structure

```text
pacusage/
  main.nf
  nextflow.config
  nextflow_schema.json
  README.md
  pyproject.toml
  workflows/              # top-level workflow composition
  subworkflows/local/     # preparation, discovery, quantification, statistics
  modules/local/          # one Nextflow process per major operation
  conf/                   # base, local, Slurm, and test profiles
  src/pacusage/           # tested Python scientific logic
  bin/                    # small Nextflow-facing Python entry points
  scripts/                # R statistical model
  envs/                   # Conda environments
  tests/                  # unit, integration, simulation, and fixtures
```

Name modules after the workflow operations above, including
`PREPARE_REFERENCE`, `PREPARE_ALIGNMENT`, `CALIBRATE_LIBRARY_PROFILE`,
`INFER_STRANDEDNESS`, `EXTRACT_3PRIME_EVIDENCE`, `CLUSTER_PACS`,
`ANNOTATE_PACS`, `QUANTIFY_PACS`, `FIT_USAGE_MODEL`, and `BUILD_REPORT`.

## Coding Standards

Readability is a primary requirement.

- Use descriptive process, function, and variable names.
- Keep one major operation in each Nextflow module.
- Keep channel construction explicit.
- Put scientific logic in tested Python or R functions, not dense shell blocks.
- Use Python type hints, dataclasses where helpful, and short docstrings.
- Comment strand, coordinate, CIGAR, and statistical logic.
- Avoid deeply nested expressions and unnecessary metaprogramming.
- Make errors identify the sample, invalid value, expected format, and remedy.
- Keep example files small and realistic.

Prefer clarity over shorter code or small speed improvements. Still avoid major
runtime costs:

- stream BAM/CRAM records;
- use name-collated streams for PE reconstruction;
- aggregate endpoints instead of storing read names;
- use columnar tables for large operations;
- parallelize by sample or chromosome;
- benchmark code that touches every fragment;
- document any optimization that makes code less obvious.

## HPC Requirements

Use Nextflow DSL2 and its Slurm executor.

Provide:

- `local`, `slurm`, `conda`, `apptainer`, and `test` profiles;
- configurable Slurm account, partition, and QoS;
- process labels for low, medium, and high resource jobs;
- CPU, memory, runtime, and temporary-storage settings for every process;
- retry and resource escalation for transient or memory failures;
- configurable queue size and submission rate;
- one log per process;
- clean cancellation;
- reliable `-resume` behavior.

Parallelize sample-level extraction, alignment preparation, quantification, and
track creation. Parallelize atlas work by chromosome and strand where useful.
Run the initial statistical model in one sufficiently resourced R job unless
testing proves that chunking preserves global precision estimation.

Use the configured scratch location for temporary files and avoid large files
in `$HOME`. Never overwrite source inputs. Publish only completed outputs.

## Testing

### Unit Tests

Cover:

- unindexed FASTA preparation;
- sorted, unindexed, and unsorted BAM/CRAM preparation;
- source-file preservation;
- Plasmidsaurus profile defaults and calibration;
- exact-boundary and proximal-tag profile behavior;
- SE read-end and PE fragment-end coordinates on both strands;
- soft clips and splice-aware CIGAR operations;
- proper, improper, orphaned, interchromosomal, duplicate, secondary,
  supplementary, QC-failed, low-MAPQ, and multimapped fragments;
- forward, reverse, and ambiguous strandedness;
- control-condition validation and multiple comparison families;
- PAC seed ranking, radius boundaries, and deterministic ties;
- deterministic proximal kernels, peak ranking, minimum resolvable separation,
  regional peak merging, and assay-resolution intervals;
- calibrated proximal-tag assignment likelihood ratios and ambiguity handling;
- protocol evidence sources and rejection of incompatible comparison families;
- prevention of motif use during primary PAC discovery;
- exclusion of rescue-only PACs from primary motif testing;
- strand-aware sequence annotation;
- motif classification on both genomic strands;
- equal-gene motif-usage summaries and arcsine-square-root model inputs;
- motif-preference tests under known simulated shifts;
- overlapping and ambiguous gene assignments;
- count conservation and PAU sums;
- model matrix construction and confounding errors;
- minimum replicate validation for control and treatment conditions;
- contrast-specific gene screening and stage-wise PAC adjustment;
- deterministic all-zero stabilization, instability flags, and parametric
  bootstrap intervals;
- separate plus- and minus-strand browser tracks;
- YAML and command-line parameter precedence.

### Integration Fixture

Create a small synthetic FASTA, annotation, and alignment collection containing:

- an unindexed FASTA;
- sorted and unsorted BAM/CRAM files;
- separate successful runs for Plasmidsaurus-like SE, exact-boundary SE, and
  exact-boundary PE samples;
- an intentionally incompatible mixed-protocol run that must fail before atlas
  construction;
- DMSO controls referenced by TreatmentA and TreatmentB;
- at least one second comparison family;
- forward- and reverse-stranded samples;
- one-, two-, and three-PAC genes;
- a replicated treatment-only PAC;
- constant PAU with changing gene abundance;
- changing PAU with unchanged original-PAC count;
- internal-priming-like sequence;
- ambiguous gene assignment;
- proper and improper pairs.

The complete `test` profile must run this fixture and compare stable tables with
expected outputs.

### Statistical Simulation

Simulate null and non-null genes across:

- low and high counts;
- low and high overdispersion;
- balanced and unbalanced replicates;
- treatment PACs absent from controls;
- treatment preferences for a known PAS motif class;
- optional batch and donor covariates.

Verify:

- reasonable null p-value calibration;
- correct delta-PAU direction and magnitude;
- greater uncertainty at lower gene totals;
- unchanged APA results when every PAC count for a gene is scaled equally;
- recovery of replicated treatment-specific PACs;
- controlled null error after contrast-specific stage-wise adjustment;
- reproducible handling and explicit instability flags for all-zero group PACs;
- recovery of simulated motif-preference shifts without sensitivity to
  proportional gene-expression changes.

## Definition of Done

The work is complete when a beginner can provide a sample sheet, deduplicated
BAM/CRAM files, an indexed or unindexed FASTA, and an annotation, then run one
Nextflow command to obtain:

- prepared and validated inputs;
- a documented protocol profile and calibration result;
- a resolved evidence source and compatibility decision for every comparison;
- strandedness and fragment-filtering QC;
- a versioned PAC atlas whose coordinate precision and unresolved resolution
  groups are reported honestly;
- raw PAC counts and within-gene PAU;
- descriptive omnibus tests within each control family;
- contrast-specific gene tests and hierarchical PAC tests;
- every non-control condition compared with its declared control;
- fitted PAU, delta PAU, uncertainty, and FDR;
- a direct table of gained, lost, switched, and redistributed PAC events;
- per-PAC motif sequences and condition-level motif-preference results;
- browser tracks and a static report;
- enough metadata to reproduce the run offline.
