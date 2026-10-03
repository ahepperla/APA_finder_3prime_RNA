# PACusage design

## Goal

PACusage is a readable, reproducible command-line pipeline that:

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
5. A referenced condition may itself have a control, so comparisons can be
   nested, for example `WT -> disease_vehicle -> disease_drug`.
6. Every root control condition (blank `control`) is referenced by at least
   one other condition.
7. A condition cannot reference itself, and control relationships must be
   acyclic.
8. Every modeled control and treatment condition has at least two biological
   replicates by default.

Create a normalized `control_condition` field:

- control rows use their own `condition`;
- non-control rows use their declared `control`.

Conditions sharing a `control_condition` form a comparison family. For example,
`DMSO`, `TreatmentA`, and `TreatmentB` form the `DMSO` family. Multiple
comparison families are allowed. A condition that is both a treatment and a
control belongs to two families: it is a treatment in its control's family and
the control of its own family.

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

- a known-PAC BED file, whose records each give a PAC at their 3′ edge in
  transcript orientation (`end` on the plus strand, `start` on the minus
  strand);
- chromosome-name aliases, mapping alignment contig names to FASTA and
  annotation names;
- a PAS motif catalog replacing the built-in one.

The FASTA must be uncompressed and does not need to be indexed. A
`PREPARE_REFERENCE` process must:

1. Link the FASTA into the work directory and index the link with
   `samtools faidx`, so the index is written in the work directory.
2. Never modify the source FASTA or its directory.
3. Emit the prepared FASTA and index together.
4. Fail clearly if the FASTA is compressed or not indexable.

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
excluded_contigs: [chrM, MT, chrMT]   # mitochondrial contigs, by default

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
# Whole numbers are sample counts; values in (0, 1) are within-condition
# fractions rounded up (for example, 0.5 requires 3 of 5 samples).
pac_min_supporting_samples: 2
known_pac_rescue_total: 5
known_pac_match_radius: 12
proximal_kernel_overlap_threshold: 0.50
proximal_assignment_likelihood_ratio: 3.0
proximal_bin_size: 25
internal_exon_end_filter: true
constitutive_readthrough_filter: true
constitutive_readthrough_min_junction_count: 2

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
min_site_usage_samples: null
min_site_usage_dropouts: 0
min_site_usage_gene_reads: 10
min_test_supporting_samples: 2

model_covariates: []
gene_fdr: 0.05
site_fdr: 0.05
min_abs_delta_pau: 0.10
apa_pattern_min_change: 0.10
dm_bootstrap_replicates: 200
dm_bootstrap_min_success_fraction: 0.80
dm_zero_sensitivity_repeats: 5
dm_zero_max_delta_pau_spread: 0.02
random_seed: 1729
event_min_treatment_pau: 0.05
active_pac_min_pau: 0.10
event_max_control_pau: 0.01
event_min_supporting_samples: 2
potential_internal_priming_withheld_calls: true
motif_preference_min_genes: 50
motif_kmer_length: 6
run_kmer_enrichment: true

save_prepared_reference: false
save_prepared_alignments: false
save_intermediates: false
```

Defaults must permit a standard run with a shorter YAML containing only
`input`, `outdir`, `assembly`, `fasta`, and `gtf`.

`nextflow_schema.json` is the single definition of every parameter: its type,
range, default, and description. `VALIDATE_INPUTS` resolves the parameters
against it once, and every later step reads the resolved file.

`insufficient_replicates_policy` accepts `error` or `warn`. The default
`error` stops before modeling. `warn` permits an explicitly exploratory run.
A condition with fewer than `min_replicates_per_condition` samples then makes
every comparison of each family that includes it exploratory, not only the
comparisons that use it. `exploratory_insufficient_replicates` marks these
comparisons in their `.genes`, `.pacs`, and `.calls` tables and in the
family's `gene_omnibus` table, and their gained and lost PACs are reported
only as candidates.

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

A reused alignment is a symbolic link to the source, not a copy, and a
reusable index is linked the same way. Hash each source once, and record its
size and modification time. Every later step that reads a linked alignment
checks both before and after its pass. The prepared FASTA is also a link;
validate or generate its index beside the link in the work directory, and
pass the FASTA and index together to every step that needs them.

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

Describe the selected source's pooled kernel in
`qc/calibration_kernel_diagnostics.tsv`. Read ends around one site give one
mode. For a proximal-tag run, warn when the kernel, smoothed with a centred
moving average as wide as its minimum resolvable separation (rounded down to
an odd width), has more than one separated mode, or when the median
per-sample central interval width exceeds six times that separation. A mode
reaches at least a quarter of the highest one, and neighbouring modes merge
unless the valley between them drops below half the lower one. A kernel of
PAC spacings fails these checks, for example from exact-boundary reads under
a proximal-tag profile or from unannotated alternative polyadenylation near
annotated ends. The warning goes to the Nextflow log and the report and does
not stop the run. Exact-boundary runs record the same values with status
`not_applicable`.

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

Use annotated exonic sequence that belongs unambiguously to one gene, with no
opposite-strand overlap, on contigs that are not excluded. Reservoir-sample
eligible alignments so coordinate order does not bias the
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
- otherwise: fail as ambiguous unless the user or the library profile supplied
  a value.

Stop the run when a sample's inferred strandedness conflicts with the value
its profile or the user requested. Samples may differ from one another; each
is resolved on its own evidence.

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

Read each alignment once for calibration and extraction. `SCAN_ALIGNMENT`
makes one pass, with one name-collate for paired-end data, and collects for
every candidate evidence source:
- aggregated observations;
- filtering counters, kept per source in first-seen order;
- direct splice continuations, when the readthrough filter can apply.

`EXTRACT_3PRIME_EVIDENCE` then writes the resolved source's evidence from
that scan.

### 5. Build the PAC Atlas

Build one atlas from all samples. Find candidates in every sample's read ends
pooled, without CPM, the treatment–control contrasts, effect sizes, or PAS
motifs. Condition labels are used only to require replicate support within
the condition that supports a candidate, and in the readthrough rule below,
which asks each condition that supports a candidate.

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
3. Bin observed endpoints and the pooled calibration kernel at the smaller of
   `proximal_bin_size` and the minimum resolvable separation. For every
   implied regional PAC coordinate \(x\), calculate the matched score
   \(L(x)=\sum_{e,s} \min(y_{es},3)K(e-x)\), where \(e\) is an observed
   endpoint, \(y_{es}\) is its count in sample \(s\), and \(K\) is the pooled
   kernel.
4. Identify local maxima of the matched score. Rank maxima by matched score,
   then genomic coordinate.
5. Apply the same total-count, per-sample-count, and supporting-sample
   requirements used for exact-boundary candidates. The supporting-sample
   threshold must be met by replicates within at least one condition; samples
   from unrelated conditions cannot be combined to admit a candidate.

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
representative coordinate, retain the merged maxima, and report the
resolution region, of width `resolution_nt` around the representative, as the
PAC's `start` and `end`. Never report unresolved
maxima as independently quantified PACs or interpret the representative as a
nucleotide-resolution cleavage site.

#### Proximal-Tag Filters

Reads that end inside a spliced transcript are not transcript ends. Two
filters reject the candidates they make, before annotation, and record each
in `rejected_candidates.tsv.gz` with its reason. Both locate a candidate's
read pile: a sharp pile of read ends peaks `shift` bins downstream of it,
where `shift` is the binned kernel's first bin plus the offset of its highest
bin. So a candidate's pile lies `shift` bins upstream of it.

- **Internal exon ends** (`internal_exon_end_filter`): reads that cross a
  splice junction, but whose short overhang into the next exon is not
  aligned, end at the donor.
  - The donors are the 3' ends, in transcript orientation, of every exon a
    transcript of the annotation splices onward from. Genes without
    transcript IDs have none.
  - A candidate within one bin of a donor's bin plus `shift` is rejected as
    `internal_exon_end`. That reaches one and a half bins from the donor, so
    donors within two bins of any refined terminal exon's 3' end on the
    strand are left out, and a real end there is never rejected.
  - The cost: a real site whose reads pile up that close to an internal
    donor is rejected too, whatever its reads show.
- **Readthrough** (`constitutive_readthrough_filter`): continuation reads are
  direct CIGAR splices out of an aligned exon block that holds the pile,
  within one bin of either end of the block, on either strand.
  - Discovery records each candidate's supporting conditions, those whose
    replicates pass the support rule, as `supporting_conditions`.
  - A candidate is rejected as `constitutive_readthrough` when, in every
    supporting condition, the continuation reads pooled over that
    condition's samples reach `constitutive_readthrough_min_junction_count`.
  - A supporting condition without them keeps the candidate, since the exon
    may end there in that condition. Conditions that do not support the
    candidate are not consulted, so their coverage cannot protect it.
  - The cost: a real site is rejected too when its read pile lies inside an
    exon that every supporting condition splices onward from, within an
    aligned block's reach of the donor, about one read length. That covers
    an internal-exon site used by a minority of transcripts, and a site just
    past the donor when the kernel peaks far downstream of the reads. Sites
    further into the intron are untouched.

For `proximal_tag`, known transcript ends may be used for calibration and
annotation, but PAS motif sequence must not be used to select or reposition
primary candidates. This prevents circular motif-preference results. PACusage
has no motif-assisted rescue pass; one added later must keep its rescued sites
out of motif-preference testing.

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
`known_rescue_only`. Exclude these sites from primary motif-preference and k-mer
analyses. Produce a sensitivity
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

- `chrom`, `start`, `end`: the PAC's interval in the BED convention, the same
  as its atlas BED record: `[coordinate, coordinate + 1)` for an exact PAC and
  the resolution region for a proximal-tag PAC; and `locus`, the same interval
  1-based for genome browsers;
- `endpoint_model`;
- `coordinate_precision`: `exact` or `estimated`;
- assay-resolution width;
- merged candidate coordinates;
- calibration profile and version;
- resolved evidence source;
- primary evidence type (`primary_evidence_type`): `exact_boundary` for an
  exact PAC and `endpoint_peak` for a proximal-tag PAC.

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
- gene assignment (`gene_id`, `gene_name`), with `ambiguous_gene_assignment`
  marking a PAC assigned to several genes, but no transcript assignment;
- gene region (`gene_region`): `last_exon`, `internal_exon`, `intron`,
  `downstream_of_gene`, or `intergenic`;
- proximal-to-distal rank;
- confidence tier (`confidence`), recorded without reasons: `low` for a PAC
  assigned to several genes or flagged for internal priming; otherwise `high`
  for a PAC that matches a known PAC, has at least two supporting samples,
  and is not `known_rescue_only`; and `moderate` for the rest.

Internal-priming evidence is a flag by default, not an automatic exclusion.
Strict filtering may be exposed as an optional parameter.

Assign genes in this order:

1. Same-strand terminal exon (`last_exon`).
2. Same-strand other exon (`internal_exon`).
3. Same-strand intron (`intron`).
4. Downstream of a transcript end within the configured distance, without
   crossing an intervening same-strand gene (`downstream_of_gene`).
5. Intergenic (`intergenic`).

**Terminal exons.** A transcript's final exon is a terminal exon unless:
- it overlaps an internal exon of another transcript of the gene, as the
  final exons of retained-intron and 3'-incomplete models do; or
- it is a single-exon model apart from every exon of the gene's multi-exon
  transcripts, such as a fragment inside an intron.

PACs in such exons are `internal_exon` PACs. Without transcript IDs, a gene
is one transcript, so only its 3'-most exon is terminal.

The rule cannot tell a composite terminal exon, an internal exon extended into
its intron, from a canonical last exon that a minor isoform splices through.
Both have a final exon containing another transcript's internal exon, so both
count as internal. In genes of the second kind, tandem 3' UTR changes
therefore read as intronic patterns. Telling them apart would need transcript
biotypes or tags, which the pipeline does not use.

The single-exon rule compares a model with every exon of the gene's multi-exon
transcripts. So a single-exon model lying 3' of all of them is excluded, and a
PAC in it counts as upstream region. An intronic single-exon model that
overlaps a retained-intron transcript's final exon is kept, and can form a
false last exon.

**Last exons.** A gene's terminal exons that overlap, sharing even a single
coordinate, form one last exon. The atlas records each terminal-exon PAC's
last exon as a 1-based locus in `last_exon_locus`. A PAC downstream of a gene gets
the gene's 3'-most last exon. Every other PAC gets none.

Retain ambiguous assignments in the atlas but exclude them from primary
within-gene testing.

Name each gene once: its gene record's name (`gene_name`, a GFF3 gene's
`Name`, or NCBI's `gene`), else the first name on its other records, else its
`gene_id`. A PAC assigned to several genes lists their names in the order of
their IDs. Every published table carries `gene_name` next to `gene_id`.

### 7. Quantify the Frozen Atlas

Recount every sample after the atlas is frozen.

For exact-boundary profiles, define nonoverlapping PAC territories by:

- capping assignment distance at `pac_cluster_radius`;
- assigning each end observation to the nearest representative, which splits
  nearby territories at their midpoint, with a tie going to the lower
  coordinate;
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
support in >= 2 samples
PAC usage >= 1 percent in as many samples as the smallest condition has,
  counting samples with >= 10 gene reads
```

Every PAC filter uses the family's samples without regard to condition. The usage
filter counts the samples in which a PAC has `min_site_usage` of the gene's
reads, from any of the family's conditions, and needs `min_site_usage_samples`
of them.
- `min_site_usage_samples` defaults to the size of the family's smallest
  condition, as in the DRIMSeq workflow (Love et al., 2018), which asks for a
  10% share in at least as many samples as the smallest group. It is never
  more than the family's samples.
- A sample counts only when the gene has `min_site_usage_gene_reads` reads
  in it.
- `min_site_usage_dropouts` (default 0) lowers the default count by that
  many samples, to allow for bad replicates. It never takes the count below
  2, or below the smallest condition if that has fewer, so it never makes
  the rule stricter; it can't be combined with `min_site_usage_samples`.
  - It departs from standard practice at small replicate numbers. The
    DRIMSeq workflow and edgeR's `filterByExpr` ask for the full smallest
    group size, and edgeR relaxes that only when every group has more than
    10 samples, keeping at least 70% of the smallest group.
  - It stays blind to condition, since it too uses only condition sizes. In
    S7's design at one dropout (2 of 11 samples, three seeds), the filter
    admitted 84% to 85% of the borderline null PACs, 5% to 7% of those had
    p <= 0.05, and 3% to 9% of the treatment's gene calls were false.
- **Why.** A share pooled over the family dilutes a site used in one
  condition, more so the more conditions share the control, and a share in
  one sample can be noise. Every condition has at least as many samples as
  the smallest, so a site used throughout any one condition passes, however
  many conditions share the control. A site that reaches the share in fewer
  samples, such as one noisy sample, does not.
- **Independence.** The filter uses the condition sizes, never which sample
  belongs to which condition, so permuting the labels leaves it unchanged.
  Under a null with exchangeable labels it is then independent of
  label-permutation statistics, and approximately of the
  Dirichlet-multinomial test: the condition under which filtering leaves
  type I error control intact (Bourgon et al., 2010). Simulation S7 checks
  this at `min_site_usage` 0.10 with
  `min_abs_delta_pau` 0.05, for null PACs at 8% usage:
  - the filter admitted 60% to 65% of them;
  - 5% to 8% of those admitted had p <= 0.05, against a nominal 5%;
  - with 100 shifted genes among 500, 2% to 9% of the treatment's gene calls
    were false;
  - null genes keeping a borderline PAC were called about as often as other
    null genes, 1.0% against 0.8%.

  Null genes that kept the borderline PAC had p <= 0.05 a little more often
  than those that dropped it, about 6.5% against 3.7%. That is composition,
  not selection: the filter keeps fewer of the deep, precise genes, whose
  null p-values are rarely small. Within each depth and precision stratum
  the two groups' rates agree.

Record, for every PAC and family, whether it was tested and why not, in
`statistics/FAMILY.statistical_filtering.tsv.gz`; PACs without a gene or with
an ambiguous gene assignment are never tested.

After the family fit, each comparison tests a gene only when its control and
its treatment each have depth: `event_min_supporting_samples` samples, or all
of a smaller group, with `min_site_usage_gene_reads` reads at the gene's
tested PACs. A nearly silent group's handful of reads would otherwise read as
a large shift in usage. The rule counts reads because the uncertainty of a
usage estimate depends on reads, not on the gene's share of the library.
`min_group_gene_cpm` (default 0, off) optionally adds an expression floor: a
sample then also needs that many reads per million of its assigned reads at
the gene. Library sizes are fixed, so the floor is still a function of the
gene totals and keeps the rule independent of the usage test.
- **Independence.** Unlike the PAC filters, this rule reads the condition
  labels, but it reads only per-sample gene totals, never how they split among
  PACs. The Dirichlet-multinomial test conditions on those totals, so under
  the null the rule is independent of the test (Bourgon et al., 2010), as far
  as the test's asymptotic null holds. The genes kept have the same p-values
  as before; only the set the comparison's corrections run over changes.
  Simulation S8 checks it.
- The family fit, its omnibus test, precision, and zero-count stabilization
  are unchanged; only the comparison's tables leave the genes out.
- The genes a comparison skips are in
  `statistics/CONDITION_vs_CONTROL.genes_without_depth.tsv.gz`. A gene is
  `turned_off` when the control has depth and fewer treatment samples than a
  group needs have any read, although at the control's mean CPM the
  treatment's libraries would have given it depth; `turned_on` is the mirror
  image. That is a description, not a test: not detected, at a depth where it
  would have shown. The rest are `too_low_in_treatment`, `too_low_in_control`,
  or `too_low_in_both`.
- Reads here are at the gene's tested PACs, the model's totals; the usage
  filter's `min_site_usage_gene_reads` counts a sample's reads at all of the
  gene's PACs, so the two can differ.

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
- raw replicate counts, and observed PAU over the gene's PACs tested in the
  comparison's family, so it can differ from `counts/observed_pau.tsv.gz`;
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
matrix first. A gene is a boundary gene when its unmodified fit is non-finite
anywhere, or when a PAC has no counts across a condition or covariate level.
Stabilize boundary genes as a model-fitting fallback, once per comparison
family:

1. In each repeat, replace the boundary genes' zero counts with draws from
   U(0, 0.1), DRIMSeq's `addUniform` rule. Seed each gene's draws from
   `random_seed`, atlas checksum, family, gene ID, and repeat number.
2. Fit precision, the full model, and every null model to the same perturbed
   counts. DRIMSeq's own `add_uniform` flag is not used, because `dmTest`
   refits the null model to the original counts.
3. Estimate precision across the whole family with the unmodified fit's common
   precision, so the moderation matches the unmodified fit.
4. Keep PAC ordering deterministic by genomic coordinate.
5. Repeat the stabilized fit `dm_zero_sensitivity_repeats` times. Write back
   medians of the repeats: proportions, precision, and likelihood ratios, with
   p-values from the median likelihood ratio. A gene needs at least 80 percent
   of repeats to succeed.
6. Record that stabilization was used (`model_status`
   `fitted_with_zero_count_stabilization`; an unstabilized fit is `fitted`).
   Never write the perturbed values to a count or PAU output.
7. Flag the PAC as `zero_boundary_unstable`, with a reason, when any of these
   holds:
   - the effect direction changes beyond ±0.001;
   - fewer than 80 percent of fits converge;
   - fitted delta PAU spans more than `dm_zero_max_delta_pau_spread`.
8. Do not assign a PAC-level p-value when stable finite fits cannot be obtained.
   Keep the detection evidence and label the inferential result unavailable
   rather than inventing a value.
9. When a condition has no counts for a gene, mark the family-wide test
   `group_without_counts`. Every comparison using that condition lacks depth
   for the gene, so it leaves the gene untested and lists it in
   `genes_without_depth`.

For PACs in genes that pass the contrast-specific screen or satisfy the
pre-statistical gained/lost event criteria, estimate delta-PAU uncertainty with
a parametric bootstrap conditional on the family-fit precision:

1. Simulate sample count vectors from the fitted Dirichlet-multinomial model at
   each sample's observed gene total and design row.
2. Refit the proportions with each gene's precision held at the family-fit
   estimate. Perturb zero groups exactly as in the family fit. Recompute fitted
   treatment-control delta PAU for every selected comparison of the family from
   the same draw.
3. Run `dm_bootstrap_replicates` replicates with deterministic seeds per gene
   and replicate, so results do not depend on batching or worker count.
4. Report percentile 95 percent intervals, bootstrap success counts, and a
   `bootstrap_status`. The intervals do not include uncertainty in the
   precision. In simulations with four replicates per group, nominal 95
   percent intervals covered the true change 80-93 percent of the time (about
   88 percent on average).
5. Leave the interval unavailable and flag the reason when the success
   fraction is below `dm_bootstrap_min_success_fraction`.
6. `dm_bootstrap_include_candidates` defaults to `true`; setting it to `false`
   retains intervals for genes passing the contrast-specific gene FDR screen
   while omitting them for non-significant pre-statistical candidates.

#### Researcher-Facing PAC Events

Create one gene-level event summary for every treatment-control comparison.
Classify PACs using both statistical evidence and detection evidence:

- **gained:** not detectably used in adequately covered controls, reproducibly
  detected in treatment, positive delta PAU, and significant PAC-level FDR;
- **lost:** the reverse of gained;
- **increased usage:** detected in both groups, with a significant delta PAU
  of at least `min_abs_delta_pau`;
- **decreased usage:** detected in both groups, with a significant delta PAU
  of at most minus `min_abs_delta_pau`;
- **dominant switch:** in a gene that passes the contrast-specific gene test,
  the highest fitted-usage PAC differs between treatment and control;
- **more or fewer active PACs:** in a gene that passes the contrast-specific
  gene test, the number of active PACs differs between the groups. An active
  PAC has fitted PAU of at least `active_pac_min_pau` in a group. The default,
  0.10, matches `min_abs_delta_pau`, so a PAC whose usage falls from active to
  none changes by enough to be called lost.

The first four are PAC calls, a PAC's `event_type`. The last two are gene
events, and only the `.genes` table records them, in `dominant_switch` and
`active_pacs_change` (`more`, `fewer`, or `none`), after `control_condition`,
followed by the APA pattern and the gene's shift:

```text
gene_id  gene_name  condition  control_condition  dominant_switch  active_pacs_change
apa_pattern  shift_direction  shift_from_pac_id  shift_from_gene_region
shift_from_event_type  shift_to_pac_id  shift_to_gene_region  shift_to_event_type
delta_intronic_share  delta_utr_distal_share  last_exon_switch
gene_fdr  gene_pvalue  gene_likelihood_ratio  gene_degrees_of_freedom  model_status
stabilization_successes  exploratory_insufficient_replicates
```

#### Shifts

Name each gene's change in usage once. PAU are shares of a gene, so usage
gained at one PAC is lost at others, and one shift gives PAC calls in both
directions. Counting calls therefore counts every shift twice. A gene with a
confirmed PAC call has a shift:

- **to:** its confirmed gain or increase with the largest change in fitted
  PAU;
- **from:** its confirmed loss or decrease with the largest fall;
- **a side without a confirmed call:** the gene's other PAC with the largest
  fitted change that way. Fitted PAU sum to 1, so a confirmed change on one
  side always has usage moving the other way. The PAC's own `event_type`,
  `none` or a candidate, shows that it has no confirmed call;
- **ties:** the lower `pac_fdr`, then the PAC ID;
- **direction:** `distal` when the to-PAC lies 3' of the from-PAC in
  transcript orientation, and `proximal` otherwise.

The shift describes the calls and makes none, so it needs no test of its
own. The figures read it to count and label each gene once.

#### APA Patterns

Classify each gene's change by where its usage moved. This follows the split
between 3' UTR APA, among tandem sites of a last exon, and upstream-region APA,
at sites in introns or internal exons (Tian and Manley, 2017).

**Regions.**
- **Upstream region:** `intron` and `internal_exon` PACs.
- **3' region:** `last_exon` and `downstream_of_gene` PACs, grouped by
  `last_exon_locus`.

**Gating calls.**
- A gating call is a confirmed call (`gained`, `increased_usage`, `lost`,
  `decreased_usage`) on a PAC that is not low confidence.
- A flagged call is a confirmed call on a PAC flagged for possible internal
  priming. With `potential_internal_priming_withheld_calls`, a flagged PAC's
  `gained_candidate` or `lost_candidate` counts too when classifying it with
  the flag lifted makes it `gained` or `lost`: the flag alone withheld it.
- Ambiguous PACs are never tested, so tested low-confidence PACs are exactly
  the flagged ones.
- The threshold is `apa_pattern_min_change`, with 1e-9 of slack for sums
  of fitted proportions. It is set apart from `min_abs_delta_pau`, so PAC
  calls can take a smaller change than patterns.

**Numbers,** reported for every gene with fitted usage:
- `delta_intronic_share`: the change in the upstream region's share of the
  gene's fitted usage.
- `last_exon_switch`: the smaller of the largest gain and the largest loss
  in any last exon's share of the gene. It needs at least 2 last exons.
- `delta_utr_distal_share`: the change in the distal PAC's share of the main
  last exon, the 3' counterpart of DaPars' PDUI.
  - The main last exon is the one with the most usage over both groups; a
    tie goes to the 3'-most.
  - The number needs 2 or more tested PACs in that exon, and at least
    `event_min_treatment_pau` of the gene's usage there in both groups.

**Patterns.** A gene is classified when it passes `gene_fdr`, has fitted
usage, and has a confirmed call or a counted withheld one.
- `intronic_gain`: `delta_intronic_share` of at least the threshold, and a
  gating increase in the upstream region. `intronic_loss` is the mirror
  image.
- `alternative_last_exon`: `last_exon_switch` of at least the threshold, a
  gating increase in the gaining last exon, and a gating decrease in the
  losing one.
- `utr_shortening`: `delta_utr_distal_share` of at most minus the threshold,
  and either a gating increase at a more proximal PAC of the main last exon
  or a gating decrease at its distal PAC. `utr_lengthening` is the mirror
  image.
- **Shifts into or out of the main last exon:** when the gene also has an
  upstream-region or last-exon shift of at least the threshold, all of the
  main last exon's PACs move one way. Only calls against that movement count
  then: increases if the exon's share falls, decreases if it rises.
- **Potential internal priming:** a pattern that holds only once flagged
  calls also count takes the suffix `_potential_internal_priming`. Every rule
  only gains patterns from more calls, so a gene has at most one form of each.
- `apa_pattern` joins the patterns found with `;`: the gated ones in this
  order, then the potential ones in the same order. It is
  `unclassified_change` for a classified gene with no pattern, and `none` for
  every other gene.

When the detection, fitted-PAU, and effect-size requirements for `gained` or
`lost` pass but any other requirement below fails, report `gained_candidate`
or `lost_candidate`. The other requirements are significance (the gene and PAC
FDRs), coverage, stability, confidence (neither low confidence nor an
internal-priming flag), and replication (a comparison that is not
exploratory). A `zero_boundary_unstable` PAC has no PAC-level p-value, so it
fails both significance and stability. These labels are detection-supported
but not statistically confirmed and must never be merged with significant
gained/lost counts.

Default gained-PAC requirements:

- the parent gene passes the relevant contrast-specific gene test;
- the control samples together have at least `min_gene_total` reads at the
  gene's PACs tested in the family (`control_gene_total`; for lost calls, the
  treatment samples' `treatment_gene_total`);
- fitted control PAU is at most `event_max_control_pau`;
- fitted treatment PAU is at least `event_min_treatment_pau`;
- treatment support meets `event_min_supporting_samples`;
- delta PAU is at least `min_abs_delta_pau`; at the defaults this is the
  stricter of the two usage rules, so fitted treatment PAU must be at least
  0.10;
- PAC-level FDR passes `site_fdr`;
- the PAC is not low confidence, `zero_boundary_unstable`, or flagged as likely
  internal priming;
- the comparison is not exploratory.

Use the label `gained` or `not detected in control`, not `de novo`, unless
independent evidence establishes biological absence in the control.

The `.pacs` table, and its `.calls` subset, start with the identity block of
every PAC-level table, then run from the answer to the evidence:

```text
pac_id  gene_id  gene_name  chrom  start  end  strand  locus
condition  control_condition  event_type
fitted_control_pau  fitted_treatment_pau  delta_pau  delta_pau_ci_low  delta_pau_ci_high
pac_fdr  gene_fdr  pac_pvalue  gene_pvalue  pac_likelihood_ratio  pac_degrees_of_freedom
gene_region  last_exon_locus  confidence  internal_priming_flag  known_pac  known_rescue_only
primary_pas_motif  primary_pas_motif_rna  primary_motif_class
control_/treatment_supporting_samples  control_/treatment_gene_total
raw_control_/treatment_counts  observed_control_/treatment_pau
effect_exceeds_threshold  dominant_pac_control/treatment
control_/treatment_active_pacs
model_status  precision  alpha_control  alpha_treatment
stabilization_successes  stabilization_delta_pau_spread
zero_boundary_unstable  zero_boundary_reason
bootstrap_status  bootstrap_successes  bootstrap_perturbed
exploratory_insufficient_replicates
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

Steps 1-5 score every sample once for the run; steps 6-7 test every
treatment-control comparison:

1. Select before looking at condition effects a fixed set of adequately covered
   genes with at least two primary-discovery PACs: genes with at least
   `min_gene_total` reads in every sample.
2. Remove `known_rescue_only` PACs, then renormalize usage
   across the remaining PACs within each sample and gene. Exclude a gene when
   fewer than two primary PACs remain or their retained total is zero.
3. For each sample and gene, sum the renormalized PAU for PACs belonging to each
   motif class, and separately for PACs with each primary hexamer.
4. Average each sum within a sample over the genes that have a PAC in that
   class, or with that hexamer, giving every gene equal weight. This mean is
   the sample's `motif_usage` in `motif_class_scores.tsv` or
   `motif_scores.tsv`. A gene without such a PAC is left out of the mean, not
   counted as zero, so each class and hexamer has its own informative genes
   (`informative_genes`) within the fixed set.
5. Transform each sample-level score with
   `asin(sqrt(score))`, which accepts exact zero and one without pseudocounts.
6. Fit a `limma` linear model using the same condition and covariate design as
   the main model: one fit per comparison, on the family's samples. Fit it
   twice, once with each primary hexamer as a row (`.preference`) and once
   with each motif class as a row (`.preference_class`). A row is tested in a
   comparison only when every sample of that comparison has at least
   `motif_preference_min_genes` informative genes for it.
7. Test the treatment-control coefficient and adjust across the rows of each
   table within the comparison using BH.

This tests whether treatment shifts transcript usage toward PACs carrying a
particular motif class without allowing highly expressed genes to dominate.
Report untransformed group means and delta motif usage for interpretation,
along with transformed-scale coefficients and p-values. Report the number of
informative genes and do not test when it is below
`motif_preference_min_genes`. Repeat the descriptive summaries, but not the
primary significance claim, with `known_rescue_only` sites included as a
sensitivity analysis.

As a secondary exploratory analysis, test presence or absence of every upstream
k-mer among gained or increased-usage PACs against the other PACs of the same
genes that were tested in that comparison. Use a Cochran-Mantel-Haenszel test stratified by gene,
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
- calibration kernel diagnostics, with a warning at the top of the report when
  a proximal-tag kernel fails them;
- end-observation and PAC counts;
- proximal-kernel resolution and merged resolution groups;
- PAC width, peakiness, annotation, motif, and internal-priming summaries;
- motif-class usage by sample and condition;
- treatment-versus-control motif-preference tests;
- enriched upstream k-mers for gained or increased PACs;
- known-PAC overlap;
- testable-gene counts;
- PAU correlation for adequately covered genes, with the PAU PCA figure;
- model convergence and p-value distributions;
- zero-boundary stabilization and bootstrap success summaries;
- searchable gained, lost, switched, and redistributed PAC events;
- top genes for every condition-versus-control test;
- per-gene plots of mean observed PAU per condition, for the genes with the
  most reads;
- the treatment-control figures below, embedded as PNGs.

The report must remain usable after being copied off the cluster.

`scripts/plot_usage_figures.R` draws the figures with ggplot2, as PDF and PNG,
into `figures/`. They show the calls in the `.pacs` and `.genes` tables and
make no calls of their own. For each comparison:

- **volcano:** each tested PAC's change in fitted PAU against its PAC-level
  p-value, with confirmed gains and losses colored and candidates open.
  `pac_fdr` is missing outside screened genes, so the raw p-value is the
  axis. Each gene is labeled once, at its shift's to-PAC, or at its from-PAC
  when only that side has a confirmed call. Up to 20 genes on each side are
  labeled, the most significant first and the larger change first among
  equal p-values.
- **distal usage:** each tested gene's distal PAC, its most 3' tested PAC in
  a last exon or downstream, with its fitted usage in the control against
  the treatment. `direction` is that PAC's own confirmed call (`distal_up`,
  `distal_down`, their candidates, or `none`), so it carries stageR's error
  control. Points are colored by the gene's APA pattern, filled when the
  distal PAC has a call, and drawn as diamonds when only flagged PACs support
  the pattern. Up to 20 genes whose distal PAC rose and 20 whose
  distal PAC fell, the largest changes, are labeled by gene name.
  `CONDITION_vs_CONTROL.distal_usage.tsv.gz` lists the genes.
- **shifts by gene region:** each gene with a shift, once, by the regions
  its usage moved from and to. Proximal shifts are drawn left of zero and
  distal ones right. The rows are the region pairs that shifts take in any
  comparison, so the comparisons' figures line up.

Across comparisons:

- **event counts:** `event_counts` shows the genes by their PAC calls and
  the gene events, then the genes each comparison left untested for depth
  (turned off, turned on, or too low to test), then the genes per APA
  pattern.
  - Each gene with a PAC call counts once by its calls: a gained and a lost
    PAC, a gained PAC, a lost PAC, changes in usage only, or candidate calls
    only.
  - The genes per APA pattern are counted in three panels: with all tested
    PACs, unflagged and flagged; for the patterns that unflagged PACs
    support; and for those that only flagged PACs support. A gene never has
    both forms of one pattern, so each count in the first is the sum of the
    other two.
  - The patterns that only flagged PACs support keep their pattern's color,
    since a lighter tint of intronic gain would match intronic loss.
- **effect against coverage:** `effect_vs_coverage` shows each PAC's change in
  PAU against the reads at its gene's tested PACs in the less-covered group.
- **shared APA patterns:** `apa_pattern_grid` draws genes with an APA pattern
  in at least two comparisons against the comparisons. A `*` marks a cell
  whose pattern only flagged PACs support.
  - At most 50 genes are drawn, most shared first, then by best gene FDR.
  - `apa_patterns_by_comparison.tsv.gz` has every tested gene.
- **concordance:** `concordance` plots, for each pair of comparisons, each
  shared PAC's change in PAU in one against the other.
  - The panels form a matrix: a pair's earlier comparison names its column
    and the later one its row, so the pairs fill the lower triangle.
  - It has at most 15 panels; beyond that, only pairs that share a control or
    a condition are drawn.
  - `concordance_matrix` shows Pearson r for every pair, in the same cells.
  - `concordance.tsv.gz` gives each pair's relation (`shared_control`,
    `chained`, or `unrelated`), shared PACs, PACs called in both, and r.
  - Comparisons that share a control correlate positively through its
    estimate. Chained ones estimate the middle condition from the same
    samples with opposite signs, so they correlate negatively. r is therefore
    read with the relation.
- **PAU PCA:** `pau_pca` places the samples on the first two principal
  components of observed PAU.
  - It uses the genes with at least `min_gene_total` reads in every sample,
    and their PACs observed in every sample, with no zero-filling.
  - Each PAC is centered across samples. Each component's sign makes its
    largest loading positive, the first PAC by ID breaking ties.
  - Samples are labeled by ID. When every ID starts with its condition's
    name and a separator, the labels leave the name out.
  - `pau_pca.tsv` has the coordinates.

These figures set each comparison's own results side by side. They make no
treatment-versus-treatment test, so every call remains a comparison with the
direct control.

PDFs are written without dates or a producer, so reruns reproduce them byte
for byte. ggrepel places the labels: its search starts from a fixed seed and
has no time limit, which would otherwise make the layout depend on the
machine's speed. Uncalled PACs are drawn as a density above 5,000 per panel,
on a log scale, which keeps real-data figures small.

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
    calibration_kernel_diagnostics.tsv
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
    FAMILY.statistical_filtering.tsv.gz
    CONDITION_vs_CONTROL.genes.tsv.gz
    CONDITION_vs_CONTROL.genes_without_depth.tsv.gz
    CONDITION_vs_CONTROL.pacs.tsv.gz
    CONDITION_vs_CONTROL.calls.tsv.gz
    fitted_pau.tsv.gz
    gene_precision.tsv.gz
  motifs/
    pac_motifs.tsv.gz
    motif_scores.tsv
    motif_scores_known_rescue_sensitivity.tsv
    motif_class_scores.tsv
    motif_class_scores_known_rescue_sensitivity.tsv
    CONDITION_vs_CONTROL.preference.tsv.gz
    CONDITION_vs_CONTROL.preference_known_rescue_sensitivity.tsv.gz
    CONDITION_vs_CONTROL.preference_class.tsv.gz
    CONDITION_vs_CONTROL.preference_class_known_rescue_sensitivity.tsv.gz
    CONDITION_vs_CONTROL.kmer_enrichment.tsv.gz
    kmer_enrichment_status.tsv
  tracks/
    SAMPLE.plus.3prime_evidence.bedGraph.gz
    SAMPLE.minus.3prime_evidence.bedGraph.gz
  figures/
    CONDITION_vs_CONTROL.volcano.pdf/.png
    CONDITION_vs_CONTROL.distal_usage.pdf/.png
    CONDITION_vs_CONTROL.distal_usage.tsv.gz
    CONDITION_vs_CONTROL.shifts_by_gene_region.pdf/.png
    event_counts.pdf/.png
    effect_vs_coverage.pdf/.png
    apa_pattern_grid.pdf/.png
    apa_patterns_by_comparison.tsv.gz
    concordance.pdf/.png
    concordance_matrix.pdf/.png
    concordance.tsv.gz
    pau_pca.pdf/.png
    pau_pca.tsv
  report/
    index.html
  pipeline_info/               # Nextflow's reports for the latest run
    execution_report.html
    execution_timeline.html
    execution_trace.txt
    pipeline_dag.html
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
  conf/                   # resources, Slurm, and test profiles
  src/pacusage/           # tested Python scientific logic
  bin/                    # the pacusage command Nextflow tasks run
  scripts/                # R statistical model and figures
  containers/             # Apptainer recipe
  envs/                   # Conda environment
  tests/                  # unit, integration, R, pipeline, and fixtures
  docs/                   # this design and the decisions log
```

Name modules after the workflow operations above, including
`PREPARE_REFERENCE`, `PREPARE_ALIGNMENT`, `RECORD_INPUT_CHECKSUMS`,
`INFER_STRANDEDNESS`, `SCAN_ALIGNMENT` (per-sample calibration and evidence),
`AGGREGATE_CALIBRATION`, `EXTRACT_3PRIME_EVIDENCE`, `CLUSTER_PACS`,
`ANNOTATE_PACS`, `QUANTIFY_PACS`, `FIT_USAGE_MODEL`, `PLOT_FIGURES`, and
`BUILD_REPORT`.

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
- parallelize preparation, calibration, extraction, and quantification by
  sample, discovery by chromosome/strand internally, and differential usage
  modeling as independent Nextflow jobs by direct comparison family;
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
Fit each comparison family's model in one sufficiently resourced R job, so that
precision is estimated across the whole family, and scatter only the bootstrap
into independent batches.

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
- calibration kernel warnings for PAC spacings, separated modes, and wide
  spreads, and their absence for single-site kernels;
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
- YAML and command-line parameter precedence;
- gene-level dominant-switch and active-PAC events in genes whose PACs have
  their own calls, and PAC rows that carry only PAC calls;
- a column guide that documents exactly the columns of every published table;
- terminal exons without retained-intron, 3'-incomplete, or stray
  single-exon ends, and last exons merged by overlap, on both strands;
- APA patterns from region shares and direction-matched calls, including
  low-confidence PACs that cannot gate and ties between last exons;
- each pattern's potential-internal-priming form, confident patterns taking
  precedence, and withheld gains counted only when significant, not
  exploratory, and enabled;
- proximal-tag filters: internal exon donors on both strands, kept clear of
  transcript ends, and the kernel shift to a candidate's read pile; the
  readthrough rule pooled over each supporting condition, and inclusive at
  the donor on both strands;
- figure data: the distal PAC on both strands, gene-region and event counts,
  and dropped p-values and coverage rows;
- figures across comparisons: the shared-pattern grid's order and 50-gene
  cap, related pairs and the 15-panel cap, concordance calls and Pearson r,
  and the PAU PCA against a direct SVD, with its coverage rule and tied
  loadings;
- figures without dates, byte-identical across processes, for empty
  comparisons too.

### Integration Fixture

Create a small synthetic FASTA, annotation, and alignment collection containing:

- an unindexed FASTA;
- sorted and unsorted BAM/CRAM files;
- separate successful runs for Plasmidsaurus-like SE, exact-boundary SE, and
  exact-boundary PE samples;
- an intentionally incompatible mixed-protocol run that must fail before atlas
  construction;
- exact-boundary reads run under the Plasmidsaurus profile, which must warn at
  calibration;
- DMSO controls referenced by TreatmentA and TreatmentB;
- at least one second comparison family;
- forward- and reverse-stranded samples;
- one-, two-, and three-PAC genes;
- a replicated treatment-only PAC;
- constant PAU with changing gene abundance;
- changing PAU with unchanged original-PAC count;
- internal-priming-like sequence;
- ambiguous gene assignment;
- proper and improper pairs;
- a gene that gains an intronic PAC and a gene that switches between two
  alternative last exons (chr3), with the intronic PAC kept more than 1,000
  nt from any annotated end so calibration is unchanged.

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
  proportional gene-expression changes;
- calibration after the usage filter: how often it admits a null PAC near
  the usage threshold, whether those PACs' p-values stay calibrated, and the
  false share of calls in a family with real changes (S7).
- the depth each comparison needs: genes nearly silent in a group are not
  tested, the null genes kept stay calibrated, real shifts at normal depth
  are still called, and genes silent in one group are classed turned off or
  on (S8).

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
