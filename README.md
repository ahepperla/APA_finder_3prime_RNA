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
  figures/     per comparison: volcano, distal-usage, and site-class figures
               (PDF and PNG) and CONDITION_vs_CONTROL.distal_usage.tsv.gz;
               across comparisons: event counts, shared APA patterns,
               concordance, effect against coverage, and the PAU PCA, with
               apa_patterns_by_comparison.tsv.gz, concordance.tsv.gz, and
               pau_pca.tsv
  report/      index.html, a self-contained report with the figures embedded
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
`.events` tables keep only the PACs with an event. The `.genes` tables have
one row per tested gene: the gene-level events (`dominant_switch` and
`complexity_change`), the APA pattern and its numbers, then the gene FDR.

| Event | When |
|---|---|
| `gained`, `lost` | A PAC is detected (reads in at least two samples) in only one group, and that group's usage passes the fitted-PAU thresholds. The group without it has at least `min_gene_total` reads at the gene. The gene and the PAC pass their FDRs. The PAC is neither ambiguous nor flagged for internal priming. |
| `gained_candidate`, `lost_candidate` | The detection rules hold, but significance, stability, coverage, or confidence does not, or the comparison is exploratory. These are not confirmed calls. |
| `increased_usage`, `decreased_usage` | The PAC is detected in both groups, passes both FDRs, and its usage changes by at least `min_abs_delta_pau`. |
| `dominant_switch` | The most-used PAC differs between the groups, in a gene that passes `gene_fdr`. The label goes on the treatment's most-used PAC only when that PAC has no call of its own. |
| `complexity_gain`, `complexity_loss` | The number of PACs with fitted PAU of at least `event_min_treatment_pau` differs, in a gene that passes `gene_fdr`. The label goes only on the gene's PACs without a call of their own. |

A PAC has one `event_type`, so the gene-level labels undercount genes: a gene
whose PACs all have their own calls carries no label. Every gene with a
dominant switch or a complexity change is flagged in the `.genes` table's
`dominant_switch` (`TRUE` or `FALSE`) and `complexity_change` (`gain`,
`loss`, or `none`).

### APA patterns

Each `.genes` table's `apa_pattern` says how a gene's usage moved, following
the usual split between 3′ UTR APA, among tandem sites of a last exon, and
upstream-region APA, at sites in introns or internal exons (Tian and Manley,
2017). A gene gets a pattern only when it passes `gene_fdr`, and each pattern
needs a confirmed PAC call (`gained`, `increased_usage`, `lost`, or
`decreased_usage`) on a PAC that is not low confidence. The threshold is
`min_abs_delta_pau`.

| `apa_pattern` | When |
|---|---|
| `intronic_gain` | The share of the gene's usage at intronic and internal-exon PACs (`delta_intronic_share`) rises by at least the threshold, with a confirmed increase at one of them. This is premature intronic polyadenylation. |
| `intronic_loss` | That share falls by at least the threshold, with a confirmed decrease at one of them. |
| `alternative_last_exon` | One last exon gains and another loses at least the threshold of the gene's usage (`last_exon_switch`), each with a confirmed call in that direction. |
| `utr_shortening`, `utr_lengthening` | In the gene's main last exon, the distal PAC's share of that exon (`delta_utr_distal_share`, the 3′-end counterpart of DaPars' PDUI) falls or rises by at least the threshold. A confirmed call must point the same way: for shortening, an increase at a more proximal PAC or a decrease at the distal one. |
| `..._potential_internal_priming` | A pattern above that holds only once calls on PACs flagged for possible internal priming count, such as `intronic_gain_potential_internal_priming`. |
| `other` | The gene passes `gene_fdr` and has calls, but no pattern applies. |
| `none` | Anything else. |

- **Several patterns:** a gene can have more than one, joined by `;`, such as
  `intronic_gain;utr_shortening`.
- **Potential internal priming:** oligo-dT can prime on a genomic A-run, so a
  flagged PAC's change may mark more RNA over the A-run rather than a new
  cleavage site. 3′ tags can't tell the two apart.
  - A pattern that only flagged PACs support takes the
    `_potential_internal_priming` suffix and comes after every pattern
    without it. A gene has at most one form of each pattern.
  - With `potential_internal_priming_withheld_calls` (the default), a
    flagged PAC's gain or loss that the flag alone withheld also counts:
    classified with the flag lifted, it would be a confirmed `gained` or
    `lost` call. So a flagged site that appears from nothing, such as a new
    intronic site, can carry the pattern. Set it to `false` to count only
    flagged PACs' confirmed calls.
- **Shifts into or out of the last exon:** when a gene's usage also moves
  into or out of its main last exon, every PAC there moves the same way. So
  only calls in the other direction can show a UTR shift: after an intronic
  gain, only increases in the last exon count.
- **The numbers:** the three columns are given for every gene with fitted
  usage, screened or not, so genes can be ranked by effect size.
  `delta_utr_distal_share` is empty when the main last exon has fewer than 2
  tested PACs, or less than `event_min_treatment_pau` of the gene's usage in
  either group.
- **The main last exon** is the one with the most usage over both groups; a
  tie goes to the 3′-most.

**Where PACs sit.**
- **Last exons:** the atlas column `last_exon` names the last exon of each
  terminal-exon PAC, as a 1-based locus. Overlapping last exons of a gene's
  transcripts merge into one. A PAC downstream of a gene belongs to the
  gene's 3′-most last exon.
- **Terminal exons:** a transcript's final exon counts as a terminal exon
  unless:
  - it overlaps an internal exon of another transcript of the gene, as the
    final exons of retained-intron and 3′-incomplete models do; or
  - it is a single-exon model apart from every exon of the gene's multi-exon
    transcripts, such as a fragment inside an intron.
  PACs in such exons are `other_exon`, so they count toward the intronic
  share.

**Caveats.**
- **Spliced-through last exons.** A final exon that another transcript of
  the gene splices through counts as internal, so its PACs are upstream
  region.
  - This rightly treats a composite terminal exon, an internal exon extended
    into its intron and ended there, as intronic polyadenylation.
  - It also demotes a canonical last exon whenever a minor isoform continues
    past it, into a downstream exon or through an intron in its 3′ UTR. In
    such genes a tandem 3′ UTR change reads as `intronic_gain` or
    `intronic_loss`.
- **Single-exon models.**
  - A single-exon model within a gene, apart from its spliced transcripts,
    isn't a terminal exon even when it lies 3′ of all of them. A PAC in it
    counts as upstream region. Without the model, the PAC would be
    downstream and join the 3′-most last exon.
  - An intronic single-exon model that overlaps a retained-intron
    transcript's final exon is kept. It then forms its own last exon, and
    can give a false `alternative_last_exon`.
- **Same-direction changes under a shift.** When usage moves into or out of
  the main last exon, a within-exon change in the same direction can't be
  confirmed by a call. It shows only in `delta_utr_distal_share`.
- With three or more PACs in a last exon, losing a middle PAC raises the
  distal share, so the gene reads as `utr_lengthening`.
- A PAC just past a last exon that is not the gene's 3′-most is `intronic`.
- Genes annotated without transcript IDs have only their 3′-most exon as a
  terminal exon, so they cannot show `alternative_last_exon`.
- In an exploratory comparison, gained and lost calls are candidates, so
  only increases and decreases can make a pattern.

### Figures

`figures/` holds, for each comparison, three figures as PDF and PNG, and six
figures across comparisons. The report embeds the PNGs. The figures show the
calls in the statistics tables; they never make calls of their own. The
figures that set comparisons side by side describe each comparison's own
results, each treatment against its direct control, and add no test between
treatments.

- **`CONDITION_vs_CONTROL.volcano`**: every tested PAC's change in fitted PAU
  against its PAC-level p-value. Confirmed gains (`gained`,
  `increased_usage`) and losses (`lost`, `decreased_usage`) are colored,
  candidates are open circles, and PACs without a PAC call are grey (a grey
  density above 5,000 PACs). Up to 20 genes in each direction are labeled
  with their gene names, each at its most significant PAC. The y-axis uses the raw p-value because
  `pac_fdr` is missing outside screened genes. PACs without a p-value, such
  as unstable zero-boundary fits, are counted in the subtitle, and p-values
  of 0 are drawn as triangles at the top.
- **`CONDITION_vs_CONTROL.distal_usage`**: for each tested gene, the fitted
  usage of its distal PAC in the control against the treatment.
  - The distal PAC is the gene's most 3′ tested PAC in a last exon or
    downstream of the annotated end.
  - Points are colored by the gene's APA pattern, the first one when it has
    several. They are filled when the distal PAC has a confirmed call, and
    drawn as diamonds when only flagged PACs support the pattern.
  - `direction` in `CONDITION_vs_CONTROL.distal_usage.tsv.gz` is that PAC's
    own call: `distal_up`, `distal_down`, their `_candidate` forms, or
    `none`. A distal PAC can fall because usage moved to a tandem site or
    into an intron; the color tells which.
  - Up to 20 genes whose distal PAC rose and 20 whose distal PAC fell, those
    with the largest changes, are labeled with their gene names.
- **`CONDITION_vs_CONTROL.site_classes`**: confirmed PAC calls by where the
  PAC lies (terminal exon, other exon, intron, downstream of the gene), with
  losses to the left of zero and gains to the right. A shift to intronic
  polyadenylation shows as intron gains.
- **`event_counts`**: PAC events, candidates in lighter shades, and gene
  events from the `.genes` tables, for every comparison. A lower panel counts
  genes per APA pattern; a gene with two patterns counts in both. Patterns
  that only flagged PACs support are counted in a panel of their own.
- **`apa_pattern_grid`**: genes with an APA pattern in two or more
  comparisons, as rows, against the comparisons.
  - Each cell shows the gene's first pattern there. A + marks two or more,
    and a * a pattern that only flagged PACs support, drawn in its pattern's
    color. White is no pattern, and grey is a gene the comparison didn't
    test.
  - At most 50 genes are drawn: those shared by the most comparisons first,
    then by best gene FDR.
  - `apa_patterns_by_comparison.tsv.gz` lists every tested gene in that
    order, with its pattern in each comparison and a count of the comparisons
    where it has one. A cell is empty where the gene wasn't tested.
- **`concordance`** and **`concordance_matrix`**: how far comparisons agree.
  - Each panel of `concordance` takes two comparisons and plots every PAC
    tested in both: its change in PAU in one against the other. The panel
    gives Pearson r and the number of PACs.
  - PACs with a confirmed call in both comparisons, or in one, are colored.
  - Beyond 15 pairs (more than six comparisons), only pairs that share a
    control or a condition get a panel.
  - `concordance_matrix` shows r for every pair.
  - `concordance.tsv.gz` has one row per pair: its relation (`shared_control`,
    `chained`, or `unrelated`), the shared PACs, the PACs called in both, and
    `pearson_r`, which is empty with fewer than three shared PACs.
  - Read r with its relation in mind:
    - Two comparisons against one control share that control's estimate, so
      their changes correlate positively without any shared biology.
    - Chained comparisons, where one's treatment is the other's control,
      estimate the middle condition from the same samples, with opposite
      signs, so they correlate negatively. A rescue that reverses a treatment
      gives a strongly negative r partly for that reason.
- **`effect_vs_coverage`**: each PAC's change in fitted PAU against the reads
  at its gene in the less-covered group, one panel per comparison, with
  `min_gene_total` marked. Calls driven by low coverage would cluster on the
  left.
- **`pau_pca`**: samples on the first two principal components of observed
  PAU, colored by condition.
  - It uses the genes with at least `min_gene_total` reads in every sample,
    and those genes' PACs observed in every sample. There is no zero-filling
    and no pseudocount.
  - Each PAC is centered across samples. Each component's sign makes its
    largest loading positive; when loadings tie, as a two-PAC gene's do, the
    first PAC by ID decides.
  - Replicates should sit together. The report shows this figure beside the
    PAU sample correlation table, and `pau_pca.tsv` has the coordinates and
    the variance each component explains.

The PDFs use the standard Helvetica font, which is not embedded, and carry no
dates, so reruns reproduce them byte for byte. PNG rendering depends on the
fonts on the machine, so PNGs match between reruns on the same machine with
the same image or environment.

## How it works

1. **Validate.** The parameters, sample sheet, references, and alignment
   contigs are checked. The R packages for the statistics and the figures are
   loaded alongside, so a missing one stops the run within minutes.
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
     site. Two filters then reject reads that end inside a spliced transcript:
     - **Internal exon ends** (`internal_exon_end_filter`): reads that cross a
       splice junction, but whose short overhang into the next exon isn't
       aligned, end at the donor. A candidate where such a pile would peak,
       within one 25-nt bin of the annotated donor plus the kernel's peak
       offset, is rejected. Donors within two bins of any transcript's end
       are left alone. A real site whose reads pile up that close to an
       internal donor is rejected too.
     - **Readthrough** (`constitutive_readthrough_filter`): a candidate whose
       reads end in an exon block that reads splice onward from is rejected
       when every condition that supports it splices onward, with
       `constitutive_readthrough_min_junction_count` reads pooled over that
       condition's replicates. A supporting condition without them keeps the
       candidate. The cost: a real site used by a minority of transcripts is
       rejected too when its reads pile up inside such an exon, within about a
       read length of the donor.
6. **Annotate.** Genes, last exons, poly(A) signals, internal-priming flags,
   and known-PAC matches are added, and the atlas is frozen and checksummed.
7. **Count.** Each sample's read ends are assigned to the frozen atlas, giving
   raw counts, gene totals, and PAU.
8. **Test.** Each comparison family (a control and the treatments that name
   it) is fitted and tested, and each gene's APA pattern is classified.
   Motif-class preference and an exploratory k-mer enrichment follow.
9. **Figures and report.** ggplot2 draws the figures in `figures/` from the
   statistics tables, and everything is summarized in `report/index.html`.

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

Gzip files and PDFs carry no timestamps and every random step is seeded, so a
fresh rerun with the same inputs, parameters, and software reproduces every
published file byte for byte. The exception is Nextflow's own reports in
`pipeline_info/`, which record the run's times.

Keep the Nextflow work directory until the analysis is final: `-resume` reuses
finished steps from it.

- **Which steps rerun.** Every step from the alignment scan onward reads the
  resolved parameters, so changing any parameter reruns them all. That
  includes `outdir`, the `save_*` flags, and the Slurm and CPU settings.
  Alignment preparation and strandedness inference are reused unless their
  own inputs change. Set those parameters before the first run.
- **After updating PACusage.** Rebuild the image first under
  `-profile apptainer`.
  - **Python package:** Nextflow does not track most of it, so start a fresh
    run: use a new work directory, or leave out `-resume`.
  - **Annotation:** the one exception. Its step stages `annotation.py` and
    `reference.py`, so `-resume` reruns annotation and every step after it
    when they change.
  - **R scripts:** tracked. Changes to the statistics script rerun the
    statistics, and changes to the figure script redraw the figures.
  - **Version 0.2.0** changes the atlas (the `last_exon` column and refined
    terminal exons). The atlas checksum seeds the statistics, so a rerun
    moves stabilized genes' p-values and bootstrap intervals slightly.

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

The R tests need R 4.5 or newer, with DRIMSeq, stageR, limma, BiocParallel,
ggplot2, and yaml. Set `R_LIBS` to a library that has them:

```bash
Rscript tests/r/test_usage_model.R scripts/fit_usage_model.R
Rscript tests/r/test_usage_model_simulation.R scripts/fit_usage_model.R
Rscript tests/r/test_usage_figures.R scripts/plot_usage_figures.R scripts/fit_usage_model.R
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
| `scripts/plot_usage_figures.R` | The figures, drawn with ggplot2 from the statistics tables |
| `bin/pacusage` | The command Nextflow tasks run |
| `containers/`, `envs/` | The Apptainer recipe and the Conda environment |
| `tests/unit/`, `tests/integration/` | pytest tests |
| `tests/r/` | R tests for the statistics and the figures |
| `tests/pipeline/` | The end-to-end script and its result checks |
| `tests/fixtures/` | Synthetic references and alignments, and their builders |
| `docs/` | The design, the decisions log, and the methods document for peer review (`docs/methods.html`) |
