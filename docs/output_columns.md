# PACusage output columns

This guide explains every table that PACusage publishes: what one row is, what each column
means, its units, and its possible values. Every table is tab-separated text with one
header line, `.gz` files are gzip-compressed, and the `.parquet` files hold the same kind
of table in Parquet format, for programs. Coordinates are zero-based and half-open, as in
BED files, except `locus`, which is 1-based for genome browsers.

## Terms used throughout

- **PAC**: a polyadenylation site cluster, one place where transcripts end: a single
  cleavage and polyadenylation site, or several sites closer together than the library can
  separate. Its `pac_id`, such as `PACv1.GRCh38.chr1.+.1234567`, holds the atlas version,
  genome assembly, chromosome, strand, and representative coordinate.
- **PAC atlas**: the frozen list of PACs (`atlas/pacs.v1.*`), found once in the read ends of
  every sample pooled, before and without any treatment-control comparison. Every sample is
  then counted at these PACs.
- **Reads**: counts are of reads for single-end data and of read pairs (fragments) for
  paired-end data, after the read filters. Each gives one read end.
- **PAU** (PAC usage): a PAC's reads divided by all reads assigned to its gene's PACs in the
  same sample, from 0 to 1, so a gene's PAUs add up to 1. No pseudocounts are added: PAU is
  empty when the gene has no reads in that sample.
- **Fitted PAU and observed PAU**: observed PAU comes straight from one sample's counts.
  Fitted PAU is the statistical model's estimate for a group of samples, the mean of the
  model's fitted proportions over the group's samples. Changes in usage (`delta_pau`) and
  calls use fitted PAU.
- **Comparison**: one treatment condition tested against the control it names in the sample
  sheet, its direct control, named `CONDITION_vs_CONTROL` (for example
  `TreatmentA_vs_DMSO`). Treatments are never tested against each other.
- **Comparison family**: a control condition and every treatment that names it, fitted
  together in one model, so each gene's precision is estimated from all of their samples.
  Family tables are named after the control (`FAMILY`). A condition that is a treatment and
  also another condition's control belongs to two families.
- **Tested PAC**: a PAC that passed its comparison family's statistical filters, in a gene
  with at least two such PACs. The `.pacs`, `.calls`, and `fitted_pau` tables hold only
  tested PACs, and the `.genes` tables only their genes; `FAMILY.statistical_filtering`
  says why the other PACs were not tested.
- **Gene FDR and PAC FDR**: `gene_fdr` is the gene-level screen, the Benjamini-Hochberg
  adjusted p-value of the comparison's test of whether the gene's reads are split
  differently among its tested PACs. `pac_fdr` says which of the gene's PACs changed; it is
  stageR's stage-wise adjusted p-value, given only in genes that pass the screen. A
  confirmed call needs `gene_fdr` of at most gene_fdr (default 0.05) and `pac_fdr` of at
  most site_fdr (default 0.05); both limit false discoveries within one comparison.
- **Call and candidate call**: a PAC's call in a comparison is its `event_type`. The
  confirmed calls, `gained`, `lost`, `increased_usage`, and `decreased_usage`, pass both
  FDRs and the effect and detection rules. The candidate calls, `gained_candidate` and
  `lost_candidate`, meet the detection and effect rules of a gain or loss but miss
  significance, stability, coverage, confidence, or replication; they are not confirmed.
- **Detected** (in a group): at least event_min_supporting_samples (default 2) of the
  group's samples have a read at the PAC. Gained and lost calls use this.
- **Active PAC**: a PAC with at least active_pac_min_pau (default 0.10) of its gene's
  fitted usage in a group. `active_pacs_change` in the `.genes` tables compares their
  number between the groups.
- **Shift**: where a gene's usage moved, named once. PAU are shares of the gene, so usage
  gained at one PAC is lost at others, and one shift gives calls in both directions. The
  `shift_` columns of the `.genes` tables name the PAC the usage moved from and the PAC it
  moved to, and say whether the new PAC is more distal or more proximal.
- **Exact-boundary and proximal-tag libraries**: in an exact-boundary library, read ends
  mark the transcript's poly(A) junction, so a PAC is placed to the nucleotide
  (`coordinate_precision` `exact`). In a proximal-tag library, such as Plasmidsaurus 3'
  tag-seq, reads end some distance upstream of the site. Calibration against annotated
  transcript ends measures that distance, and each PAC is reported as an estimated
  position with a resolution region (`estimated`), not as a cleavage site.
- **Internal priming flag**: oligo-dT can prime on an A-rich stretch of the genome instead
  of on a poly(A) tail, which makes a false read end. A PAC is flagged
  (`internal_priming_flag`) when the genome just downstream of it is A-rich. Flagged PACs
  are low confidence: they are counted and tested, and can increase or decrease in usage,
  but are never called gained or lost.
- **Exploratory comparison**: a condition with fewer than min_replicates_per_condition
  (default 2) samples stops the run, unless insufficient_replicates_policy is `warn`
  (default `error`). The run then goes on, every comparison in a family that includes
  such a condition is marked `exploratory_insufficient_replicates`, and its gained and
  lost PACs are reported only as candidates.
- **Upstream and downstream**: always in transcript orientation, 5' to 3'. On the minus
  strand, upstream means higher coordinates.

## Which table answers which question

- Which PACs changed in a comparison? `statistics/CONDITION_vs_CONTROL.calls.tsv.gz`. Every
  tested PAC, called or not, is in `statistics/CONDITION_vs_CONTROL.pacs.tsv.gz`.
- Which genes changed, and how (the APA pattern)?
  `statistics/CONDITION_vs_CONTROL.genes.tsv.gz`.
- Where are the PACs, and what is known about each? `atlas/pacs.v1.metadata.tsv.gz`.
- Why did a candidate site not enter the atlas? `atlas/rejected_candidates.tsv.gz`.
- How many reads does each PAC have in each sample, and what is its usage there?
  `counts/pac_counts.tsv.gz` and `counts/observed_pau.tsv.gz`.
- Why was a PAC or gene not tested? `statistics/FAMILY.statistical_filtering.tsv.gz`.
- Did a treatment shift usage toward sites with a particular poly(A) signal?
  `motifs/CONDITION_vs_CONTROL.preference.tsv.gz` and
  `motifs/CONDITION_vs_CONTROL.preference_class.tsv.gz`.
- Do the comparisons agree with each other? `figures/concordance.tsv.gz` and
  `figures/apa_patterns_by_comparison.tsv.gz`.
- Did the library and the samples behave as expected? `qc/library_calibration.tsv`,
  `qc/SAMPLE.strandedness.tsv`, `qc/SAMPLE.fragment_filtering.tsv`, and
  `figures/pau_pca.tsv`.
- How exactly was the run set up? `manifest/normalized_samples.tsv`,
  `manifest/resolved_params.yaml`, and `manifest/software_versions.tsv`.

## Columns shared by many tables

### Identity columns

Every PAC-level table starts with these eight columns, except the per-sample count tables
(`counts/per_sample/`), and gene-level tables start with `gene_id` and `gene_name`, which
there name one gene. Where another table has one of these columns, such as `strand`, it
means the same, and the sections below do not repeat it.

| Column | Meaning |
|---|---|
| `pac_id` | The PAC's ID, `PACv1.ASSEMBLY.CHROM.STRAND.COORDINATE`: atlas version, genome assembly (the assembly parameter), chromosome, strand, and representative coordinate. |
| `gene_id` | The gene the PAC is assigned to, from the annotation; several IDs separated by commas when it is assigned to more than one gene; empty for an intergenic PAC. |
| `gene_name` | The gene's name from the annotation (`gene_name`, a GFF3 gene's `Name`, or `gene`), or its `gene_id` when it has none; comma-separated in the same order as `gene_id`. |
| `chrom` | Chromosome (contig), named as in the FASTA and annotation. |
| `start` | Start of the PAC's interval, zero-based: the coordinate itself for an exact-boundary PAC, the start of the resolution region for a proximal-tag PAC. |
| `end` | End of the interval, exclusive: the coordinate plus 1 for an exact-boundary PAC, the end of the resolution region for a proximal-tag PAC. The same interval is in `atlas/pacs.v1.bed.gz`. |
| `strand` | The transcript strand, `+` or `-`. |
| `locus` | The same interval, 1-based and inclusive, as `chrom:start-end` (for example `chr1:1234568-1234568`), ready to paste into IGV or the UCSC browser. |

## Tables

Parameter names such as min_gene_total refer to `nextflow_schema.json`, with their defaults
in parentheses; a run's own values are in `manifest/resolved_params.yaml`. Yes-or-no columns
hold `True` or `False`, except those that the statistics compute in `statistics/`, which
hold `TRUE` or `FALSE`; atlas columns copied into those tables keep `True` or `False`. An
empty cell, or `nan` in a few tables, means the value does not apply or could not be
computed.

*`atlas/`: the PAC atlas, and the candidate sites it left out.*

### `atlas/pacs.v1.metadata.tsv.gz`

One row per PAC in the frozen atlas: its position, gene, discovery evidence, and sequence
annotation. The PACs of every other table come from here.

| Column | Meaning |
|---|---|
| `coordinate` | The PAC's representative position, zero-based, counting the boundaries between bases; `pac_id` ends with it. For an exact-boundary PAC it is the read-end position with the strongest support; for a proximal-tag PAC, the strongest peak of the calibrated read-end signal, an estimate rather than a cleavage site. |
| `coordinate_precision` | `exact` for a nucleotide-resolution transcript end (exact-boundary runs), `estimated` for a position estimated from proximal tags. |
| `resolution_nt` | Width in nt of the PAC's interval (`start` to `end`): `1` for exact-boundary PACs; for proximal-tag PACs, the minimum resolvable separation from calibration, the smallest distance at which two PACs can be told apart. |
| `merged_candidate_coordinates` | Positions merged into the PAC during discovery, comma-separated and zero-based: the read-end positions within pac_cluster_radius (default 12 nt) clustered into an exact-boundary PAC, or the signal peaks too close to separate that were merged into a proximal-tag PAC. |
| `gene_region` | Where the PAC lies, on its strand: `last_exon` (a terminal exon: a transcript's final exon, unless it overlaps an internal exon of another transcript of the gene or is a single-exon model apart from the gene's spliced transcripts), `internal_exon` (any other exon), `intron` (in the gene but in no exon), `downstream_of_gene` (past the 3' end of the nearest gene, within max_downstream_distance, default 5000 nt), or `intergenic`. The first that applies wins, and only the genes where it applies are assigned. |
| `last_exon_locus` | For a `last_exon` PAC, its last exon as a 1-based locus (`chrom:start-end`), with the gene's overlapping terminal exons merged into one; for a `downstream_of_gene` PAC, the gene's 3'-most last exon. Empty for other PACs and for PACs assigned to several genes. |
| `ambiguous_gene_assignment` | `True` when the PAC is assigned to two or more genes, which `gene_id` lists; such PACs are counted but never tested. |
| `proximal_distal_rank` | The PAC's rank among all of its gene's atlas PACs in transcript orientation: `1` is the most upstream (proximal), and the largest number the most downstream (distal). Empty for intergenic PACs and PACs assigned to several genes. |
| `proximal_distal_label` | `proximal` for rank 1, `distal` for the last rank, and `middle` for the others; a gene's only PAC is `proximal`. Empty where the rank is empty. |
| `total_count` | Reads at the PAC during discovery, pooled over all samples. The final per-sample counts are in `counts/`. |
| `supporting_samples` | Samples with at least pac_min_sample_count (default 2) reads at the PAC, counted within its best-supported condition. |
| `total_supporting_samples` | Samples with at least pac_min_sample_count reads at the PAC, over all conditions. |
| `best_supporting_condition` | The condition whose replicates support the PAC best: one that meets pac_min_supporting_samples (default 2 samples) comes first, then the one with the most supporting samples (or the largest share of its samples, when that parameter is a fraction). |
| `candidate_status` | `primary`: passed discovery, with at least pac_min_total_count reads (default 10) and replicate support within one condition. `known_rescue_only`: an exact-boundary PAC near a known PAC that had too few reads for `primary` but at least known_pac_rescue_total (default 5), with the same replicate support. |
| `confidence` | `low` when the PAC is assigned to several genes or flagged for internal priming; otherwise `high` when it matches a known PAC and has at least 2 supporting samples, and `moderate` for the rest, including known_rescue_only PACs. Without a known_pacs file, no PAC is `high`. |
| `known_pac` | `True` when a PAC in the known_pacs BED file lies on the same strand within known_pac_match_radius (default 12 nt) of the coordinate; always `False` without that file. |
| `known_pac_coordinates` | Those known PACs' positions, comma-separated and zero-based (the 3' edge of each BED record); empty when none. |
| `known_rescue_only` | `True` for a PAC kept only by the relaxed known-PAC rule (candidate_status `known_rescue_only`); the primary motif and k-mer analyses leave these PACs out. |
| `internal_priming_flag` | `True` when the genome just downstream could have primed oligo-dT: in the internal_priming_window (default 20 nt), a run of A longer than internal_priming_max_a_run (default 6) or an A share above internal_priming_max_a_fraction (default 0.6). Flagged PACs are low confidence and are never called gained or lost. |
| `downstream_a_fraction` | Share of A bases in that downstream window (`downstream_sequence`), from 0 to 1. |
| `downstream_longest_a_run` | Longest run of consecutive A bases in that window, in nt. |
| `primary_pas_motif` | The chosen poly(A) signal (PAS) hexamer upstream of the PAC, in DNA letters (for example `AATAAA`); empty when none was found. It is chosen by catalog priority, then by a position in the core window, then by closeness to 20 nt upstream. |
| `primary_pas_motif_rna` | The same hexamer in RNA letters (for example `AAUAAA`); `none` when no motif was found. |
| `primary_motif_class` | The hexamer's class: `canonical` (AATAAA), `common_variant` (ATTAAA), `other_variant` (TATAAA, AGTAAA, AAGAAA, AATACA, CATAAA, GATAAA, AATATA, or AATAGA), or `no_recognized_motif`. A pas_motif_catalog file replaces these hexamers and classes. |
| `primary_motif_position` | Distance in nt upstream from the PAC coordinate to the hexamer's first base; empty when none was found. |
| `primary_motif_in_core` | `True` when that distance lies in the core window where a PAS is expected, pas_core_upstream_near to pas_core_upstream_far (default 10 to 35 nt). |
| `all_pas_motifs` | Every catalog hexamer found in the scan window, separated by `;`, each as `motif:position:class:in_core` (in_core `1` or `0`), the primary one first. |
| `upstream_sequence` | The PAS scan window: genomic sequence from pas_scan_upstream_far to pas_scan_upstream_near (default 50 to 5 nt) upstream of the PAC, read 5' to 3' on the PAC's strand. |
| `downstream_sequence` | The internal_priming_window (default 20 nt) of genomic sequence just downstream of the PAC, read 5' to 3' on the PAC's strand. |
| `fraction_within_2nt` | Share of the PAC's discovery reads that end within 2 nt of its coordinate, near 1 for a sharp site. Exact-boundary PACs only; `0` for proximal-tag PACs. |
| `width_90` | Width in nt of the central 90% of the PAC's read ends (5th to 95th percentile position), `0` when all end at one position. Exact-boundary PACs only; `0` for proximal-tag PACs. |
| `local_strand_enrichment` | The PAC's reads divided by the same-strand read ends in the flanking windows from pac_cluster_radius to twice that on each side (at least 1 in the denominator), so higher means a peak that stands out. Exact-boundary PACs only; `0` for proximal-tag PACs. |
| `poly_a_clip_fraction` | Share of the PAC's discovery reads whose 3' end carries a soft-clipped poly(A)-like tail (at least 6 nt, at least 80% A), direct evidence of the poly(A) junction. Exact-boundary PACs only; `0` for proximal-tag PACs. |
| `endpoint_model` | How the run interpreted read ends: `exact_boundary` or `proximal_tag`, chosen by calibration or set with endpoint_model. |
| `resolved_evidence_source` | The read end the run used: `read_3p` (the transcript-oriented 3' edge of the read's alignment; of read 1, for pairs), `read_5p` (the transcript-oriented 5' edge), `fragment_3p` (the outer 3' edge of a read pair), or `polyA_junction` (the 3' edge of reads that end in a poly(A)-like soft clip). |
| `primary_evidence_type` | `exact_boundary` for PACs placed at observed transcript ends, `endpoint_peak` for PACs placed at a peak of the calibrated proximal-tag signal. |
| `calibration_profile` | Version label of the calibration method, `PACusage-0.1`. |

### `atlas/rejected_candidates.tsv.gz`

One row per candidate site that discovery considered and rejected, with the reason. Its
`chrom` and `strand` columns are as in the identity columns.

| Column | Meaning |
|---|---|
| `coordinate` | The candidate's representative position, zero-based, as for atlas PACs. |
| `total_count` | Reads at the candidate, pooled over all samples. |
| `supporting_samples` | Samples with at least pac_min_sample_count (default 2) reads at the candidate, within its best-supported condition. |
| `capped_support` | Reads at the candidate summed over samples, counting at most 3 per sample, so no single sample dominates. |
| `member_coordinates` | Positions merged into the candidate, comma-separated and zero-based: its read-end positions (exact-boundary) or its merged signal peaks (proximal-tag). |
| `status` | Always `rejected`. |
| `rejection_reason` | Why: `total_count<N` (fewer than pac_min_total_count reads, default 10); `supporting_samples<N` or `supporting_sample_fraction<F` (no condition met pac_min_supporting_samples, default 2); or `internal_exon_end` or `constitutive_readthrough` (the proximal-tag filters, below). Several reasons are joined by `;`. |
| `fraction_within_2nt` | Share of its reads that end within 2 nt of its coordinate (exact-boundary; `0` for proximal-tag). |
| `width_90` | Width in nt of the central 90% of its read ends (exact-boundary; `0` for proximal-tag). |
| `poly_a_clip_fraction` | Share of its reads with a poly(A)-like soft clip at the 3' end (exact-boundary; `0` for proximal-tag). |
| `local_enrichment` | Its reads divided by the read ends in the flanking windows just outside its cluster, as `local_strand_enrichment` in the atlas (exact-boundary; `0` for proximal-tag). |
| `region_start` | Start of the candidate's region, zero-based: its first read-end position (exact-boundary), or the start of its resolution region (proximal-tag). |
| `region_end` | End of the region, exclusive: one past its last read-end position, or the end of its resolution region. |
| `resolution_nt` | `1` for exact-boundary candidates; for proximal-tag candidates, the calibrated minimum resolvable separation in nt. |
| `total_supporting_samples` | Samples with at least pac_min_sample_count reads at the candidate, over all conditions. |
| `best_supporting_condition` | The condition whose replicates support it best, as in the atlas; empty when no sample does. |
| `supporting_conditions` | Every condition whose replicates meet pac_min_supporting_samples, separated by `;`; empty when none. The readthrough filter consults only these. |

The two filters apply only to proximal-tag runs, and reject candidates made by reads that
end inside a spliced transcript:
- `internal_exon_end` (internal_exon_end_filter, default true): reads that cross a splice
  junction, but whose short overhang into the next exon is not aligned, end at the exon's
  3' end. A candidate is rejected when it lies within one discovery bin (at most
  proximal_bin_size, default 25 nt) of where such reads, at an annotated internal exon's 3'
  end, would place a peak. Exon ends within two bins of a terminal exon's 3' end are
  exempt, and a real site that close to an internal exon end is rejected too.
- `constitutive_readthrough` (constitutive_readthrough_filter, default true): a candidate
  is rejected when its reads end in an aligned block from which reads splice directly
  onward, and every condition that supports it has at least
  constitutive_readthrough_min_junction_count (default 2) such spliced reads, pooled over
  the condition's samples. A supporting condition without them keeps the candidate.

*`counts/`: raw counts at the frozen atlas, and observed usage.*

### `counts/pac_counts.tsv.gz`

One row per atlas PAC, in genomic order, with one column of raw counts per sample. Every
atlas PAC is included, also intergenic ones and those assigned to several genes.

| Column | Meaning |
|---|---|
| `SAMPLE` | One column per sample, named by its sample_id, in sorted order: the raw number of reads assigned to the PAC in that sample (an integer, not normalized). |

### `counts/pac_counts.long.parquet`

The same counts as `counts/pac_counts.tsv.gz`, in long form: one row per PAC and sample.

| Column | Meaning |
|---|---|
| `sample_id` | The sample. |
| `count` | Raw reads assigned to the PAC in this sample. |

### `counts/gene_totals.tsv.gz`

One row per gene and sample: the reads assigned to the gene's PACs, the denominator of PAU.
PACs assigned to several genes, and intergenic PACs, are left out.

| Column | Meaning |
|---|---|
| `sample_id` | The sample. |
| `gene_total` | Reads assigned to the gene's PACs in this sample, the sum of their counts. |

### `counts/observed_pau.tsv.gz`

One row per PAC and sample, for PACs assigned to exactly one gene: the raw count and the
observed PAU.

| Column | Meaning |
|---|---|
| `sample_id` | The sample. |
| `count` | Raw reads assigned to the PAC in this sample. |
| `gene_total` | Reads assigned to all of the gene's PACs in this sample. |
| `pau` | Observed PAU: count divided by gene_total, from 0 to 1; empty when gene_total is 0. No pseudocount is added. |

### `counts/per_sample/SAMPLE.pac_counts.tsv.gz`

Published only with save_intermediates (default false): one sample's raw count at every
atlas PAC, one row per PAC, before the samples are merged. Unlike the other PAC tables it
has just three columns, `sample_id`, `pac_id`, and `count`, in that order.

| Column | Meaning |
|---|---|
| `sample_id` | The sample. |
| `count` | Raw reads assigned to the PAC in this sample. |

*`evidence/`: each sample's read ends, before any PAC is called.*

### `evidence/SAMPLE.3prime_evidence.tsv.gz`, `evidence/SAMPLE.3prime_evidence.parquet`

One row per position where the sample's reads end, with the number of reads: the evidence
that discovery and counting use. The two files hold the same rows, and `strand` is the
transcript strand, set from the library's strandedness.

| Column | Meaning |
|---|---|
| `sample_id` | The sample. |
| `contig` | Chromosome, named as in the FASTA and annotation (after chromosome_aliases). |
| `coordinate` | Where the reads end, zero-based: the transcript-oriented edge of the alignment that evidence_source names. |
| `count` | Reads that passed the read filters and end here. |
| `poly_a_clip_count` | Of those, reads whose 3' end carries a soft-clipped poly(A)-like tail (at least 6 nt, at least 80% A; read 1's, for pairs). |
| `evidence_source` | Which read end this is: `read_3p`, `read_5p`, `fragment_3p`, or `polyA_junction`, as `resolved_evidence_source` in the atlas. |

### `evidence/SAMPLE.splice_continuations.tsv.gz`

One row per splice in the sample's reads: an aligned block joined directly to the next by a
gap in the alignment (a CIGAR N), which the readthrough filter uses. Rows are written only
in proximal-tag runs with constitutive_readthrough_filter on (default true), and otherwise
the file has only its header; `strand` is the transcript strand, so on the minus strand the
upstream block has the higher coordinates.

| Column | Meaning |
|---|---|
| `sample_id` | The sample. |
| `contig` | Chromosome, named as in the FASTA and annotation. |
| `upstream_start` | Start, zero-based, of the aligned block before the splice, in transcript orientation. |
| `upstream_end` | End, exclusive, of that block. |
| `downstream_start` | Start, zero-based, of the aligned block after the splice. |
| `downstream_end` | End, exclusive, of the block after the splice. |
| `count` | Reads with this splice. |

*`figures/`: the tables behind the figures, which are listed under Other files.*

### `figures/CONDITION_vs_CONTROL.distal_usage.tsv.gz`

One row per tested gene that has a tested PAC in a last exon or downstream of the gene: the
usage of its distal PAC, as drawn in the distal-usage figure. Rows are sorted by gene FDR;
this distal PAC can differ from the one behind `delta_utr_distal_share` in the `.genes`
table, which looks only at the main last exon.

| Column | Meaning |
|---|---|
| `condition` | The treatment condition. |
| `control_condition` | Its direct control. |
| `direction` | The distal PAC's own call: `distal_up` (gained or increased usage), `distal_down` (lost or decreased usage), `distal_up_candidate`, `distal_down_candidate`, or `none`. |
| `apa_pattern` | The gene's APA pattern from the `.genes` table; `none` when it has none. |
| `distal_pac_id` | The distal PAC's pac_id: the gene's most 3' tested PAC in any last exon or downstream of the gene. |
| `distal_locus` | That PAC's locus, 1-based. |
| `distal_gene_region` | That PAC's gene_region: `last_exon` or `downstream_of_gene`. |
| `fitted_control_distal_pau` | The distal PAC's fitted PAU in the control, as a share of the whole gene. |
| `fitted_treatment_distal_pau` | Its fitted PAU in the treatment. |
| `delta_distal_pau` | Treatment minus control. |
| `distal_event_type` | The distal PAC's event_type from the `.pacs` table. |
| `gene_fdr` | The gene's FDR in this comparison. |
| `distal_pac_fdr` | The distal PAC's pac_fdr; empty outside genes that pass the screen. |
| `tested_pacs` | Number of the gene's tested PACs, in all regions. |

### `figures/apa_patterns_by_comparison.tsv.gz`

One row per gene tested in at least one comparison, with its APA pattern in each: genes with
patterns in the most comparisons come first, then by their best gene FDR. The
`apa_pattern_grid` figure draws the first 50 genes with a pattern in two or more
comparisons.

| Column | Meaning |
|---|---|
| `patterned_comparisons` | Number of comparisons in which the gene's apa_pattern is not `none` (`unclassified_change` counts). |
| `CONDITION_vs_CONTROL` | One column per comparison: the gene's apa_pattern there; empty when that comparison did not test the gene. |

### `figures/concordance.tsv.gz`

One row per pair of comparisons: how well their changes in PAU agree at the PACs both
tested. It describes the comparisons side by side and tests nothing between treatments.

| Column | Meaning |
|---|---|
| `comparison_a` | The first comparison of the pair, CONDITION_vs_CONTROL. |
| `comparison_b` | The second comparison. |
| `relation` | `shared_control` (both use the same control), `chained` (one's treatment is the other's control), or `unrelated`. A shared control makes the changes correlate positively, and a chain negatively, even without shared biology. |
| `shared_pacs` | PACs tested in both comparisons, with a fitted change in both. |
| `called_in_both` | Of those, PACs with a confirmed call (gained, lost, increased or decreased usage) in both. |
| `pearson_r` | Pearson correlation of the two comparisons' delta_pau over the shared PACs, from -1 to 1; empty with fewer than 3 shared PACs or no spread. |

### `figures/pau_pca.tsv`

One row per sample: its place on the first two principal components of observed PAU, as
drawn in the `pau_pca` figure. It uses the genes with at least min_gene_total (default 20)
reads in every sample, each PAC centered across samples, with no zero-filling or
pseudocount.

| Column | Meaning |
|---|---|
| `sample_id` | The sample. |
| `condition` | Its condition. |
| `PC1` | The sample's score on the first principal component, in PAU units; the sign makes the largest loading positive. |
| `PC2` | Its score on the second principal component. |
| `pc1_variance_fraction` | Share of the total variance that PC1 explains, from 0 to 1; the same on every row. |
| `pc2_variance_fraction` | Share of the total variance that PC2 explains. |

*`manifest/`: what the run used, for reproducing it.*

### `manifest/normalized_samples.tsv`

One row per sample: the sample sheet as validated, with defaults filled in. Any other
sample-sheet columns, such as model covariates, follow in alphabetical order.

| Column | Meaning |
|---|---|
| `sample_id` | The sample's ID. |
| `alignment` | Absolute path of its BAM or CRAM file. |
| `condition` | Its condition. |
| `control` | The control condition named in the sample sheet; empty for a control condition. |
| `control_condition` | The condition's direct control; a control condition lists itself. Conditions that share one form a comparison family. |
| `replicate` | Replicate label from the sample sheet, if given; descriptive only. |
| `batch` | Batch label, if given; it enters the model only when named in model_covariates. |
| `donor` | Donor label, if given; it too enters the model only through model_covariates. |
| `layout` | Requested layout, `SE`, `PE`, or `auto`: the sample sheet's value, else the layout parameter (default `auto`). The layout used is in `qc/SAMPLE.strandedness.tsv`. |
| `strandedness` | Requested strandedness, `forward`, `reverse`, or `auto` (default `auto`). The value used is in `qc/SAMPLE.strandedness.tsv`. |
| `library_profile` | Protocol profile: `generic_3prime` (default), `plasmidsaurus_3prime`, or `exact_boundary`. |
| `evidence_source` | Requested read end: `auto` (default), `read_3p`, `read_5p`, `fragment_3p`, or `polyA_junction`. The one the run used is the atlas's `resolved_evidence_source`. |

### `manifest/input_checksums.tsv`

One row per input file, with its SHA-256 checksum, to confirm later that the same inputs
were used.

| Column | Meaning |
|---|---|
| `role` | What the file is: `sample_sheet`, `fasta`, `annotation`, or `alignment` (one row per sample, in sample-sheet order). |
| `path` | The file's absolute path. |
| `sha256` | SHA-256 checksum of the file's contents. |

### `manifest/software_versions.tsv`

One row per program or package the run used.

| Column | Meaning |
|---|---|
| `software` | The program: `pacusage`, `python`, `pysam`, `samtools (in pysam)`, `htslib (in pysam)`, `R`, `DRIMSeq`, `stageR`, `limma`, `BiocParallel`, `ggplot2`, or `ggrepel`. |
| `version` | Its version. |

### `manifest/calibration_kernel.tsv`

The pooled calibration kernel: how far read ends fall from annotated transcript ends,
smoothed and pooled with equal weight per sample. Proximal-tag discovery and counting use it
to place reads at PACs; one row per offset with a nonzero weight.

| Column | Meaning |
|---|---|
| `offset` | Distance in nt from the annotated transcript end to the read end, in transcript orientation; positive means the read ends upstream of the transcript end. It runs from minus to plus calibration_max_distance (default 1000). |
| `weight` | Share of the pooled, smoothed read ends at this offset; the weights sum to 1. |

*`motifs/`: poly(A) signals at each PAC, and whether treatments shift usage among them.*

### `motifs/pac_motifs.tsv.gz`

One row per atlas PAC: its poly(A) signal annotation, copied from the atlas.

| Column | Meaning |
|---|---|
| `primary_pas_motif` | As in the atlas: the chosen PAS hexamer in DNA letters; empty when none. |
| `primary_pas_motif_rna` | As in the atlas: the hexamer in RNA letters; `none` when none. |
| `primary_motif_class` | As in the atlas: `canonical`, `common_variant`, `other_variant`, or `no_recognized_motif`. |
| `primary_motif_position` | As in the atlas: distance in nt upstream from the PAC to the hexamer's first base. |
| `primary_motif_in_core` | As in the atlas: `True` when the hexamer lies in the core window (default 10 to 35 nt upstream). |
| `all_pas_motifs` | As in the atlas: every hexamer found, as `motif:position:class:in_core`, separated by `;`. |
| `upstream_sequence` | As in the atlas: the PAS scan window's sequence, 5' to 3' on the PAC's strand. |
| `known_rescue_only` | As in the atlas: `True` for PACs kept only by the relaxed known-PAC rule, which the primary motif analyses leave out. |

### `motifs/motif_scores.tsv`, `motifs/motif_scores_known_rescue_sensitivity.tsv`

One row per sample and primary PAS hexamer: how much of a typical gene's usage goes to PACs
with that hexamer, the input to the `.preference` tests. It uses genes with at least
min_gene_total (default 20) reads in every sample and at least two PACs; known-rescue-only
PACs are removed and each gene's usage renormalized over the rest, except in the
`known_rescue_sensitivity` file, which keeps them.

| Column | Meaning |
|---|---|
| `sample_id` | The sample. |
| `primary_pas_motif_rna` | The primary PAS hexamer in RNA letters; `none` for PACs without a recognized motif. |
| `primary_motif_class` | The hexamer's class. |
| `motif_usage` | Mean, over the informative genes, of the share of each gene's reads at its PACs with this hexamer, each gene counting once; from 0 to 1. Genes without such a PAC are not averaged in. |
| `transformed_motif_usage` | arcsin(sqrt(motif_usage)), in radians: the scale the preference test uses, which handles 0 and 1 without pseudocounts. |
| `informative_genes` | Genes averaged: covered genes with at least one PAC whose primary PAS is this hexamer. |

### `motifs/motif_class_scores.tsv`, `motifs/motif_class_scores_known_rescue_sensitivity.tsv`

One row per sample and motif class: the same scores as `motifs/motif_scores.tsv`, with the
hexamers of a class pooled within each gene; the input to the `.preference_class` tests.

| Column | Meaning |
|---|---|
| `sample_id` | The sample. |
| `primary_motif_class` | The class: `canonical`, `common_variant`, `other_variant`, or `no_recognized_motif`. |
| `motif_usage` | Mean, over the informative genes, of the share of each gene's reads at its PACs of this class, each gene counting once; from 0 to 1. |
| `transformed_motif_usage` | arcsin(sqrt(motif_usage)), in radians, the scale the test uses. |
| `informative_genes` | Genes averaged: covered genes with at least one PAC of this class. |

### `motifs/CONDITION_vs_CONTROL.preference.tsv.gz`, `motifs/CONDITION_vs_CONTROL.preference_known_rescue_sensitivity.tsv.gz`

One row per primary PAS hexamer: whether the treatment shifted usage toward PACs with that
signal, tested with limma on the motif scores of the comparison family's samples, using the
main model's design. A row is tested only when every sample of the comparison has at least
motif_preference_min_genes (default 50) informative genes, and a comparison with no such row
gets no file; the `known_rescue_sensitivity` file is a check with known-rescue-only PACs
kept, not a primary result.

| Column | Meaning |
|---|---|
| `primary_pas_motif_rna` | The primary PAS hexamer in RNA letters; `none` for PACs without a recognized motif. |
| `primary_motif_class` | The hexamer's class. |
| `condition` | The treatment condition. |
| `control_condition` | Its direct control. |
| `delta_motif_usage` | treatment_mean minus control_mean: how much of a typical gene's usage moved toward (positive) or away from (negative) PACs with this hexamer. |
| `fdr` | Benjamini-Hochberg adjusted p-value across the rows of this table. |
| `pvalue` | Raw p-value of limma's moderated t-test of the treatment effect, on the arcsine-square-root scale. |
| `control_mean` | Mean motif_usage over the control samples. |
| `treatment_mean` | Mean motif_usage over the treatment samples. |
| `transformed_coefficient` | The treatment effect on the arcsine-square-root scale that the test uses, adjusted for model_covariates. |
| `informative_genes` | The fewest informative genes in any sample of the comparison. |

### `motifs/CONDITION_vs_CONTROL.preference_class.tsv.gz`, `motifs/CONDITION_vs_CONTROL.preference_class_known_rescue_sensitivity.tsv.gz`

One row per motif class: the same test as the `.preference` tables, on the class scores
(canonical, common variant, other variants, and no recognized motif).

| Column | Meaning |
|---|---|
| `primary_motif_class` | The class: `canonical`, `common_variant`, `other_variant`, or `no_recognized_motif`. |
| `condition` | The treatment condition. |
| `control_condition` | Its direct control. |
| `delta_motif_usage` | treatment_mean minus control_mean: how much of a typical gene's usage moved toward (positive) or away from (negative) PACs of this class. |
| `fdr` | Benjamini-Hochberg adjusted p-value across the classes of this table. |
| `pvalue` | Raw p-value of limma's moderated t-test of the treatment effect, on the arcsine-square-root scale. |
| `control_mean` | Mean motif_usage over the control samples. |
| `treatment_mean` | Mean motif_usage over the treatment samples. |
| `transformed_coefficient` | The treatment effect on the arcsine-square-root scale, adjusted for model_covariates. |
| `informative_genes` | The fewest informative genes in any sample of the comparison. |

### `motifs/CONDITION_vs_CONTROL.kmer_enrichment.tsv.gz`

Exploratory: one row per upstream k-mer, testing whether a comparison's gained and
increased-usage PACs carry it more often than the other tested PACs of the same genes (a
Cochran-Mantel-Haenszel test stratified by gene), with known-rescue-only PACs and PACs
assigned to several genes left out. Sequence, PAC position, and the choice of calls can be
confounded, so treat it as a source of hypotheses.

| Column | Meaning |
|---|---|
| `kmer` | The k-mer, motif_kmer_length (default 6) DNA letters, found in the PACs' upstream_sequence (the PAS scan window). |
| `common_odds_ratio` | Mantel-Haenszel odds ratio, pooled over genes, for a gained or increased PAC to carry the k-mer compared with the gene's other tested PACs; above 1 means enriched at gained or increased PACs. It can be `0` or `inf` at the extremes. |
| `ci_low` | Lower end of the odds ratio's 95% confidence interval; `nan` when it cannot be computed. |
| `ci_high` | Upper end of that interval. |
| `fdr` | Benjamini-Hochberg adjusted p-value across this comparison's k-mers. |
| `pvalue` | Raw p-value of the Cochran-Mantel-Haenszel test (chi-square, 1 degree of freedom). |
| `event_pacs` | Gained or increased-usage PACs in the genes where some tested PAC carries the k-mer. |
| `background_pacs` | The other tested PACs of those genes. |
| `informative_genes` | Genes that inform the test: they have both kinds of PAC, and the k-mer is in some but not all of their PACs. |
| `analysis_label` | Always `exploratory`. |

### `motifs/kmer_enrichment_status.tsv`

One row per comparison: how many k-mers the exploratory enrichment tested. With
run_kmer_enrichment set to false (default true), the file instead has one column, `status`,
with the value `not_requested`, and no k-mer tables are written.

| Column | Meaning |
|---|---|
| `comparison` | The comparison, CONDITION_vs_CONTROL. |
| `tested_kmers` | Number of k-mers with at least one informative gene: the rows of its k-mer table. |

Columns that only some runs have:

| Column | Meaning |
|---|---|
| `status` | Only with run_kmer_enrichment false, in place of the other columns: `not_requested`. |

*`qc/`: checks and summaries from every step, from input validation to counting.*

### `qc/input_validation.tsv`

One row per check that passed: each sample's row of the sample sheet, and each alignment. A
failed check stops the run with a message instead of a row.

| Column | Meaning |
|---|---|
| `sample_id` | The sample checked. |
| `check` | `sample_sheet` (the sample's row is valid) or `alignment_open` (its alignment opens, and its contigs match the FASTA and annotation). |
| `status` | `PASS`. |
| `detail` | `validated`, or the alignment's layout and sort order, for example `SE; sort=coordinate`. |
| `exploratory_insufficient_replicates` | On sample_sheet rows, `true` when the sample's condition has fewer than min_replicates_per_condition (default 2) samples, which only insufficient_replicates_policy `warn` allows, and `false` otherwise; empty on alignment_open rows. |

### `qc/control_mapping.tsv`

One row per condition: its direct control, and its role in the comparisons.

| Column | Meaning |
|---|---|
| `condition` | The condition. |
| `control_condition` | Its direct control; a control condition lists itself. |
| `role` | `control` (a control that names no control of its own), `treatment` (compared with its control), or `treatment_and_control` (compared with its control, and also the control of another condition). |

### `qc/reference_preparation.tsv`

One row: how the genome FASTA was prepared. It is linked, never copied or changed, and
indexed in the work directory.

| Column | Meaning |
|---|---|
| `source_fasta` | Absolute path of the FASTA. |
| `prepared_fasta` | File name of the link in the work directory, `genome.fa`. |
| `prepared_fai` | File name of the index built for it, `genome.fa.fai`. |
| `action` | `generated_index`: the index is always built afresh. |
| `fasta_sha256` | SHA-256 checksum of the FASTA. |
| `fai_sha256` | SHA-256 checksum of the index. |

### `qc/SAMPLE.alignment_preparation.tsv`

One row: how the sample's alignment was prepared. Sources are never changed: a sorted
alignment is linked, and only an unsorted one is copied, sorted, into the work directory.

| Column | Meaning |
|---|---|
| `source_alignment` | Absolute path of the BAM or CRAM from the sample sheet. |
| `prepared_alignment` | File name of the prepared alignment in the work directory. |
| `prepared_index` | File name of its index (`.bai` or `.crai`). |
| `format` | `BAM` or `CRAM`. |
| `action` | `reused_alignment_and_index` (the sorted source and its index linked), `indexed` (the sorted source linked and a new index built), or `sorted_and_indexed` (the unsorted source sorted into a new file and indexed). |
| `layout_detected` | Layout seen in the prepared file's first 10,000 records: `SE` or `PE` (at least 95% unpaired or paired), `mixed`, or `unknown` (no records). |
| `source_sha256` | SHA-256 checksum of the source file. |
| `prepared_sha256` | SHA-256 checksum of the prepared file, the same as the source's when it is a link. |
| `source_size` | Size of the source file in bytes; later steps check that it has not changed. |
| `source_mtime_ns` | The source file's modification time, in nanoseconds since 1970; later steps check it too. |
| `sample_id` | The sample. |

### `qc/SAMPLE.strandedness.tsv`

One row: the read orientation inferred for the sample from reads over annotated exons, and
the settings resolved for it. An inferred strandedness that conflicts with a requested one
stops the run.

| Column | Meaning |
|---|---|
| `sample_id` | The sample. |
| `informative_fragments` | Reads sampled for the test: primary alignments (read 1, for pairs) over exons of exactly one gene, drawn at random from the whole file, at most strand_max_sampled_fragments (default 200,000). |
| `forward_count` | Sampled reads aligned to the same strand as their gene. |
| `reverse_count` | Sampled reads aligned to the opposite strand. |
| `forward_fraction` | forward_count divided by informative_fragments. |
| `reverse_fraction` | reverse_count divided by informative_fragments. |
| `inferred_strandedness` | `forward` or `reverse` when that share reaches strand_decision_fraction (default 0.8), from at least strand_min_informative_fragments (default 10,000) reads; otherwise `ambiguous`. |
| `reason` | The inference's reason, in words. |
| `library_profile` | The sample's protocol profile. |
| `layout` | Layout used, `SE` or `PE`. |
| `strandedness` | Strandedness used: the requested value, or the inferred one when `auto` was requested. |
| `evidence_source` | Requested evidence source, after the profile's default: `auto` or a named source. The one the run used is in the atlas. |
| `endpoint_model` | Requested endpoint model, after the profile's default: `auto`, `exact_boundary`, or `proximal_tag`. Calibration makes the final choice. |
| `alignment_basename` | File name of the prepared alignment. |

### `qc/library_calibration.tsv`

One row per sample: where its read ends fall relative to annotated transcript ends, for the
run's evidence source. It decides the endpoint model, and every sample must classify alike.

| Column | Meaning |
|---|---|
| `sample_id` | The sample. |
| `evidence_source` | The evidence source measured, the one the run uses. |
| `calibration_genes` | Genes whose transcripts all share one annotated end and that have at least pac_min_sample_count (default 2) read ends within calibration_max_distance (default 1000 nt) of it; at least calibration_min_genes (default 100) are needed. |
| `observations` | Read ends measured near those genes' ends. |
| `median_offset` | Median distance in nt from the annotated end to the read end; positive means reads end upstream, as in `manifest/calibration_kernel.tsv`. |
| `modal_offset` | The most common offset, in nt. |
| `central_low` | Offset at the calibration_quantile_low quantile (default 0.05), in nt. |
| `central_high` | Offset at the calibration_quantile_high quantile (default 0.95), in nt; with central_low it bounds the central interval of read ends. |
| `adjacent_boundary_fraction` | Share of read ends within 1 nt of the annotated end. |
| `poly_a_clip_fraction` | Share of read ends with a poly(A)-like soft clip at the 3' end. |
| `reproducibility` | Pearson correlation between this sample's smoothed offset profile and the average profile of the other samples (`1` with a single sample); at least calibration_min_kernel_correlation (default 0.8) is needed. |
| `classification` | `exact` (absolute median offset at most exact_max_median_abs_offset, 2 nt; central interval at most exact_max_central_width, 12 nt; adjacent fraction at least exact_min_boundary_fraction, 0.5), `proximal` (median offset at least proximal_min_median_upstream_offset, 3 nt, upstream), or `ambiguous`. |
| `reason` | The classification's reason, in words. |

### `qc/calibration_kernel_diagnostics.tsv`

One row describing the pooled calibration kernel, `manifest/calibration_kernel.tsv`. In a
proximal-tag run it warns, without stopping the run, when the kernel does not look like read
ends scattered around one site; the warning also appears at the top of the report.

| Column | Meaning |
|---|---|
| `endpoint_model` | The run's endpoint model. |
| `evidence_source` | The run's evidence source. |
| `kernel_modes` | Number of separate peaks in the kernel after smoothing it at its own resolution; read ends around one site give `1`. |
| `minimum_resolvable_separation` | The smallest shift, in nt, at which the kernel overlaps a shifted copy of itself by at most proximal_kernel_overlap_threshold (default 0.5): PACs closer than this cannot be told apart. |
| `median_central_interval_width` | Median over samples of central_high minus central_low in `qc/library_calibration.tsv`, in nt. |
| `spread_to_resolution_ratio` | median_central_interval_width divided by minimum_resolvable_separation; above 6 warns. |
| `status` | `ok`, `warning` (a proximal-tag kernel with more than one peak or a ratio above 6), or `not_applicable` (exact-boundary runs, which do not place reads through the kernel). |
| `reason` | The warning's reasons, separated by `;`; empty otherwise. |

### `qc/SAMPLE.fragment_filtering.tsv`

One row: how many of the sample's reads passed the read filters and gave a read end for the
run's evidence source, and how many spliced reads the readthrough filter recorded. Columns
appear in the order their counts were first met, so the order can differ between samples.

| Column | Meaning |
|---|---|
| `accepted_fragments` | Reads (read pairs, for paired-end samples) that passed every filter and gave a read end; for `polyA_junction`, only reads with a poly(A)-like soft clip. |
| `unique_observations` | Distinct read-end positions (chromosome, strand, coordinate): the rows of the sample's evidence table. |
| `sample_id` | The sample. |
| `layout` | Layout used, `SE` or `PE`. |
| `strandedness` | Strandedness used, `forward` or `reverse`. |
| `evidence_source` | The run's evidence source. |
| `splice_filter_enabled` | `True` when spliced reads were recorded for the readthrough filter (proximal-tag runs with constitutive_readthrough_filter true), else `False`. |
| `splice_accepted_fragments` | Reads or pairs that passed the filters and were checked for splices; `0` when the readthrough filter is off. |
| `splice_direct_edges` | Splices counted in those reads, each a jump from one aligned block directly to the next; `0` when it is off. |
| `splice_unique_continuations` | Distinct splices: the rows of the sample's splice_continuations table; `0` when it is off. |

Columns that only some runs have, depending on the layout, the filters, and which reasons
removed reads; a reason's count appears only when it removed at least one read:

| Column | Meaning |
|---|---|
| `records_examined` | Single-end samples: alignment records read from the file, including those then removed. |
| `query_groups_examined` | Paired-end samples, in place of records_examined: read names examined, each a pair with any other alignments of its reads. |
| `splice_records_examined` | Records examined for splices, in single-end samples with the readthrough filter on; `0` whenever it is off. |
| `splice_query_groups_examined` | Read names examined for splices, in paired-end samples with the readthrough filter on. |
| `splice_layout`, `splice_strandedness` | With the readthrough filter on: the splice scan's layout and strandedness, the same as layout and strandedness. |
| `unmapped`, `secondary`, `supplementary`, `qc_failed` | Records removed as unmapped, as secondary or supplementary alignments, or as marked failing quality checks. |
| `duplicate` | Records marked as duplicates, removed with exclude_duplicates (default true). |
| `low_mapq` | Records with mapping quality below min_mapq (default 20). |
| `excluded_contig` | Records on a contig in excluded_contigs (default `chrM`, `MT`, and `chrMT`). |
| `multimapped` | Records with an NH tag other than 1, removed with require_unique (default true). |
| `no_poly_a_clip` | With the `polyA_junction` evidence source: reads that passed the filters but have no poly(A)-like soft clip, so they give no read end. |
| `orphan_or_multiple_primary` | Paired-end samples: read names without exactly two primary alignments. |
| `query_name_mismatch`, `mate_designation`, `interchromosomal` | Paired-end samples: pairs removed because the mates' names differ, their read 1 and read 2 flags don't pair up, or they align to different chromosomes. A pair with a mate removed for one of the reasons above counts under that reason. |
| `improper_pair` | Paired-end samples: pairs not flagged as properly paired, removed with require_proper_pair (default true). |
| `splice_unmapped`, `splice_secondary`, `splice_supplementary`, `splice_qc_failed`, `splice_duplicate`, `splice_low_mapq`, `splice_excluded_contig`, `splice_multimapped`, `splice_orphan_or_multiple_primary`, `splice_query_name_mismatch`, `splice_mate_designation`, `splice_interchromosomal`, `splice_improper_pair` | With the readthrough filter on: the same counts for the splice scan. |

### `qc/pac_discovery.tsv`

One row summarizing how the atlas was built: candidates accepted and rejected, and what each
proximal-tag filter removed.

| Column | Meaning |
|---|---|
| `endpoint_model` | The run's endpoint model. |
| `accepted_pacs` | PACs accepted into the atlas. |
| `rejected_candidates` | Candidates rejected, each listed in `atlas/rejected_candidates.tsv.gz`. |
| `accepted_before_internal_exon_end_filter` | Candidates accepted before the internal-exon-end filter ran. |
| `internal_exon_end_filter_enabled` | `True` when that filter ran (proximal-tag runs with internal_exon_end_filter true), else `False`. |
| `internal_exon_end_rejected` | Candidates that filter rejected. |
| `accepted_before_constitutive_readthrough_filter` | Candidates accepted before the readthrough filter ran. |
| `constitutive_readthrough_filter_enabled` | `True` when the readthrough filter ran (proximal-tag runs with constitutive_readthrough_filter true), else `False`. |
| `constitutive_readthrough_min_junction_count` | The spliced-read threshold the readthrough filter used (default 2); empty when it did not run. |
| `constitutive_readthrough_rejected` | Candidates the readthrough filter rejected. |
| `minimum_resolvable_separation` | The atlas's resolution in nt: pac_cluster_radius (default 12) in exact-boundary runs, the calibrated minimum resolvable separation in proximal-tag runs. |
| `support_requirement` | The replicate-support rule in words, from pac_min_supporting_samples, for example `2 samples within one condition`. |

### `qc/SAMPLE.quantification.tsv`

One row: how the sample's read ends were assigned to the frozen atlas. An exact-boundary
read end goes to the nearest PAC within pac_cluster_radius (default 12 nt; a tie goes to the
lower coordinate), and a proximal-tag read end to the PAC the calibration kernel makes most
likely, when it is at least proximal_assignment_likelihood_ratio (default 3) times likelier
than the next.

| Column | Meaning |
|---|---|
| `accepted_fragments` | Read ends in the sample's evidence table: assigned plus unassigned. |
| `assigned_fragments` | Read ends assigned to a PAC; they make up the counts. |
| `unassigned_fragments` | Read ends left out of the counts: no PAC within reach, or ambiguous. |
| `ambiguous_fragments` | Proximal-tag read ends that two PACs explain about equally well, included in unassigned_fragments; always `0` in exact-boundary runs. |
| `sample_id` | The sample. |

*`statistics/`: the differential usage results.*

### `statistics/CONDITION_vs_CONTROL.pacs.tsv.gz`

The main result of a comparison: one row per tested PAC, with its call, the change in its
usage, the tests, its annotation, the data behind the call, and model diagnostics. Rows are
grouped by gene.

| Column | Meaning |
|---|---|
| `condition` | The treatment condition. |
| `control_condition` | Its direct control, the reference of the comparison. |
| `event_type` | The PAC's call: `gained`, `lost`, `increased_usage`, or `decreased_usage` (confirmed calls), `gained_candidate` or `lost_candidate` (candidate calls), or `none`. The rules are listed below the table. |
| `fitted_control_pau` | The model's PAU for this PAC in the control group, the mean of the control samples' fitted proportions, from 0 to 1; empty when either group has no reads at the gene. |
| `fitted_treatment_pau` | The same for the treatment group. |
| `delta_pau` | The change in usage, fitted_treatment_pau minus fitted_control_pau, from -1 to 1. |
| `delta_pau_ci_low` | Lower end of a 95% interval for delta_pau: the 2.5th percentile over dm_bootstrap_replicates (default 200) parametric bootstrap refits. Given only in genes that pass gene_fdr or, with dm_bootstrap_include_candidates (default true), have a possible gain or loss (a PAC meeting its fitted-PAU, support, and effect rules); bootstrap_status says why it is empty. |
| `delta_pau_ci_high` | Upper end of the interval, the 97.5th percentile. The interval holds the gene's precision fixed, so it is approximate: in simulations with four replicates per group, these 95% intervals held the true change about 88% of the time. |
| `pac_fdr` | Stage-wise adjusted p-value from stageR, compared with site_fdr: did this PAC change, given that its gene did? Empty when the gene fails the screen (gene_fdr above site_fdr) or the PAC has no p-value; in a gene with two tested PACs that passes, both are `0`. |
| `gene_fdr` | Benjamini-Hochberg adjusted p-value of the comparison's gene-level test, across its tested genes: did the gene's PAC usage change at all? |
| `pac_pvalue` | Raw p-value of the PAC-level test, DRIMSeq's beta-binomial test of this PAC's share against the rest of its gene; empty when no stable fit gives one. |
| `gene_pvalue` | Raw p-value of the gene-level test. |
| `pac_likelihood_ratio` | The PAC-level test's likelihood-ratio statistic; larger means stronger evidence of a change. |
| `pac_degrees_of_freedom` | Its degrees of freedom, 1 for one treatment against its control. |
| `gene_region` | As in the atlas: `last_exon`, `internal_exon`, `intron`, or `downstream_of_gene`. |
| `last_exon_locus` | As in the atlas: the PAC's last exon, 1-based; empty when it has none. |
| `confidence` | As in the atlas: `high`, `moderate`, or `low`. |
| `internal_priming_flag` | As in the atlas: `True` when the PAC is flagged for possible internal priming. |
| `known_pac` | As in the atlas: `True` when the PAC matches a known PAC. |
| `known_rescue_only` | As in the atlas: `True` when the PAC was kept only by the relaxed known-PAC rule. |
| `primary_pas_motif` | As in the atlas: the primary PAS hexamer in DNA letters; empty when none. |
| `primary_pas_motif_rna` | As in the atlas: the hexamer in RNA letters; `none` when none. |
| `primary_motif_class` | As in the atlas: the hexamer's class. |
| `control_supporting_samples`, `treatment_supporting_samples` | Samples in each group with at least one read at this PAC; with at least event_min_supporting_samples (default 2), the PAC is detected in that group. |
| `control_gene_total`, `treatment_gene_total` | Reads at the gene's tested PACs, summed over each group's samples. A gained call needs control_gene_total, and a lost call treatment_gene_total, of at least min_gene_total (default 20). |
| `raw_control_counts`, `raw_treatment_counts` | Each sample's raw count at this PAC, as `sample=count` pairs separated by commas. |
| `observed_control_pau`, `observed_treatment_pau` | Each sample's observed PAU, as `sample=value` pairs: the PAC's reads over the reads at the gene's tested PACs, so it can differ from `counts/observed_pau.tsv.gz`; `NA` for a sample without reads there. |
| `effect_exceeds_threshold` | `TRUE` when the absolute delta_pau is at least min_abs_delta_pau (default 0.1). |
| `dominant_pac_control`, `dominant_pac_treatment` | pac_id of the gene's tested PAC with the highest fitted PAU in each group; the same on every row of the gene. |
| `control_active_pacs`, `treatment_active_pacs` | Number of the gene's tested PACs that are active in each group, with fitted PAU of at least active_pac_min_pau (default 0.10). |
| `model_status` | How the gene was fitted: `fitted`; `fitted_with_zero_count_stabilization` (a PAC had no reads in a whole group, or the plain fit failed, so the results are medians of refits with small seeded values in place of zero counts, and the precision and p-values depend on random_seed); `fit_unavailable` (no usable fit); or `group_without_counts` (the control or treatment group has no reads at the gene, so it is not tested in this comparison). |
| `precision` | The gene's Dirichlet-multinomial precision, shared by the conditions of the family: higher means the replicates vary less around the fitted usage. |
| `alpha_control`, `alpha_treatment` | The fitted Dirichlet parameters for this PAC in each group, fitted PAU times precision. |
| `stabilization_successes` | For a stabilized gene, how many of the dm_zero_sensitivity_repeats refits (default 5) gave a usable result for this PAC; `0` otherwise. |
| `stabilization_delta_pau_spread` | For a stabilized gene, the range of delta_pau across those refits, largest minus smallest; empty otherwise. |
| `zero_boundary_unstable` | `TRUE` when the refits disagree: fewer than 80% succeed, the change flips direction (beyond 0.001), or delta_pau spreads by more than dm_zero_max_delta_pau_spread (default 0.02). Such a PAC gets no p-value, so no confirmed call. |
| `zero_boundary_reason` | Which of those: `insufficient_successes`, `direction_change`, or `delta_spread`, joined by `;`; empty when stable. |
| `bootstrap_status` | `ok` (interval given), `insufficient_successes` (fewer than dm_bootstrap_min_success_fraction, default 0.8, of the refits succeeded), `fit_unavailable` (the gene's fit could not seed a bootstrap), `not_selected` (the gene was not bootstrapped), or `disabled` (dm_bootstrap_replicates is 0). |
| `bootstrap_successes` | Bootstrap refits that succeeded; `0` when none ran. |
| `bootstrap_perturbed` | Of those, refits whose simulated counts left a group without reads, and which were stabilized the same way as the family fit. |
| `exploratory_insufficient_replicates` | `TRUE` when a condition in this comparison's family has fewer than min_replicates_per_condition (default 2) samples (insufficient_replicates_policy `warn`); gained and lost calls are then only candidates. |

How `event_type` is decided, with T the min_abs_delta_pau (default 0.1):
- `gained`: detected in the treatment but not in the control; fitted_control_pau at most
  event_max_control_pau (default 0.01) and fitted_treatment_pau at least
  event_min_treatment_pau (default 0.05); delta_pau at least T; gene_fdr at most gene_fdr
  and pac_fdr at most site_fdr (both default 0.05); control_gene_total at least
  min_gene_total (default 20), so the control had the reads to show the PAC; and the PAC is
  not zero_boundary_unstable, not low confidence, not flagged for internal priming, and not
  in an exploratory comparison. It means not detected in the control, not necessarily
  absent there. At the defaults the delta_pau rule is the stricter one, so
  fitted_treatment_pau must be at least 0.10.
- `lost`: the mirror image, detected in the control but not in the treatment, with
  treatment_gene_total at least min_gene_total.
- `gained_candidate`, `lost_candidate`: the detection, fitted-PAU, and delta_pau rules of a
  gain or loss hold, but another requirement does not. They are not statistically
  confirmed, so do not count them as gains or losses.
- `increased_usage`, `decreased_usage`: detected in both groups, delta_pau at least T (or
  at most -T), and both FDRs pass.
- `none`: anything else, including significant changes smaller than T.

### `statistics/CONDITION_vs_CONTROL.calls.tsv.gz`

Same columns as `statistics/CONDITION_vs_CONTROL.pacs.tsv.gz`. It holds the rows whose
`event_type` is not `none`: every PAC with a confirmed or candidate call in the comparison.

### `statistics/CONDITION_vs_CONTROL.genes.tsv.gz`

One row per tested gene in a comparison: its gene-level test, its gene-level events, its
APA pattern, the shift in its usage, and the numbers behind the pattern.

| Column | Meaning |
|---|---|
| `condition` | The treatment condition. |
| `control_condition` | Its direct control. |
| `dominant_switch` | `TRUE` when the gene passes gene_fdr (default 0.05) and its most-used tested PAC, the one with the highest fitted PAU, differs between control and treatment. |
| `active_pacs_change` | `more` or `fewer` when the gene has more or fewer active PACs (fitted PAU of at least active_pac_min_pau, default 0.10) in the treatment than in the control; `none` when the number is the same, the gene does not pass gene_fdr, or it has no fitted usage. |
| `apa_pattern` | How the gene's usage moved: `intronic_gain`, `intronic_loss`, `alternative_last_exon`, `utr_shortening`, or `utr_lengthening`, each possibly with the suffix `_potential_internal_priming`; `unclassified_change`; or `none`. Several patterns are joined by `;`, and the rules are listed below the table. |
| `shift_direction` | For a gene with a confirmed PAC call: `distal` when the PAC its usage moved to lies 3' of the PAC it moved from, along the transcript, and `proximal` when it lies 5'. Empty, as are the other `shift_` columns, for a gene without a confirmed call. |
| `shift_from_pac_id` | The PAC the gene's usage moved from: its confirmed `lost` or `decreased_usage` PAC with the largest fall in fitted PAU. Without one, its PAC with the largest fitted fall, which has no confirmed call of its own. Ties go to the lower `pac_fdr`, then the PAC ID. |
| `shift_from_gene_region` | That PAC's `gene_region`. |
| `shift_from_event_type` | That PAC's `event_type`: `lost` or `decreased_usage`, or, in a gene without a confirmed loss, `none` or a candidate call. |
| `shift_to_pac_id` | The PAC the gene's usage moved to: its confirmed `gained` or `increased_usage` PAC with the largest rise in fitted PAU, or without one, its PAC with the largest fitted rise. Ties go as for `shift_from_pac_id`. |
| `shift_to_gene_region` | That PAC's `gene_region`. |
| `shift_to_event_type` | That PAC's `event_type`: `gained` or `increased_usage`, or, in a gene without a confirmed gain, `none` or a candidate call. |
| `delta_intronic_share` | Change, treatment minus control, in the share of the gene's fitted usage at its tested PACs in introns and internal exons (the upstream region), from -1 to 1; `0` when it has no tested PAC there, and empty without fitted usage. |
| `delta_utr_distal_share` | Change in the distal PAC's share of the main last exon's usage, the 3' counterpart of DaPars' PDUI: negative means a shorter 3' UTR. Empty unless that exon has 2 or more tested PACs and at least event_min_treatment_pau (default 0.05) of the gene's usage in both groups. |
| `last_exon_switch` | How much usage moved between last exons: the smaller of the largest gain and the largest loss in any last exon's share of the gene (negative when every last exon moves the same way). Empty when fewer than 2 last exons have tested PACs. |
| `gene_fdr` | Benjamini-Hochberg adjusted p-value of the gene-level test, across this comparison's tested genes: the screen for every call. |
| `gene_pvalue` | Raw p-value of the gene-level test: does the treatment change how the gene's reads are split among its tested PACs? |
| `gene_likelihood_ratio` | The test's likelihood-ratio statistic; larger means stronger evidence. |
| `gene_degrees_of_freedom` | Its degrees of freedom: the gene's tested PACs minus 1. |
| `model_status` | How the gene was fitted: `fitted`, `fitted_with_zero_count_stabilization`, `fit_unavailable`, or `group_without_counts`, as in the `.pacs` table. |
| `stabilization_successes` | For a stabilized gene, how many refits gave a usable gene-level test; `0` otherwise. |
| `exploratory_insufficient_replicates` | `TRUE` when the comparison is exploratory, as in the `.pacs` table. |

How `apa_pattern` is decided, with T the apa_pattern_min_change (default 0.1):
- A gene can have a pattern only when it passes gene_fdr, has fitted usage in both groups,
  and has at least one confirmed call, or a withheld call as described below; other genes
  are `none`. A call supports a pattern when it is confirmed and its PAC is not low
  confidence.
- The upstream region is the gene's tested `intron` and `internal_exon` PACs. Its
  `last_exon` and `downstream_of_gene` PACs are grouped into last exons by
  `last_exon_locus`; the main last exon has the most fitted usage over both groups (a tie
  goes to the 3'-most), and its distal PAC is its most 3' tested PAC.
- `intronic_gain`: delta_intronic_share at least T, with a supporting increase (`gained` or
  `increased_usage`) in the upstream region, as in intronic polyadenylation.
  `intronic_loss` is the mirror image.
- `alternative_last_exon`: last_exon_switch at least T, with a supporting increase in the
  last exon that gained most and a supporting decrease in the one that lost most.
- `utr_shortening`: delta_utr_distal_share at most -T, with a supporting increase at a more
  proximal PAC of the main last exon or a supporting decrease at its distal PAC.
  `utr_lengthening` is the mirror image. When usage also moves into or out of the main last
  exon (delta_intronic_share or last_exon_switch at least T in size), only calls against
  that movement count.
- `_potential_internal_priming`, added to any of these: the pattern holds only once calls
  on PACs flagged for possible internal priming count too. With
  potential_internal_priming_withheld_calls (default true), a flagged PAC's
  `gained_candidate` or `lost_candidate` also counts when the flag alone withheld a
  confirmed call, except in exploratory comparisons.
- `unclassified_change`: the gene qualifies, but no pattern applies.
- Several patterns are listed in the order above, the supported ones first, then those
  with the `_potential_internal_priming` suffix.

### `statistics/FAMILY.gene_omnibus.tsv.gz`

One row per tested gene of a comparison family: a family-wide test of whether the gene's
usage differs among any of the family's conditions, every treatment against the control at
once. It is descriptive only; calls use each comparison's own gene test.

| Column | Meaning |
|---|---|
| `family` | The comparison family, named after its control condition. |
| `gene_fdr` | Benjamini-Hochberg adjusted p-value of the family-wide test, across the family's tested genes. |
| `gene_pvalue` | Raw p-value of the family-wide test. |
| `gene_likelihood_ratio` | Its likelihood-ratio statistic; larger means stronger evidence. |
| `gene_degrees_of_freedom` | Its degrees of freedom: the gene's tested PACs minus 1, times the number of treatments. |
| `model_status` | As in the `.genes` table; `group_without_counts` when any condition of the family has no reads at the gene. |
| `stabilization_successes` | For a stabilized gene, how many refits gave a usable test; `0` otherwise. |
| `exploratory_insufficient_replicates` | `TRUE` when a condition of the family has fewer than min_replicates_per_condition (default 2) samples. |

### `statistics/FAMILY.statistical_filtering.tsv.gz`

One row per atlas PAC for each comparison family: whether it was tested there, and if not,
why. The filters use the family's samples without regard to condition.

| Column | Meaning |
|---|---|
| `family` | The comparison family, named after its control condition. |
| `tested` | `TRUE` when the PAC entered the family's model, so it appears in each of the family's `.pacs` tables; `FALSE` otherwise. |
| `reason` | Why not, joined by `;`: `no_gene_assignment` (intergenic), `ambiguous_gene_assignment` (several genes), `site_count<N` (fewer than min_site_count reads, default 5), `supporting_samples<N` (reads in fewer than min_test_supporting_samples samples, default 2), or `sample_usage<F` (at least min_site_usage of the gene's reads, default 0.01, in fewer than min_site_usage_samples of the family's samples, by default as many as its smallest condition has less min_site_usage_dropouts (default 0, never below 2, or below the smallest condition if it has fewer), counting only samples with min_site_usage_gene_reads reads at the gene, default 10). A PAC that passes these can still fail with its gene, as `gene_total<N` (fewer than min_gene_total reads at the gene, default 20) or `fewer_than_2_testable_pacs`; empty when tested. |

### `statistics/fitted_pau.tsv.gz`

One row per tested PAC and comparison, for every comparison of the run in one file: the
fitted usage and model parameters of the `.pacs` tables.

| Column | Meaning |
|---|---|
| `condition` | The treatment condition. |
| `control_condition` | Its direct control. |
| `fitted_control_pau` | The PAC's fitted PAU in the control group, as in the `.pacs` table. |
| `fitted_treatment_pau` | Its fitted PAU in the treatment group. |
| `delta_pau` | fitted_treatment_pau minus fitted_control_pau. |
| `model_status` | How the gene was fitted, as in the `.pacs` table. |
| `precision` | The gene's Dirichlet-multinomial precision, as in the `.pacs` table. |
| `alpha_control`, `alpha_treatment` | The fitted Dirichlet parameters for this PAC in each group, fitted PAU times precision. |

### `statistics/gene_precision.tsv.gz`

One row per tested gene and comparison family: the gene's Dirichlet-multinomial precision,
which says how closely replicates agree.

| Column | Meaning |
|---|---|
| `family` | The comparison family, named after its control condition. |
| `precision` | The gene's precision, shared by the family's conditions: higher means the replicates vary less around the fitted usage. |
| `model_status` | How the gene was fitted in the family: `fitted`, `fitted_with_zero_count_stabilization`, `fit_unavailable`, or `group_without_counts` (a condition of the family has no reads at the gene). |

## Other files

- `atlas/pacs.v1.bed.gz`: the atlas as BED6 (chrom, start, end, pac_id, score, strand), for
  genome browsers; the score is total_count, capped at 1000.
- `atlas/pacs.v1.sha256`: the SHA-256 checksum of `atlas/pacs.v1.metadata.tsv.gz`, which
  identifies the frozen atlas.
- `manifest/library_resolution.json`: the run's protocol decision: library profile, evidence
  source, endpoint model, the kernel's offset range, the calibration version, and each
  sample's layout and strandedness.
- `manifest/resolved_params.yaml`: every parameter's value for the run: the defaults,
  overridden by analysis.yaml, then by the command line.
- `manifest/run_manifest.json`: the PACusage version, the number of samples, the
  conditions, and the checksum of `manifest/resolved_params.yaml`.
- `figures/CONDITION_vs_CONTROL.volcano.pdf` and `.png`: each tested PAC's change in fitted
  PAU against its PAC-level p-value, with confirmed calls colored and candidates open. Each
  gene is labeled once, at `shift_to_pac_id`, or at `shift_from_pac_id` when only that
  side has a confirmed call.
- `figures/CONDITION_vs_CONTROL.distal_usage.pdf` and `.png`: each tested gene's distal PAC
  usage, control against treatment, colored by the gene's APA pattern.
- `figures/CONDITION_vs_CONTROL.shifts_by_gene_region.pdf` and `.png`: each gene with a
  confirmed call, once, by `shift_from_gene_region` and `shift_to_gene_region`, with
  proximal shifts left of zero and distal shifts right.
- `figures/event_counts.pdf` and `.png`: genes by their PAC calls and gene events per
  comparison, and the genes per APA pattern.
- `figures/apa_pattern_grid.pdf` and `.png`: up to 50 genes with an APA pattern in two or
  more comparisons, against the comparisons.
- `figures/concordance.pdf` and `.png`: for pairs of comparisons (up to 15 panels), each
  shared PAC's change in PAU in one against the other.
- `figures/concordance_matrix.pdf` and `.png`: the Pearson r of every pair of comparisons.
- `figures/effect_vs_coverage.pdf` and `.png`: each PAC's change in fitted PAU against the
  reads at its gene in the less-covered group, with min_gene_total marked.
- `figures/pau_pca.pdf` and `.png`: the samples on the first two principal components of
  observed PAU, colored by condition.
- `tracks/SAMPLE.plus.3prime_evidence.bedGraph.gz`: the sample's plus-strand read ends as a
  bedGraph, the value being the number of reads ending at each position, for genome
  browsers.
- `tracks/SAMPLE.minus.3prime_evidence.bedGraph.gz`: the same for the minus strand; values
  stay positive, and the strand is in the file name.
- `report/index.html`: a self-contained HTML report of the run's checks, summaries, top
  results, and figures, which stays readable when copied elsewhere.
- `pipeline_info/`: Nextflow's execution report, timeline, trace, and DAG for the latest
  run, replaced by each run.
- `prepared_reference/`: only with save_prepared_reference (default false), the prepared
  FASTA link, `genome.fa`, and its index.
- `prepared_alignments/`: only with save_prepared_alignments (default false), each sample's
  prepared BAM or CRAM and its index.
