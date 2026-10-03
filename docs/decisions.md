# Decisions

Decisions about how PACusage is built and maintained, newest first. The
project lead accepts decisions. `fable-overseer` may add entries with status
`proposed`, which stay proposed until the project lead accepts or rejects them.
To change an accepted decision, add a new entry that supersedes it, and only
with the project lead's approval.

Entry format: a dated heading, a status line, the decision, and the reason.

## 2026-10-03: A comparison tests only genes with depth in both groups

Status: accepted (project lead)

In one treatment, many of the lead's calls came from genes nearly silent
after treatment. Two PACs with one read each read as an even split against
the control's well-supported 10/90. The family filter, blind to condition,
admits such genes on the control's reads, and only gained and lost calls
checked coverage. The lead kept the remedy to PAU analysis: no
gene-expression test and no shrunken effect sizes.

**Decision** (version 0.9.0)
- **The rule:** a comparison tests a gene only when its control and its
  treatment each have `event_min_supporting_samples` (2) samples with
  `min_site_usage_gene_reads` (10) reads at the gene's tested PACs. A
  smaller group needs all of its samples, so a one-replicate exploratory
  group needs its one sample deep. No parameters were added.
- **What changes:** genes failing the rule leave the comparison before BH and
  stageR, and leave its `.pacs`, `.genes`, `.calls`, and `fitted_pau` rows.
  The family fit, its omnibus test, precision, stabilization, the family
  filter, and every call rule are unchanged, as is the gained/lost
  `min_gene_total` coverage rule, which the new rule implies at the
  defaults in groups of two or more samples.
- **The new table:** `CONDITION_vs_CONTROL.genes_without_depth.tsv.gz` lists
  the skipped genes with each group's reads, samples, and mean CPM. Its
  statuses are `turned_off` and `turned_on` (one group has reads in fewer
  samples than depth needs, where the other's mean CPM would have given it
  depth), and otherwise `too_low_in_treatment`, `too_low_in_control`, or
  `too_low_in_both`. `event_counts` counts them in a third plot, and the
  report shows each comparison's table.

**Why**
- **Independence:** the family filter reads PAC shares, so it must be blind
  to the labels (2026-10-01). This rule reads the labels but only per-sample
  gene totals, never how they split among PACs. The Dirichlet-multinomial
  test conditions on those totals, so under the null the rule is independent
  of the test (Bourgon et al., 2010), as far as the test's asymptotic null
  holds. It is a comparison-level step after the family fit, beside the
  label-blind filter, and supersedes none of the 2026-10-01 entries.
- **Simulation S8, three seeds:**
  - 71% to 74% of genes whose treatment had a few reads at an even split
    had p <= 0.05, and none is tested now;
  - the null genes kept had p <= 0.05 in 2.7% to 5.1% of cases;
  - every shift at normal depth was still called.
- **Turned on and off:** the guard against shallow libraries makes these a
  description, not a test: not detected, at a depth where the gene would
  have shown.
- **Dropped:** a gene-expression test with edgeR, which would be a separate
  analysis, and the shrunken effect sizes. The depth rule removes the
  low-read calls directly.

## 2026-10-03: Figure text kept on the page and labels kept apart

Status: accepted (project lead; supersedes the 10,000-iteration detail of
"Figure labels placed by ggrepel, and the figures reworked")

The lead found the PCA legend cut off at the right edge, its last condition
missing, and asked for an audit of every figure.

**Decision** (version 0.8.2)
- The PCA's legend sits beside the panel, in columns of at most 20
  conditions. This applies the 2026-10-01 rule that one-panel figures put
  their legends to the right; the PCA was its one exception. The figure is
  5.5 inches tall, and 5.5 inches wide plus 0.4 + 0.075 inch per character
  of the longest condition name for each legend column.
- Every figure's right margin is 14 pt, which holds half an axis label.
- ggrepel places labels with force 32 and at most 100,000 iterations, from
  the same seed and with no time limit. This replaces the 10,000 iterations
  in the 2026-10-01 entry.
- Comparison titles over 30 characters take two or more lines, on facet
  strips and on the count axes. A condition name over 30 characters breaks
  after an underscore, hyphen, or period.

**Why**

The audit drew every figure from synthetic tables in the lead's design (13
conditions, 12 comparisons, 10,000 genes, three seeds), from the lead's real
PCA coordinates, and from a set of long, wide condition names. It checked
for three faults: anything drawn past a figure's edges (each figure drawn on
a larger page), text that a panel or strip cut off (the text drawn alone,
with and without clipping), and overlapping ggrepel labels.
- The PCA's 13-condition legend sat under the panel, wider than the
  5.5-inch figure.
- The last x-axis label is centered on the panel's right edge. With 439
  genes each way, the shifts figure's 600 crossed the page's edge.
- At ggrepel's default force and 10,000 iterations, 28 to 40 label pairs
  overlapped on each set's 24 volcano and distal plots, one label's halo
  hiding part of the other. With the new settings, overlapping pairs fell
  from 103 to 2 among 2,880 labels, and the 2 left nearly touch but stay
  legible. The median leader line grew from about 0.16 to 0.23 inch, while
  the longest lines, and the number over an inch, fell.
- Strip titles of up to 36 characters stayed on one line, but a third of
  the figure holds about 32 capitals. A long comparison title on the event
  counts' axis also narrowed the pattern panels until a facet title was cut.
- Limit: when both names in a comparison exceed 30 characters, its
  count-axis label takes up to four lines and can touch its neighbours.
  Only the stress names reach this.

## 2026-10-03: The pattern counts include all tested PACs

Status: accepted (project lead)

The lead asked for a panel of all tested PACs, unflagged and flagged, beside
the two that split the APA patterns by whether flagged PACs support them.

**Decision** (version 0.8.1)
- `event_counts` counts genes per APA pattern in three panels: all tested
  PACs, the patterns that unflagged PACs support, and those that only
  flagged PACs support.
- Under all tested PACs, a pattern counts whether or not it needs flagged
  PACs. A gene never has both forms of one pattern, so each count there is
  the sum of the other two.
- It regroups the published patterns, as the entry below allows, and
  changes no table.

## 2026-10-02: Figures may regroup published calls, but apply no rule of their own

Status: accepted (project lead; proposed by fable-overseer at the 0.8.0
review)

How to read the AGENTS.md invariant "PLOT_FIGURES draws the calls in the
`.pacs` and `.genes` tables and never makes its own":

**Rule**
- A figure may count, group, or intersect the calls the tables publish. In
  0.8.0, `gene_call_counts` in `scripts/plot_usage_figures.R` counts each gene
  once by which `event_type` values its PACs carry (a gain and a loss, a
  gain, a loss, changed usage only, candidates only). A reader can reproduce
  such a grouping from the tables with no parameter or threshold, so it is a
  legend, not a call.
- A figure may not apply a threshold, parameter, test, or positional rule to
  produce a category, direction, or label the tables do not carry. Those
  belong to `fit_usage_model.R`, as the shift does in 0.8.0.

**Why.** The invariant keeps figure-only results from passing for calls. A
lossless regrouping of published calls creates none; a new rule would.

## 2026-10-02: Each gene's shift is named once

Status: accepted (project lead, option 1 of three)

The lead noticed that almost every `increased_usage` call has a matching
`decreased_usage` call at another PAC of the same gene, as the biology
predicts: PAU are shares of a gene, so usage gained at one PAC is lost at
others. The figures counted and labeled calls, so they looked balanced
between increases and decreases even when the treatment moved usage one
way.

**Decision** (version 0.8.0)
- The `.genes` tables name each gene's shift once, in seven `shift_` columns
  after `apa_pattern`.
  - The to-PAC is the confirmed gain or increase with the largest change in
    fitted PAU. The from-PAC is the confirmed loss or decrease with the
    largest fall.
  - A side without a confirmed call names the gene's PAC with the largest
    fitted change that way, and its event type shows it has no confirmed
    call.
  - Ties go to the lower `pac_fdr`, then the PAC ID.
  - `shift_direction` is `distal` or `proximal`.
- `calls_by_gene_region` becomes `shifts_by_gene_region`. It counts genes by
  where their usage moved from and to, with proximal shifts left of zero and
  distal ones right, and every comparison's figure has the same rows.
- `event_counts` counts genes, each once, by their PAC calls.
- The volcano stays one point per PAC and labels each gene once, at the PAC
  its usage moved to, or at the one it moved from when only that side has a
  confirmed call.
- No call, pattern, or existing column changes.

**Not chosen**
- Showing only the gaining side in the figures. It is simpler, but hides
  treatments whose main effect is shutting a site off: in a gene with three
  or more PACs, the loss can be the event, with the gain spread too thin to
  call.
- Keeping the figures and adding a note to each that every shift appears as
  one increase and one decrease.

## 2026-10-01: Active PACs default to the minimum change

Status: accepted (project lead, option 2 of three)

The lead noticed that with the defaults an active PAC could never be called
lost. `active_pac_min_pau` was 0.05, and a `gained` or `lost` call needs a
change of `min_abs_delta_pau`, 0.10, besides its usage rules. A PAC between
0.05 and 0.10 that disappeared counted as one fewer active PAC, but got no
PAC call.

**Decision** (version 0.7.0)
- The `active_pac_min_pau` default rises to 0.10, the same as
  `min_abs_delta_pau`. A PAC whose usage falls from active to none now
  changes by enough to be called.
- One edge remains: a PAC at 0.10 in one group and up to the 0.01 detection
  bar in the other changes by as little as 0.09.
- `event_min_treatment_pau` stays at 0.05. At the defaults the change rule is
  stricter, so a gained or lost PAC needs at least 0.10 in its group. The
  docs now say so.
- `active_pacs_change` results change at the defaults; 0.05 keeps the old
  counts.

**Not chosen**
- Dropping the change requirement from gained and lost calls, so that a site
  switching on or off is called at any size above the presence bar.
- Leaving the defaults and only documenting how they combine.

## 2026-10-01: A dropout allowance for the usage filter

Status: accepted (project lead)

The lead asked whether requiring the smallest condition's size less one, to
allow for a bad replicate, is fair, and then asked for a setting.

**Rule** (version 0.6.1)
- `min_site_usage_dropouts` (integer, default 0) lowers the default count,
  the smallest condition's size, by that many samples. With 4 replicates per
  condition and 1 dropout, a PAC needs 3 qualifying samples.
- It never takes the count below 2, or below the smallest condition when
  that has fewer, so an allowance never makes the rule stricter.
- An explicit `min_site_usage_samples` is used as it is. VALIDATE_INPUTS
  stops when one is set beside a non-zero allowance, so neither is silently
  ignored.
- Results are unchanged at the default.

**Why it is fair.** The count still depends only on condition sizes, so the
filter stays blind to the comparisons. In S7's design with one dropout (2 of
11 samples, three seeds): 84-85% of the borderline null PACs were admitted,
5-7% of those had p <= 0.05, and 3-9% of the treatment's gene calls were
false.

**Cost.** It admits more borderline PACs, and it cannot tell 3 of 4 samples
in one condition from 1 + 1 + 1 across three. It is more lenient than
standard practice at small replicate numbers. The DRIMSeq workflow and
edgeR's `filterByExpr` (4.8.2) ask for the full smallest group, and edgeR
relaxes that only when every group has more than 10 samples, keeping at
least 70% of the smallest group.

## 2026-10-01: The usage filter counts samples across the family

Status: accepted (project lead)

This supersedes "The usage filter asks for consistent usage within a
condition", below. After 0.5.0 shipped that rule, the lead returned to the
option first recommended: the DRIMSeq workflow's rule (Love et al., 2018),
which asks for a 10% share in at least as many samples as the smallest group.

**Rule** (version 0.6.0)
- A sample qualifies when the gene has `min_site_usage_gene_reads` reads in
  it (default 10) and the PAC has `min_site_usage` of them.
- A PAC is tested when `min_site_usage_samples` of the family's samples
  qualify, from any of its conditions. By default (empty) that is the size of
  the family's smallest condition. It is never more than the family's
  samples, and there is no floor of 2.
- `min_site_usage_replicates` is gone, and is ignored if still set. The
  filtering reason is `sample_usage<F`.
- The read-count, supporting-sample, and gene-total filters are unchanged,
  and so is `apa_pattern_min_change`.

**Why.** The rule uses the condition sizes but not which sample belongs to
which condition, so the filter is blind to the comparisons. That keeps the
independence argument for filtering (Bourgon et al., 2010) while still
keeping a site used throughout one condition, however many conditions share
the control. At the lead's 0.10, the two false Tacc2 PACs go:
130338300's best sample is 0.036, and only one sample of 130353825 reaches
0.10. The PHA alternative last exon stays, at 0.35 or more in all four PHA
samples.

**S7** now checks calibration, guarded at S1's tolerances (0.12 of PAC tests
and 0.15 of genes with p <= 0.05), with S4's call-level guards. The
accepted-cost ceiling rule below still stands, but it has no cost to cap
here. At `min_site_usage` 0.10 and `min_abs_delta_pau` 0.05, over five seeds:
- 60-65% of the borderline null PACs at 8% were admitted;
- 5-8% of those admitted had p <= 0.05;
- the genes keeping one had gene p <= 0.05 in 4-9% of cases, against 1-6% for
  genes that dropped it, in line with S1's 7.5% for unfiltered null genes.
  The gap is composition: the filter keeps fewer of the deep, precise genes,
  whose null p-values are rarely small. Within each depth and precision
  stratum the two groups agree (Mantel-Haenszel common odds ratio 0.75-1.60,
  p 0.55-0.97, seeds 11-13);
- 2-9% of T1's gene calls were false;
- null genes keeping a borderline PAC were called in 12 of 1,236 cases
  (1.0%), and other null genes in 6 of 764 (0.8%);
- calls at a change of 0.10 had at most one null gene fewer than calls at
  0.05.

## 2026-10-01: The usage filter asks for consistent usage within a condition

Status: accepted (project lead)

The lead's Tacc2 review found intronic PACs tested at 2-4% of their gene's
reads. Raising `min_site_usage` would not have removed them, because the
share was pooled over the comparison family. Pooling also dilutes a
treatment-specific site, more so the more conditions share the control.

**Rule**
- A PAC is tested when, in at least one condition of the family, at least
  `min_site_usage_replicates` of the replicates give it `min_site_usage` of
  the gene's reads.
- The setting is a fraction of the condition's replicates, rounded up, or a
  whole number capped at the condition's size. It is never below 2 unless the
  condition has one replicate. The default, 0.75, asks for 2 of 2, 3 of 3,
  3 of 4, 4 of 5, and 5 of 6.
- A replicate counts only when the gene has `min_site_usage_gene_reads`
  reads in it, by default 10, as in the DRIMSeq workflow.
- The read count, supporting-sample, and gene-total filters stay pooled.
- The filtering reason is `replicate_usage<F`.

**Pattern threshold.** APA patterns take `apa_pattern_min_change` (default
0.10). The lead plans 0.10 inclusion with PAC calls at a change of 0.05
(`min_abs_delta_pau`), and patterns keep the larger change.

**Cost.** The filter now uses condition labels, so a PAC near the threshold
that it admits has an optimistic p-value. The lead accepted this, the same
kind of choice as discovery's within-condition support rule, with the cost
measured first. S7 runs null genes with a borderline PAC at 8% beside 100
genes shifted in T1, at `min_site_usage` 0.10, replicates 0.75, and
`min_abs_delta_pau` 0.05. Over five seeds:
- 14-19% of the borderline null PACs were admitted;
- 10-21% of those admitted had p <= 0.05, at a nominal 5%;
- the genes keeping one had gene p <= 0.05 in 6-18% of cases, against 2-6%
  for genes that dropped it;
- 4-9% of T1's gene calls were false;
- null genes keeping a borderline PAC were called in 8 of 332 cases (2.4%),
  and other null genes in 14 of 1,668 (0.8%);
- calls at a change of 0.10 had at most one null gene fewer than calls at
  0.05.

S7 guards the call-level error as S4 does, a false share of gene calls up to
0.25. It also caps the accepted cost at 0.30, as the next entry describes.

**Confirmation.** The lead confirmed the rule and its defaults from these
numbers, and accepted S7's measurement at the planned settings without a
second scenario at the defaults. A scratch run at the defaults
(`min_site_usage` 0.01 with the borderline PAC at 0.8%, `min_abs_delta_pau`
0.10, three seeds) found the cost near nominal:
- 7-9% of the borderline null PACs were admitted;
- 0-13% of those had p <= 0.05, about 5% on average;
- 5-9% of T1's gene calls were false;
- no admitted borderline PAC was called.

**Not chosen.** The lead set aside the alternative of a share in most of the
family's samples regardless of condition, which is label-blind but admits
fewer treatment-specific sites.

## 2026-10-01: A measured cost gets an explicit ceiling in the simulation suite

Status: accepted (project lead), with the ceiling at 0.30

S7 was planned to guard the admitted borderline PACs' p <= 0.05 rate at the
suite's 2-3x nominal tolerances, as S1 guards unselected nulls. Over five
seeds that rate was 10-21%, so the guard was dropped and the rate recorded
only. The overseer's view: a 2-3x calibration guard is the wrong instrument
for a subset selected by construction, but a metric with no guard protects
nothing. A later change to the filter that worsens the selection effect (the
floor of 2, the gene-read floor, or the rounding) should trip a test.

**Proposed rule.** When the suite measures an accepted cost rather than a
calibration claim, it guards the metric at an explicit ceiling, named in the
test as the accepted cost, not as calibration. For S7: the admitted
borderline PACs' p <= 0.05 rate, and the gene p <= 0.05 rate of null genes
keeping one, each at most 0.30 (six times nominal; the five-seed maximum
was 0.21 on 55-75 PACs, whose binomial 99% upper bound is about 0.27). The
call-level guards stay as they are. The test's title should name what it
guards; "keeps calls near nominal error" overstates a 0.25 false-share
ceiling.

## 2026-10-01: Figure labels placed by ggrepel, and the figures reworked

Status: accepted (project lead)

The lead found the figures unprofessional at a real run's scale: overlapping
and missing labels, a clipped legend, unformatted numbers, and a colorbar with
one tick. They were reworked against synthetic tables at that scale: 12,000
genes, 60,610 PACs, and six comparisons with long condition names.

**Dependency**
- The lead approved ggrepel, declared in `envs/pacusage.yml`
  (`r-ggrepel>=0.9.7`). conda-forge has 0.9.6 and 0.9.8, so the pin resolves
  to 0.9.8, the version tested.
- ggrepel stops searching for label positions after half a second by
  default, so the layout would depend on the machine's speed. The figures set
  no time limit and a fixed seed; the search then ends when no labels overlap
  or after 10,000 iterations, the same way on every run.
- `software_versions.tsv` records ggrepel beside ggplot2.

**Changes**
- Volcano, distal-usage, and PCA labels are repelled from each other and
  their points. Each starts nudged off its point, so a line always joins the
  two; one touching its point could sit between two points. On the
  distal-usage plot, labels start toward the open middle of their half of the
  plot, since the labeled genes sit on the edges.
- Volcano labels go on each gene's most significant PAC in each direction,
  so a gene whose usage moves between its PACs can be labeled twice. Among
  equal p-values, such as several of 0, the larger changes get the labels.
- One-panel figures put their legends to the right; faceted figures stack
  theirs under the panels.
- The uncalled-PAC density uses a log scale.
- Counts have thousands separators.
- Captions explain the dashed lines and marks.
- `concordance` became a matrix. Each pair's earlier comparison names its
  column and the later one its row, so the panels fill the lower triangle.
  Two-line strip titles replace per-panel titles, which were clipped.
  `concordance_matrix` draws each pair once, in the same cells.
- PCA labels drop the condition's name from sample IDs when every ID starts
  with it, and take a darker shade of the condition's color.

## 2026-10-01: Merged donor peaks and internal priming stay as they are

Status: accepted (project lead)

The lead judged two intronic Tacc2 PACs in their Plasmidsaurus data false.
- PACv1.GRCm39.chr7.+.130353825 merges five read-end peaks. Two of them sit
  at an internal exon's donor (130353640), but the PAC's representative peak
  lies 185 nt into the intron. The internal-exon-end filter tests only a PAC's
  representative coordinate, so the PAC passed.
- PACv1.GRCm39.chr7.+.130338300 sits on seven A's and is flagged for internal
  priming.

The lead's decisions:
- Discovery is unchanged: the exon-end filter keeps testing only the
  representative coordinate. Testing every merged peak would also reject real
  last-exon PACs whose region reaches a short last exon's upstream donor. The
  lead also declined a variant protected by annotated transcript ends.
- Internal priming stays a flag. Flagged PACs are tested and reported, but
  never called gained or lost and never support an APA pattern.

## 2026-09-30: Plain-language output names and a column guide

Status: accepted (project lead)

The lead asked for outputs that read without knowing the pipeline: one
shouldn't need its internals to know what a column or label means. The
version rises to 0.4.0, since scripts written against the old names break.
The atlas columns change, so its checksum, and with it the seeded statistics,
change too.

**Calls and gene events**
- `CONDITION_vs_CONTROL.events.tsv.gz` becomes `.calls.tsv.gz`, the `.pacs`
  rows with a call.
- `event_type` holds only PAC calls. The gene labels `dominant_switch`,
  `complexity_gain`, and `complexity_loss` no longer go on PAC rows. They
  marked PACs that had not changed, and undercounted genes whose PACs all had
  calls. The `.genes` tables already record every gene's events.
- `complexity_change` (`gain`, `loss`) becomes `active_pacs_change` (`more`,
  `fewer`). An active PAC has at least `active_pac_min_pau` (default 0.05) of
  its gene's fitted usage in a group.
  - The new parameter replaces `event_min_treatment_pau` for this count only.
    That parameter keeps the gained and lost calls and the UTR share, so the
    lead can raise the active-PAC threshold without changing calls.
  - The lead chose "active" over "present". "Detected" keeps its meaning,
    reads in at least `event_min_supporting_samples` samples of a group, for
    gained and lost calls. The per-PAC counts become `control_active_pacs`
    and `treatment_active_pacs`.

**Other renames**
- `assignment_class` becomes `gene_region`, with values `last_exon`,
  `internal_exon`, `intron`, `downstream_of_gene`, and `intergenic` (were
  `terminal_exon`, `other_exon`, `intronic`, `downstream`). In the figure
  tables, `distal_assignment_class` becomes `distal_gene_region`.
- `last_exon`, a locus, becomes `last_exon_locus`.
- `supporting_condition` becomes `best_supporting_condition`.
- Test columns name their level:
  - `.pacs` gets `pac_pvalue`, `gene_pvalue`, `pac_likelihood_ratio`, and
    `pac_degrees_of_freedom`.
  - `.genes` and `FAMILY.gene_omnibus` get `gene_pvalue`,
    `gene_likelihood_ratio`, and `gene_degrees_of_freedom`.
  - The motif tables' `p_value` becomes `pvalue`, like `fdr` beside it.
- `apa_pattern` `other` becomes `unclassified_change`.
- `model_status` `drimseq` becomes `fitted`, and `drimseq_add_uniform`
  becomes `fitted_with_zero_count_stabilization`.
- The `site_classes` figures become `calls_by_gene_region`, titled "PAC calls
  by gene region". The lead chose to rename them in this release rather than
  break scripts a second time later; the other figure names stay.
- Yes/no text stays as each writer prints it: `True`/`False` in the tables
  Python writes, `TRUE`/`FALSE` in the statistics tables, and `true`/`false`
  in the sample check. The guide says which is which, and every reader in the
  pipeline compares them case-insensitively. The lead chose to keep the mix.

**The guide and the report**
- `docs/output_columns.md` explains every published table and column: a
  glossary, then each column's meaning, units, and values. The README links it
  from a "Reading the results" section, and the report from a glossary line
  at its top.
- `tests/pipeline/column_guide.py` checks, in both fixture runs, that the
  guide documents exactly the columns each published table has.
- The report's headings say what each table holds, such as "TreatmentA vs
  DMSO: PACs with a call", with a one-line description under each.

## 2026-09-30: Exon-end peaks, pooled readthrough, and potential internal priming

Status: accepted (project lead)

Tacc2 in the lead's Plasmidsaurus data prompted this. Eight candidates at
internal exon 3' ends reached its atlas, and a real CR8 intronic gain at an
internal-priming-flagged site made the gene `other`.

**Internal exon ends** (`internal_exon_end_filter`, default true)
- Proximal-tag candidates where a pile of read ends at an annotated internal
  exon's donor would peak are rejected as `internal_exon_end`: within one
  25-nt bin of the donor's bin plus the kernel's peak offset.
- Donors within two bins of any transcript's 3' end are left out, since the
  filter reaches one and a half bins from a donor; a real end there is never
  rejected. Genes without transcript IDs have no donors.
- A real site whose reads pile up that close to an internal donor is rejected
  too, whatever its reads show.
- The filter uses only the annotation and the kernel, so the atlas is still
  built without contrasts.

**Readthrough**
- The CIGAR rule now pools each condition's replicates, and only conditions
  that support the candidate decide. A supporting condition without
  continuation reads keeps it.
- `constitutive_readthrough_min_replicate_support` is removed.
  `constitutive_readthrough_min_junction_count` (default 2) is now pooled
  over a supporting condition's samples.
- The block test uses the candidate's read pile, one bin either way, the same
  on both strands.
- Testing the candidate's whole resolution window was considered and not
  chosen, because it would reject real sites within about 150 nt of an exon.
- **Cost, accepted by the lead (option a):** a real site is rejected too when
  its reads pile up inside an exon that every supporting condition splices
  onward from, within about one read length of the donor. That covers
  minority-isoform internal-exon sites. README and design say so.
- Two alternatives were raised in review and not chosen: a threshold
  relative to the candidate's own reads, and exempting piles that outnumber
  the continuation reads.

**Potential internal priming**
- Every pattern gets a form with the suffix `_potential_internal_priming`,
  for a pattern that holds only once calls on flagged PACs count. The forms
  follow the confident patterns, and a gene has at most one form of each.
- `potential_internal_priming_withheld_calls` (default true) lets a flagged
  PAC's gain or loss count when only the flag withheld it. The lead asked for
  this to be a parameter.
- **Figures:** these patterns keep their pattern's color.
  - They are counted in a panel of their own, marked `*` on the grid, and
    drawn as diamonds on the distal plot.
  - Lighter tints were checked with the dataviz palette validator and
    failed: intronic gain's tint met intronic loss's tan at a normal-vision
    ΔE of 7.6, against a floor of 15.

**Left alone:** searching the whole resolution window for a PAS.

The atlas changes for proximal-tag runs, which reseeds their statistics.
Parameters are not checked against a list of names, so the removed
parameter, if still set, is ignored. The version rises to 0.3.0.

## 2026-09-29: Figures across comparisons and the PAU PCA figure

Status: accepted (project lead)

**Shared APA patterns**
- `figures/apa_pattern_grid` draws the genes with a pattern other than `none`
  in at least two comparisons, at most 50. They are ordered by the number of
  comparisons with a pattern, then by best gene FDR.
- `apa_patterns_by_comparison.tsv.gz` has every tested gene, with an empty
  cell where a comparison didn't test it.

**Concordance**
- `figures/concordance` has one panel per pair of comparisons, up to 15. Each
  panel plots the change in PAU of every PAC tested in both comparisons.
- Beyond 15 pairs, only pairs that share a control or chain through a
  condition are drawn. `concordance_matrix` shows Pearson r for every pair.
- `concordance.tsv.gz` records each pair's relation (`shared_control`,
  `chained`, or `unrelated`). A shared control's estimate correlates its
  comparisons positively. A chain estimates its middle condition from the
  same samples in both comparisons, with opposite signs, which correlates
  them negatively. r is therefore read with the relation.

**PAU PCA**
- `figures/pau_pca` draws the report's former PCA rule as a figure. The rule:
  - genes with at least `min_gene_total` reads in every sample;
  - their PACs observed in every sample, with no zero-filling;
  - each PAC centered across samples.
- The report shows the figure beside the PAU correlation table, and no longer
  tabulates the components. `pau_pca.tsv` has the coordinates.
- Each component's sign makes its largest loading positive. A two-PAC gene's
  loadings tie exactly, because its PAU sum to 1, so among loadings within a
  relative 1e-6 of the largest, the first PAC by ID decides.

These figures set each comparison's own results side by side. They add no
treatment-versus-treatment test, so every call remains a comparison with the
direct control.

## 2026-09-29: Gene-level APA patterns and refined terminal exons

Status: accepted (project lead)

**APA patterns**
- Each `.genes` table gains `apa_pattern`: `intronic_gain`, `intronic_loss`,
  `alternative_last_exon`, `utr_shortening`, `utr_lengthening`, `other`, or
  `none`. A gene can carry several, joined by `;`.
- Three numbers come with it: `delta_intronic_share`,
  `delta_utr_distal_share`, and `last_exon_switch`.
- Patterns come from each region's share of the fitted usage, gated by
  confirmed PAC calls on PACs that are not low confidence. Internal-priming
  sites therefore never make a pattern.
- A UTR pattern needs a call pointing the same way. When usage also moves
  into or out of the main last exon, only calls against that movement count.
- The alternative was a positional rule: every confirmed gain upstream of
  every confirmed loss. It was not chosen because it leaves one-sided genes
  unclassified.

**Within-last-exon metric**
- The UTR patterns use the distal PAC's share of the main last exon, the
  counterpart of DaPars' PDUI.
- LABRAT's rank-based ψ was considered. It moves by ΔPAU/(n − 1), so genes
  with more PACs would need larger shifts, and its value depends on which
  minor PACs the atlas found.

**Alternative last exons**
- The atlas gains `last_exon`: each terminal-exon PAC's last exon, made of a
  gene's overlapping terminal exons. A downstream PAC gets the gene's 3′-most
  last exon.
- The atlas bytes change, so every R seed changes. Stabilized p-values and
  bootstrap intervals move slightly on rerun; in the simulation, S6 coverage
  went from 83% to 85%.

**Terminal exons, everywhere**
- A transcript's final exon is a terminal exon unless:
  - it overlaps an internal exon of another transcript of the gene; or
  - it is a single-exon model apart from the gene's multi-exon transcripts.
- This removes the false last exons that retained-intron, 3′-incomplete, and
  intronic-fragment models create in Ensembl and GENCODE annotations.
- `assignment_class` changes with it, so the site-class figure, the distal
  PAC, and the patterns agree.

**Distal-usage wording**
- The distal-usage table and figure call the distal PAC's own call
  `distal_up` and `distal_down`. This supersedes "lengthened" and
  "shortened" from the 2026-09-28 figures entry.
- Shortening and lengthening now mean only the UTR patterns, and the figure
  colors each gene by its pattern.

**Rerunning**
- ANNOTATE_PACS stages `annotation.py` and `reference.py`, so `-resume`
  reruns annotation and every step after it when they change.
- A stale atlas stops the statistics at once, with a message that says what
  to do.
- The version rises to 0.2.0.

**Fixture**
- chr3 gains `ipa01`, which gains an intronic PAC in TreatmentA, and `ale01`,
  which switches between two last exons.

**Answers to the overseer's review.** All three rules stay as they are, and
their consequences are documented in the README, the design, and the methods.
1. **The single-exon rule.** A single-exon model lying 3′ of the gene's
   spliced transcripts is excluded, and its PACs count as upstream region.
   An intronic model that overlaps a retained-intron final exon is kept, and
   can form a false last exon.
2. **Spliced-through last exons.** A canonical last exon that a minor
   isoform splices through counts as internal. Its tandem UTR changes then
   read as intronic patterns. Telling it from a composite terminal exon would
   need transcript tags, which is out of scope for now.
3. **Same-direction changes under a region shift.** A within-exon change in
   the same direction as the shift stays unclassified, and appears only in
   `delta_utr_distal_share`.

## 2026-09-28: Treatment–control figures, drawn with ggplot2

Status: accepted (project lead)

**Figures**
- Each comparison gets three: a volcano plot, distal PAC usage, and PAC calls
  by site class.
- Two figures cover all comparisons together: event counts, and effect
  against coverage.
- The per-gene panels with fitted PAU, and a figure of motif-class shifts,
  were offered and not chosen.

**How they are made**
- `scripts/plot_usage_figures.R` draws them with ggplot2, as PDF and PNG, in
  PLOT_FIGURES after the statistics.
- They are published to `figures/`, and the report embeds the PNGs.
- The alternative was hand-drawn SVG in the report, which needs no new
  dependency. It was not chosen.

**Dependency**
- ggplot2 is declared directly in `envs/pacusage.yml` (`r-ggplot2>=4.0`), and
  `r-base` rises to `>=4.5`. Both already came in with DRIMSeq 1.38.0.
- R 4.5 is needed for `pdf(timestamp = FALSE, producer = FALSE)`, which
  leaves the dates and producer out of PDFs, so reruns are byte-identical.

**Distal usage**
- The distal PAC is a gene's most 3′ tested PAC in the terminal exon or
  downstream of the gene.
- A gene is `lengthened` or `shortened` only when that PAC has a confirmed
  call of its own: `gained` or `increased_usage`, or `lost` or
  `decreased_usage`.
- The alternative was "the gene passes `gene_fdr` and the distal PAU changes
  by at least `min_abs_delta_pau`". It is more sensitive: it would call
  fixture gene bg12, whose distal PAC rises by 0.17 with a PAC FDR of 0.10.
  But it leaves the direction untested, so it was not chosen.
- The figure only puts existing calls side by side, so its table,
  `CONDITION_vs_CONTROL.distal_usage.tsv.gz`, lives in `figures/`.

**Gene labels**
- The volcano and distal-usage plots label genes by `gene_name`, falling back
  to the ID, with at most 20 genes in each direction.
- The volcano labels a gene's most significant confirmed PAC. The distal
  plot labels the largest changes of the lengthened and shortened genes.

**Gene-level events**
- The `.genes` tables gain `dominant_switch` and `complexity_change`, and the
  event-count figure reads them.
- The `.events` labels go only on PACs without a call of their own, so they
  undercount these events. In the fixture, TreatmentA_vs_DMSO labels 1 of 7
  dominant switches, and Rescue_vs_TreatmentA 0 of 5.
- No statistic changes.

## 2026-09-28: Peer-review follow-up: event rules, motif classes, mitochondria, names, and columns

Status: accepted (project lead, answering the overseer's review of the methods
report and two later requests)

**Discovery**
- The support rule stays as it is. Candidates come from every sample's read
  ends pooled, and a PAC needs replicate support within the condition that
  supports it.
- The readthrough filter stays per condition.
- The docs now say "built without the treatment–control contrasts" instead of
  "condition-blind", and name both places condition labels are used.
- The end-to-end effect of this selection on false positives is unmeasured,
  and the methods report lists it as an open question.

**Event rules**
- `dominant_switch` and `complexity_gain`/`complexity_loss` are given only in
  genes that pass the comparison's gene-level screen.
- A `gained` call needs the control samples together to have at least
  `min_gene_total` reads at the gene, and a `lost` call the same of the
  treatment samples. Otherwise the call is `*_candidate`.

**Motif and k-mer analyses**
- The k-mer background is the other PACs of an event's genes that were
  tested in the same comparison.
- Motif preference is tested per primary hexamer, as before, and also per
  motif class, in `.preference_class` tables. This supersedes the "motif
  classes as rows" wording of the 2026-09-27 cleanup.

**Parameters**
- `dm_bootstrap_replicates` stays at 200 for the HPC run, to be revisited.
- `excluded_contigs` defaults to `[chrM, MT, chrMT]`, matched by alignment
  name or alias target. Strandedness inference honours it too.

**Output tables**
- Every gene-keyed output carries `gene_name` from one gene-level map: the
  gene record's name, else the first name on its other records, else the
  `gene_id`.
- Every PAC-level table starts with `pac_id, gene_id, gene_name, chrom, start,
  end, strand, locus`. `start`/`end` follow the atlas BED, and `locus` is
  1-based. The rest of each table runs from the answer to the evidence.
- Duplicate columns are dropped: `feature_id` is published as `pac_id` only,
  `site_class` duplicated `assignment_class`, and `family` in `.pacs`
  duplicated `control_condition`. From the atlas, `resolution_group` and
  `region_start`/`region_end` are dropped, and `contig` becomes `chrom`.
- Count rows follow the atlas's genomic order, and BED scores are capped at
  1000.

**Design document**
- Proximal peaks are ranked by matched score on the pooled binned kernel.
- Strandedness has no terminal-exon priority, and a conflict with the
  requested value stops the run.
- The exact-quantification tie rule is stated.

Reason: the overseer found these differences between the design and the code,
or a peer reviewer's objection, and the project lead chose each answer. The
names, columns, and mitochondrial default came from the lead's own review of
the HPC results. The atlas bytes change, so the seeded statistics shift once;
the next HPC run starts fresh.

## 2026-09-27: Stabilized genes' p-values are seed-dependent, and documented as such

Status: accepted (project lead, option (a) of the overseer's question)

Genes fitted with zero-count stabilization (`model_status`
`drimseq_add_uniform`) keep the stabilization and instability rules of
"Statistics repair details". The README now says that their precision and
p-values depend on the seed, and should be read alongside `model_status`.

Reason: on the fixture, a new seed moved one stabilized gene's precision 7x
and one PAC's `pac_fdr` from 0.18 to 0.001, while ΔPAU moved by under 0.001.
The rejected alternatives were a precision or p-value spread criterion in the
instability flag, more stabilization repeats, and a spread column. Any of
them would change the statistics just before the final HPC run.

## 2026-09-27: Resume reruns everything from the alignment scan after any parameter change

Status: accepted (project lead)

Every step from SCAN_ALIGNMENT on reads `resolved_params.yaml`. So any
parameter change reruns them on `-resume`, including `outdir`, the `save_*`
flags, and the Slurm and CPU settings; only preparation and strandedness
inference are reused. This stays as it is, and the README says so.

Reason: VALIDATE_INPUTS reruns on any change, and its outputs get new paths,
so downstream caches miss. Avoiding that would need content-hashed inputs,
which re-read every alignment on resume, or validation split across more
steps. That is complexity the final cleanup was meant to remove.

## 2026-09-27: Final cleanup before the HPC rerun

Status: accepted (project lead, by approving the cleanup plan and answering
its questions)

- **One implementation of each step.**
  - `scan.py` is the only evidence extraction. `evidence.py`'s
    `extract_evidence` and `extract_splice_continuations`, and the serial
    `pacusage calibrate`, are deleted. Tests assert the documented rules, plus
    calibration outputs recorded from the removed path while it existed.
  - Statistical filtering reasons are written per comparison family by
    `fit_usage_model.R`, from the filter that selects the tested PACs, to
    `statistics/FAMILY.statistical_filtering.tsv.gz`. `qc/statistical_filtering.tsv`
    and Python's `filter_testable_features` are removed.
  - `nextflow_schema.json` defines every parameter. VALIDATE_INPUTS stages it
    and validates once. Later steps read `resolved_params.yaml` as it is,
    except INFER_STRANDEDNESS, which takes its parameters on the command line
    so that other parameter changes don't rerun it. Python's `DEFAULTS` and
    its copies of the schema's rules are removed.
- **Logic fixes.**
  - Contig aliases apply to `excluded_contigs` (raw or aliased name) and to
    strandedness inference. This closes the deferral in "Read each alignment
    once".
  - The k-mer CMH interval uses the Robins-Breslow-Greenland variance, as R's
    `mantelhaen.test` does.
  - Motif preference follows design section 9:
    - it fits the family's main design, covariates included;
    - one limma fit per comparison, with motif classes as rows;
    - a class is tested only when every sample of the comparison has enough
      informative genes.

    Motif scores and the report's PAU QC use genes with at least
    `min_gene_total` reads in every sample.
  - A compressed FASTA is rejected, and the `.fai` is always generated in the
    task directory. This supersedes the sibling-index reuse in
    "Implementation choices in the BAM pass cleanup".
  - `samtools sort` runs with `--no-PG`.
  - File parameters are resolved against the launch directory.
  - The `save_*` publication flags are compared as text, since Nextflow 25.10
    and later pass command-line values as strings.
  - RECORD_SOFTWARE_VERSIONS adds R and the statistics packages to
    `software_versions.tsv`. Loading them right after validation stops a run
    that lacks one within minutes. It is a separate step, so that an edit to
    the R script reruns only the statistics on `-resume`.
  - The report stages the sample sheet, so its gene plots appear.
  - Motif preference tables are published only to `motifs/`.
  - Gene and known-PAC lookups in annotation and clustering use indexes
    instead of scanning every feature.
- **Removed.**
  - The scan and bootstrap-batch version guards.
  - The VALIDATE_INPUTS cache-busting comment, which supersedes that item of
    "Implementation choices in the BAM pass cleanup".
  - Fallbacks for older output formats.
  - Three never-set atlas columns and the never-produced
    `motif_assisted_rescue` status.
  - The Dockerfile, `envs/test.yml`, the bootstrap channel test, and the
    benchmark script.
- **Layout.**
  - The design and this log moved to `docs/design.md` and
    `docs/decisions.md`; this supersedes the file name in "Decisions are
    recorded in DECISIONS.md".
  - The R tests are in `tests/r/`, and the end-to-end script and its checks in
    `tests/pipeline/`.
- **Deployment.**
  - The Apptainer image stays self-contained, with the Python package
    installed, and is rebuilt after every pull. `bin/pacusage` prefers the
    installed package and otherwise runs `src/`.
  - With `scratch` set, tasks run in the work directory.
  - CI pins Nextflow 25.04.7, the cluster's version.

Reason: the lead asked for one last simplification pass, with no backwards
compatibility, before a fresh run on the cluster. Each item removes duplicate
or dead logic, or fixes a mismatch with the design. The atlas file lost three
columns, and its checksum seeds the statistics, so the seeded results change
once. On the fixture this moved ΔPAU by under 0.001. In the one gene with
zero-count stabilization, it moved p-values by up to two orders of magnitude
and changed one descriptive event call.

## 2026-09-27: The calibration kernel diagnostic warns and never stops a run

Status: accepted (project lead, by approving the calibration kernel warning
plan)

This settles the details that the next entry left open.
- **Warn only.** A proximal-tag run whose pooled kernel fails a check goes on.
  - The warning appears in the Nextflow log, in
    `qc/calibration_kernel_diagnostics.tsv`, and at the top of the report.
  - The empty-atlas stop at ANNOTATE_PACS still catches the worst case.
- **Two checks,** on the selected source's pooled kernel exactly as discovery
  reads it back from `calibration_kernel.tsv`. Summing the untrimmed array can
  round an overlap tie differently and shift the resolution by 1 nt.
  - **Separated modes.** The kernel is smoothed with a centred moving average
    as wide as its minimum resolvable separation, rounded down to an odd
    width. A peak counts as a mode when it reaches 25% of the highest.
    Neighbours merge when the valley between them stays at or above 50% of the
    lower one. More than one mode warns.
  - **Spread.** A median per-sample central interval width more than 6 times
    the minimum resolvable separation warns.
- **Constants, not parameters.** The thresholds are
  `KERNEL_MODE_MINIMUM_HEIGHT`, `KERNEL_MODE_VALLEY`, and
  `KERNEL_MAXIMUM_SPREAD_RATIO` in `calibration.py`. They can become
  parameters if real libraries call for tuning.
- **Exact-boundary runs** record the same values with status
  `not_applicable`, because exact discovery does not assign reads through the
  kernel.

Reason: the kernel of PAC spacings has 3 modes and a spread ratio of 50.
Single-site kernels measured 1 mode and ratios of 1.8 to 4.2, and two sites
250 nt apart gave 2 modes. Counting raw local maxima instead would split the
healthy fixture's flat top into twin peaks. Warning rather than failing keeps
a run going when an unusual but genuine kernel trips a check.

## 2026-09-27: Calibration kernels get a diagnostic, not new gene selection

Status: accepted (project lead, option b of the proximal-fixture review)

- **Gene selection stays as it is.** Calibration genes remain
  annotation-defined: genes whose transcripts share one annotated end.
- **A diagnostic is added.** It will warn, or fail, when a proximal-tag
  profile's pooled kernel is multi-modal, or when its minimum resolvable
  separation is far below the calibrated central interval width.
- **Details are still open.** Its thresholds, and whether it warns or fails,
  are planned separately.

Reason: calibrating the exact-boundary fixture as Plasmidsaurus reported a
proximal model at 0.999 reproducibility from a kernel of PAC spacings, with a
2-nt resolution, and nothing warned. Unannotated alternative polyadenylation
near annotated ends could broaden real kernels in the same way.

## 2026-09-27: A Plasmidsaurus-like fixture with its own reference

Status: accepted (project lead, by asking for it to be committed)

- **Why the old sheet failed.** `samples_plasmidsaurus.tsv` reused the
  exact-boundary alignments, whose reads end on PACs spaced 50 nt apart.
  - Calibration measures reads against annotated distal ends, so it learned
    that spacing: a kernel with spikes at 0, 50, and 100 nt.
  - Every endpoint then matched several peaks about equally, so proximal
    discovery accepted no PAC.
  - Discovery itself recovers PACs from genuine proximal-tag reads. The sheet
    is removed.
- **The replacement.** `tests/fixtures/plasmidsaurus/` has its own FASTA and
  annotation, built by `build_plasmidsaurus_fixture.py` through
  `build_fixture.py`.
  - 24 single-PAC genes calibrate the kernel.
  - 18 multi-PAC genes annotate one transcript per PAC, so calibration skips
    them.
  - Reads end 20-280 nt upstream of their PAC, at deterministic triangular
    quantiles.
  - PACs sit 400 nt apart on the 25-nt discovery grid, so every read belongs
    to exactly one PAC.
  - Three genes shift usage in TreatmentA, and two are exact 3:2 nulls.
- **Where it runs.** `tests/run_nextflow.sh` runs the fixture and checks it
  with `tests/verify_plasmidsaurus.py`.
- **Empty atlas.** An empty accepted atlas stops the run at ANNOTATE_PACS,
  with an error naming `qc/pac_discovery.tsv` and
  `atlas/rejected_candidates.tsv.gz`. It no longer fails later in
  MOTIF_SCORES.

Reason: the plan's integration fixture requires a successful
Plasmidsaurus-like SE run, which exact-boundary reads cannot provide.

## 2026-09-27: Implementation choices in the BAM pass cleanup

Status: accepted (project lead)

- **Milestone.** The project lead accepted the milestone with the proximal
  path covered by the CLI-level equivalence tests (option a). The
  Plasmidsaurus-like fixture, recorded separately, has since covered that path
  end to end.
- **Source fingerprint.** It is size plus modification time, not the planned
  size, modification time, and inode (option a). Nextflow's cache ignores
  inodes, so an inode-only change (for example a file restored by `rsync -a`)
  would fail every later step, and `-resume` could never clear it.
- **Merge order.** MERGE_USAGE_MODELS orders `gene_precision.tsv.gz` and
  `fitted_pau.tsv.gz` shards by family. Their directories are numbered in
  task-creation order, which changed between two otherwise identical runs, so
  the published row order was not reproducible.
- **VALIDATE cache.** A comment inside VALIDATE_INPUTS's script changes its
  task hash. Otherwise `-resume` would reuse an older version's
  `input_checksums.tsv`, which still has alignment rows.
- **Reference inputs.** As planned, PREPARE_REFERENCE takes the FASTA as a
  declared input instead of reading its path from `resolved_params.yaml`. So
  changing an unrelated parameter no longer reruns reference or alignment
  preparation on `-resume`.
  - The deviation: the source's sibling `.fai` is found through the FASTA's
    real path, not declared as an input.
  - This is safe because a reused `.fai` is validated against the FASTA, and
    an index is a deterministic function of an unchanged FASTA, so a cached
    task can't hold a wrong one.
- **Sorted-copy detection.** `check_prepared_source` recognizes a sorted copy
  by its recorded `action` (`sorted_and_indexed`), not by comparing paths. So
  a linked source whose path resolves differently inside a container is still
  checked.

Reason: each keeps the approved plan's targets, reproducible reruns and
reliable `-resume`, where a literal reading of the plan would not.

## 2026-09-27: Read each alignment once for calibration and evidence

Status: accepted (project lead, by approving the BAM pass cleanup plan)

- Python gzip writers set a zero header timestamp. Reruns then reproduce the
  atlas checksum and the statistics seeds within one software environment.
  Statistics change once when this lands.
- Strandedness inference keeps its full-pass reservoir sample.
- `manifest/input_checksums.tsv` stays complete.
  - VALIDATE hashes the sample sheet, FASTA, and GTF, but no longer the
    alignments.
  - PREPARE_ALIGNMENT hashes each alignment once, and a new
    RECORD_INPUT_CHECKSUMS step adds those rows.
- Prepared alignments, their indexes, and the prepared FASTA are symlinks to
  the sources. Source files must stay in place, unmodified, and visible inside
  containers until an analysis is accepted.
- SCAN_ALIGNMENT (`pacusage scan-alignment`, label `medium`) replaces
  CALIBRATE_SAMPLE. One pass over each alignment collects every candidate
  evidence source, splice continuations, and filtering counters.
- EXTRACT_3PRIME_EVIDENCE writes its outputs from the scan without reading the
  alignment, with label `serial_medium`.
- Two existing behaviors stay for exact equivalence, to be fixed separately:
  - excluded contigs match the alignment's raw contig names;
  - strandedness inference ignores contig aliases.

Reason: each alignment was hashed 4 times, copied up to 3 times, and read in
full 4-7 times, on 8-CPU reservations for single-threaded work. Evidence,
atlas, counts, and statistics stay identical apart from the one-time seed
change.

## 2026-09-27: Accept the measured bootstrap coverage and pin DRIMSeq

Status: accepted (project lead)

Bootstrap intervals conditional on the family-fit precision cover the true
change in PAU in 80-93% of simulated cases, about 88% on average at a nominal
95% with four replicates per group. That is the documented behavior of
version 1, as the README and plan section 9 describe. The simulation check S6
pools all 60 PACs of its coverage family and requires coverage of at least
0.75. This was option (b) of the design review.

The Conda environments (`envs/pacusage.yml` and `envs/test.yml`, and through
them the Apptainer and Docker images) pin `bioconductor-drimseq=1.38.0`, the
release the statistics were validated with and the newest on Bioconda. stageR
is not pinned; `stage_adjust` stops on an unexpected stageR result.

Reason: the coverage shortfall is not specific to stabilized genes and has no
cheap, correct fix. A pinned DRIMSeq stops a release with different numerics
from silently changing results.

## 2026-09-27: Statistics repair details

Status: accepted (project lead, by approving the statistics repair plan)

- Zero-count stabilization perturbs the counts themselves, following DRIMSeq's
  `addUniform` rule (each zero becomes `runif(0, 0.1)`). Precision, the full
  model, and the null model are all fitted to the same perturbed data. It runs
  once per comparison family, with seeds drawn from the family, gene, and
  repeat.
- Stabilized genes need at least 80% successful repeats. Their written-back
  values are the medians across repeats; each p-value comes from the median
  likelihood ratio.
- A direction change counts as instability only beyond ±0.001 change in PAU.
- A condition with no counts for a gene gets the status
  `group_without_counts`.
- Bootstrap intervals are conditional on the family-fit precision. One family
  draw serves every comparison, and `dm_bootstrap_replicates` counts attempts.
- Output tables gain `stabilization_delta_pau_spread`, `zero_boundary_reason`,
  `bootstrap_perturbed`, `bootstrap_status`, and gene-level `model_status`.
  DRIMSeq's `adj_pvalue` is dropped.
- Walltime scales with `task.attempt` for every label: statistics 12 h,
  bootstrap 6 h, low 2 h, medium 12 h. The default bootstrap batch size is
  500.
- The statistics script reads IDs as text and treats covariates as
  categorical. It orders PACs by coordinate.
- Simulation-test thresholds may sit at about 2-3× nominal error rates,
  because DRIMSeq is liberal for overdispersed genes at 3-4 replicates.
- The integration fixture BAMs grow to roughly 0.3-0.4 MB.
- CLAUDE.md, AGENTS.md, the README, and plan section 9 are updated to match.

Reason: these are the choices in the approved repair plan
(`~/.claude/plans/snazzy-dancing-hopper.md`). They follow from the path A
decision below.

## 2026-09-27: Repair the statistics in place and keep DRIMSeq

Status: accepted (project lead)

Path A from the runtime audit. Keep the current pipeline and the DRIMSeq
design in plan section 9. Fix the stageR lookup, the zero-count fallback, and
the incomplete write-back of stabilized results. Stabilize boundary genes once
per comparison family with seeded whole-family `add_uniform` fits. Base
bootstrap intervals on the family fit's precision. Scale walltime with each
retry. Add a realistic fixture and the plan's statistical simulation test.
Work on the upstream BAM passes and resource labels comes after that. The
rejected alternatives were a new statistics engine (satuRn or edgeR), an
existing tool (MAAPER, or PolyASite with DEXSeq), and a full rewrite.

Reason: the audit traced both the multi-day runtime and the missing
significant events to `scripts/fit_usage_model.R` and the retry settings. The
upstream stages measured at minutes per sample.

## 2026-09-27: A treatment can be another treatment's control

Status: accepted (project lead)

A condition can be tested against its own control and also be the declared
control of another condition, as in `WT -> disease_vehicle -> disease_drug`.
Control relationships must be acyclic. This is how `samples.py` and the README
already behaved since `ddf8fe4`; validation rules 5 to 7 in
`pacusage_hpc_implementation_plan.md` were updated to match.

## 2026-09-27: Agent worktrees branch from local HEAD and are gitignored

Status: accepted (project lead)

`.claude/settings.json` sets `"worktree": {"baseRef": "head"}`, and
`.gitignore` lists `.claude/worktrees/`. Worktrees merge back into the local
`main`, so they start from it, including commits that haven't been pushed.
Ignoring the directory keeps agent worktrees out of `git status` and the code
index.

## 2026-09-27: Decisions are recorded in DECISIONS.md

Status: accepted (project lead)

This file is the decisions log that `CLAUDE.md` and the agents refer to.
`fable-overseer` has no Write or Edit tools, so it adds proposed entries
through Bash.

## 2026-09-27: haiku-implementer has no effort setting

Status: accepted (project lead)

Haiku 4.5 has no effort levels, so `effort: max` would have had no effect.
The line was removed from `.claude/agents/haiku-implementer.md`, and the
implementer runs at the model's default.

## 2026-09-27: The main session's stored effort is xhigh

Status: accepted (project lead)

Settings files accept effort only up to `xhigh`, so `.claude/settings.json`
stores `"effortLevel": "xhigh"`. Choose `max` per session, in the desktop
app's model menu or with `claude --effort max`.

## 2026-09-27: Three-model Claude Code workflow

Status: accepted (project lead)

Opus 5.5 orchestrates, Haiku 4.5 implements file-scoped tasks, and Fable 5.1
reviews design questions and signs off milestones. `CLAUDE.md` describes the
roles and the review checklist, and `.claude/agents/` defines the two
subagents. `CLAUDE.md` imports `AGENTS.md`, because Claude Code stops reading
AGENTS.md on its own once a CLAUDE.md exists.
