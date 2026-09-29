# Decisions

Decisions about how PACusage is built and maintained, newest first. The
project lead accepts decisions. `fable-overseer` may add entries with status
`proposed`, which stay proposed until the project lead accepts or rejects them.
To change an accepted decision, add a new entry that supersedes it, and only
with the project lead's approval.

Entry format: a dated heading, a status line, the decision, and the reason.

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
