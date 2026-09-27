# Decisions

Decisions about how PACusage is built and maintained, newest first. The
project lead accepts decisions. `fable-overseer` may add entries with status
`proposed`, which stay proposed until the project lead accepts or rejects them.
To change an accepted decision, add a new entry that supersedes it, and only
with the project lead's approval.

Entry format: a dated heading, a status line, the decision, and the reason.

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
