@AGENTS.md

Keep the import above. Claude Code reads AGENTS.md on its own only when a
project has no CLAUDE.md.

## Model roles

| Role | Model | Effort |
|---|---|---|
| Orchestrator (main session) | Opus 5.5, `claude-opus-5-5` | max, chosen per session |
| `haiku-implementer` subagent | Haiku 4.5, `claude-haiku-4-5-20251001` | model default (Haiku has no effort levels) |
| `fable-overseer` subagent | Fable 5.1, `claude-fable-5-1` | xhigh |

The subagents take their model and effort from `.claude/agents/`.
`.claude/settings.json` sets the main session's model and stores `xhigh`, the
highest effort a settings file accepts. Choose `max` in the desktop app's model
menu, or with `claude --effort max` in a terminal. The desktop app passes its
own model and effort flags, and they override `.claude/settings.json`.

Point `fable-overseer` to the project documents: goals and design in
`pacusage_hpc_implementation_plan.md`, conventions and invariants in AGENTS.md,
and decisions in `DECISIONS.md`.

### Orchestrating

- Split each request into tasks scoped to one module or a few named files, with
  the behavior, the tests to add, and the commands to run spelled out. Give
  each task to its own `haiku-implementer`.
- Run independent tasks in parallel with `isolation: "worktree"`, and merge
  worktrees one at a time.
- Keep architecture, ambiguous work, prompt engineering, and anything
  security-, privacy-, concurrency-, or invariant-sensitive for yourself. In
  this repository that includes the AGENTS.md invariants (condition-blind
  atlas, exact versus proximal semantics, interbase coordinates, raw counts and
  PAU without pseudocounts, direct-control comparisons) and the statistics in
  `scripts/fit_usage_model.R`.
- Never accept "done" on an implementer's word. Read the diff, run the full
  lint and test suite yourself, and fix or redo what's wrong.
- When an implementer reports `UNCERTAIN:`, resolve it yourself if the docs
  settle it; otherwise ask `fable-overseer` or the project lead.
- Ask `fable-overseer` about design or scope calls the docs don't settle, and
  always before treating a milestone or large feature as done. Pass its
  `FOR THE PROJECT LEAD:` section to the project lead unchanged, and don't act
  on it until they answer.
- Use plan mode for anything touching more than one module, and share the plan
  before editing.
- When requirements are ambiguous, ask the project lead instead of guessing.
  Offer options with your recommendation first, and record the answer in
  `DECISIONS.md`.
- Before calling a task done, run the full lint and test suite and report the
  results as they are, failures included.
- Commits, pushes, and new dependencies wait for the project lead's request.

### Reviewing delegated work

- Weak tests: each new test must fail when the behavior breaks. Look for
  assertions like `> 0`, `any(...)`, or `a in x or b in x`, slicing of
  unordered results, and tests loosened until they pass. For critical tests,
  break the code, confirm the test fails, and restore it.
- Hidden changes: diff against the state before the task, including files that
  were deleted or reset to HEAD.
- "Pre-existing" failures: check the claim by running the tests against HEAD's
  version of the file.
- Lint claims: rerun the exact project command, `ruff check src tests`.
- Out-of-scope edits: look for edits outside the task's files, stray files
  (`.bak`, extra lockfiles), uncommitted work, and changed function signatures
  the task didn't ask for.
- Versions and formats: check any pinned library version against the registry,
  and check parsers against the documented format (BAM/CRAM, GTF/GFF3, BED, the
  sample-sheet TSV), not the implementer's own sample.
- Test summary: ask for the exact final test summary line, and rerun the suite
  yourself after every merge.

### Checks inside a worktree

Worktrees branch from local HEAD, so they don't include uncommitted changes in
the main checkout. A worktree also has no `.venv` or `.Rlib`, which are
gitignored, and the shared `.venv` installs `pacusage` in editable mode from
the main checkout's `src/`. Without `PYTHONPATH`, tests run in a worktree
silently test the main checkout. Run checks from the worktree root:

```bash
MAIN="$(cd "$(git rev-parse --git-common-dir)/.." && pwd)"
PYTHONPATH="$PWD/src" "$MAIN/.venv/bin/python" -m pytest
"$MAIN/.venv/bin/ruff" check src tests
R_LIBS="$MAIN/.Rlib" Rscript tests/test_usage_model.R scripts/fit_usage_model.R
```

For `tests/run_nextflow.sh`, also put `$MAIN/.venv/bin` first on `PATH`, keep
the same `PYTHONPATH`, and set `R_LIBS="$MAIN/.Rlib"` so the pipeline's R steps
find DRIMSeq and stageR.
