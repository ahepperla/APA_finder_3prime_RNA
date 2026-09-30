#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${project_dir}"
# Stale outputs from an earlier run must not satisfy the checks below.
rm -rf results-test

scratch_root="$(mktemp -d)"
trap 'rm -rf "${scratch_root}"' EXIT

Rscript tests/r/test_usage_model.R scripts/fit_usage_model.R
Rscript tests/r/test_usage_model_simulation.R scripts/fit_usage_model.R
Rscript tests/r/test_usage_figures.R scripts/plot_usage_figures.R scripts/fit_usage_model.R
python tests/fixtures/build_fixture.py

# Prepared alignments and the FASTA link to the fixture files, so every run
# must leave them byte-identical and must write nothing beside them, such as
# a FASTA index next to genome.fa.
fixture_state() {
  (cd tests/fixtures && find . -type f -not -path '*/__pycache__/*' | LC_ALL=C sort \
    | xargs shasum -a 256)
}
fixtures_before="$(fixture_state)"

nextflow run . -profile test,local -resume

test -s results-test/report/index.html
test -s results-test/atlas/pacs.v1.metadata.tsv.gz
test -s results-test/counts/pac_counts.tsv.gz
test -s results-test/statistics/TreatmentA_vs_DMSO.pacs.tsv.gz
test -s results-test/statistics/TreatmentA_vs_DMSO.calls.tsv.gz
test -s results-test/statistics/TreatmentB_vs_Vehicle.pacs.tsv.gz
test -s results-test/statistics/Rescue_vs_TreatmentA.pacs.tsv.gz
test -s results-test/statistics/gene_precision.tsv.gz
test -s results-test/statistics/fitted_pau.tsv.gz
test -s results-test/figures/TreatmentA_vs_DMSO.volcano.pdf
test -s results-test/figures/event_counts.png
test -s results-test/pipeline_info/execution_report.html
test -s results-test/pipeline_info/execution_timeline.html
test -s results-test/pipeline_info/pipeline_dag.html
python tests/pipeline/verify_integration.py
test "$(fixture_state)" = "${fixtures_before}"

# The runs below launch from scratch directories, so they never become the
# session that a later -resume from the repository continues.

# A fresh run in a new work directory must reproduce every published file.
rerun_root="${scratch_root}/rerun"
mkdir -p "${rerun_root}"
(cd "${rerun_root}" && nextflow run "${project_dir}" -profile test,local \
  --outdir "${rerun_root}/results" \
  -work-dir "${rerun_root}/work")
python tests/pipeline/compare_results.py results-test "${rerun_root}/results"
# Nextflow's reports go to the results, never to the launch directory.
test -s "${rerun_root}/results/pipeline_info/execution_trace.txt"
test -z "$(find "${rerun_root}" -maxdepth 1 -type f \( -name '*.html' -o -name '*.txt' \))"
test "$(fixture_state)" = "${fixtures_before}"

# Samples from two protocols must fail at calibration, before discovery.
incompatible_root="${scratch_root}/incompatible"
mkdir -p "${incompatible_root}"
if (cd "${incompatible_root}" && nextflow run "${project_dir}" -profile test,local \
  --input "${project_dir}/tests/fixtures/samples_incompatible.tsv" \
  --outdir "${incompatible_root}/results" \
  -work-dir "${incompatible_root}/work") > "${incompatible_root}/log.txt" 2>&1; then
  echo "The incompatible-profile sample sheet did not fail." >&2
  exit 1
fi
grep -q "Samples use incompatible library profiles" "${incompatible_root}/log.txt"
# The trace must show the step that stopped the run, or the check after it
# would pass on a missing trace.
incompatible_trace="${incompatible_root}/results/pipeline_info/execution_trace.txt"
grep -q "AGGREGATE_CALIBRATION" "${incompatible_trace}"
if grep -q "CLUSTER_PACS" "${incompatible_trace}"; then
  echo "The incompatible-profile run reached PAC discovery." >&2
  exit 1
fi
test "$(fixture_state)" = "${fixtures_before}"

# The Plasmidsaurus-like fixture runs the proximal-tag path end to end, with
# its own reference and annotation. It also publishes the per-sample count
# tables, so the column guide's section for them is checked.
plasmidsaurus_fixture="${project_dir}/tests/fixtures/plasmidsaurus"
plasmidsaurus_root="${scratch_root}/plasmidsaurus"
mkdir -p "${plasmidsaurus_root}"
(cd "${plasmidsaurus_root}" && nextflow run "${project_dir}" -profile test,local \
  --input "${plasmidsaurus_fixture}/samples.tsv" \
  --fasta "${plasmidsaurus_fixture}/genome.fa" \
  --gtf "${plasmidsaurus_fixture}/genes.gtf" \
  --endpoint_model auto \
  --save_intermediates true \
  --outdir "${plasmidsaurus_root}/results" \
  -work-dir "${plasmidsaurus_root}/work")
python tests/pipeline/verify_plasmidsaurus.py "${plasmidsaurus_root}/results"
test "$(fixture_state)" = "${fixtures_before}"

# Exact-boundary reads run as Plasmidsaurus tags calibrate to a kernel of PAC
# spacings. Calibration warns and the run goes on, until discovery finds no PAC.
mislabelled_root="${scratch_root}/mislabelled"
mkdir -p "${mislabelled_root}"
awk -F '\t' -v OFS='\t' -v fixtures="${project_dir}/tests/fixtures" '
  NR == 1 { for (i = 1; i <= NF; i++) if ($i == "alignment") column = i }
  NR > 1 { $column = fixtures "/" $column }
  { print }
' "${plasmidsaurus_fixture}/samples.tsv" > "${mislabelled_root}/samples.tsv"
if (cd "${mislabelled_root}" && nextflow run "${project_dir}" -profile test,local \
  --input "${mislabelled_root}/samples.tsv" \
  --endpoint_model auto \
  --outdir "${mislabelled_root}/results" \
  -work-dir "${mislabelled_root}/work") > "${mislabelled_root}/log.txt" 2>&1; then
  echo "The mislabelled sample sheet did not fail." >&2
  exit 1
fi
grep -q "Calibration kernel: the pooled kernel has 3 separated modes" "${mislabelled_root}/log.txt"
grep -q "No PAC candidates passed discovery" "${mislabelled_root}/log.txt"
grep -q $'\twarning\t' "${mislabelled_root}/results/qc/calibration_kernel_diagnostics.tsv"
mislabelled_trace="${mislabelled_root}/results/pipeline_info/execution_trace.txt"
grep -q "CLUSTER_PACS" "${mislabelled_trace}"
if grep -q "QUANTIFY_PACS" "${mislabelled_trace}"; then
  echo "The mislabelled run reached quantification." >&2
  exit 1
fi
test "$(fixture_state)" = "${fixtures_before}"
