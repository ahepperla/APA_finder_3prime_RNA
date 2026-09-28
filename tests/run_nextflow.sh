#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${project_dir}"
# Stale outputs from an earlier run must not satisfy the checks below.
rm -rf results-test

scratch_root="$(mktemp -d)"
trap 'rm -rf "${scratch_root}"' EXIT

channel_test_root="${scratch_root}/channel"
mkdir -p "${channel_test_root}"
touch "${channel_test_root}/batch-001.rds"
touch "${channel_test_root}/batch-002.rds"
touch "${channel_test_root}/batch-003.rds"
nextflow run tests/test_bootstrap_batch_channel.nf \
  --batch_root "${channel_test_root}" \
  -work-dir "${channel_test_root}/work" \
  -with-trace "${channel_test_root}/trace.txt"
test "$(grep -c "ASSERT_SINGLE_BOOTSTRAP_BATCH" "${channel_test_root}/trace.txt")" -eq 3

Rscript tests/test_usage_model.R scripts/fit_usage_model.R
Rscript tests/test_usage_model_simulation.R scripts/fit_usage_model.R
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
test -s results-test/statistics/TreatmentA_vs_DMSO.events.tsv.gz
test -s results-test/statistics/TreatmentB_vs_Vehicle.pacs.tsv.gz
test -s results-test/statistics/Rescue_vs_TreatmentA.pacs.tsv.gz
test -s results-test/statistics/gene_precision.tsv.gz
test -s results-test/statistics/fitted_pau.tsv.gz
python tests/verify_integration.py
test "$(fixture_state)" = "${fixtures_before}"

# The runs below launch from scratch directories, so they never become the
# session that a later -resume from the repository continues.

# A fresh run in a new work directory must reproduce every published file.
rerun_root="${scratch_root}/rerun"
mkdir -p "${rerun_root}"
(cd "${rerun_root}" && nextflow run "${project_dir}" -profile test,local \
  --outdir "${rerun_root}/results" \
  -work-dir "${rerun_root}/work" \
  -with-trace "${rerun_root}/trace.txt")
python tests/compare_results.py results-test "${rerun_root}/results"
test "$(fixture_state)" = "${fixtures_before}"

# Samples from two protocols must fail at calibration, before discovery.
incompatible_root="${scratch_root}/incompatible"
mkdir -p "${incompatible_root}"
if (cd "${incompatible_root}" && nextflow run "${project_dir}" -profile test,local \
  --input "${project_dir}/tests/fixtures/samples_incompatible.tsv" \
  --outdir "${incompatible_root}/results" \
  -work-dir "${incompatible_root}/work" \
  -with-trace "${incompatible_root}/trace.txt") > "${incompatible_root}/log.txt" 2>&1; then
  echo "The incompatible-profile sample sheet did not fail." >&2
  exit 1
fi
grep -q "Samples use incompatible library profiles" "${incompatible_root}/log.txt"
if grep -q "CLUSTER_PACS" "${incompatible_root}/trace.txt"; then
  echo "The incompatible-profile run reached PAC discovery." >&2
  exit 1
fi
test "$(fixture_state)" = "${fixtures_before}"

# The Plasmidsaurus-like fixture runs the proximal-tag path end to end, with
# its own reference and annotation.
plasmidsaurus_fixture="${project_dir}/tests/fixtures/plasmidsaurus"
plasmidsaurus_root="${scratch_root}/plasmidsaurus"
mkdir -p "${plasmidsaurus_root}"
(cd "${plasmidsaurus_root}" && nextflow run "${project_dir}" -profile test,local \
  --input "${plasmidsaurus_fixture}/samples.tsv" \
  --fasta "${plasmidsaurus_fixture}/genome.fa" \
  --gtf "${plasmidsaurus_fixture}/genes.gtf" \
  --endpoint_model auto \
  --outdir "${plasmidsaurus_root}/results" \
  -work-dir "${plasmidsaurus_root}/work" \
  -with-trace "${plasmidsaurus_root}/trace.txt")
python tests/verify_plasmidsaurus.py "${plasmidsaurus_root}/results"
test "$(fixture_state)" = "${fixtures_before}"
