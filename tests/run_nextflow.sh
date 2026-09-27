#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${project_dir}"
# Stale outputs from an earlier run must not satisfy the checks below.
rm -rf results-test

channel_test_root="$(mktemp -d)"
trap 'rm -rf "${channel_test_root}"' EXIT
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
