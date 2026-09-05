#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${project_dir}"

python tests/fixtures/build_fixture.py
nextflow run . -profile test,local -resume

test -s results-test/report/index.html
test -s results-test/atlas/pacs.v1.metadata.tsv.gz
test -s results-test/counts/pac_counts.tsv.gz
test -s results-test/statistics/TreatmentA_vs_DMSO.pacs.tsv.gz
test -s results-test/statistics/TreatmentB_vs_Vehicle.pacs.tsv.gz
test -s results-test/statistics/Rescue_vs_TreatmentA.pacs.tsv.gz
test -s results-test/statistics/gene_precision.tsv.gz
test -s results-test/statistics/fitted_pau.tsv.gz
python tests/verify_integration.py
