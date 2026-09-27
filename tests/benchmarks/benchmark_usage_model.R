# Time the statistics modes of a fit_usage_model.R version on a synthetic
# comparison family. Not part of CI; used to compare versions.
#
# Usage:
#   Rscript tests/benchmarks/benchmark_usage_model.R --script scripts/fit_usage_model.R \
#     --label after --replicates-per-group 4 --treatments 2 --genes 300 \
#     [--bootstrap-replicates 50] [--workers 4] [--max-bootstrap-genes 40] [--output results.tsv]
#
# --max-bootstrap-genes truncates each batch so that slow versions finish;
# the bootstrap cost is reported per bootstrapped gene.

options(warn = 1)
arguments <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  index <- match(paste0("--", name), arguments)
  if (is.na(index)) default else arguments[[index + 1L]]
}
script <- normalizePath(option("script", "scripts/fit_usage_model.R"), mustWork = TRUE)
label <- option("label", "run")
per_group <- as.integer(option("replicates-per-group", "4"))
treatments <- as.integer(option("treatments", "2"))
gene_count <- as.integer(option("genes", "300"))
replicates <- as.integer(option("bootstrap-replicates", "50"))
workers <- as.integer(option("workers", "4"))
max_bootstrap_genes <- as.integer(option("max-bootstrap-genes", "40"))
output <- option("output", "")
seed <- as.integer(option("seed", "20260927"))

benchmark_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
source(file.path(dirname(normalizePath(benchmark_file)), "..", "r", "usage_model_helpers.R"))
model <- load_usage_model(script)

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(seed)
conditions <- c("C", sprintf("T%d", seq_len(treatments)))
sample_condition <- rep(conditions, each = per_group)
sample_ids <- paste0(sample_condition, "_", rep(seq_len(per_group), length(conditions)))
genes <- lapply(seq_len(gene_count), function(index) {
  base <- list(c(0.5, 0.3, 0.2), c(0.4, 0.35, 0.25), c(0.6, 0.25, 0.15))[[index %% 3L + 1L]]
  proportions <- stats::setNames(rep(list(base), length(conditions)), conditions)
  # One gene in ten shifts usage in T1; one in twenty gains a PAC in T1.
  if (index %% 10L == 0L) proportions$T1 <- c(0.3, 0.3, 0.4)
  if (index %% 20L == 5L) {
    for (condition in conditions) proportions[[condition]] <- c(0.65, 0.35, 0)
    proportions$T1 <- c(0.45, 0.25, 0.30)
  }
  counts <- simulate_dm_counts(proportions, 60, c(60, 150, 300)[index %% 3L + 1L], sample_condition)
  list(
    gene_id = sprintf("g%04d", index),
    contig = "chrB",
    coordinates = 1000L * index + c(0L, 50L, 100L),
    counts = counts
  )
})
dataset <- assemble_usage_dataset(genes, sample_ids, sample_condition, "C", statistics_params(replicates))
work <- tempfile("pacusage-benchmark-")
paths <- write_usage_inputs(dataset, file.path(work, "inputs"))
fit_directory <- file.path(work, "fit")

fit_seconds <- system.time(suppressMessages(model$run_fit_mode(list(
  family = "C",
  samples = paths$samples,
  counts = paths$counts,
  atlas = paths$atlas,
  params = paths$params,
  model_workers = as.character(workers),
  bootstrap_batch_size = "500",
  output_dir = fit_directory
))))[["elapsed"]]

batches <- list_batches(fit_directory)
bootstrapped <- 0L
bootstrap_seconds <- 0
for (batch_path in batches) {
  batch <- readRDS(batch_path)
  if (isTRUE(batch$empty)) next
  if (identical(batch$schema_version, 2L)) {
    keep <- utils::head(batch$genes$gene_id, max_bootstrap_genes)
    batch$genes <- batch$genes[batch$genes$gene_id %in% keep, , drop = FALSE]
    batch$counts <- batch$counts[batch$counts$gene_id %in% keep, , drop = FALSE]
    batch$fitted <- batch$fitted[batch$fitted$gene_id %in% keep, , drop = FALSE]
    batch$selections <- batch$selections[batch$selections$gene_id %in% keep, , drop = FALSE]
    genes_in_batch <- length(keep)
  } else {
    batch$genes <- utils::head(batch$genes, max_bootstrap_genes)
    genes_in_batch <- length(batch$genes)
  }
  truncated <- file.path(work, basename(batch_path))
  saveRDS(batch, truncated)
  bootstrapped <- bootstrapped + genes_in_batch
  bootstrap_seconds <- bootstrap_seconds + system.time(suppressMessages(model$run_bootstrap_batch(
    truncated, file.path(work, paste0(basename(batch_path), ".intervals.tsv.gz")), workers
  )))[["elapsed"]]
}

# Old versions scatter one batch per comparison; the repaired version covers
# every comparison of a gene in one batch. Report cost per gene and comparison.
comparisons <- treatments
per_gene_comparison <- if (identical(readRDS(batches[[1]])$schema_version, 2L)) {
  bootstrap_seconds / max(1L, bootstrapped) / comparisons
} else {
  bootstrap_seconds / max(1L, bootstrapped)
}
row <- data.frame(
  label = label,
  samples = length(sample_ids),
  treatments = treatments,
  genes = gene_count,
  workers = workers,
  bootstrap_replicates = replicates,
  fit_seconds = round(fit_seconds, 1),
  bootstrap_batches = length(batches),
  bootstrapped_gene_batches = bootstrapped,
  bootstrap_seconds = round(bootstrap_seconds, 1),
  bootstrap_seconds_per_gene_comparison = signif(per_gene_comparison, 3),
  stringsAsFactors = FALSE
)
print(row, row.names = FALSE)
if (nzchar(output)) {
  utils::write.table(row, output, sep = "\t", quote = FALSE, row.names = FALSE,
    append = file.exists(output), col.names = !file.exists(output))
}
