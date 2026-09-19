args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript tests/test_bootstrap_total.R scripts/fit_usage_model.R")
}

lines <- readLines(args[[1]])
start <- grep("^empty_bootstrap_intervals <-", lines)[1]
end <- grep("^stabilize_boundary_gene <-", lines)[1] - 1L
if (is.na(start) || is.na(end) || start > end) {
  stop("Could not locate bootstrap helpers in ", args[[1]])
}
eval(parse(text = lines[start:end]))

stopifnot(identical(bootstrap_total(c(1, 2), "gene-1", "sample-1"), 3L))
stopifnot(identical(bootstrap_total(c(0, 0), "gene-1", "sample-1"), 0L))

expect_invalid_total <- function(values) {
  message <- tryCatch(
    {
      bootstrap_total(values, "gene-1", "sample-1")
      ""
    },
    error = conditionMessage
  )
  stopifnot(grepl(
    "Invalid bootstrap total for gene gene-1, sample sample-1",
    message,
    fixed = TRUE
  ))
}

expect_invalid_total(c(1, NA_real_))
expect_invalid_total(c(0.5, 0.6))
expect_invalid_total(c(.Machine$integer.max, 1))

stopifnot(identical(
  bootstrap_gene_ids(
    c("gene-1", "gene-2", NA_character_, "gene-4"),
    c(0.01, NA_real_, 0.01, 0.50),
    c(FALSE, NA, FALSE, TRUE),
    0.05
  ),
  c("gene-1", "gene-4")
))
stopifnot(identical(
  bootstrap_gene_ids(
    c("gene-1", "gene-2"),
    c(0.01, NA_real_),
    c(FALSE, TRUE),
    0.05,
    include_candidates = NULL
  ),
  c("gene-1", "gene-2")
))
stopifnot(identical(
  bootstrap_gene_ids(
    c("gene-1", "gene-2", "gene-3"),
    c(0.01, NA_real_, 0.50),
    c(FALSE, TRUE, TRUE),
    0.05,
    include_candidates = FALSE
  ),
  "gene-1"
))

bootstrap_probe <- function(repeat_number) {
  set.seed(1000L + repeat_number)
  c(repeat_number, runif(1))
}

serial_runs <- bootstrap_apply(1:4, 1L, bootstrap_probe)
if (.Platform$OS.type != "windows") {
  parallel_runs <- bootstrap_apply(1:4, 2L, bootstrap_probe)
  stopifnot(identical(serial_runs, parallel_runs))
}

stopifnot(identical(
  dominant_pac(c("pac-1", "pac-2"), c(NA_real_, NA_real_)),
  NA_character_
))
stopifnot(identical(
  dominant_pac(c("pac-1", "pac-2"), c(NA_real_, 0.25)),
  "pac-2"
))
stopifnot(identical(
  stable_feature_matches(
    c(1L, NA_integer_, 3L, 4L),
    c(1L, 2L, NA_integer_, 4L),
    c(FALSE, FALSE, FALSE, NA)
  ),
  c(TRUE, FALSE, FALSE, FALSE)
))

stabilize_start <- grep("^stabilize_boundary_gene <-", lines)[1]
stabilize_end <- grep("^bootstrap_gene <-", lines)[1] - 1L
if (is.na(stabilize_start) || is.na(stabilize_end) || stabilize_start > stabilize_end) {
  stop("Could not locate boundary stabilization implementation in ", args[[1]])
}
eval(parse(text = lines[stabilize_start:stabilize_end]))

dmDSdata <- function(counts, samples) list(counts = counts, samples = samples)
dmPrecision <- function(value, design, verbose) value
dmFit <- function(value, design, add_uniform = FALSE, verbose) value
dmTest <- function(value, coef, verbose) value
results <- function(value, level) {
  data.frame(
    feature_id = value$counts$feature_id,
    pvalue = c(0.02, 0.03),
    stringsAsFactors = FALSE
  )
}
proportions <- function(value) {
  data.frame(
    gene_id = value$counts$gene_id,
    feature_id = value$counts$feature_id,
    control_1 = c(0.1, 0.9),
    treatment_1 = c(0.7, 0.3),
    stringsAsFactors = FALSE
  )
}
boundary_inputs <- list(
  gene_id = "gene-1",
  counts = data.frame(
    gene_id = c("gene-1", "gene-1"),
    pac_id = c("pac-1", "pac-2"),
    control_1 = c(1L, 9L),
    treatment_1 = c(7L, 3L)
  ),
  sample_rows = data.frame(
    sample_id = c("control_1", "treatment_1"),
    condition = c("control", "treatment")
  ),
  design = matrix(c(1, 1), ncol = 1),
  coefficient = 1L,
  control = "control",
  treatment = "treatment",
  params = list(
    model_covariates = list(),
    dm_zero_sensitivity_repeats = 4L,
    dm_bootstrap_min_success_fraction = 0.8,
    dm_zero_max_delta_pau_spread = 0.02,
    random_seed = 1729
  ),
  atlas_checksum = "atlas",
  comparison = "treatment_vs_control"
)
serial_boundary <- do.call(
  stabilize_boundary_gene,
  c(boundary_inputs, list(bootstrap_workers = 1L))
)
if (.Platform$OS.type != "windows") {
  parallel_boundary <- do.call(
    stabilize_boundary_gene,
    c(boundary_inputs, list(bootstrap_workers = 2L))
  )
  stopifnot(identical(serial_boundary, parallel_boundary))
}

bootstrap_start <- grep("^bootstrap_gene <-", lines)[1]
bootstrap_end <- grep("^classify_event <-", lines)[1] - 1L
if (is.na(bootstrap_start) || is.na(bootstrap_end) || bootstrap_start > bootstrap_end) {
  stop("Could not locate bootstrap implementation in ", args[[1]])
}
eval(parse(text = lines[bootstrap_start:bootstrap_end]))

batch_payloads <- bootstrap_batches(
  comparison = "treatment_vs_control",
  bootstrap_genes = c("gene-1", "gene-2", "gene-3"),
  counts = data.frame(
    gene_id = rep(c("gene-1", "gene-2", "gene-3"), each = 2L),
    pac_id = paste0("pac-", seq_len(6L)),
    control_1 = c(3L, 1L, 2L, 2L, 1L, 3L),
    treatment_1 = c(1L, 3L, 2L, 2L, 3L, 1L)
  ),
  fitted = data.frame(
    gene_id = rep(c("gene-1", "gene-2", "gene-3"), each = 2L),
    feature_id = paste0("pac-", seq_len(6L)),
    control_1 = rep(c(0.75, 0.25), 3L),
    treatment_1 = rep(c(0.25, 0.75), 3L)
  ),
  precision = data.frame(
    gene_id = c("gene-1", "gene-2", "gene-3"),
    precision = c(10, 20, 30)
  ),
  sample_rows = data.frame(
    sample_id = c("control_1", "treatment_1"),
    condition = c("control", "treatment")
  ),
  design = matrix(c(1, 1), ncol = 1),
  control = "control",
  treatment = "treatment",
  params = list(),
  atlas_checksum = "atlas",
  batch_size = 2L
)
stopifnot(length(batch_payloads) == 2L)
stopifnot(identical(
  vapply(batch_payloads, function(batch) length(batch$genes), integer(1)),
  c(2L, 1L)
))
stopifnot(identical(
  vapply(
    batch_payloads[[1]]$genes,
    function(gene) gene$gene_id,
    character(1)
  ),
  c("gene-1", "gene-2")
))

finalize_start <- grep("^read_bootstrap_intervals <-", lines)[1]
finalize_end <- grep("^fit_motif_preferences <-", lines)[1] - 1L
if (is.na(finalize_start) || is.na(finalize_end) || finalize_start > finalize_end) {
  stop("Could not locate bootstrap finalization helpers in ", args[[1]])
}
eval(parse(text = lines[finalize_start:finalize_end]))
read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}
expect_error <- function(expression, pattern) {
  message <- tryCatch(
    {
      force(expression)
      ""
    },
    error = conditionMessage
  )
  stopifnot(grepl(pattern, message, fixed = TRUE))
}

valid_interval_path <- tempfile(fileext = ".tsv")
write.table(
  data.frame(
    comparison = "treatment_vs_control",
    gene_id = "gene-1",
    feature_id = "pac-1",
    delta_pau_ci_low = NA_real_,
    delta_pau_ci_high = NA_real_,
    bootstrap_successes = 15L
  ),
  valid_interval_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
valid_intervals <- read_bootstrap_intervals(valid_interval_path)
stopifnot(is.na(valid_intervals$delta_pau_ci_low[[1]]))
stopifnot(identical(valid_intervals$bootstrap_successes, 15L))

expect_error(
  read_bootstrap_intervals(character()),
  "No bootstrap interval files were provided"
)

missing_column_path <- tempfile(fileext = ".tsv")
write.table(
  data.frame(comparison = "treatment_vs_control", gene_id = "gene-1"),
  missing_column_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
expect_error(
  read_bootstrap_intervals(missing_column_path),
  "missing required columns"
)

duplicate_interval_path <- tempfile(fileext = ".tsv")
write.table(
  rbind(valid_intervals, valid_intervals),
  duplicate_interval_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
expect_error(
  read_bootstrap_intervals(duplicate_interval_path),
  "duplicate comparison/gene/PAC rows"
)

null_batch_path <- tempfile(fileext = ".rds")
saveRDS(NULL, null_batch_path)
expect_error(
  run_bootstrap_batch(null_batch_path, tempfile(fileext = ".tsv.gz"), 1L),
  "Bootstrap batch payload must be a list"
)

invalid_precision_intervals <- suppressWarnings(bootstrap_gene(
  gene_id = "gene-1",
  counts = data.frame(
    gene_id = c("gene-1", "gene-1"),
    pac_id = c("pac-1", "pac-2")
  ),
  sample_rows = NULL,
  design = NULL,
  fitted_model = NULL,
  precision = NA_real_,
  control = NULL,
  treatment = NULL,
  params = NULL,
  atlas_checksum = NULL,
  comparison = NULL,
  bootstrap_workers = NULL
))
stopifnot(identical(
  invalid_precision_intervals$feature_id,
  c("pac-1", "pac-2")
))
stopifnot(all(is.na(invalid_precision_intervals$delta_pau_ci_low)))
stopifnot(all(is.na(invalid_precision_intervals$delta_pau_ci_high)))
stopifnot(identical(
  invalid_precision_intervals$bootstrap_successes,
  c(0L, 0L)
))

disabled_bootstrap_intervals <- bootstrap_gene(
  gene_id = "gene-1",
  counts = data.frame(
    gene_id = c("gene-1", "gene-1"),
    pac_id = c("pac-1", "pac-2")
  ),
  sample_rows = NULL,
  design = NULL,
  fitted_model = NULL,
  precision = 10,
  control = NULL,
  treatment = NULL,
  params = list(
    dm_bootstrap_replicates = 0L,
    dm_bootstrap_min_success_fraction = 0.8
  ),
  atlas_checksum = NULL,
  comparison = NULL,
  bootstrap_workers = NULL
)
stopifnot(all(is.na(disabled_bootstrap_intervals$delta_pau_ci_low)))
stopifnot(all(is.na(disabled_bootstrap_intervals$delta_pau_ci_high)))
stopifnot(identical(
  disabled_bootstrap_intervals$bootstrap_successes,
  c(0L, 0L)
))

proportions <- function(value) {
  data.frame(
    gene_id = c("gene-1", "gene-1"),
    feature_id = c("pac-1", "pac-2"),
    control_1 = c(0.75, 0.25),
    treatment_1 = c(0.25, 0.75)
  )
}
bootstrap_apply <- function(repeat_numbers, workers, worker) {
  lapply(
    repeat_numbers,
    function(repeat_number) structure(
      "Error in FUN(X[[i]], ...) : simulated worker failure\n",
      class = "try-error"
    )
  )
}
worker_failure_intervals <- suppressWarnings(bootstrap_gene(
  gene_id = "gene-1",
  counts = data.frame(
    gene_id = c("gene-1", "gene-1"),
    pac_id = c("pac-1", "pac-2"),
    control_1 = c(6L, 2L),
    treatment_1 = c(2L, 6L)
  ),
  sample_rows = data.frame(
    sample_id = c("control_1", "treatment_1"),
    condition = c("control", "treatment")
  ),
  design = NULL,
  fitted_model = NULL,
  precision = 10,
  control = "control",
  treatment = "treatment",
  params = list(
    model_covariates = list(),
    dm_bootstrap_replicates = 2L,
    dm_bootstrap_min_success_fraction = 0.8,
    random_seed = 1729
  ),
  atlas_checksum = "atlas",
  comparison = "treatment_vs_control",
  bootstrap_workers = 2L
))
stopifnot(all(is.na(worker_failure_intervals$delta_pau_ci_low)))
stopifnot(all(is.na(worker_failure_intervals$delta_pau_ci_high)))
stopifnot(identical(
  worker_failure_intervals$bootstrap_successes,
  c(0L, 0L)
))

classify_start <- grep("^classify_event <-", lines)[1]
classify_end <- grep("^fit_family <-", lines)[1] - 1L
if (is.na(classify_start) || is.na(classify_end) || classify_start > classify_end) {
  stop("Could not locate event classifier in ", args[[1]])
}
eval(parse(text = lines[classify_start:classify_end]))

event_params <- list(
  event_min_supporting_samples = 2,
  min_abs_delta_pau = 0.10,
  gene_fdr = 0.05,
  site_fdr = 0.05,
  event_max_control_pau = 0.01,
  event_min_treatment_pau = 0.05
)
gained_row <- list(
  control_supporting_samples = 0,
  treatment_supporting_samples = 2,
  delta_pau = 0.20,
  gene_fdr = 0.01,
  pac_fdr = 0.01,
  zero_boundary_unstable = FALSE,
  confidence = "high",
  internal_priming_flag = FALSE,
  exploratory_insufficient_replicates = FALSE,
  fitted_control_pau = 0,
  fitted_treatment_pau = 0.20
)
stopifnot(identical(classify_event(gained_row, event_params), "gained"))

missing_fitted_row <- gained_row
missing_fitted_row$fitted_control_pau <- NA_real_
stopifnot(identical(classify_event(missing_fitted_row, event_params), "none"))

missing_support_row <- gained_row
missing_support_row$control_supporting_samples <- NA_real_
stopifnot(identical(classify_event(missing_support_row, event_params), "none"))
