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

bootstrap_start <- grep("^bootstrap_gene <-", lines)[1]
bootstrap_end <- grep("^classify_event <-", lines)[1] - 1L
if (is.na(bootstrap_start) || is.na(bootstrap_end) || bootstrap_start > bootstrap_end) {
  stop("Could not locate bootstrap implementation in ", args[[1]])
}
eval(parse(text = lines[bootstrap_start:bootstrap_end]))

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
