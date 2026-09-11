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
