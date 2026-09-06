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

stopifnot(identical(bootstrap_total(c(1, 2)), 3L))
stopifnot(identical(bootstrap_total(c(0, 0)), 0L))
stopifnot(is.na(bootstrap_total(c(1, NA_real_))))
stopifnot(is.na(bootstrap_total(c(0.5, 0.6))))
stopifnot(is.na(bootstrap_total(c(.Machine$integer.max, 1))))
