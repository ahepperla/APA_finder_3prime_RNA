# Statistical simulation checks for scripts/fit_usage_model.R (design section
# "Statistical Simulation"): null calibration, effect direction and size,
# recovery of gained and lost PACs, stage-wise error, count-scaling
# invariance, and bootstrap interval coverage.
#
# Usage: Rscript tests/r/test_usage_model_simulation.R scripts/fit_usage_model.R [seed]
#
# Thresholds sit at about 2-3x nominal error rates because DRIMSeq is liberal
# for overdispersed genes at 3-4 replicates (docs/decisions.md, 2026-09-27). Interval
# coverage is about 88% at a nominal 95%, which the project lead accepted; S6
# guards the coverage pooled over all 60 PACs at 0.75.

arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments) %in% c(1L, 2L)) {
  stop("Usage: Rscript tests/r/test_usage_model_simulation.R scripts/fit_usage_model.R [seed]")
}
script <- normalizePath(arguments[[1]], mustWork = TRUE)
seed <- if (length(arguments) == 2L) as.integer(arguments[[2]]) else 11L
test_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
source(file.path(dirname(normalizePath(test_file)), "usage_model_helpers.R"))

model <- load_usage_model(script)
work <- tempfile("pacusage-simulation-")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
metrics <- list()
record <- function(name, value) {
  metrics[[name]] <<- value
  value
}

design_samples <- function(replicates) {
  conditions <- rep(names(replicates), times = replicates)
  ids <- paste0(conditions, "_", stats::ave(seq_along(conditions), conditions, FUN = seq_along))
  list(sample_ids = ids, sample_condition = conditions)
}

gene_entry <- function(gene_id, index, counts) {
  list(
    gene_id = gene_id,
    contig = "chrS",
    coordinates = 1000L * index + 50L * (seq_len(nrow(counts)) - 1L),
    counts = counts
  )
}

run_family <- function(genes, layout, control, name, replicates = 0L, workers = 2L) {
  dataset <- assemble_usage_dataset(
    genes, layout$sample_ids, layout$sample_condition, control, statistics_params(replicates)
  )
  paths <- write_usage_inputs(dataset, file.path(work, name, "inputs"))
  run <- run_statistics_in_process(
    model, paths, file.path(work, name), family = control,
    model_workers = workers, batch_size = 500L, bootstrap_workers = workers
  )
  list(dataset = dataset, run = run)
}

called <- function(pacs) {
  as_number(pacs$pac_fdr) <= 0.05 & as_number(pacs$gene_fdr) <= 0.05
}

# ---- Main family: C (4), T1 (4), T2 (3), 3 PACs per gene -------------------

set.seed(seed)
main_layout <- design_samples(c(C = 4L, T1 = 4L, T2 = 3L))
main_genes <- list()
classes <- character()
add <- function(gene_id, class, proportions, precision, depth) {
  index <- length(main_genes) + 1L
  counts <- simulate_dm_counts(proportions, precision, depth, main_layout$sample_condition)
  main_genes[[index]] <<- gene_entry(gene_id, index, counts)
  classes[[gene_id]] <<- class
}
base <- c(0.5, 0.3, 0.2)
shifted <- c(0.3, 0.3, 0.4)
for (depth in c(40, 400)) for (precision in c(20, 200)) {
  for (index in 1:30) {
    add(sprintf("null_d%d_p%d_%02d", depth, precision, index), "null",
      list(C = base, T1 = base, T2 = base), precision, depth)
  }
  for (index in 1:6) {
    add(sprintf("shift_d%d_p%d_%02d", depth, precision, index), "shift_T1",
      list(C = base, T1 = shifted, T2 = base), precision, depth)
  }
}
without <- c(0.65, 0.35, 0)
with_site <- c(0.45, 0.25, 0.30)
for (depth in c(60, 300)) for (index in 1:3) {
  add(sprintf("gained_T1_d%d_%d", depth, index), "gained_T1",
    list(C = without, T1 = with_site, T2 = without), 100, depth)
  add(sprintf("gained_T2_d%d_%d", depth, index), "gained_T2",
    list(C = without, T1 = without, T2 = with_site), 100, depth)
  add(sprintf("lost_T1_d%d_%d", depth, index), "lost_T1",
    list(C = with_site, T1 = without, T2 = with_site), 100, depth)
}
main <- tryCatch(
  run_family(main_genes, main_layout, "C", "main"),
  error = function(error) error
)
require_main <- function() {
  if (inherits(main, "error")) stop("Main simulation failed: ", conditionMessage(main))
  main$run
}
main_table <- function(comparison) {
  pacs <- read_result(require_main()$final_directory, paste0(comparison, ".pacs.tsv.gz"))
  pacs$class <- unname(classes[pacs$gene_id])
  pacs$position <- stats::ave(seq_len(nrow(pacs)), pacs$gene_id, FUN = seq_along)
  pacs
}
gene_rows <- function(pacs) pacs[!duplicated(pacs$gene_id), , drop = FALSE]

test_case("S1", "null genes are calibrated", {
  for (comparison in c("T1_vs_C", "T2_vs_C")) {
    pacs <- main_table(comparison)
    nulls <- gene_rows(pacs[pacs$class == "null", , drop = FALSE])
    p <- as_number(nulls$pvalue_gene)
    check(length(p) == 120L && all(is.finite(p)), comparison, ": null gene p-values are missing.")
    at_05 <- record(paste0(comparison, " null p<=0.05"), mean(p <= 0.05))
    at_01 <- record(paste0(comparison, " null p<=0.01"), mean(p <= 0.01))
    middle <- record(paste0(comparison, " null median p"), stats::median(p))
    check(at_05 <= 0.15, comparison, ": ", format(at_05), " of null genes have p <= 0.05.")
    check(at_01 <= 0.07, comparison, ": ", format(at_01), " of null genes have p <= 0.01.")
    check(middle >= 0.35 && middle <= 0.75, comparison, ": null median p is ", format(middle), ".")
    pac_p <- as_number(pacs$pvalue_pac[pacs$class == "null"])
    pac_05 <- record(paste0(comparison, " null PAC p<=0.05"), mean(pac_p[is.finite(pac_p)] <= 0.05))
    check(pac_05 <= 0.12, comparison, ": ", format(pac_05), " of null PAC tests have p <= 0.05.")
  }
})

test_case("S2", "shifted genes have the right direction and size", {
  pacs <- main_table("T1_vs_C")
  shift <- pacs[pacs$class == "shift_T1", , drop = FALSE]
  third <- as_number(shift$delta_pau[shift$position == 3L])
  first <- as_number(shift$delta_pau[shift$position == 1L])
  check(sum(third > 0) >= 22L, "PAC3 increased in ", sum(third > 0), " of 24 genes.")
  check(sum(first < 0) >= 22L, "PAC1 decreased in ", sum(first < 0), " of 24 genes.")
  check(abs(record("shift mean PAC3 delta", mean(third)) - 0.2) <= 0.05, "mean PAC3 delta ", mean(third))
  check(abs(record("shift mean PAC1 delta", mean(first)) + 0.2) <= 0.05, "mean PAC1 delta ", mean(first))
  precise <- grepl("^shift_d400_p200", shift$gene_id) & shift$position == 3L
  check(
    all(abs(as_number(shift$delta_pau[precise]) - 0.2) <= 0.08),
    "a deep, precise gene misses the true PAC3 delta by more than 0.08."
  )
})

test_case("S3", "gained and lost PACs are recovered without cross-talk", {
  t1 <- main_table("T1_vs_C")
  t2 <- main_table("T2_vs_C")
  third <- function(pacs, class) pacs[pacs$class == class & pacs$position == 3L, , drop = FALSE]
  recovered <- c(
    gained_T1 = sum(third(t1, "gained_T1")$event_type == "gained"),
    gained_T2 = sum(third(t2, "gained_T2")$event_type == "gained"),
    lost_T1 = sum(third(t1, "lost_T1")$event_type == "lost")
  )
  record("recovered gained_T1/gained_T2/lost_T1", paste(recovered, collapse = "/"))
  check(all(recovered >= 5L), "recovered only ", paste(names(recovered), recovered, sep = "=", collapse = ", "))
  unrelated <- rbind(
    t1[t1$class == "gained_T2", , drop = FALSE],
    t2[t2$class %in% c("gained_T1", "lost_T1"), , drop = FALSE]
  )
  labels <- c("gained", "lost", "increased_usage", "decreased_usage", "gained_candidate", "lost_candidate")
  check(!any(unrelated$event_type %in% labels), "events in comparisons without an effect.")
  unrelated_genes <- gene_rows(unrelated)
  check(nrow(unrelated_genes) == 18L, "expected 18 unaffected boundary gene comparisons.")
  check(all(is.finite(as_number(unrelated_genes$pvalue_gene))), "unaffected boundary genes lack p.")
  check(all(as_number(unrelated_genes$gene_fdr) > 0.05), "an unaffected boundary gene was called.")
})

test_case("S4", "stage-wise PAC calls control false discoveries", {
  t1 <- main_table("T1_vs_C")
  calls <- called(t1)
  # True PAC effects in T1: PAC1 and PAC3 of shifted genes, and every PAC of
  # genes gained or lost in T1 (the other PACs shift to compensate).
  true_effect <- (t1$class == "shift_T1" & t1$position %in% c(1L, 3L)) |
    t1$class %in% c("gained_T1", "lost_T1")
  total <- record("T1 PAC calls", sum(calls, na.rm = TRUE))
  false_rate <- record("T1 false discovery proportion", sum(calls & !true_effect, na.rm = TRUE) / max(1, total))
  check(total >= 40L, "only ", total, " PAC calls in T1_vs_C.")
  check(false_rate <= 0.25, "false discovery proportion is ", format(false_rate), ".")
  null_genes_t1 <- length(unique(t1$gene_id[calls & t1$class == "null"]))
  check(record("T1 null genes called", null_genes_t1) <= 8L, null_genes_t1, " null genes called in T1.")
  t2 <- main_table("T2_vs_C")
  null_genes_t2 <- length(unique(t2$gene_id[called(t2) & t2$class == "null"]))
  check(record("T2 null genes called", null_genes_t2) <= 5L, null_genes_t2, " null genes called in T2.")
})

# ---- Scaling: multiplying every count of a gene leaves usage unchanged -----

set.seed(seed + 1L)
scale_layout <- design_samples(c(C = 4L, T1 = 4L))
scale_genes <- list()
scale_classes <- character()
for (index in 1:40) {
  shift <- index <= 20L
  gene_id <- sprintf("%s_%02d", if (shift) "shift" else "null", index)
  counts <- simulate_dm_counts(
    list(C = base, T1 = if (shift) c(0.3, 0.3, 0.4) else base), 100, 200,
    scale_layout$sample_condition
  )
  scale_genes[[index]] <- gene_entry(gene_id, index, counts)
  scale_classes[[gene_id]] <- if (shift) "shift" else "null"
}
scaled_family <- function(multiplier_for) {
  lapply(scale_genes, function(gene) {
    gene$counts <- gene$counts * multiplier_for(gene)
    storage.mode(gene$counts) <- "integer"
    gene
  })
}
t1_columns <- scale_layout$sample_condition == "T1"
scale_runs <- tryCatch(
  list(
    original = run_family(scale_genes, scale_layout, "C", "scale_original"),
    all = run_family(scaled_family(function(gene) 3L), scale_layout, "C", "scale_all"),
    t1_nulls = run_family(scaled_family(function(gene) {
      if (scale_classes[[gene$gene_id]] == "null") {
        matrix(rep(ifelse(t1_columns, 3L, 1L), each = nrow(gene$counts)), nrow = nrow(gene$counts))
      } else {
        1L
      }
    }), scale_layout, "C", "scale_t1_nulls")
  ),
  error = function(error) error
)
scale_table <- function(name) {
  if (inherits(scale_runs, "error")) stop("Scaling runs failed: ", conditionMessage(scale_runs))
  pacs <- read_result(scale_runs[[name]]$run$final_directory, "T1_vs_C.pacs.tsv.gz")
  pacs[order(pacs$pac_id), , drop = FALSE]
}

test_case("S5", "scaling a gene's counts leaves its usage and calls unchanged", {
  original <- scale_table("original")
  everything <- scale_table("all")
  t1_nulls <- scale_table("t1_nulls")
  change <- function(a, b, column) max(abs(as_number(a[[column]]) - as_number(b[[column]])))
  check(record("scale all max PAU change", change(original, everything, "fitted_control_pau")) <= 0.04,
    "fitted control PAU moved by more than 0.04.")
  check(change(original, everything, "delta_pau") <= 0.04, "delta PAU moved by more than 0.04.")
  shift_rows <- scale_classes[original$gene_id] == "shift"
  agreement <- sum(called(original)[shift_rows] == called(everything)[shift_rows], na.rm = TRUE) / 3
  check(record("scale all shift call agreement (genes)", agreement) >= 19, "shift calls agree for only ", agreement, " genes.")
  null_rows <- !shift_rows
  null_difference <- abs(sum(called(original)[null_rows], na.rm = TRUE) - sum(called(everything)[null_rows], na.rm = TRUE))
  check(null_difference <= 6L, "null PAC calls changed by ", null_difference, ".")
  check(record("scale T1 nulls max delta change", change(original[null_rows, ], t1_nulls[null_rows, ], "delta_pau")) <= 0.04,
    "scaling one group's counts moved null deltas by more than 0.04.")
  t1_difference <- abs(sum(called(original)[null_rows], na.rm = TRUE) - sum(called(t1_nulls)[null_rows], na.rm = TRUE))
  check(t1_difference <= 3L, "null PAC calls changed by ", t1_difference, " after scaling T1.")
})

# ---- Bootstrap interval coverage -------------------------------------------

set.seed(seed + 2L)
coverage_layout <- design_samples(c(C = 4L, T1 = 4L))
coverage_genes <- list()
coverage_truth <- numeric()
for (index in 1:60) {
  gained <- index <= 30L
  gene_id <- sprintf("%s_%02d", if (gained) "gained" else "shift", index)
  control <- if (gained) without else base
  treatment <- if (gained) with_site else shifted
  counts <- simulate_dm_counts(list(C = control, T1 = treatment), 100, 200, coverage_layout$sample_condition)
  coverage_genes[[index]] <- gene_entry(gene_id, index, counts)
  coverage_truth[[gene_id]] <- treatment[[3]] - control[[3]]
}
coverage <- tryCatch(
  run_family(coverage_genes, coverage_layout, "C", "coverage", replicates = 100L, workers = 4L),
  error = function(error) error
)

test_case("S6", "bootstrap intervals cover the true change in usage", {
  if (inherits(coverage, "error")) stop("Coverage run failed: ", conditionMessage(coverage))
  pacs <- read_result(coverage$run$final_directory, "T1_vs_C.pacs.tsv.gz")
  pacs$position <- stats::ave(seq_len(nrow(pacs)), pacs$gene_id, FUN = seq_along)
  third <- pacs[pacs$position == 3L, , drop = FALSE]
  low <- as_number(third$delta_pau_ci_low)
  high <- as_number(third$delta_pau_ci_high)
  truth <- unname(coverage_truth[third$gene_id])
  gained <- grepl("^gained", third$gene_id)
  check(all(is.finite(low) & is.finite(high)), sum(!is.finite(low)), " PACs have no interval.")
  covered <- low <= truth & truth <= high
  record("coverage, gained (stabilized) PACs", mean(covered[gained]))
  record("coverage, shifted PACs", mean(covered[!gained]))
  pooled_coverage <- record("coverage, all PACs", mean(covered))
  check(pooled_coverage >= 0.75, "pooled coverage is ", format(pooled_coverage), ".")
})

cat("\nMetrics (seed ", seed, "):\n", sep = "")
for (name in names(metrics)) cat(sprintf("  %-42s %s\n", name, format(metrics[[name]], digits = 4)))
finish_tests()
