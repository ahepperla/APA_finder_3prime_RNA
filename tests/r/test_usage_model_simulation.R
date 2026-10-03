# Statistical simulation checks for scripts/fit_usage_model.R (design section
# "Statistical Simulation"): null calibration, effect direction and size,
# recovery of gained and lost PACs, stage-wise error, count-scaling
# invariance, bootstrap interval coverage, and calibration after the usage
# filter.
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

run_family <- function(genes, layout, control, name, replicates = 0L, workers = 2L,
                       params = statistics_params(replicates)) {
  dataset <- assemble_usage_dataset(
    genes, layout$sample_ids, layout$sample_condition, control, params
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
    p <- as_number(nulls$gene_pvalue)
    check(length(p) == 120L && all(is.finite(p)), comparison, ": null gene p-values are missing.")
    at_05 <- record(paste0(comparison, " null p<=0.05"), mean(p <= 0.05))
    at_01 <- record(paste0(comparison, " null p<=0.01"), mean(p <= 0.01))
    middle <- record(paste0(comparison, " null median p"), stats::median(p))
    check(at_05 <= 0.15, comparison, ": ", format(at_05), " of null genes have p <= 0.05.")
    check(at_01 <= 0.07, comparison, ": ", format(at_01), " of null genes have p <= 0.01.")
    check(middle >= 0.35 && middle <= 0.75, comparison, ": null median p is ", format(middle), ".")
    pac_p <- as_number(pacs$pac_pvalue[pacs$class == "null"])
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
  check(all(is.finite(as_number(unrelated_genes$gene_pvalue))), "unaffected boundary genes lack p.")
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

# ---- Usage filter -----------------------------------------------------------

# Null genes, alike in every condition, with a borderline PAC at 8% between a
# dominant one and a minor one, beside genes whose usage shifts in T1, so the
# false discovery rate has real discoveries to work with. With min_site_usage
# at 0.10, the borderline PAC is tested only when chance lifts it to 10% in
# as many samples as the smallest condition has (3 here), from any
# conditions. The filter never sees which samples belong to which condition,
# so the PACs it admits should stay calibrated. S7 checks that at the lead's
# planned settings: a PAC-call threshold of 0.05, and calls at 0.10 taken as
# those with a change of at least 0.10.
set.seed(seed + 3L)
filter_layout <- design_samples(c(C = 4L, T1 = 4L, T2 = 3L))
borderline <- c(0.60, 0.08, 0.32)
moved <- c(0.40, 0.08, 0.52)
filter_genes <- list()
filter_classes <- character()
add_filter_gene <- function(gene_id, class, proportions, precision, depth) {
  counts <- simulate_dm_counts(proportions, precision, depth, filter_layout$sample_condition)
  index <- length(filter_genes) + 1L
  filter_genes[[index]] <<- gene_entry(gene_id, index, counts)
  filter_classes[[gene_id]] <<- class
}
for (depth in c(40, 400)) for (precision in c(20, 200)) {
  for (index in 1:100) {
    add_filter_gene(sprintf("null_d%d_p%d_%03d", depth, precision, index), "null",
      list(C = borderline, T1 = borderline, T2 = borderline), precision, depth)
  }
  for (index in 1:25) {
    add_filter_gene(sprintf("shift_d%d_p%d_%02d", depth, precision, index), "shift_T1",
      list(C = borderline, T1 = moved, T2 = borderline), precision, depth)
  }
}
filter_params <- utils::modifyList(statistics_params(0L), list(
  min_site_usage = 0.10, min_site_usage_gene_reads = 10L, min_abs_delta_pau = 0.05
))
filter_run <- tryCatch(
  run_family(filter_genes, filter_layout, "C", "filter", params = filter_params),
  error = function(error) error
)

test_case("S7", "PACs the usage filter admits near its threshold stay calibrated", {
  if (inherits(filter_run, "error")) stop("Filter run failed: ", conditionMessage(filter_run))
  directory <- filter_run$run$final_directory
  filtering <- read_result(directory, "C.statistical_filtering.tsv.gz")
  filtering <- filtering[order(filtering$gene_id, as_number(filtering$start)), , drop = FALSE]
  filtering$position <- stats::ave(seq_len(nrow(filtering)), filtering$gene_id, FUN = seq_along)
  check(all(as_flag(filtering$tested[filtering$position != 2L])), "a dominant or minor PAC was filtered.")
  border <- filtering[filtering$position == 2L, , drop = FALSE]
  kept <- as_flag(border$tested)
  check(all(border$reason[!kept] == "sample_usage<0.1"), "borderline PACs filtered for another reason: ",
    paste(unique(border$reason[!kept]), collapse = ", "))
  nulls <- names(filter_classes)[filter_classes == "null"]
  keeping <- intersect(border$gene_id[kept], nulls)
  share <- record("S7 null borderline PACs kept", length(keeping) / length(nulls))
  check(share > 0.1 && share < 0.9, "the filter kept ", format(share), " of borderline PACs; S7 needs both kinds.")
  # Every metric is recorded before any guard is checked.
  rates <- list()
  for (comparison in c("T1_vs_C", "T2_vs_C")) {
    pacs <- read_result(directory, paste0(comparison, ".pacs.tsv.gz"))
    label <- function(text) paste("S7", comparison, text)
    genes <- pacs[!duplicated(pacs$gene_id) & pacs$gene_id %in% nulls, , drop = FALSE]
    gene_p <- as_number(genes$gene_pvalue)
    with_border <- genes$gene_id %in% keeping
    border_p <- as_number(pacs$pac_pvalue[pacs$pac_id %in% border$pac_id[kept] & pacs$gene_id %in% nulls])
    calls <- pacs$event_type %in% c("gained", "lost", "increased_usage", "decreased_usage")
    large <- calls & abs(as_number(pacs$delta_pau)) >= 0.10
    called <- unique(pacs$gene_id[calls])
    rates[[comparison]] <- list(
      gene_with = record(label("null gene p<=0.05, borderline kept"), mean(gene_p[with_border] <= 0.05)),
      gene_without = record(label("null gene p<=0.05, borderline filtered"), mean(gene_p[!with_border] <= 0.05)),
      border = record(label("kept borderline PAC p<=0.05"), mean(border_p[is.finite(border_p)] <= 0.05)),
      called = record(label("genes called at 0.05"), length(called)),
      null_05 = record(label("null genes called at 0.05"), sum(called %in% nulls)),
      null_10 = record(label("null genes called at 0.10"), length(intersect(unique(pacs$gene_id[large]), nulls))),
      null_kept = record(label("null genes called at 0.05, borderline kept"), sum(called %in% keeping))
    )
  }
  t1 <- rates[["T1_vs_C"]]
  false_share <- record("S7 T1_vs_C false share of genes called at 0.05", t1$null_05 / max(1, t1$called))
  check(t1$called >= 50L, "only ", t1$called, " of 100 shifted genes called in T1_vs_C.")
  # As S4: a false discovery proportion of at most 0.25, and at most 8 null
  # genes called per 120.
  check(false_share <= 0.25, "false share of T1_vs_C gene calls is ", format(false_share), ".")
  for (comparison in names(rates)) {
    rate <- rates[[comparison]]
    check(rate$null_05 <= 27L, comparison, ": ", rate$null_05, " null genes called.")
    # As S1 for unfiltered nulls: at most 0.12 of PAC tests, and 0.15 of
    # genes, with p <= 0.05.
    check(rate$border <= 0.12, comparison, ": ", format(rate$border),
      " of admitted borderline PACs have p <= 0.05.")
    check(rate$gene_with <= 0.15, comparison, ": ", format(rate$gene_with),
      " of null genes keeping a borderline PAC have p <= 0.05.")
  }
})

# ---- Depth in each group ----------------------------------------------------

# A two-condition family, so the family test is the comparison's test, with
# the control at depth 400 and the treatment's depth set per gene:
# - lead's case: the control at 10/90; the treatment at 50/50 but nearly
#   silent, at about 0.3 reads ("off") or 3 reads ("low") per sample, so a
#   handful of reads read as an even split;
# - nulls alike in both groups, the treatment at depths from 2 to 400;
# - genes turned on: the control nearly silent, the treatment at depth 400;
# - shifts at depth 400 in both groups.
set.seed(seed + 4L)
depth_layout <- design_samples(c(C = 4L, T = 4L))
depth_genes <- list()
depth_classes <- character()
add_depth_gene <- function(gene_id, class, proportions, depths, precision = 100) {
  counts <- do.call(cbind, lapply(names(depths), function(condition) {
    columns <- depth_layout$sample_condition == condition
    simulate_dm_counts(proportions[condition], precision, depths[[condition]],
      depth_layout$sample_condition[columns])
  }))
  index <- length(depth_genes) + 1L
  depth_genes[[index]] <<- gene_entry(gene_id, index, counts)
  depth_classes[[gene_id]] <<- class
}
for (index in 1:20) {
  for (level in c(off = 0.3, low = 3)) {
    name <- if (level < 1) "off" else "low"
    add_depth_gene(sprintf("lead_%s_%02d", name, index), paste0("lead_", name),
      list(C = c(0.1, 0.9), T = c(0.5, 0.5)), c(C = 400, T = level))
  }
}
for (depth in c(2, 5, 10, 20, 40, 400)) for (index in 1:20) {
  add_depth_gene(sprintf("null_t%d_%02d", depth, index), "null",
    list(C = c(0.5, 0.3, 0.2), T = c(0.5, 0.3, 0.2)), c(C = 400, T = depth))
}
for (index in 1:15) {
  add_depth_gene(sprintf("on_%02d", index), "on", list(C = c(0.5, 0.5), T = c(0.5, 0.5)),
    c(C = 0.2, T = 400))
}
for (index in 1:12) {
  add_depth_gene(sprintf("shift_%02d", index), "shift",
    list(C = c(0.5, 0.3, 0.2), T = c(0.3, 0.3, 0.4)), c(C = 400, T = 400))
}
depth_run <- tryCatch(
  run_family(depth_genes, depth_layout, "C", "depth"),
  error = function(error) error
)

test_case("S8", "a comparison tests only genes with depth in both groups", {
  if (inherits(depth_run, "error")) stop("Depth run failed: ", conditionMessage(depth_run))
  directory <- depth_run$run$final_directory
  pacs <- read_result(directory, "T_vs_C.pacs.tsv.gz")
  genes <- read_result(directory, "T_vs_C.genes.tsv.gz")
  skipped <- read_result(directory, "T_vs_C.genes_without_depth.tsv.gz")
  omnibus <- read_result(directory, "C.gene_omnibus.tsv.gz")
  filtering <- read_result(directory, "C.statistical_filtering.tsv.gz")
  status <- stats::setNames(skipped$depth_status, skipped$gene_id)
  class_of <- function(ids) unname(depth_classes[ids])

  # The gate, recomputed from the counts at the family's tested PACs.
  counts <- depth_run$dataset$counts
  counts <- counts[counts$pac_id %in% filtering$pac_id[as_flag(filtering$tested)], , drop = FALSE]
  totals <- rowsum(as.matrix(counts[, depth_layout$sample_ids]), counts$gene_id)
  deep <- function(condition) {
    rowSums(totals[, depth_layout$sample_condition == condition, drop = FALSE] >= 10) >= 2
  }
  expected <- rownames(totals)[deep("C") & deep("T")]
  check(setequal(genes$gene_id, expected), "tested genes differ from the rule: ",
    paste(utils::head(setdiff(union(genes$gene_id, expected), intersect(genes$gene_id, expected))), collapse = ", "))
  check(setequal(c(genes$gene_id, skipped$gene_id), omnibus$gene_id) && !any(skipped$gene_id %in% genes$gene_id),
    "the comparison's genes and its genes without depth do not partition the family's.")
  check(!any(pacs$gene_id %in% skipped$gene_id), "a gene without depth has PAC rows.")

  # Lead's case: the family test, which still sees these genes, often calls
  # them; the comparison tests none of them.
  lead <- names(depth_classes)[startsWith(depth_classes, "lead_")]
  lead_p <- as_number(omnibus$gene_pvalue[omnibus$gene_id %in% lead])
  lead_share <- record("S8 lead's-case genes with family p<=0.05", mean(lead_p[is.finite(lead_p)] <= 0.05))
  off <- names(depth_classes)[depth_classes == "lead_off"]
  low <- names(depth_classes)[depth_classes == "lead_low"]
  off_share <- record("S8 nearly silent genes turned off", mean(status[off] %in% "turned_off"))
  low_share <- record("S8 low genes too low in the treatment", mean(status[low] %in% "too_low_in_treatment"))
  check(sum(is.finite(lead_p)) >= 15L && lead_share >= 0.2, "only ", format(lead_share), " of ",
    sum(is.finite(lead_p)), " lead's-case genes had family p <= 0.05; S8 needs them to.")
  check(all(lead %in% skipped$gene_id), "lead's-case genes were tested: ", paste(setdiff(lead, skipped$gene_id), collapse = ", "))
  check(all(status[lead] %in% c("turned_off", "too_low_in_treatment")), "lead's-case statuses: ",
    paste(unique(status[lead]), collapse = ", "))
  check(off_share >= 0.5 && low_share >= 0.5, "turned off ", format(off_share), ", too low ", format(low_share), ".")

  # Nulls the comparison keeps stay calibrated, as in S1.
  nulls <- genes[class_of(genes$gene_id) == "null", , drop = FALSE]
  p <- as_number(nulls$gene_pvalue)
  kept <- record("S8 null genes tested", length(p))
  at_05 <- record("S8 tested null p<=0.05", mean(p <= 0.05))
  at_01 <- record("S8 tested null p<=0.01", mean(p <= 0.01))
  middle <- record("S8 tested null median p", stats::median(p))
  check(kept >= 60L && kept < 120L, kept, " of 120 null genes were tested; S8 needs both kinds.")
  check(at_05 <= 0.15 && at_01 <= 0.07, "tested nulls: ", format(at_05), " at 0.05, ", format(at_01), " at 0.01.")
  check(middle >= 0.35 && middle <= 0.75, "tested null median p is ", format(middle), ".")
  check(isTRUE(all.equal(as_number(genes$gene_fdr), model$bh(as_number(genes$gene_pvalue)))),
    "gene_fdr is not BH over the tested genes.")

  # Genes turned on, and shifts at full depth.
  on <- names(depth_classes)[depth_classes == "on"]
  on_share <- record("S8 nearly silent controls turned on", mean(status[on] %in% "turned_on"))
  check(all(status[on] %in% c("turned_on", "too_low_in_control")) && on_share >= 0.6,
    "turned on: ", paste(status[on], collapse = ", "))
  shifts <- names(depth_classes)[depth_classes == "shift"]
  calls <- pacs$event_type %in% c("gained", "lost", "increased_usage", "decreased_usage")
  shifted <- record("S8 shifts called", length(intersect(unique(pacs$gene_id[calls]), shifts)))
  check(all(shifts %in% genes$gene_id) && shifted >= 10L, shifted, " of 12 shifts called.")
})

cat("\nMetrics (seed ", seed, "):\n", sep = "")
for (name in names(metrics)) cat(sprintf("  %-42s %s\n", name, format(metrics[[name]], digits = 4)))
finish_tests()
