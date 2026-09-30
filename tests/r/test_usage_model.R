# Tests for scripts/fit_usage_model.R with real DRIMSeq and stageR.
#
# Usage: Rscript tests/r/test_usage_model.R scripts/fit_usage_model.R
#
# Cases marked M cover the model and its outputs; cases marked H cover helpers.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: Rscript tests/r/test_usage_model.R scripts/fit_usage_model.R")
}
script <- normalizePath(arguments[[1]], mustWork = TRUE)
test_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
source(file.path(dirname(normalizePath(test_file)), "usage_model_helpers.R"))

model <- load_usage_model(script)
install_drimseq_counters()

replicates <- 40L
dataset <- build_usage_dataset(replicates)
work <- tempfile("pacusage-usage-model-")
paths <- write_usage_inputs(dataset, file.path(work, "inputs"))

# Run A: in-process, serial, one gene per bootstrap batch, DRIMSeq calls counted.
# The ambient seed is set to show that results do not depend on it.
set.seed(1)
run_a <- tryCatch(
  run_statistics_in_process(
    model, paths, file.path(work, "run_a"), batch_size = 1L, count = TRUE
  ),
  error = function(error) error
)
# Run B: separate processes, two workers, one large batch.
run_b <- tryCatch(
  run_statistics_subprocess(script, paths, file.path(work, "run_b")),
  error = function(error) error
)

require_run <- function(run) {
  if (inherits(run, "error")) stop("Statistics run failed: ", conditionMessage(run))
  run
}

pac <- function(gene_index, position) {
  sprintf("PACv1.test.chrT.+.%d", 1000L * gene_index + 50L * (position - 1L))
}
gain_pac3 <- pac(1L, 3L)
unrelated_pac1 <- pac(2L, 1L)
unrelated_pac2 <- pac(2L, 2L)
unrelated_pac3 <- pac(2L, 3L)
boundary_genes <- c("gene_gain", "gene_unrelated_zero")

# Independent stageR computation on a table's own p-values, with simple keys.
stager_oracle <- function(pacs, alpha) {
  genes <- unique(pacs[, c("gene_id", "gene_pvalue")])
  check(!anyDuplicated(genes$gene_id), "Gene p-values differ within a gene.")
  gene_keys <- paste0("g", seq_len(nrow(genes)))
  screen <- as_number(genes$gene_pvalue)
  screen[!is.finite(screen)] <- NA_real_
  names(screen) <- gene_keys
  if (!any(is.finite(screen))) return(rep(NA_real_, nrow(pacs)))
  feature_keys <- paste0("f", seq_len(nrow(pacs)))
  confirmation <- as_number(pacs$pac_pvalue)
  missing <- !is.finite(confirmation)
  confirmation[missing] <- 1
  confirmation <- matrix(confirmation, ncol = 1L, dimnames = list(feature_keys, "contrast"))
  tx2gene <- data.frame(
    txID = feature_keys,
    geneID = gene_keys[match(pacs$gene_id, genes$gene_id)]
  )
  object <- suppressMessages(stageR::stageRTx(
    pScreen = screen,
    pConfirmation = confirmation,
    pScreenAdjusted = FALSE,
    tx2gene = tx2gene
  ))
  object <- suppressMessages(stageR::stageWiseAdjustment(
    object, method = "dtu", alpha = alpha, allowNA = TRUE
  ))
  adjusted <- suppressMessages(stageR::getAdjustedPValues(
    object, order = FALSE, onlySignificantGenes = FALSE
  ))
  result <- adjusted$transcript[match(feature_keys, adjusted$txID)]
  result[missing] <- NA_real_
  result
}

test_case("M-01", "pac_fdr is finite for screened genes and matches stageR", {
  run <- require_run(run_a)
  screened <- list(
    T1_vs_C = c(
      sprintf("shift_T1_%02d", 1:6), "shift_both_01", "shift_both_02",
      "gene:colon", "gene_gain", "gene_unrelated_zero"
    ),
    T2_vs_C = c("shift_both_01", "shift_both_02", "gene_unrelated_zero")
  )
  for (comparison in names(screened)) {
    pacs <- read_result(run$final_directory, paste0(comparison, ".pacs.tsv.gz"))
    rows <- pacs[pacs$gene_id %in% screened[[comparison]], , drop = FALSE]
    check(
      setequal(unique(rows$gene_id), screened[[comparison]]),
      comparison, ": screened genes are missing from the table."
    )
    fdr <- as_number(rows$pac_fdr)
    check(
      all(is.finite(fdr)),
      comparison, ": ", sum(!is.finite(fdr)), " of ", length(fdr),
      " PACs of screened genes lack pac_fdr."
    )
    expected <- stager_oracle(pacs, alpha = 0.05)
    observed <- as_number(pacs$pac_fdr)
    check(
      identical(is.na(observed), is.na(expected)),
      comparison, ": the NA pattern of pac_fdr differs from stageR."
    )
    # stageR's DTU procedure confirms both PACs of a significant two-PAC gene
    # and reports 0 for them, so pac_fdr can be below the raw p-value.
    finite <- !is.na(observed)
    check(
      max(abs(observed[finite] - expected[finite])) < 1e-10,
      comparison, ": pac_fdr differs from stageR."
    )
  }
})

test_case("M-02", "a PAC silent in the control is called gained", {
  run <- require_run(run_a)
  pacs <- read_result(run$final_directory, "T1_vs_C.pacs.tsv.gz")
  row <- one_row(pacs, pac_id = gain_pac3)
  check(identical(row$event_type, "gained"), "event_type is ", row$event_type, ".")
  check(as_number(row$gene_fdr) <= 0.05, "gene_fdr is ", row$gene_fdr, ".")
  check(as_number(row$pac_fdr) <= 0.05, "pac_fdr is ", row$pac_fdr, ".")
  check(as_number(row$fitted_control_pau) <= 0.01, "fitted control PAU is ", row$fitted_control_pau)
  check(as_number(row$fitted_treatment_pau) >= 0.25, "fitted treatment PAU is ", row$fitted_treatment_pau)
  check(as_number(row$delta_pau) >= 0.25, "delta PAU is ", row$delta_pau, ".")
  check(!as_flag(row$zero_boundary_unstable), "the PAC is flagged unstable.")
  check(
    identical(row$model_status, "fitted_with_zero_count_stabilization"),
    "model_status is ", row$model_status
  )
  check(as_number(row$stabilization_successes) == 5, "stabilization successes: ", row$stabilization_successes)
  check(identical(row$raw_control_counts, "C_1=0,C_2=0,C_3=0"), "raw control counts changed.")
  check(identical(row$observed_control_pau, "C_1=0,C_2=0,C_3=0"), "observed control PAU changed.")
  calls <- read_result(run$final_directory, "T1_vs_C.calls.tsv.gz")
  check(
    sum(calls$pac_id == gain_pac3 & calls$event_type == "gained") == 1L,
    "the calls table does not contain exactly one gained row."
  )
})

test_case("M-03", "a PAC silent only in another treatment leaves this comparison testable", {
  run <- require_run(run_a)
  t1 <- read_result(run$final_directory, "T1_vs_C.pacs.tsv.gz")
  rows <- t1[t1$gene_id == "gene_unrelated_zero", , drop = FALSE]
  check(nrow(rows) == 3L, "expected 3 rows, found ", nrow(rows), ".")
  for (column in c("gene_pvalue", "pac_pvalue", "fitted_control_pau", "fitted_treatment_pau")) {
    check(all(is.finite(as_number(rows[[column]]))), "non-finite ", column, " in T1_vs_C.")
  }
  check(!any(as_flag(rows$zero_boundary_unstable)), "rows are flagged unstable in T1_vs_C.")
  check(all(as_number(rows$gene_fdr) <= 0.05), "gene_fdr exceeds 0.05 in T1_vs_C.")
  increased <- one_row(t1, pac_id = unrelated_pac2)
  check(identical(increased$event_type, "increased_usage"), "PAC2 is ", increased$event_type, ".")
  check(as_number(increased$pac_fdr) <= 0.05, "PAC2 pac_fdr is ", increased$pac_fdr, ".")
  decreased <- one_row(t1, pac_id = unrelated_pac1)
  check(identical(decreased$event_type, "decreased_usage"), "PAC1 is ", decreased$event_type, ".")
  check(identical(one_row(t1, pac_id = unrelated_pac3)$event_type, "none"), "PAC3 is not none.")
  t2 <- read_result(run$final_directory, "T2_vs_C.pacs.tsv.gz")
  lost <- one_row(t2, pac_id = unrelated_pac3)
  check(identical(lost$event_type, "lost"), "PAC3 in T2_vs_C is ", lost$event_type, ".")
  gain_rows <- t2[t2$gene_id == "gene_gain", , drop = FALSE]
  check(all(is.finite(as_number(gain_rows$gene_pvalue))), "gene_gain lacks a gene p in T2_vs_C.")
  labelled <- c(
    "gained", "lost", "increased_usage", "decreased_usage", "gained_candidate", "lost_candidate"
  )
  check(!any(gain_rows$event_type %in% labelled), "gene_gain has an event in T2_vs_C.")
  check(
    !as_flag(one_row(t2, pac_id = gain_pac3)$zero_boundary_unstable),
    "a PAC silent in both groups is flagged unstable."
  )
})

test_case("M-04", "stabilized genes get precision, proportions, and gene p-values", {
  run <- require_run(run_a)
  for (comparison in c("T1_vs_C", "T2_vs_C")) {
    pacs <- read_result(run$final_directory, paste0(comparison, ".pacs.tsv.gz"))
    for (gene in boundary_genes) {
      rows <- pacs[pacs$gene_id == gene, , drop = FALSE]
      check(nrow(rows) == 3L, comparison, "/", gene, ": expected 3 rows.")
      precision <- as_number(rows$precision)
      check(all(is.finite(precision) & precision > 0), comparison, "/", gene, ": no precision.")
      stable <- !as_flag(rows$zero_boundary_unstable)
      for (side in c("control", "treatment")) {
        pau <- as_number(rows[[paste0("fitted_", side, "_pau")]])
        alpha <- as_number(rows[[paste0("alpha_", side)]])
        check(abs(sum(pau) - 1) < 1e-6, comparison, "/", gene, ": ", side, " PAU does not sum to 1.")
        expected <- pau[stable] * precision[stable]
        check(
          all(abs(alpha[stable] - expected) <= 1e-8 * pmax(1, abs(expected))),
          comparison, "/", gene, ": alpha is not PAU times precision."
        )
      }
    }
    genes <- read_result(run$final_directory, paste0(comparison, ".genes.tsv.gz"))
    check(
      all(is.finite(as_number(genes$gene_pvalue[genes$gene_id %in% boundary_genes]))),
      comparison, ": the genes table lacks p-values for stabilized genes."
    )
  }
  precision <- read_result(run$final_directory, "gene_precision.tsv.gz")
  values <- as_number(precision$precision[precision$gene_id %in% boundary_genes])
  check(length(values) == 2L && all(is.finite(values)), "gene_precision lacks stabilized genes.")
  fitted <- read_result(run$final_directory, "fitted_pau.tsv.gz")
  rows <- fitted[fitted$gene_id %in% boundary_genes, , drop = FALSE]
  check(nrow(rows) == 12L, "fitted_pau has ", nrow(rows), " stabilized rows; expected 12.")
  check(all(is.finite(as_number(rows$fitted_control_pau))), "fitted_pau has missing values.")
})

test_case("M-05", "bootstrap intervals exist for gained, lost, and shifted PACs", {
  for (run in list(require_run(run_a), require_run(run_b))) {
    t1 <- read_result(run$final_directory, "T1_vs_C.pacs.tsv.gz")
    t2 <- read_result(run$final_directory, "T2_vs_C.pacs.tsv.gz")
    gained <- one_row(t1, pac_id = gain_pac3)
    low <- as_number(gained$delta_pau_ci_low)
    high <- as_number(gained$delta_pau_ci_high)
    delta <- as_number(gained$delta_pau)
    check(is.finite(low) && is.finite(high), "the gained PAC has no interval.")
    check(
      as_number(gained$bootstrap_successes) >= ceiling(0.8 * replicates),
      "the gained PAC has ", gained$bootstrap_successes, " bootstrap successes."
    )
    check(low > 0 && low <= delta && delta <= high, "the gained interval does not contain delta.")
    targets <- rbind(
      one_row(t1, pac_id = unrelated_pac2),
      one_row(t2, pac_id = unrelated_pac3),
      t1[t1$gene_id %in% c("shift_both_01", "shift_both_02"), , drop = FALSE],
      t2[t2$gene_id %in% c("shift_both_01", "shift_both_02"), , drop = FALSE]
    )
    check(nrow(targets) == 14L, "expected 14 interval rows, found ", nrow(targets), ".")
    check(
      all(is.finite(as_number(targets$delta_pau_ci_low)) &
        is.finite(as_number(targets$delta_pau_ci_high))),
      "some selected PACs lack intervals."
    )
  }
})

test_case("M-06", "results do not depend on workers, batch size, or the ambient seed", {
  a <- require_run(run_a)
  b <- require_run(run_b)
  check(length(a$batches) >= 3L, "run A wrote ", length(a$batches), " batches.")
  check(length(b$batches) == 1L, "run B wrote ", length(b$batches), " batches; expected 1.")
  tables_a <- final_tables(a$final_directory)
  tables_b <- final_tables(b$final_directory)
  check(identical(names(tables_a), names(tables_b)), "the runs wrote different tables.")
  for (name in names(tables_a)) {
    check(identical(tables_a[[name]], tables_b[[name]]), name, " differs between runs.")
  }
})

test_case("M-07", "the fit calls DRIMSeq once plus once per stabilization repeat", {
  calls <- require_run(run_a)$fit_calls
  check(calls[["dmPrecision"]] == 6L, "dmPrecision was called ", calls[["dmPrecision"]], " times.")
  check(calls[["dmFit"]] == 6L, "dmFit was called ", calls[["dmFit"]], " times.")
})

test_case("M-08", "the bootstrap keeps precision fixed and refits each batch once per draw", {
  run <- require_run(run_a)
  calls <- run$bootstrap_calls
  batches <- length(run$batches)
  check(calls[["dmPrecision"]] == 0L, "dmPrecision was called ", calls[["dmPrecision"]], " times.")
  check(
    calls[["dmFit"]] >= replicates * batches && calls[["dmFit"]] <= 2L * replicates * batches,
    "dmFit was called ", calls[["dmFit"]], " times for ", batches, " batches."
  )
})

test_case("M-09", "sourcing the script does not run the command line", {
  environment <- new.env()
  outcome <- tryCatch(
    {
      source(script, local = environment)
      "ok"
    },
    error = function(error) conditionMessage(error)
  )
  check(identical(outcome, "ok"), "source() failed: ", outcome)
  check(exists("run_fit_mode", envir = environment, inherits = FALSE), "run_fit_mode is missing.")
})

test_case("M-10", "bootstrap batches default to 500 genes", {
  check(identical(model$parse_batch_size(NULL), 500L), "default batch size is not 500.")
})

test_case("H-01", "bootstrap totals must be finite non-negative integers", {
  check(identical(model$bootstrap_total(c(1, 2), "gene-1", "sample-1"), 3L), "total of 1 and 2.")
  check(identical(model$bootstrap_total(c(0, 0), "gene-1", "sample-1"), 0L), "total of zeros.")
  for (values in list(c(1, NA_real_), c(0.5, 0.6), c(.Machine$integer.max, 1))) {
    message_text <- tryCatch(
      {
        model$bootstrap_total(values, "gene-1", "sample-1")
        ""
      },
      error = conditionMessage
    )
    check(
      grepl("Invalid bootstrap total for gene gene-1, sample sample-1", message_text, fixed = TRUE),
      "invalid total accepted: ", paste(values, collapse = ",")
    )
  }
})

test_case("H-02", "bootstrap gene selection follows gene FDR and candidates", {
  check(identical(
    model$bootstrap_gene_ids(
      c("gene-1", "gene-2", NA_character_, "gene-4"),
      c(0.01, NA_real_, 0.01, 0.50),
      c(FALSE, NA, FALSE, TRUE),
      0.05
    ),
    c("gene-1", "gene-4")
  ), "default selection.")
  check(identical(
    model$bootstrap_gene_ids(
      c("gene-1", "gene-2", "gene-3"), c(0.01, NA_real_, 0.50), c(FALSE, TRUE, TRUE), 0.05,
      include_candidates = FALSE
    ),
    "gene-1"
  ), "selection without candidates.")
})

test_case("H-03", "parallel bootstrap application matches serial application", {
  probe <- function(repeat_number) {
    set.seed(1000L + repeat_number)
    c(repeat_number, stats::runif(1))
  }
  check(
    identical(model$bootstrap_apply(1:4, 1L, probe), model$bootstrap_apply(1:4, 2L, probe)),
    "serial and parallel results differ."
  )
})

test_case("H-04", "dominant PAC ignores missing values", {
  check(identical(model$dominant_pac(c("p1", "p2"), c(NA_real_, NA_real_)), NA_character_), "all NA.")
  check(identical(model$dominant_pac(c("p1", "p2"), c(NA_real_, 0.25)), "p2"), "one NA.")
})

test_case("H-05", "event classification", {
  params <- list(
    event_min_supporting_samples = 2, min_abs_delta_pau = 0.10, gene_fdr = 0.05,
    site_fdr = 0.05, event_max_control_pau = 0.01, event_min_treatment_pau = 0.05,
    min_gene_total = 20, potential_internal_priming_withheld_calls = TRUE
  )
  gained <- list(
    control_supporting_samples = 0, treatment_supporting_samples = 2, delta_pau = 0.20,
    gene_fdr = 0.01, pac_fdr = 0.01, zero_boundary_unstable = FALSE, confidence = "high",
    internal_priming_flag = FALSE, exploratory_insufficient_replicates = FALSE,
    fitted_control_pau = 0, fitted_treatment_pau = 0.20,
    control_gene_total = 40, treatment_gene_total = 40
  )
  modified <- function(row, ...) utils::modifyList(row, list(...))
  lost <- modified(gained,
    control_supporting_samples = 2, treatment_supporting_samples = 0, delta_pau = -0.20,
    fitted_control_pau = 0.20, fitted_treatment_pau = 0
  )
  expectations <- list(
    list(gained, "gained"),
    list(modified(gained, pac_fdr = NA_real_), "gained_candidate"),
    list(modified(gained, confidence = "low"), "gained_candidate"),
    list(modified(gained, internal_priming_flag = TRUE), "gained_candidate"),
    list(modified(gained, zero_boundary_unstable = TRUE), "gained_candidate"),
    list(modified(gained, fitted_control_pau = NA_real_), "none"),
    list(modified(gained, control_supporting_samples = NA_real_), "none"),
    # Absence from the control counts only when the control covers the gene.
    list(modified(gained, control_gene_total = 20), "gained"),
    list(modified(gained, control_gene_total = 19), "gained_candidate"),
    list(modified(gained, control_gene_total = NA_real_), "gained_candidate"),
    list(modified(gained, treatment_gene_total = 5), "gained"),
    list(lost, "lost"),
    list(modified(lost, treatment_gene_total = 20), "lost"),
    list(modified(lost, treatment_gene_total = 19), "lost_candidate"),
    list(modified(lost, control_gene_total = 5), "lost"),
    list(modified(gained, control_supporting_samples = 2, fitted_control_pau = 0.30,
      fitted_treatment_pau = 0.50), "increased_usage"),
    list(modified(gained, control_supporting_samples = 2, delta_pau = -0.20,
      fitted_control_pau = 0.50, fitted_treatment_pau = 0.30), "decreased_usage")
  )
  for (index in seq_along(expectations)) {
    observed <- model$classify_event(expectations[[index]][[1]], params)
    check(
      identical(observed, expectations[[index]][[2]]),
      "case ", index, " gave ", observed, "; expected ", expectations[[index]][[2]], "."
    )
  }
})

test_case("H-07", "gene-level events need a gene that passes the screen, and stay off PAC rows", {
  params <- list(
    event_min_supporting_samples = 2, min_abs_delta_pau = 0.10, gene_fdr = 0.05,
    site_fdr = 0.05, event_max_control_pau = 0.01, event_min_treatment_pau = 0.05,
    active_pac_min_pau = 0.05, min_gene_total = 20,
    potential_internal_priming_withheld_calls = TRUE
  )
  # Every PAC is detected in both groups with a PAC FDR of 0.5, so no PAC has
  # a call and only the gene-level events can apply.
  gene <- function(gene_id, gene_fdr, control, treatment) {
    data.frame(
      feature_id = paste0(gene_id, "_p", seq_along(control)),
      gene_id = gene_id,
      fitted_control_pau = control,
      fitted_treatment_pau = treatment,
      delta_pau = treatment - control,
      gene_fdr = gene_fdr,
      pac_fdr = 0.5,
      control_supporting_samples = 2,
      treatment_supporting_samples = 2,
      control_gene_total = 100,
      treatment_gene_total = 100,
      zero_boundary_unstable = FALSE,
      confidence = "moderate",
      internal_priming_flag = FALSE,
      exploratory_insufficient_replicates = FALSE,
      stringsAsFactors = FALSE
    )
  }
  switched <- list(c(0.7, 0.3), c(0.3, 0.7))
  more_active <- list(c(1.0, 0.0), c(0.8, 0.2))
  pacs <- rbind(
    gene("switch_screened", 0.01, switched[[1]], switched[[2]]),
    gene("switch_unscreened", 0.5, switched[[1]], switched[[2]]),
    gene("switch_untested", NA_real_, switched[[1]], switched[[2]]),
    gene("more_screened", 0.01, more_active[[1]], more_active[[2]]),
    gene("more_unscreened", 0.5, more_active[[1]], more_active[[2]])
  )
  labelled <- model$assign_events(pacs, params)
  check(
    all(labelled$event_type == "none"),
    "PAC rows carry ", paste(unique(labelled$event_type), collapse = ", ")
  )
  genes <- model$gene_level_events(labelled, params)
  expected <- list(
    gene_id = c(
      "switch_screened", "switch_unscreened", "switch_untested", "more_screened",
      "more_unscreened"
    ),
    dominant_switch = c(TRUE, FALSE, FALSE, FALSE, FALSE),
    active_pacs_change = c("none", "none", "none", "more", "none")
  )
  check(
    identical(as.list(genes), expected),
    "gene events were ",
    paste(genes$gene_id, genes$dominant_switch, genes$active_pacs_change, collapse = "; ")
  )
})

test_case("H-08", "active_pac_min_pau sets the active-PAC count, not gained or lost calls", {
  params <- list(
    event_min_supporting_samples = 2, min_abs_delta_pau = 0.10, gene_fdr = 0.05,
    site_fdr = 0.05, event_max_control_pau = 0.01, event_min_treatment_pau = 0.05,
    active_pac_min_pau = 0.05, min_gene_total = 20,
    potential_internal_priming_withheld_calls = TRUE
  )
  # p1 falls, p2 barely moves, and p3, absent from the control, is gained
  # with 0.12 of the gene in the treatment.
  pacs <- data.frame(
    feature_id = c("p1", "p2", "p3"),
    gene_id = "g",
    fitted_control_pau = c(0.94, 0.06, 0.00),
    fitted_treatment_pau = c(0.80, 0.08, 0.12),
    delta_pau = c(-0.14, 0.02, 0.12),
    gene_fdr = 0.01,
    pac_fdr = c(0.001, 0.5, 0.001),
    control_supporting_samples = c(2, 2, 0),
    treatment_supporting_samples = 2,
    control_gene_total = 100,
    treatment_gene_total = 100,
    zero_boundary_unstable = FALSE,
    confidence = "moderate",
    internal_priming_flag = FALSE,
    exploratory_insufficient_replicates = FALSE,
    stringsAsFactors = FALSE
  )
  outcome <- function(changes) {
    labelled <- model$assign_events(pacs, utils::modifyList(params, changes))
    genes <- model$gene_level_events(labelled, utils::modifyList(params, changes))
    list(
      events = labelled$event_type,
      active = c(labelled$control_active_pacs[[1]], labelled$treatment_active_pacs[[1]]),
      change = genes$active_pacs_change
    )
  }
  calls <- c("decreased_usage", "none", "gained")
  describe <- function(result) {
    paste0(
      paste(result$events, collapse = ","), " active ", paste(result$active, collapse = "/"),
      " ", result$change
    )
  }
  default <- outcome(list())
  check(
    identical(default, list(events = calls, active = c(2L, 3L), change = "more")),
    "defaults gave ", describe(default)
  )
  # Only PACs with at least 0.13 of the gene count as active; the calls stay.
  raised <- outcome(list(active_pac_min_pau = 0.13))
  check(
    identical(raised, list(events = calls, active = c(1L, 1L), change = "none")),
    "active_pac_min_pau = 0.13 gave ", describe(raised)
  )
  # A gained call now needs 0.13 in the treatment; the active count stays.
  stricter <- outcome(list(event_min_treatment_pau = 0.13))
  check(
    identical(
      stricter,
      list(events = c("decreased_usage", "none", "none"), active = c(2L, 3L), change = "more")
    ),
    "event_min_treatment_pau = 0.13 gave ", describe(stricter)
  )
})

test_case("H-06", "bootstrap interval files are validated", {
  expect_error_text <- function(code, pattern) {
    message_text <- tryCatch(
      {
        force(code)
        ""
      },
      error = conditionMessage
    )
    check(grepl(pattern, message_text, fixed = TRUE), "expected an error containing '", pattern, "'.")
  }
  expect_error_text(model$read_bootstrap_intervals(character()), "No bootstrap interval files")
  missing_path <- tempfile(fileext = ".tsv")
  utils::write.table(
    data.frame(comparison = "a_vs_b", gene_id = "g1"), missing_path,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  expect_error_text(model$read_bootstrap_intervals(missing_path), "missing required columns")
  duplicate_path <- tempfile(fileext = ".tsv")
  row <- data.frame(
    comparison = "a_vs_b", gene_id = "g1", feature_id = "p1",
    delta_pau_ci_low = NA_real_, delta_pau_ci_high = NA_real_, bootstrap_successes = 15L,
    bootstrap_perturbed = 0L, bootstrap_status = "insufficient_successes"
  )
  utils::write.table(
    rbind(row, row), duplicate_path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
  expect_error_text(model$read_bootstrap_intervals(duplicate_path), "duplicate comparison/gene/PAC rows")
  null_batch <- tempfile(fileext = ".rds")
  saveRDS(NULL, null_batch)
  expect_error_text(
    model$run_bootstrap_batch(null_batch, tempfile(fileext = ".tsv.gz"), 1L),
    "Bootstrap batch payload must be a list"
  )
})

expect_error_containing <- function(code, pattern) {
  message_text <- tryCatch(
    {
      force(code)
      ""
    },
    error = conditionMessage
  )
  check(grepl(pattern, message_text, fixed = TRUE), "expected an error containing '", pattern, "'.")
}

test_case("M-11", "PACs are ordered by genomic coordinate within each gene", {
  counts <- data.frame(
    gene_id = c("g1", "g1", "g2", "g2"),
    pac_id = c("PACv1.t.chr1.+.1000", "PACv1.t.chr1.+.200", "PACv1.t.chr1.-.90", "PACv1.t.chr1.-.80"),
    s1 = 1:4,
    stringsAsFactors = FALSE
  )
  atlas <- data.frame(pac_id = counts$pac_id, coordinate = c(1000, 200, 90, 80))
  ordered <- model$order_family_counts(counts, atlas)
  check(
    identical(ordered$pac_id, c("PACv1.t.chr1.+.200", "PACv1.t.chr1.+.1000", "PACv1.t.chr1.-.80", "PACv1.t.chr1.-.90")),
    "unexpected order: ", paste(ordered$pac_id, collapse = ", ")
  )
  expect_error_containing(model$order_family_counts(counts, atlas[-1, ]), "lacks coordinates")
  expect_error_containing(
    model$order_family_counts(counts[-1, ], atlas),
    "fewer than 2 testable PACs"
  )
})

small_layout <- function() {
  samples <- data.frame(
    sample_id = c("a1", "a2", "b1", "b2"),
    condition = c("A", "A", "B", "B"),
    control_condition = "A",
    batch = c("x", "y", "x", "y"),
    stringsAsFactors = FALSE
  )
  model$family_layout(
    "A", samples,
    list(model_covariates = list("batch"), min_replicates_per_condition = 2)
  )
}

test_case("M-12", "zero groups are found across condition and covariate levels", {
  layout <- small_layout()
  check(is.factor(layout$samples$batch), "covariates are not categorical.")
  counts <- data.frame(
    gene_id = c("condition_zero", "condition_zero", "batch_zero", "batch_zero", "none", "none"),
    pac_id = paste0("p", 1:6),
    a1 = c(0, 5, 0, 5, 1, 5), a2 = c(0, 5, 3, 5, 2, 5),
    b1 = c(4, 5, 0, 5, 0, 5), b2 = c(4, 5, 4, 5, 3, 5),
    stringsAsFactors = FALSE
  )
  check(
    identical(sort(model$zero_group_genes(counts, layout)), c("batch_zero", "condition_zero")),
    "unexpected zero-group genes: ", paste(model$zero_group_genes(counts, layout), collapse = ", ")
  )
})

test_case("M-13", "zero perturbation changes only zeros, within (0, 0.1), reproducibly", {
  values <- matrix(c(0, 3, 0, 7, 0, 1), nrow = 2)
  first <- model$perturb_zero_cells(values, 42L)
  again <- model$perturb_zero_cells(values, 42L)
  other <- model$perturb_zero_cells(values, 43L)
  zeros <- values == 0
  check(identical(first[!zeros], values[!zeros]), "non-zero counts changed.")
  check(all(first[zeros] > 0 & first[zeros] < 0.1), "perturbed values are outside (0, 0.1).")
  check(identical(first, again), "the same seed gave different values.")
  check(!identical(first, other), "different seeds gave identical values.")
})

test_case("M-14", "event classification reads typed flags in a mixed table", {
  params <- list(
    event_min_supporting_samples = 2, min_abs_delta_pau = 0.10, gene_fdr = 0.05,
    site_fdr = 0.05, event_max_control_pau = 0.01, event_min_treatment_pau = 0.05,
    min_gene_total = 20, potential_internal_priming_withheld_calls = TRUE
  )
  table <- data.frame(
    control_supporting_samples = c(0, 0, 0, 0),
    treatment_supporting_samples = c(3, 3, 3, 3),
    delta_pau = c(0.3, 0.3, 0.3, 0.3),
    gene_fdr = c(0.001, 0.001, 0.001, 0.001),
    pac_fdr = c(0.001, 0.001, 0.001, 0.001),
    zero_boundary_unstable = c(TRUE, FALSE, FALSE, FALSE),
    confidence = c("high", "high", "high", "high"),
    internal_priming_flag = c(FALSE, TRUE, FALSE, FALSE),
    exploratory_insufficient_replicates = c(FALSE, FALSE, FALSE, FALSE),
    fitted_control_pau = c(0, 0, 0, 0),
    fitted_treatment_pau = c(0.3, 0.3, 0.3, 0.3),
    control_gene_total = c(40, 40, 40, 10),
    treatment_gene_total = c(40, 40, 40, 40),
    stringsAsFactors = FALSE
  )
  events <- model$classify_table_events(table, params)
  check(
    identical(events, c("gained_candidate", "gained_candidate", "gained", "gained_candidate")),
    "events were ", paste(events, collapse = ", ")
  )
})

test_case("M-15", "bootstrap batches are validated before use", {
  run <- require_run(run_a)
  batch <- readRDS(run$batches[[1]])
  model$validate_bootstrap_batch(batch)
  # Replace fields directly: modifyList would merge data frames column-wise.
  with_field <- function(name, value) {
    changed <- batch
    changed[[name]] <- value
    changed
  }
  broken <- list(
    list(with_field("counts", batch$counts[, -3L, drop = FALSE]), "counts columns that do not match"),
    list(with_field("selections", data.frame(
      gene_id = "unknown", comparison = batch$comparisons$comparison[[1]]
    )), "invalid gene selections"),
    list(with_field("genes", batch$genes[0, , drop = FALSE]), "lists genes that do not match")
  )
  for (case in broken) expect_error_containing(model$validate_bootstrap_batch(case[[1]]), case[[2]])
})

test_case("M-19", "bootstrap seeds are distinct across replicates of a gene", {
  seeds <- vapply(1:1000, function(replicate_number) {
    model$stable_seed(1729, "0123456789abcdef0123456789abcdef", "DMSO", "gene-1", "bootstrap", replicate_number)
  }, integer(1))
  check(length(unique(seeds)) == 1000L, "only ", length(unique(seeds)), " distinct seeds for 1000 replicates.")
  check(seeds[[120]] != seeds[[201]], "replicates 120 and 201 share a seed.")
})

test_case("M-16", "finalize stops when a selected PAC receives no interval", {
  run <- require_run(run_a)
  empty_intervals <- tempfile(fileext = ".tsv.gz")
  connection <- gzfile(empty_intervals, "wt")
  writeLines(paste(model$INTERVAL_COLUMNS, collapse = "\t"), connection)
  close(connection)
  expect_error_containing(
    model$run_finalize_mode(list(
      preliminary_dir = file.path(run$fit_directory, "preliminary"),
      intervals = empty_intervals,
      output_dir = tempfile("final-")
    )),
    "received no bootstrap intervals"
  )
})

test_case("M-20", "a stabilized gene whose first repeat failed keeps its df and p-values", {
  sample_ids <- c("C_1", "C_2", "T1_1", "T1_2")
  layout <- list(
    family = "C", control = "C", sample_ids = sample_ids,
    groups = list(C = c("C_1", "C_2"), T1 = c("T1_1", "T1_2")),
    comparisons = data.frame(comparison = "T1_vs_C", treatment = "T1", stringsAsFactors = FALSE)
  )
  genes <- c("gA", "gB")
  features <- data.frame(
    gene_id = c("gA", "gA", "gB", "gB"), feature_id = c("a1", "a2", "b1", "b2"),
    stringsAsFactors = FALSE
  )
  repeat_fit <- function(failed_gene = NULL, shift = 0) {
    proportions <- rbind(c(0.3, 0.3, 0.6, 0.6), c(0.7, 0.7, 0.4, 0.4), c(0.5, 0.5, 0.2, 0.2), c(0.5, 0.5, 0.8, 0.8))
    dimnames(proportions) <- list(features$feature_id, sample_ids)
    precision <- c(gA = 100, gB = 120)
    gene_tests <- data.frame(gene_id = genes, lr = c(10, 12) + shift, df = 1, pvalue = 0.001)
    feature_tests <- data.frame(
      gene_id = features$gene_id, feature_id = features$feature_id,
      lr = c(10, 10, 12, 12) + shift, df = 1, pvalue = 0.001
    )
    if (!is.null(failed_gene)) {
      rows <- features$gene_id == failed_gene
      proportions[rows, ] <- NA_real_
      precision[[failed_gene]] <- NA_real_
      gene_tests[gene_tests$gene_id == failed_gene, c("lr", "df", "pvalue")] <- NA_real_
      feature_tests[rows, c("lr", "df", "pvalue")] <- NA_real_
    }
    list(
      gene_ids = genes, features = features, proportions = proportions, precision = precision,
      omnibus = gene_tests,
      contrasts = list(T1_vs_C = list(genes = gene_tests, features = feature_tests))
    )
  }
  # Repeat 1 fails for gA only; the remaining four succeed.
  repeats <- list(
    repeat_fit("gA"), repeat_fit(shift = 0.1), repeat_fit(shift = -0.1), repeat_fit(), repeat_fit(shift = 0.2)
  )
  summary <- model$summarize_zero_sensitivity(
    repeats, genes, layout,
    list(dm_zero_sensitivity_repeats = 5L, dm_zero_max_delta_pau_spread = 0.02)
  )
  omnibus <- summary$omnibus[summary$omnibus$gene_id == "gA", , drop = FALSE]
  check(omnibus$stabilization_successes == 4L, "gA has ", omnibus$stabilization_successes, " successes.")
  check(is.finite(omnibus$df) && is.finite(omnibus$pvalue), "gA lost its gene-level df or p-value.")
  pacs <- summary$contrasts$T1_vs_C$features
  pacs <- pacs[pacs$gene_id == "gA", , drop = FALSE]
  check(all(is.finite(pacs$df) & is.finite(pacs$pvalue)), "gA lost its PAC-level df or p-values.")
})

# A variant with bootstrap intervals disabled and one gene with no counts in T2.
variant <- build_usage_dataset(0L)
t2_samples <- c("T2_1", "T2_2", "T2_3")
variant$counts[variant$counts$gene_id == "null_01", t2_samples] <- 0L
variant_paths <- write_usage_inputs(variant, file.path(work, "variant_inputs"))
run_variant <- tryCatch(
  run_statistics_in_process(model, variant_paths, file.path(work, "variant")),
  error = function(error) error
)

test_case("M-17", "with zero bootstrap replicates, rows are marked disabled", {
  run <- require_run(run_variant)
  check(length(run$batches) == 1L, "expected one placeholder batch, found ", length(run$batches), ".")
  check(isTRUE(readRDS(run$batches[[1]])$empty), "the only batch is not the empty placeholder.")
  for (comparison in c("T1_vs_C", "T2_vs_C")) {
    pacs <- read_result(run$final_directory, paste0(comparison, ".pacs.tsv.gz"))
    check(all(pacs$bootstrap_status == "disabled"), comparison, ": not every row is disabled.")
  }
})

test_case("M-18", "a condition with no counts leaves only its comparisons untested", {
  run <- require_run(run_variant)
  t2 <- read_result(run$final_directory, "T2_vs_C.pacs.tsv.gz")
  rows <- t2[t2$gene_id == "null_01", , drop = FALSE]
  check(nrow(rows) == 3L, "expected 3 T2 rows for null_01.")
  check(all(rows$model_status == "group_without_counts"), "T2 rows are not group_without_counts.")
  check(all(is.na(as_number(rows$pac_pvalue)) & is.na(as_number(rows$gene_pvalue))), "T2 has p-values.")
  check(all(rows$event_type == "none"), "T2 rows have events: ", paste(rows$event_type, collapse = ", "))
  t1 <- read_result(run$final_directory, "T1_vs_C.pacs.tsv.gz")
  rows <- t1[t1$gene_id == "null_01", , drop = FALSE]
  check(all(is.finite(as_number(rows$gene_pvalue))), "T1_vs_C lost its gene p-value.")
  check(!any(rows$model_status == "group_without_counts"), "T1_vs_C was masked.")
  omnibus <- read_result(run$final_directory, "C.gene_omnibus.tsv.gz")
  check(one_row(omnibus, gene_id = "null_01")$model_status == "group_without_counts", "omnibus not masked.")
  precision <- read_result(run$final_directory, "gene_precision.tsv.gz")
  check(one_row(precision, gene_id = "null_01")$model_status == "group_without_counts", "precision status.")
})

test_case("M-21", "statistical filtering gives every PAC one row: tested or its reasons", {
  counts <- data.frame(
    gene_id = c(
      "gA", "gA", "gA", "gB", "gB", NA, "gC,gD", "gE", "gE", "gF", "gF", "gG", "gG"
    ),
    gene_name = c(
      "GeneA", "GeneA", "GeneA", "GeneB", "GeneB", NA, "GeneC,GeneD", "GeneE", "GeneE",
      "GeneF", "GeneF", "GeneG", "GeneG"
    ),
    pac_id = sprintf("p%02d", 1:13),
    chrom = "chr1",
    start = seq(100L, 1300L, by = 100L),
    end = seq(101L, 1301L, by = 100L),
    strand = "+",
    locus = sprintf("chr1:%d-%d", seq(101L, 1301L, by = 100L), seq(101L, 1301L, by = 100L)),
    s1 = c(30, 10, 0, 3, 4, 50, 50, 10, 10, 15, 1, 3, 1),
    s2 = c(30, 10, 1, 3, 4, 50, 50, 10, 10, 15, 0, 3, 0),
    outside = rep(99, 13),
    stringsAsFactors = FALSE
  )
  params <- list(
    min_site_count = 5L, min_test_supporting_samples = 2L,
    min_site_usage = 0.01, min_gene_total = 20L
  )
  result <- model$family_filter(counts, c("s1", "s2"), params, "F")
  expected_reasons <- c(
    "", "", "site_count<5;supporting_samples<2", "gene_total<20", "gene_total<20",
    "no_gene_assignment", "ambiguous_gene_assignment", "", "",
    "fewer_than_2_testable_pacs", "site_count<5;supporting_samples<2",
    "gene_total<20;fewer_than_2_testable_pacs", "site_count<5;supporting_samples<2"
  )
  reasons <- result$reasons
  check(identical(names(reasons), model$FILTERING_COLUMNS), "columns: ", paste(names(reasons), collapse = ", "))
  check(identical(reasons$pac_id, counts$pac_id), "rows are not in input order.")
  check(
    identical(reasons$gene_name, ifelse(is.na(counts$gene_name), "", counts$gene_name)),
    "gene names were ", paste(reasons$gene_name, collapse = ", ")
  )
  check(all(reasons$family == "F"), "family column is wrong.")
  check(identical(reasons$reason, expected_reasons), "reasons were ", paste(reasons$reason, collapse = " | "))
  check(identical(reasons$tested, expected_reasons == ""), "tested flags do not match the reasons.")
  check(identical(result$counts$pac_id, c("p01", "p02", "p08", "p09")), "tested PACs: ",
    paste(result$counts$pac_id, collapse = ", "))
  check(identical(names(result$counts), c("gene_id", "pac_id", "s1", "s2")), "samples outside the family leaked.")
  check(identical(reasons$locus, counts$locus), "the filtering table lost the PAC locations.")

  run <- require_run(run_a)
  written <- read_result(run$final_directory, "C.statistical_filtering.tsv.gz")
  check(nrow(written) == nrow(dataset$counts), "the family table does not list every PAC once.")
  tested <- unique(read_result(run$final_directory, "T1_vs_C.pacs.tsv.gz")$pac_id)
  check(setequal(written$pac_id[written$tested == "TRUE"], tested), "tested flags differ from the tested PACs.")
})

test_case("M-22", "motif preference uses the family design, covariates included", {
  samples <- data.frame(
    sample_id = c("a1", "a2", "a3", "b1", "b2", "b3"),
    condition = c("A", "A", "A", "B", "B", "B"),
    control_condition = "A",
    batch = c("x", "x", "y", "x", "y", "y"),
    stringsAsFactors = FALSE
  )
  layout <- model$family_layout(
    "A", samples, list(model_covariates = list("batch"), min_replicates_per_condition = 2)
  )
  noise <- c(0.01, -0.02, 0.015, -0.01, 0.02, -0.005)
  canonical <- 0.5 + 0.1 * (samples$condition == "B") + 0.2 * (samples$batch == "y") + noise
  scores <- rbind(
    data.frame(
      sample_id = samples$sample_id, primary_pas_motif_rna = "AAUAAA",
      primary_motif_class = "canonical", motif_usage = sin(canonical)^2,
      transformed_motif_usage = canonical, informative_genes = 60L, stringsAsFactors = FALSE
    ),
    data.frame(
      sample_id = samples$sample_id, primary_pas_motif_rna = "AUUAAA",
      primary_motif_class = "common_variant", motif_usage = 0.2,
      transformed_motif_usage = asin(sqrt(0.2)) + noise,
      informative_genes = c(60L, 60L, 60L, 10L, 60L, 60L), stringsAsFactors = FALSE
    )
  )
  directory <- tempfile("motif-")
  dir.create(directory)
  model$fit_motif_preferences(scores, layout, list(motif_preference_min_genes = 50L), directory, "preference")
  output <- read_result(directory, "B_vs_A.preference.tsv.gz")
  check(identical(names(output), model$MOTIF_PREFERENCE_COLUMNS), "unexpected columns.")
  check(identical(output$primary_pas_motif_rna, "AAUAAA"), "a class below motif_preference_min_genes in one sample was tested.")
  adjusted <- unname(stats::coef(stats::lm(canonical ~ batch + condition, data = samples))[["conditionB"]])
  unadjusted <- unname(stats::coef(stats::lm(canonical ~ condition, data = samples))[["conditionB"]])
  observed <- as_number(output$transformed_coefficient)
  check(abs(observed - adjusted) < 1e-10, "coefficient ", observed, " is not the batch-adjusted ", adjusted, ".")
  check(abs(adjusted - unadjusted) > 1e-3, "the fixture does not distinguish adjusted from unadjusted fits.")
  check(abs(as_number(output$delta_motif_usage) -
    (mean(sin(canonical[4:6])^2) - mean(sin(canonical[1:3])^2))) < 1e-10, "delta motif usage is wrong.")
  check(as_number(output$informative_genes) == 60, "informative genes should be the comparison minimum.")

  # The class-level table uses the same design, with one row per class.
  class_scores <- rbind(
    scores[, setdiff(names(scores), "primary_pas_motif_rna")],
    data.frame(
      sample_id = samples$sample_id, primary_motif_class = "other_variant", motif_usage = 0.1,
      transformed_motif_usage = asin(sqrt(0.1)) - noise, informative_genes = 55L,
      stringsAsFactors = FALSE
    )
  )
  model$fit_motif_preferences(
    class_scores, layout, list(motif_preference_min_genes = 50L), directory, "preference_class",
    keys = "primary_motif_class"
  )
  classes <- read_result(directory, "B_vs_A.preference_class.tsv.gz")
  class_columns <- c(
    "primary_motif_class", "condition", "control_condition", "delta_motif_usage", "fdr",
    "pvalue", "control_mean", "treatment_mean", "transformed_coefficient", "informative_genes"
  )
  check(identical(names(classes), class_columns), "unexpected class columns.")
  check(
    identical(classes$primary_motif_class, c("canonical", "other_variant")),
    "class rows were ", paste(classes$primary_motif_class, collapse = ", "), "."
  )
  class_coefficient <- as_number(classes$transformed_coefficient[[1]])
  check(abs(class_coefficient - adjusted) < 1e-10, "the class coefficient is not batch-adjusted.")
})

test_case("M-23", "versions mode writes a table, then appends to it", {
  path <- tempfile(fileext = ".tsv")
  versions <- c(DRIMSeq = "1.38.0", stageR = "1.2.3")
  model$run_versions_mode(list(output = path), versions)
  model$run_versions_mode(list(output = path), versions)
  table <- utils::read.delim(path, colClasses = "character")
  check(identical(names(table), c("software", "version")), "header: ", paste(names(table), collapse = ", "))
  check(identical(table$software, rep(c("R", "DRIMSeq", "stageR"), 2)), "rows: ", paste(table$software, collapse = ", "))
  check(identical(table$version[2:3], c("1.38.0", "1.2.3")), "package versions were not written.")
})

test_case("M-24", "the genes table records every gene-level event, with or without PAC calls", {
  params <- list(
    event_min_supporting_samples = 2, min_abs_delta_pau = 0.10, gene_fdr = 0.05,
    site_fdr = 0.05, event_max_control_pau = 0.01, event_min_treatment_pau = 0.05,
    active_pac_min_pau = 0.05, min_gene_total = 20,
    potential_internal_priming_withheld_calls = TRUE
  )
  gene <- function(gene_id, gene_fdr, pac_fdr, control, treatment) {
    data.frame(
      feature_id = paste0(gene_id, "_p", seq_along(control)),
      gene_id = gene_id,
      fitted_control_pau = control,
      fitted_treatment_pau = treatment,
      delta_pau = treatment - control,
      gene_fdr = gene_fdr,
      pac_fdr = pac_fdr,
      control_supporting_samples = 2,
      treatment_supporting_samples = 2,
      control_gene_total = 100,
      treatment_gene_total = 100,
      zero_boundary_unstable = FALSE,
      confidence = "moderate",
      internal_priming_flag = FALSE,
      exploratory_insufficient_replicates = FALSE,
      stringsAsFactors = FALSE
    )
  }
  # In switch_called and gain_called both PACs have calls; in loss none does.
  pacs <- rbind(
    gene("switch_called", 0.01, 0.001, c(0.7, 0.3), c(0.3, 0.7)),
    gene("switch_unscreened", 0.5, 0.001, c(0.7, 0.3), c(0.3, 0.7)),
    gene("gain_called", 0.01, 0.001, c(1.0, 0.0), c(0.8, 0.2)),
    gene("loss", 0.01, 0.5, c(0.5, 0.3, 0.2), c(0.6, 0.38, 0.02))
  )
  labelled <- model$assign_events(pacs, params)
  events <- stats::setNames(labelled$event_type, labelled$feature_id)
  expected_events <- c(
    switch_called_p1 = "decreased_usage", switch_called_p2 = "increased_usage",
    switch_unscreened_p1 = "none", switch_unscreened_p2 = "none",
    gain_called_p1 = "decreased_usage", gain_called_p2 = "increased_usage",
    loss_p1 = "none", loss_p2 = "none", loss_p3 = "none"
  )
  check(
    identical(events, expected_events),
    "events were ", paste(names(events), events, sep = "=", collapse = ", ")
  )
  genes <- model$gene_level_events(labelled, params)
  expected <- list(
    gene_id = c("switch_called", "switch_unscreened", "gain_called", "loss"),
    dominant_switch = c(TRUE, FALSE, FALSE, FALSE),
    active_pacs_change = c("none", "none", "more", "fewer")
  )
  check(
    identical(as.list(genes), expected),
    "gene events were ",
    paste(genes$gene_id, genes$dominant_switch, genes$active_pacs_change, collapse = "; ")
  )

  # The fitted tables agree: each gene's flags follow from its PAC rows, and
  # the PAC rows carry only PAC calls.
  run <- require_run(run_a)
  for (comparison in c("T1_vs_C", "T2_vs_C")) {
    genes <- read_result(run$final_directory, paste0(comparison, ".genes.tsv.gz"))
    check(
      identical(names(genes), model$GENE_COLUMNS),
      comparison, ": genes columns were ", paste(names(genes), collapse = ", ")
    )
    pacs <- read_result(run$final_directory, paste0(comparison, ".pacs.tsv.gz"))
    first <- pacs[!duplicated(pacs$gene_id), , drop = FALSE]
    screened <- !is.na(as_number(first$gene_fdr)) & as_number(first$gene_fdr) <= 0.05
    fitted <- tapply(
      is.finite(as_number(pacs$fitted_control_pau)) & is.finite(as_number(pacs$fitted_treatment_pau)),
      pacs$gene_id,
      all
    )[first$gene_id]
    switched <- screened & !is.na(first$dominant_pac_control) &
      !is.na(first$dominant_pac_treatment) &
      first$dominant_pac_control != first$dominant_pac_treatment
    change <- as_number(first$treatment_active_pacs) - as_number(first$control_active_pacs)
    change[!screened | !fitted] <- 0
    index <- match(first$gene_id, genes$gene_id)
    check(!anyNA(index) && nrow(genes) == nrow(first), comparison, ": genes and PACs differ.")
    check(
      identical(as_flag(genes$dominant_switch[index]), unname(switched)),
      comparison, ": dominant_switch disagrees with the PAC rows."
    )
    check(
      identical(
        genes$active_pacs_change[index],
        ifelse(change > 0, "more", ifelse(change < 0, "fewer", "none"))
      ),
      comparison, ": active_pacs_change disagrees with the PAC rows."
    )
    pac_calls <- c(
      "gained", "lost", "increased_usage", "decreased_usage", "gained_candidate",
      "lost_candidate", "none"
    )
    check(
      all(pacs$event_type %in% pac_calls),
      comparison, ": event_type holds ", paste(setdiff(pacs$event_type, pac_calls), collapse = ", ")
    )
  }
})

test_case("M-25", "APA patterns follow region shares and direction-matched calls", {
  params <- list(
    min_abs_delta_pau = 0.10, event_min_treatment_pau = 0.05, gene_fdr = 0.05,
    site_fdr = 0.05, event_max_control_pau = 0.01, event_min_supporting_samples = 2,
    min_gene_total = 20, potential_internal_priming_withheld_calls = TRUE
  )
  # One row per PAC: coordinate, class, last exon, fitted PAU in each group,
  # event, and confidence. A low-confidence PAC is flagged for possible
  # internal priming; the other columns let a withheld call be classified.
  gene <- function(gene_id, strand, coordinate, region, last_exon, control, treatment, event,
                   confidence = "high", gene_fdr = 0.01, control_samples = 2, pac_fdr = 0.001,
                   exploratory = FALSE) {
    data.frame(
      gene_id = gene_id, strand = strand, coordinate = coordinate, gene_region = region,
      last_exon_locus = last_exon, fitted_control_pau = control, fitted_treatment_pau = treatment,
      event_type = event, confidence = confidence, gene_fdr = gene_fdr,
      internal_priming_flag = confidence == "low", control_supporting_samples = control_samples,
      treatment_supporting_samples = 2, delta_pau = treatment - control, pac_fdr = pac_fdr,
      zero_boundary_unstable = FALSE, exploratory_insufficient_replicates = exploratory,
      control_gene_total = 100, treatment_gene_total = 100, stringsAsFactors = FALSE
    )
  }
  intron <- "intron"
  last <- "last_exon"
  cases <- list(
    # Tanc2 in AS_NT_PHA vs AS_NT_DMSO: usage moves from the last exon into
    # introns. Every call in the last exon is a decrease, which the intronic
    # shift alone explains, so there is no UTR class.
    tanc2 = gene(
      "tanc2", "+",
      c(105481650, 105488375, 105613600, 105748725, 105752550, 105754075, 105776150,
        105814650, 105815100, 105817300, 105819650, 105820125),
      c(rep(intron, 7), rep(last, 5)), c(rep(NA, 7), rep("chr11:105812001-105820300", 5)),
      c(0.026083873, 0.012972031, 0.023594885, 0.000386844, 0.004493291, 0.032067591,
        0.020962887, 0.087370443, 0.364917433, 0.043324496, 0.018143946, 0.36568228),
      c(0.060759901, 0.028172593, 0.061173605, 0.189594627, 0.213368603, 0.029396933,
        0.025721693, 0.041423914, 0.176342297, 0.050776227, 0.011048109, 0.112221497),
      c("none", "none", "none", "gained", "increased_usage", "none", "none", "none",
        "decreased_usage", "none", "none", "decreased_usage"),
      c("low", "moderate", "moderate", "moderate", "moderate", "moderate", "low", "moderate",
        "moderate", "low", "moderate", "moderate"),
      gene_fdr = 0
    ),
    # The intronic rise is on an internal-priming site, which cannot gate a
    # pattern, only a potential one.
    primed = gene(
      "primed", "+", c(100, 500), c(intron, last), c(NA, "chr1:401-500"),
      c(0.1, 0.9), c(0.4, 0.6), c("increased_usage", "decreased_usage"), c("low", "high")
    ),
    intronic_loss = gene(
      "intronic_loss", "+", c(100, 500), c(intron, last), c(NA, "chr1:401-500"),
      c(0.5, 0.5), c(0.1, 0.9), c("decreased_usage", "increased_usage")
    ),
    # Two last exons of equal usage: the main one is the 3'-most.
    ale = gene(
      "ale", "+", c(1100, 2100), c(last, last), c("chr1:1001-1100", "chr1:2001-2100"),
      c(0.8, 0.2), c(0.2, 0.8), c("decreased_usage", "increased_usage")
    ),
    # Minus strand: the distal PAC has the smaller coordinate.
    shortening = gene(
      "shortening", "-", c(5100, 5000), c(last, last), rep("chr1:4901-5100", 2),
      c(0.3, 0.7), c(0.6, 0.4), c("increased_usage", "decreased_usage")
    ),
    lengthening = gene(
      "lengthening", "+", c(5000, 5100), c(last, last), rep("chr1:4901-5100", 2),
      c(0.6, 0.4), c(0.3, 0.7), c("decreased_usage", "increased_usage")
    ),
    # An intronic gain whose last exon also shifts toward its proximal PAC:
    # the proximal increase is a call the intronic shift cannot make.
    two_classes = gene(
      "two_classes", "+", c(100, 400, 500), c(intron, last, last), c(NA, rep("chr1:301-500", 2)),
      c(0.1, 0.2, 0.7), c(0.4, 0.4, 0.2), c("increased_usage", "increased_usage", "decreased_usage")
    ),
    candidate_only = gene(
      "candidate_only", "+", c(400, 500), c(last, last), rep("chr1:301-500", 2),
      c(0.0, 1.0), c(0.3, 0.7), c("gained_candidate", "none")
    ),
    unscreened = gene(
      "unscreened", "+", c(100, 500), c(intron, last), c(NA, "chr1:401-500"),
      c(0.1, 0.9), c(0.6, 0.4), c("increased_usage", "decreased_usage"), gene_fdr = 0.5
    ),
    unfitted = gene(
      "unfitted", "+", c(100, 500), c(intron, last), c(NA, "chr1:401-500"),
      c(NA, NA), c(NA, NA), c("none", "none")
    ),
    # Equal usage in two last exons: the 3'-most is the main one, and its
    # shift toward the distal PAC makes a UTR call.
    tie = gene(
      "tie", "+", c(1100, 2050, 2100), rep(last, 3),
      c("chr1:1001-1100", "chr1:2001-2100", "chr1:2001-2100"),
      c(0.5, 0.3, 0.2), c(0.5, 0.1, 0.4), c("none", "decreased_usage", "increased_usage")
    ),
    # An intronic gain that every last exon loses to is no last-exon switch.
    all_lose = gene(
      "all_lose", "+", c(100, 1100, 2100), c(intron, last, last),
      c(NA, "chr1:1001-1100", "chr1:2001-2100"),
      c(0.1, 0.45, 0.45), c(0.5, 0.25, 0.25),
      c("increased_usage", "decreased_usage", "decreased_usage")
    ),
    # A last exon holding under 5% of the gene has no UTR share.
    faint_exon = gene(
      "faint_exon", "+", c(100, 400, 500), c(intron, last, last), c(NA, rep("chr1:301-500", 2)),
      c(0.97, 0.02, 0.01), c(0.97, 0.01, 0.02), c("none", "none", "none")
    ),
    # Each pattern carried only by a flagged PAC's call.
    flagged_loss = gene(
      "flagged_loss", "+", c(100, 500), c(intron, last), c(NA, "chr1:401-500"),
      c(0.5, 0.5), c(0.1, 0.9), c("decreased_usage", "increased_usage"), c("low", "high")
    ),
    flagged_ale = gene(
      "flagged_ale", "+", c(1100, 2100), c(last, last), c("chr1:1001-1100", "chr1:2001-2100"),
      c(0.8, 0.2), c(0.2, 0.8), c("decreased_usage", "increased_usage"), c("high", "low")
    ),
    flagged_shortening = gene(
      "flagged_shortening", "+", c(400, 500), c(last, last), rep("chr1:301-500", 2),
      c(0.3, 0.7), c(0.6, 0.4), c("increased_usage", "none"), c("low", "high")
    ),
    flagged_lengthening = gene(
      "flagged_lengthening", "+", c(400, 500), c(last, last), rep("chr1:301-500", 2),
      c(0.6, 0.4), c(0.3, 0.7), c("none", "increased_usage"), c("high", "low")
    ),
    # A confirmed call at an unflagged site makes the pattern whatever the
    # flagged sites do.
    both_sites = gene(
      "both_sites", "+", c(100, 200, 500), c(intron, intron, last), c(NA, NA, "chr1:401-500"),
      c(0.05, 0.05, 0.9), c(0.25, 0.25, 0.5),
      c("increased_usage", "increased_usage", "decreased_usage"), c("low", "high", "high")
    ),
    # A confirmed intronic gain, and a UTR shift carried only by a flagged PAC.
    mixed = gene(
      "mixed", "+", c(100, 400, 500), c(intron, last, last), c(NA, rep("chr1:301-500", 2)),
      c(0.1, 0.2, 0.7), c(0.4, 0.4, 0.2), c("increased_usage", "increased_usage", "decreased_usage"),
      c("high", "low", "high")
    ),
    # A new flagged intronic site: the flag alone withheld its gain.
    withheld = gene(
      "withheld", "+", c(100, 500), c(intron, last), c(NA, "chr1:401-500"),
      c(0.0, 1.0), c(0.3, 0.7), c("gained_candidate", "decreased_usage"), c("low", "high"),
      control_samples = c(0, 2)
    ),
    only_withheld = gene(
      "only_withheld", "+", c(100, 500), c(intron, last), c(NA, "chr1:401-500"),
      c(0.0, 1.0), c(0.3, 0.7), c("gained_candidate", "none"), c("low", "high"),
      control_samples = c(0, 2)
    ),
    # Withheld also for lack of significance, or in an exploratory comparison.
    weak_withheld = gene(
      "weak_withheld", "+", c(100, 500), c(intron, last), c(NA, "chr1:401-500"),
      c(0.0, 1.0), c(0.3, 0.7), c("gained_candidate", "decreased_usage"), c("low", "high"),
      control_samples = c(0, 2), pac_fdr = c(0.5, 0.001)
    ),
    exploratory_withheld = gene(
      "exploratory_withheld", "+", c(100, 500), c(intron, last), c(NA, "chr1:401-500"),
      c(0.0, 1.0), c(0.3, 0.7), c("gained_candidate", "decreased_usage"), c("low", "high"),
      control_samples = c(0, 2), exploratory = TRUE
    )
  )
  pacs <- do.call(rbind, unname(cases))
  patterns <- model$apa_patterns(pacs, pacs$coordinate, params)
  observed <- stats::setNames(patterns$apa_pattern, patterns$gene_id)
  potential <- function(pattern) paste0(pattern, "_potential_internal_priming")
  expected <- c(
    tanc2 = "intronic_gain", primed = potential("intronic_gain"),
    intronic_loss = "intronic_loss", ale = "alternative_last_exon",
    shortening = "utr_shortening", lengthening = "utr_lengthening",
    two_classes = "intronic_gain;utr_shortening", candidate_only = "none",
    unscreened = "none", unfitted = "none", tie = "utr_lengthening",
    all_lose = "intronic_gain", faint_exon = "none",
    flagged_loss = potential("intronic_loss"), flagged_ale = potential("alternative_last_exon"),
    flagged_shortening = potential("utr_shortening"),
    flagged_lengthening = potential("utr_lengthening"), both_sites = "intronic_gain",
    mixed = paste0("intronic_gain;", potential("utr_shortening")),
    withheld = potential("intronic_gain"), only_withheld = potential("intronic_gain"),
    weak_withheld = "unclassified_change", exploratory_withheld = "unclassified_change"
  )
  check(
    identical(observed, expected),
    "patterns were ", paste(names(observed), observed, sep = "=", collapse = ", ")
  )
  # Without withheld calls, a flagged gain that the flag alone withheld
  # supports nothing, but a flagged PAC's confirmed call still does.
  off <- model$apa_patterns(
    pacs, pacs$coordinate, utils::modifyList(params, list(potential_internal_priming_withheld_calls = FALSE))
  )
  off <- stats::setNames(off$apa_pattern, off$gene_id)
  check(
    identical(off[c("withheld", "only_withheld", "primed")],
      c(withheld = "unclassified_change", only_withheld = "none",
        primed = potential("intronic_gain"))),
    "without withheld calls: ", paste(names(off), off, sep = "=", collapse = ", ")
  )
  unset <- tryCatch(
    model$apa_patterns(pacs, pacs$coordinate, params[names(params) != "potential_internal_priming_withheld_calls"]),
    error = function(error) conditionMessage(error)
  )
  check(identical(unset, "potential_internal_priming_withheld_calls must be true or false."), unset)
  metric <- function(gene_id, column) patterns[[column]][patterns$gene_id == gene_id]
  near <- function(value, target) isTRUE(abs(value - target) < 1e-9)
  tanc2 <- cases$tanc2
  exon <- tanc2$gene_region == last
  check(
    near(metric("tanc2", "delta_intronic_share"),
      sum(tanc2$fitted_treatment_pau[!exon]) - sum(tanc2$fitted_control_pau[!exon])),
    "Tanc2 intronic share: ", metric("tanc2", "delta_intronic_share")
  )
  check(
    near(metric("tanc2", "delta_utr_distal_share"),
      0.112221497 / sum(tanc2$fitted_treatment_pau[exon]) -
        0.36568228 / sum(tanc2$fitted_control_pau[exon])),
    "Tanc2 distal share: ", metric("tanc2", "delta_utr_distal_share")
  )
  check(is.na(metric("tanc2", "last_exon_switch")), "one last exon has no switch.")
  check(near(metric("ale", "last_exon_switch"), 0.6), "ALE switch: ", metric("ale", "last_exon_switch"))
  check(is.na(metric("ale", "delta_utr_distal_share")), "a one-PAC main last exon has a UTR share.")
  check(near(metric("shortening", "delta_utr_distal_share"), -0.3), "minus-strand distal share.")
  check(near(metric("unscreened", "delta_intronic_share"), 0.5), "unscreened genes lack numbers.")
  check(
    all(is.na(unlist(patterns[patterns$gene_id == "unfitted", model$APA_METRIC_COLUMNS]))),
    "an unfitted gene has numbers."
  )
  check(is.na(metric("faint_exon", "delta_utr_distal_share")), "a faint last exon has a UTR share.")
  check(near(metric("tie", "delta_utr_distal_share"), 0.4), "the tie took the wrong last exon.")
  check(near(metric("all_lose", "last_exon_switch"), -0.2), "switch: ", metric("all_lose", "last_exon_switch"))

  # The fitted run's genes tables carry one valid pattern per gene.
  run <- require_run(run_a)
  classes <- c(
    model$APA_PATTERN_CLASSES, paste0(model$APA_PATTERN_CLASSES, model$POTENTIAL_INTERNAL_PRIMING),
    "unclassified_change", "none"
  )
  for (comparison in c("T1_vs_C", "T2_vs_C")) {
    genes <- read_result(run$final_directory, paste0(comparison, ".genes.tsv.gz"))
    parts <- strsplit(genes$apa_pattern, ";", fixed = TRUE)
    check(
      all(vapply(parts, function(part) length(part) && all(part %in% classes), logical(1))),
      comparison, ": invalid patterns ", paste(unique(genes$apa_pattern), collapse = ", ")
    )
  }
})

finish_tests()
