#!/usr/bin/env Rscript

# Differential PAC usage within one comparison family, fitted with DRIMSeq and
# adjusted with stageR. The pipeline runs this file in four modes:
#   versions   append R and package versions to a software-versions table;
#   fit        statistical filtering, family model, zero-count stabilization,
#              motif preference, preliminary tables, and bootstrap batches;
#   bootstrap  percentile intervals for one batch of genes;
#   finalize   intervals merged into the preliminary tables, then event tables.

# DRIMSeq's addUniform rule: each zero count becomes a draw from U(0, 0.1).
ZERO_PERTURBATION_MAX <- 0.1
# Design section 9: a stabilized fit needs at least 80 percent of its repeats.
ZERO_MIN_SUCCESS_FRACTION <- 0.8
# Changes in fitted PAU smaller than this count as zero when checking whether
# the direction of an effect flips between stabilization repeats.
DIRECTION_TOLERANCE <- 0.001
DEFAULT_BOOTSTRAP_BATCH_SIZE <- 500L
REQUIRED_PRECISION_SLOTS <- c(
  "mean_expression", "common_precision", "genewise_precision",
  "design_precision", "counts", "samples"
)
TEXT_COLUMNS <- c(
  "gene_id", "gene_name", "chrom", "locus", "pac_id", "feature_id", "sample_id", "condition", "control",
  "control_condition", "comparison", "bootstrap_status", "model_status"
)
BOOTSTRAP_STATUSES <- c(
  "ok", "insufficient_successes", "fit_unavailable", "not_selected", "disabled"
)
INTERVAL_COLUMNS <- c(
  "comparison", "gene_id", "feature_id", "delta_pau_ci_low", "delta_pau_ci_high",
  "bootstrap_successes", "bootstrap_perturbed", "bootstrap_status"
)
# Every PAC-level table starts with the identity block; gene-level tables
# start with the gene's ID and name. The other columns run from the answer to
# the evidence: comparison, call, effect, significance, annotation, data, and
# model diagnostics.
LOCATION_COLUMNS <- c("chrom", "start", "end", "strand", "locus")
IDENTITY_COLUMNS <- c("pac_id", "gene_id", "gene_name", LOCATION_COLUMNS)
GENE_PRECISION_COLUMNS <- c("gene_id", "gene_name", "family", "precision", "model_status")
FILTERING_COLUMNS <- c(IDENTITY_COLUMNS, "family", "tested", "reason")
MOTIF_PREFERENCE_COLUMNS <- c(
  "primary_pas_motif_rna", "primary_motif_class", "condition", "control_condition",
  "delta_motif_usage", "fdr", "p_value", "control_mean", "treatment_mean",
  "transformed_coefficient", "informative_genes"
)
# Motif preference is tested per primary hexamer and, in a second table, per
# motif class.
MOTIF_PREFERENCE_KEYS <- c("primary_pas_motif_rna", "primary_motif_class")
FITTED_PAU_COLUMNS <- c(
  IDENTITY_COLUMNS, "condition", "control_condition", "fitted_control_pau",
  "fitted_treatment_pau", "delta_pau", "model_status", "precision", "alpha_control",
  "alpha_treatment"
)
OMNIBUS_COLUMNS <- c(
  "gene_id", "gene_name", "family", "gene_fdr", "pvalue", "lr", "df", "model_status",
  "stabilization_successes", "exploratory_insufficient_replicates"
)
GENE_COLUMNS <- c(
  "gene_id", "gene_name", "condition", "control_condition", "gene_fdr", "pvalue", "lr",
  "df", "model_status", "stabilization_successes", "exploratory_insufficient_replicates"
)
ATLAS_ANNOTATION_COLUMNS <- c(
  "assignment_class", "confidence", "internal_priming_flag", "known_pac",
  "known_rescue_only", "primary_pas_motif", "primary_pas_motif_rna",
  "primary_motif_class"
)
PAC_COLUMNS <- c(
  IDENTITY_COLUMNS,
  "condition", "control_condition", "event_type",
  "fitted_control_pau", "fitted_treatment_pau", "delta_pau", "delta_pau_ci_low",
  "delta_pau_ci_high",
  "pac_fdr", "gene_fdr", "pvalue_pac", "pvalue_gene", "lr", "df",
  ATLAS_ANNOTATION_COLUMNS,
  "control_supporting_samples", "treatment_supporting_samples", "control_gene_total",
  "treatment_gene_total", "raw_control_counts", "raw_treatment_counts",
  "observed_control_pau", "observed_treatment_pau",
  "effect_exceeds_threshold", "dominant_pac_control", "dominant_pac_treatment",
  "control_detected_complexity", "treatment_detected_complexity",
  "model_status", "precision", "alpha_control", "alpha_treatment",
  "stabilization_successes", "stabilization_delta_pau_spread", "zero_boundary_unstable",
  "zero_boundary_reason",
  "bootstrap_status", "bootstrap_successes", "bootstrap_perturbed",
  "exploratory_insufficient_replicates"
)

load_statistics_packages <- function() {
  required <- c("BiocParallel", "DRIMSeq", "stageR", "limma", "yaml")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". Use the PACusage Conda environment or Apptainer image."
    )
  }
  suppressPackageStartupMessages({
    library(DRIMSeq)
    library(stageR)
    library(limma)
  })
  # Fixed-precision refits build DRIMSeq precision objects directly, so a
  # DRIMSeq release with different slots must stop here, not mid-run.
  absent <- setdiff(REQUIRED_PRECISION_SLOTS, methods::slotNames("dmDSprecision"))
  if (length(absent)) {
    stop(
      "DRIMSeq ", as.character(utils::packageVersion("DRIMSeq")),
      " lacks dmDSprecision slots: ", paste(absent, collapse = ", "),
      ". PACusage was validated with DRIMSeq 1.38."
    )
  }
  versions <- vapply(
    c("DRIMSeq", "stageR", "limma", "BiocParallel"),
    function(package) as.character(utils::packageVersion(package)),
    character(1)
  )
  message(
    "PACusage statistics with ",
    paste(names(versions), versions, collapse = ", ")
  )
  invisible(versions)
}

parse_args <- function(arguments) {
  result <- list()
  index <- 1
  while (index <= length(arguments)) {
    key <- sub("^--", "", arguments[[index]])
    if (index == length(arguments)) stop("Missing value for --", key)
    result[[gsub("-", "_", key)]] <- arguments[[index + 1]]
    index <- index + 2
  }
  if (is.null(result$mode)) result$mode <- "fit"
  result
}

require_args <- function(arguments, required) {
  missing <- required[!vapply(
    required,
    function(key) !is.null(arguments[[key]]) && arguments[[key]] != "",
    logical(1)
  )]
  if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))
}

parse_worker_count <- function(value, option_name) {
  if (is.null(value)) return(1L)
  workers <- suppressWarnings(as.numeric(value))
  if (
    length(workers) != 1L ||
      !is.finite(workers) ||
      workers < 1 ||
      workers != floor(workers)
  ) {
    stop("--", option_name, " must be a positive integer.")
  }
  as.integer(workers)
}

parse_bootstrap_workers <- function(value) {
  parse_worker_count(value, "bootstrap-workers")
}

parse_model_workers <- function(value) {
  parse_worker_count(value, "model-workers")
}

parse_batch_size <- function(value) {
  if (is.null(value)) return(DEFAULT_BOOTSTRAP_BATCH_SIZE)
  size <- suppressWarnings(as.numeric(value))
  if (
    length(size) != 1L ||
      !is.finite(size) ||
      size < 1 ||
      size != floor(size)
  ) {
    stop("--bootstrap-batch-size must be a positive integer.")
  }
  as.integer(size)
}

drimseq_bpparam <- function(workers, seed) {
  if (workers == 1L || .Platform$OS.type == "windows") {
    return(BiocParallel::SerialParam())
  }
  BiocParallel::MulticoreParam(
    workers = workers,
    RNGseed = seed,
    progressbar = FALSE,
    stop.on.error = TRUE
  )
}

read_header <- function(path) {
  connection <- file(path, "rt")
  on.exit(close(connection), add = TRUE)
  line <- readLines(connection, n = 1L, warn = FALSE)
  if (!length(line)) return(character())
  strsplit(line, "\t", fixed = TRUE)[[1]]
}

# Identifiers stay text: numeric-looking gene or sample IDs would otherwise
# become integers, which DRIMSeq rejects.
read_tsv <- function(path) {
  text_columns <- intersect(TEXT_COLUMNS, read_header(path))
  utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA"),
    colClasses = stats::setNames(rep("character", length(text_columns)), text_columns)
  )
}

write_gzip_tsv <- function(value, path) {
  connection <- gzfile(path, "wt")
  on.exit(close(connection), add = TRUE)
  write.table(value, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

select_columns <- function(table, columns, label) {
  missing <- setdiff(columns, names(table))
  if (length(missing)) {
    stop(label, " lacks columns: ", paste(missing, collapse = ", "), ".")
  }
  table[, columns, drop = FALSE]
}

empty_table <- function(columns) {
  as.data.frame(
    stats::setNames(replicate(length(columns), character(), simplify = FALSE), columns),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

bh <- function(values) {
  result <- rep(NA_real_, length(values))
  finite <- is.finite(values)
  result[finite] <- p.adjust(values[finite], method = "BH")
  result
}

chi_square_p <- function(lr, df) {
  result <- rep(NA_real_, length(lr))
  finite <- is.finite(lr) & is.finite(df)
  result[finite] <- stats::pchisq(lr[finite], df[finite], lower.tail = FALSE)
  result
}

# Statistical filtering within one comparison family (design section 8), on
# the raw counts of the family's samples only. Returns the tested PACs and one
# row per PAC of the count table: tested, or the reasons it was not. A PAC
# without a gene, or assigned to more than one gene, is never tested.
family_filter <- function(counts, sample_ids, params, family) {
  gene_ids <- counts$gene_id
  reason <- rep("", nrow(counts))
  reason[is.na(gene_ids) | gene_ids == ""] <- "no_gene_assignment"
  reason[reason == "" & grepl(",", gene_ids, fixed = TRUE)] <- "ambiguous_gene_assignment"
  assigned <- reason == ""

  values <- count_matrix(counts, sample_ids)
  site_total <- rowSums(values)
  supporting <- rowSums(values > 0)
  gene_total <- rep(NA_real_, nrow(counts))
  gene_total[assigned] <- ave(site_total[assigned], gene_ids[assigned], FUN = sum)
  site_usage <- ifelse(assigned & gene_total > 0, site_total / gene_total, 0)
  reason[assigned] <- joined_reasons(
    site_total < params$min_site_count,
    paste0("site_count<", params$min_site_count),
    supporting < params$min_test_supporting_samples,
    paste0("supporting_samples<", params$min_test_supporting_samples),
    site_usage < params$min_site_usage,
    paste0("site_usage<", params$min_site_usage)
  )[assigned]

  site_ok <- assigned & reason == ""
  eligible <- rep(0, nrow(counts))
  eligible[assigned] <- ave(as.numeric(site_ok[assigned]), gene_ids[assigned], FUN = sum)
  gene_reason <- joined_reasons(
    gene_total < params$min_gene_total,
    paste0("gene_total<", params$min_gene_total),
    eligible < 2,
    "fewer_than_2_testable_pacs"
  )
  reason[site_ok] <- gene_reason[site_ok]
  tested <- site_ok & reason == ""
  list(
    counts = counts[tested, c("gene_id", "pac_id", sample_ids), drop = FALSE],
    reasons = data.frame(
      pac_id = counts$pac_id,
      gene_id = ifelse(is.na(gene_ids), "", gene_ids),
      gene_name = ifelse(is.na(counts$gene_name), "", counts$gene_name),
      counts[, LOCATION_COLUMNS, drop = FALSE],
      family = rep(family, nrow(counts)),
      tested = tested,
      reason = reason,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
}

# Arguments alternate between a logical vector and its label. Returns, per
# element, the labels of the conditions that hold, joined by ";".
joined_reasons <- function(...) {
  arguments <- list(...)
  conditions <- arguments[c(TRUE, FALSE)]
  labels <- arguments[c(FALSE, TRUE)]
  if (!length(conditions[[1]])) return(character())
  parts <- mapply(function(condition, label) ifelse(condition, label, ""), conditions, labels)
  if (is.null(dim(parts))) parts <- matrix(parts, nrow = length(conditions[[1]]))
  apply(parts, 1, function(row) paste(row[nzchar(row)], collapse = ";"))
}

supporting_samples <- function(counts, sample_rows, group_name) {
  ids <- sample_rows$sample_id[sample_rows$condition == group_name]
  rowSums(counts[, ids, drop = FALSE] > 0)
}

group_gene_totals <- function(counts, sample_rows, group_name) {
  ids <- sample_rows$sample_id[sample_rows$condition == group_name]
  values <- rowSums(counts[, ids, drop = FALSE])
  as.numeric(ave(values, counts$gene_id, FUN = sum))
}

raw_count_text <- function(counts, sample_rows, group_name) {
  sample_ids <- sample_rows$sample_id[sample_rows$condition == group_name]
  apply(
    counts[, sample_ids, drop = FALSE],
    1,
    function(values) paste(paste(sample_ids, values, sep = "="), collapse = ",")
  )
}

observed_pau_text <- function(counts, sample_rows, group_name) {
  sample_ids <- sample_rows$sample_id[sample_rows$condition == group_name]
  usage <- sapply(sample_ids, function(sample_id) {
    totals <- ave(counts[[sample_id]], counts$gene_id, FUN = sum)
    ifelse(totals > 0, counts[[sample_id]] / totals, NA_real_)
  })
  if (is.null(dim(usage))) usage <- matrix(usage, ncol = 1)
  apply(
    usage,
    1,
    function(values) {
      paste(
        paste(sample_ids, format(values, digits = 6, trim = TRUE), sep = "="),
        collapse = ","
      )
    }
  )
}

# Deterministic seeds from a polynomial hash of the text modulo the prime
# 2^31 - 1. A position-weighted character sum would give, for example,
# bootstrap replicates 120 and 201 of a gene the same seed. The fields are
# joined with spaces, paste's default separator.
SEED_MODULUS <- 2147483647
SEED_POWERS <- local({
  powers <- numeric(4096L)
  value <- 1
  for (index in seq_along(powers)) {
    powers[[index]] <- value
    value <- (value * 131) %% SEED_MODULUS
  }
  powers
})

stable_seed <- function(base_seed, ...) {
  values <- utf8ToInt(paste(..., collapse = "|"))
  if (length(values) > length(SEED_POWERS)) stop("Seed text is too long.")
  # Horner's hash, sum(value_i * 131^(n - i)), with every product below 2^53.
  terms <- (values * SEED_POWERS[rev(seq_along(values))]) %% SEED_MODULUS
  offset <- sum(terms) %% SEED_MODULUS
  as.integer((as.numeric(base_seed) + offset) %% .Machine$integer.max)
}

bootstrap_total <- function(values, gene_id, sample_id) {
  total <- sum(values)
  if (
    length(total) != 1L ||
      !is.finite(total) ||
      total < 0 ||
      total > .Machine$integer.max ||
      abs(total - round(total)) > sqrt(.Machine$double.eps) * max(1, abs(total))
  ) {
    stop(
      "Invalid bootstrap total for gene ", gene_id,
      ", sample ", sample_id,
      ": ", format(total),
      ". Expected a finite non-negative integer no greater than ",
      .Machine$integer.max, "."
    )
  }
  as.integer(round(total))
}

bootstrap_gene_ids <- function(
  gene_ids,
  gene_fdr,
  detection_candidate,
  gene_fdr_cutoff,
  include_candidates = TRUE
) {
  selected <- !is.na(gene_fdr) & gene_fdr <= gene_fdr_cutoff
  if (is.null(include_candidates) || isTRUE(include_candidates)) {
    selected <- selected | (!is.na(detection_candidate) & detection_candidate)
  }
  unique(gene_ids[selected & !is.na(gene_ids) & gene_ids != ""])
}

bootstrap_apply <- function(repeat_numbers, workers, worker) {
  if (!length(repeat_numbers)) return(list())
  workers <- as.integer(workers)
  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("Bootstrap worker count must be a positive integer.")
  }
  workers <- min(workers, length(repeat_numbers))
  if (workers == 1L || .Platform$OS.type == "windows") {
    return(lapply(repeat_numbers, worker))
  }
  parallel::mclapply(
    repeat_numbers,
    worker,
    mc.cores = workers,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
}

bootstrap_settings <- function(params) {
  bootstrap_replicates_value <- suppressWarnings(
    as.numeric(params$dm_bootstrap_replicates)
  )
  if (
    length(bootstrap_replicates_value) != 1L ||
      !is.finite(bootstrap_replicates_value) ||
      bootstrap_replicates_value < 0 ||
      bootstrap_replicates_value != floor(bootstrap_replicates_value)
  ) {
    stop("dm_bootstrap_replicates must be a non-negative integer.")
  }
  bootstrap_min_success_fraction <- suppressWarnings(
    as.numeric(params$dm_bootstrap_min_success_fraction)
  )
  if (
    length(bootstrap_min_success_fraction) != 1L ||
      !is.finite(bootstrap_min_success_fraction) ||
      bootstrap_min_success_fraction < 0 ||
      bootstrap_min_success_fraction > 1
  ) {
    stop("dm_bootstrap_min_success_fraction must be a finite value in [0, 1].")
  }
  replicates <- as.integer(bootstrap_replicates_value)
  list(
    replicates = replicates,
    min_success_fraction = bootstrap_min_success_fraction,
    min_successes = max(1L, as.integer(ceiling(
      replicates * bootstrap_min_success_fraction - 1e-9
    )))
  )
}

zero_sensitivity_settings <- function(params) {
  repeats <- suppressWarnings(as.numeric(params$dm_zero_sensitivity_repeats))
  if (length(repeats) != 1L || !is.finite(repeats) || repeats < 1 || repeats != floor(repeats)) {
    stop("dm_zero_sensitivity_repeats must be a positive integer.")
  }
  spread <- suppressWarnings(as.numeric(params$dm_zero_max_delta_pau_spread))
  if (length(spread) != 1L || !is.finite(spread) || spread < 0) {
    stop("dm_zero_max_delta_pau_spread must be a finite non-negative value.")
  }
  repeats <- as.integer(repeats)
  list(
    repeats = repeats,
    max_spread = spread,
    min_successes = as.integer(ceiling(ZERO_MIN_SUCCESS_FRACTION * repeats - 1e-9))
  )
}

dominant_pac <- function(feature_ids, fitted_pau) {
  finite_indices <- which(is.finite(fitted_pau))
  if (!length(finite_indices)) return(NA_character_)
  feature_ids[finite_indices[[which.max(fitted_pau[finite_indices])]]]
}

safe_file_component <- function(value) {
  value <- gsub("[^A-Za-z0-9._-]+", "_", as.character(value))
  if (!nzchar(value)) value <- "comparison"
  value
}

# ---- Family layout and counts --------------------------------------------

family_layout <- function(family, sample_rows, params) {
  treatments <- sort(
    unique(sample_rows$condition[
      sample_rows$control_condition == family & sample_rows$condition != family
    ]),
    method = "radix"
  )
  if (!length(treatments)) return(NULL)
  samples <- sample_rows[sample_rows$condition %in% c(family, treatments), , drop = FALSE]
  rownames(samples) <- NULL
  samples$condition <- factor(samples$condition, levels = c(family, treatments))
  covariates <- as.character(unlist(params$model_covariates))
  for (covariate in covariates) {
    if (!covariate %in% names(samples)) {
      stop("Model covariate ", covariate, " is missing from the sample sheet.")
    }
    values <- as.character(samples[[covariate]])
    if (anyNA(values) || any(values == "")) {
      stop("Model covariate ", covariate, " is incomplete in comparison family ", family, ".")
    }
    # Covariates are categorical.
    samples[[covariate]] <- factor(values)
  }
  formula_text <- paste("~", paste(c(covariates, "condition"), collapse = " + "))
  design <- model.matrix(as.formula(formula_text), data = samples)
  if (qr(design)$rank < ncol(design)) {
    stop(
      "Comparison family ", family,
      " has a rank-deficient design. Check condition/covariate confounding."
    )
  }
  coefficients <- match(paste0("condition", treatments), colnames(design))
  if (anyNA(coefficients)) {
    stop(
      "No model coefficient found for ",
      paste0(treatments[is.na(coefficients)], "_vs_", family, collapse = ", ")
    )
  }
  conditions <- c(family, treatments)
  groups <- lapply(conditions, function(condition) {
    samples$sample_id[samples$condition == condition]
  })
  names(groups) <- conditions
  zero_levels <- c(
    list(condition = samples$condition),
    stats::setNames(lapply(covariates, function(covariate) samples[[covariate]]), covariates)
  )
  list(
    family = family,
    control = family,
    treatments = treatments,
    comparisons = data.frame(
      comparison = paste0(treatments, "_vs_", family),
      treatment = treatments,
      coefficient = coefficients,
      stringsAsFactors = FALSE
    ),
    samples = samples,
    sample_ids = samples$sample_id,
    dm_samples = samples[, c("sample_id", "condition", covariates), drop = FALSE],
    design = design,
    covariates = covariates,
    condition_coefficients = coefficients,
    groups = groups,
    zero_levels = zero_levels,
    exploratory = any(table(samples$condition) < params$min_replicates_per_condition)
  )
}

# DRIMSeq treats the last PAC of a gene as the reference, so the order must
# be deterministic: genomic coordinate within each gene (design section 9).
order_family_counts <- function(counts, atlas) {
  coordinate <- suppressWarnings(as.numeric(atlas$coordinate[match(counts$pac_id, atlas$pac_id)]))
  if (anyNA(coordinate)) {
    stop(
      "The atlas lacks coordinates for PACs: ",
      paste(utils::head(counts$pac_id[is.na(coordinate)], 5L), collapse = ", "), "."
    )
  }
  ordered <- counts[order(counts$gene_id, coordinate, counts$pac_id, method = "radix"), , drop = FALSE]
  rownames(ordered) <- NULL
  per_gene <- table(ordered$gene_id)
  if (any(per_gene < 2L)) {
    stop(
      "Genes with fewer than 2 testable PACs reached the model: ",
      paste(utils::head(names(per_gene)[per_gene < 2L], 5L), collapse = ", "), "."
    )
  }
  ordered
}

count_matrix <- function(counts, sample_ids) {
  values <- as.matrix(counts[, sample_ids, drop = FALSE])
  storage.mode(values) <- "double"
  values
}

# Genes with a PAC that has no counts in any sample of a condition level or a
# covariate level. DRIMSeq's fit for such genes is at a boundary, even when it
# silently returns a fitted proportion of zero.
zero_group_genes <- function(counts, layout) {
  values <- count_matrix(counts, layout$sample_ids)
  flagged <- rep(FALSE, nrow(values))
  for (levels in layout$zero_levels) {
    labels <- as.character(levels)
    for (level in unique(labels)) {
      flagged <- flagged | rowSums(values[, labels == level, drop = FALSE]) == 0
    }
  }
  unique(counts$gene_id[flagged])
}

# Gene by condition: TRUE when a gene has no counts at all in that condition.
empty_condition_groups <- function(counts, layout) {
  totals <- rowsum(count_matrix(counts, layout$sample_ids), counts$gene_id, reorder = FALSE)
  empty <- vapply(
    layout$groups,
    function(sample_ids) rowSums(totals[, sample_ids, drop = FALSE]) == 0,
    logical(nrow(totals))
  )
  if (is.null(dim(empty))) empty <- matrix(empty, nrow = 1L, dimnames = list(rownames(totals), names(layout$groups)))
  rownames(empty) <- rownames(totals)
  empty
}

perturb_zero_cells <- function(values, seed) {
  set.seed(seed)
  zeros <- values == 0
  values[zeros] <- stats::runif(sum(zeros), 0, ZERO_PERTURBATION_MAX)
  values
}

# Returns a copy in which each listed gene's zero counts are perturbed with a
# seed specific to that gene. The input table is never modified.
perturb_gene_rows <- function(counts, gene_ids, sample_ids, seed_for_gene) {
  values <- count_matrix(counts, sample_ids)
  rows_by_gene <- split(seq_len(nrow(counts)), counts$gene_id)
  for (gene_id in gene_ids) {
    rows <- rows_by_gene[[gene_id]]
    values[rows, ] <- perturb_zero_cells(values[rows, , drop = FALSE], seed_for_gene(gene_id))
  }
  perturbed <- counts
  for (index in seq_along(sample_ids)) perturbed[[sample_ids[[index]]]] <- values[, index]
  perturbed
}

# ---- DRIMSeq wrappers -----------------------------------------------------

dm_data <- function(counts, sample_ids, dm_samples) {
  table <- counts[, c("gene_id", "pac_id", sample_ids), drop = FALSE]
  names(table)[names(table) == "pac_id"] <- "feature_id"
  DRIMSeq::dmDSdata(counts = table, samples = dm_samples)
}

estimate_family_precision <- function(counts, layout, bpparam, prec_init = NULL) {
  data <- dm_data(counts, layout$sample_ids, layout$dm_samples)
  if (is.null(prec_init)) {
    return(DRIMSeq::dmPrecision(data, design = layout$design, verbose = 0, BPPARAM = bpparam))
  }
  # With the common precision fixed, the precision grid matches the
  # unmodified fit exactly and no random subset of genes is drawn.
  DRIMSeq::dmPrecision(
    data,
    design = layout$design,
    common_precision = FALSE,
    prec_init = prec_init,
    verbose = 0,
    BPPARAM = bpparam
  )
}

fixed_precision_object <- function(data, design, precision, common_precision, mean_expression) {
  genes <- names(data@counts)
  values <- precision[genes]
  if (anyNA(values) || any(!is.finite(values)) || any(values <= 0)) {
    stop("Fixed precision is missing or invalid for genes: ", paste(utils::head(genes, 5L), collapse = ", "))
  }
  names(values) <- genes
  expression <- mean_expression[genes]
  names(expression) <- genes
  methods::new(
    "dmDSprecision",
    mean_expression = as.numeric(expression),
    common_precision = as.numeric(common_precision),
    genewise_precision = values,
    design_precision = design,
    counts = data@counts,
    samples = data@samples
  )
}

gene_test_table <- function(result, gene_ids) {
  index <- match(gene_ids, as.character(result$gene_id))
  data.frame(
    gene_id = gene_ids,
    lr = as.numeric(result$lr[index]),
    df = as.numeric(result$df[index]),
    pvalue = as.numeric(result$pvalue[index]),
    stringsAsFactors = FALSE
  )
}

feature_test_table <- function(result, features) {
  index <- match(
    paste(features$gene_id, features$feature_id, sep = "\r"),
    paste(result$gene_id, result$feature_id, sep = "\r")
  )
  data.frame(
    gene_id = features$gene_id,
    feature_id = features$feature_id,
    lr = as.numeric(result$lr[index]),
    df = as.numeric(result$df[index]),
    pvalue = as.numeric(result$pvalue[index]),
    stringsAsFactors = FALSE
  )
}

# Fit and test one precision object, then return plain tables in the order of
# `features`, so later steps never depend on DRIMSeq's internal ordering. With
# allow_missing, features absent from the object are returned as NA.
fit_and_test <- function(precision_object, layout, bpparam, features, allow_missing = FALSE) {
  fit <- DRIMSeq::dmFit(precision_object, design = layout$design, verbose = 0, BPPARAM = bpparam)
  omnibus <- DRIMSeq::dmTest(
    fit,
    coef = layout$condition_coefficients,
    bb_model = FALSE,
    verbose = 0,
    BPPARAM = bpparam
  )
  contrasts <- lapply(layout$comparisons$coefficient, function(coefficient) {
    DRIMSeq::dmTest(fit, coef = coefficient, verbose = 0, BPPARAM = bpparam)
  })
  names(contrasts) <- layout$comparisons$comparison
  proportions_table <- DRIMSeq::proportions(fit)
  index <- match(
    paste(features$gene_id, features$feature_id, sep = "\r"),
    paste(proportions_table$gene_id, proportions_table$feature_id, sep = "\r")
  )
  if (anyNA(index) && !allow_missing) {
    stop("DRIMSeq returned proportions for a different set of PACs.")
  }
  proportions <- matrix(
    NA_real_, nrow(features), length(layout$sample_ids),
    dimnames = list(features$feature_id, layout$sample_ids)
  )
  found <- !is.na(index)
  proportions[found, ] <- as.matrix(proportions_table[index[found], layout$sample_ids, drop = FALSE])
  gene_ids <- unique(features$gene_id)
  precision_table <- DRIMSeq::genewise_precision(fit)
  precision <- as.numeric(precision_table$genewise_precision[match(gene_ids, precision_table$gene_id)])
  names(precision) <- gene_ids
  list(
    gene_ids = gene_ids,
    features = features,
    proportions = proportions,
    precision = precision,
    omnibus = gene_test_table(DRIMSeq::results(omnibus, level = "gene"), gene_ids),
    contrasts = lapply(contrasts, function(test) {
      list(
        genes = gene_test_table(DRIMSeq::results(test, level = "gene"), gene_ids),
        features = feature_test_table(DRIMSeq::results(test, level = "feature"), features)
      )
    })
  )
}

# ---- Zero-count stabilization ---------------------------------------------

boundary_gene_sets <- function(fit, counts, layout) {
  gene_ids <- fit$gene_ids
  rows <- split(seq_len(nrow(fit$features)), factor(fit$features$gene_id, levels = gene_ids))
  nonfinite_rows <- rowSums(!is.finite(fit$proportions)) > 0
  for (contrast in fit$contrasts) {
    nonfinite_rows <- nonfinite_rows | !is.finite(contrast$features$pvalue)
  }
  nonfinite <- vapply(rows, function(index) any(nonfinite_rows[index]), logical(1))
  nonfinite <- nonfinite | !(is.finite(fit$precision) & fit$precision > 0) |
    !is.finite(fit$omnibus$pvalue)
  for (contrast in fit$contrasts) nonfinite <- nonfinite | !is.finite(contrast$genes$pvalue)
  boundary <- gene_ids[nonfinite | gene_ids %in% zero_group_genes(counts, layout)]
  values <- count_matrix(counts, layout$sample_ids)
  has_zero <- unique(counts$gene_id[rowSums(values == 0) > 0])
  perturbable <- boundary[boundary %in% has_zero]
  list(
    boundary = boundary,
    perturbable = perturbable,
    unavailable = setdiff(boundary, perturbable)
  )
}

run_zero_sensitivity <- function(counts, genes, layout, common_precision, params, atlas_checksum, bpparam) {
  settings <- zero_sensitivity_settings(params)
  rows <- counts$gene_id %in% genes
  features <- data.frame(
    gene_id = counts$gene_id[rows],
    feature_id = counts$pac_id[rows],
    stringsAsFactors = FALSE
  )
  lapply(seq_len(settings$repeats), function(repeat_number) {
    perturbed <- perturb_gene_rows(
      counts,
      genes,
      layout$sample_ids,
      function(gene_id) {
        stable_seed(
          params$random_seed, atlas_checksum, layout$family, gene_id,
          "zero_sensitivity", repeat_number
        )
      }
    )
    tryCatch(
      {
        # Whole-family precision keeps the trended moderation of the unmodified
        # fit; the fits and tests then only need the boundary genes.
        precision <- estimate_family_precision(perturbed, layout, bpparam, common_precision)
        values <- precision@genewise_precision[genes]
        # A gene without a usable precision fails only its own repeat.
        fitted_genes <- genes[is.finite(values) & values > 0]
        if (!length(fitted_genes)) stop("No stabilized gene has a finite precision.")
        data <- dm_data(
          perturbed[perturbed$gene_id %in% fitted_genes, , drop = FALSE],
          layout$sample_ids,
          layout$dm_samples
        )
        object <- fixed_precision_object(
          data,
          layout$design,
          precision@genewise_precision,
          common_precision,
          precision@mean_expression
        )
        fit_and_test(object, layout, bpparam, features, allow_missing = TRUE)
      },
      error = function(error) error
    )
  })
}

first_finite_by_row <- function(values, rows) {
  if (is.null(dim(values))) values <- matrix(values, nrow = rows)
  apply(values, 1, function(row) {
    finite <- row[is.finite(row)]
    if (length(finite)) finite[[1]] else NA_real_
  })
}

median_or_na <- function(values) {
  values <- values[is.finite(values)]
  if (!length(values)) return(NA_real_)
  stats::median(values)
}

# Combine the stabilization repeats of each boundary gene into one set of
# values: medians of proportions, precision, and likelihood ratios, with the
# stability rules of design section 9 applied per PAC and comparison.
summarize_zero_sensitivity <- function(repeats, genes, layout, params) {
  settings <- zero_sensitivity_settings(params)
  failed <- vapply(repeats, inherits, logical(1), what = "error")
  if (all(failed)) {
    stop("Every zero-count stabilization repeat failed: ", conditionMessage(repeats[[1]]))
  }
  if (any(failed)) {
    warning(
      sum(failed), " of ", length(repeats), " zero-count stabilization repeats failed: ",
      conditionMessage(repeats[[which(failed)[[1]]]])
    )
  }
  fits <- repeats[!failed]
  features <- fits[[1]]$features
  rows_by_gene <- split(seq_len(nrow(features)), factor(features$gene_id, levels = genes))
  control_ids <- layout$groups[[layout$control]]
  gene_success <- vapply(fits, function(fit) {
    vapply(genes, function(gene_id) {
      precision <- fit$precision[[gene_id]]
      is.finite(precision) && precision > 0 &&
        all(is.finite(fit$proportions[rows_by_gene[[gene_id]], , drop = FALSE]))
    }, logical(1))
  }, logical(length(genes)))
  if (is.null(dim(gene_success))) {
    gene_success <- matrix(gene_success, nrow = length(genes), dimnames = list(genes, NULL))
  }
  rownames(gene_success) <- genes
  enough <- rowSums(gene_success) >= settings$min_successes

  precision <- stats::setNames(rep(NA_real_, length(genes)), genes)
  proportions <- matrix(
    NA_real_, nrow(features), length(layout$sample_ids),
    dimnames = list(features$feature_id, layout$sample_ids)
  )
  for (gene_id in genes[enough]) {
    successful <- fits[gene_success[gene_id, ]]
    precision[[gene_id]] <- stats::median(vapply(successful, function(fit) fit$precision[[gene_id]], numeric(1)))
    rows <- rows_by_gene[[gene_id]]
    draws <- simplify2array(lapply(successful, function(fit) fit$proportions[rows, , drop = FALSE]))
    median_values <- apply(draws, c(1, 2), stats::median)
    if (is.null(dim(median_values))) median_values <- matrix(median_values, nrow = length(rows))
    proportions[rows, ] <- sweep(median_values, 2, colSums(median_values), "/")
  }

  gene_summary <- function(table_of) {
    lr <- vapply(fits, function(fit) {
      table <- table_of(fit)
      table$lr[match(genes, table$gene_id)]
    }, numeric(length(genes)))
    if (is.null(dim(lr))) lr <- matrix(lr, nrow = length(genes))
    # df is constant per gene and test; take it from any repeat that has it.
    df <- first_finite_by_row(vapply(fits, function(fit) {
      table <- table_of(fit)
      table$df[match(genes, table$gene_id)]
    }, numeric(length(genes))), length(genes))
    usable <- is.finite(lr) & gene_success
    successes <- rowSums(usable)
    median_lr <- vapply(seq_along(genes), function(index) {
      if (successes[[index]] < settings$min_successes) return(NA_real_)
      stats::median(lr[index, usable[index, ]])
    }, numeric(1))
    data.frame(
      gene_id = genes,
      lr = median_lr,
      df = df,
      pvalue = chi_square_p(median_lr, df),
      stabilization_successes = as.integer(successes),
      stringsAsFactors = FALSE
    )
  }

  omnibus <- gene_summary(function(fit) fit$omnibus)
  contrasts <- lapply(seq_len(nrow(layout$comparisons)), function(index) {
    comparison <- layout$comparisons$comparison[[index]]
    treatment_ids <- layout$groups[[layout$comparisons$treatment[[index]]]]
    genes_table <- gene_summary(function(fit) fit$contrasts[[comparison]]$genes)
    lr <- vapply(fits, function(fit) fit$contrasts[[comparison]]$features$lr, numeric(nrow(features)))
    delta <- vapply(fits, function(fit) {
      rowMeans(fit$proportions[, treatment_ids, drop = FALSE]) -
        rowMeans(fit$proportions[, control_ids, drop = FALSE])
    }, numeric(nrow(features)))
    if (is.null(dim(lr))) {
      lr <- matrix(lr, nrow = nrow(features))
      delta <- matrix(delta, nrow = nrow(features))
    }
    feature_gene_success <- gene_success[features$gene_id, , drop = FALSE]
    usable <- is.finite(lr) & is.finite(delta) & feature_gene_success
    successes <- rowSums(usable)
    reasons <- character(nrow(features))
    median_lr <- rep(NA_real_, nrow(features))
    spread <- rep(NA_real_, nrow(features))
    for (row in seq_len(nrow(features))) {
      values <- delta[row, usable[row, ]]
      if (length(values)) spread[[row]] <- max(values) - min(values)
      problems <- character()
      if (successes[[row]] < settings$min_successes) problems <- c(problems, "insufficient_successes")
      if (any(values > DIRECTION_TOLERANCE) && any(values < -DIRECTION_TOLERANCE)) {
        problems <- c(problems, "direction_change")
      }
      if (length(values) && spread[[row]] > settings$max_spread) problems <- c(problems, "delta_spread")
      reasons[[row]] <- paste(problems, collapse = ";")
      if (!length(problems)) median_lr[[row]] <- stats::median(lr[row, usable[row, ]])
    }
    df <- first_finite_by_row(vapply(fits, function(fit) {
      fit$contrasts[[comparison]]$features$df
    }, numeric(nrow(features))), nrow(features))
    list(
      genes = genes_table,
      features = data.frame(
        gene_id = features$gene_id,
        feature_id = features$feature_id,
        lr = median_lr,
        df = df,
        pvalue = chi_square_p(median_lr, df),
        stabilization_successes = as.integer(successes),
        stabilization_delta_pau_spread = spread,
        zero_boundary_unstable = nzchar(reasons),
        zero_boundary_reason = reasons,
        stringsAsFactors = FALSE
      )
    )
  })
  names(contrasts) <- layout$comparisons$comparison
  list(precision = precision, proportions = proportions, omnibus = omnibus, contrasts = contrasts)
}

initialize_status_columns <- function(fit) {
  fit$omnibus$stabilization_successes <- 0L
  fit$omnibus$model_status <- "drimseq"
  for (comparison in names(fit$contrasts)) {
    genes <- fit$contrasts[[comparison]]$genes
    genes$stabilization_successes <- 0L
    genes$model_status <- "drimseq"
    features <- fit$contrasts[[comparison]]$features
    features$stabilization_successes <- 0L
    features$stabilization_delta_pau_spread <- NA_real_
    features$zero_boundary_unstable <- FALSE
    features$zero_boundary_reason <- ""
    features$model_status <- "drimseq"
    fit$contrasts[[comparison]] <- list(genes = genes, features = features)
  }
  fit$gene_status <- stats::setNames(rep("drimseq", length(fit$gene_ids)), fit$gene_ids)
  fit
}

replace_rows <- function(target, source, key_columns, value_columns) {
  target_keys <- do.call(paste, c(target[key_columns], sep = "\r"))
  source_keys <- do.call(paste, c(source[key_columns], sep = "\r"))
  index <- match(source_keys, target_keys)
  if (anyNA(index)) stop("Stabilized results refer to rows absent from the family fit.")
  for (column in value_columns) target[[column]][index] <- source[[column]]
  target
}

# Stabilized genes take every value from the repeats, so each gene has one
# consistent source; perturbed counts themselves never reach an output.
write_back_boundary <- function(fit, summary, sets) {
  genes <- sets$perturbable
  if (length(genes)) {
    rows <- fit$features$gene_id %in% genes
    fit$proportions[rows, ] <- summary$proportions[fit$features$feature_id[rows], , drop = FALSE]
    fit$precision[genes] <- summary$precision[genes]
    gene_columns <- c("lr", "df", "pvalue", "stabilization_successes")
    summary$omnibus$model_status <- "drimseq_add_uniform"
    fit$omnibus <- replace_rows(fit$omnibus, summary$omnibus, "gene_id", c(gene_columns, "model_status"))
    for (comparison in names(fit$contrasts)) {
      stabilized <- summary$contrasts[[comparison]]
      stabilized$genes$model_status <- "drimseq_add_uniform"
      stabilized$features$model_status <- "drimseq_add_uniform"
      fit$contrasts[[comparison]]$genes <- replace_rows(
        fit$contrasts[[comparison]]$genes, stabilized$genes, "gene_id",
        c(gene_columns, "model_status")
      )
      fit$contrasts[[comparison]]$features <- replace_rows(
        fit$contrasts[[comparison]]$features, stabilized$features, c("gene_id", "feature_id"),
        c(
          "lr", "df", "pvalue", "stabilization_successes", "stabilization_delta_pau_spread",
          "zero_boundary_unstable", "zero_boundary_reason", "model_status"
        )
      )
    }
    fit$gene_status[genes] <- "drimseq_add_uniform"
  }
  unavailable <- sets$unavailable
  if (length(unavailable)) {
    fit$omnibus$model_status[fit$omnibus$gene_id %in% unavailable] <- "fit_unavailable"
    for (comparison in names(fit$contrasts)) {
      contrast <- fit$contrasts[[comparison]]
      contrast$genes$model_status[contrast$genes$gene_id %in% unavailable] <- "fit_unavailable"
      contrast$features$model_status[contrast$features$gene_id %in% unavailable] <- "fit_unavailable"
      fit$contrasts[[comparison]] <- contrast
    }
    fit$gene_status[unavailable] <- "fit_unavailable"
  }
  fit
}

# A condition with no counts for a gene has no usage to estimate, so every
# comparison involving it is left untested instead of fitted to noise.
mask_groups_without_counts <- function(fit, empty, layout) {
  empty <- empty[fit$gene_ids, , drop = FALSE]
  any_empty <- fit$gene_ids[rowSums(empty) > 0]
  if (!length(any_empty)) return(fit)
  omnibus_rows <- fit$omnibus$gene_id %in% any_empty
  fit$omnibus[omnibus_rows, c("lr", "pvalue")] <- NA_real_
  fit$omnibus$model_status[omnibus_rows] <- "group_without_counts"
  for (index in seq_len(nrow(layout$comparisons))) {
    comparison <- layout$comparisons$comparison[[index]]
    treatment <- layout$comparisons$treatment[[index]]
    masked <- fit$gene_ids[empty[, layout$control] | empty[, treatment]]
    if (!length(masked)) next
    contrast <- fit$contrasts[[comparison]]
    gene_rows <- contrast$genes$gene_id %in% masked
    contrast$genes[gene_rows, c("lr", "pvalue")] <- NA_real_
    contrast$genes$model_status[gene_rows] <- "group_without_counts"
    feature_rows <- contrast$features$gene_id %in% masked
    contrast$features[feature_rows, c("lr", "pvalue")] <- NA_real_
    contrast$features$model_status[feature_rows] <- "group_without_counts"
    fit$contrasts[[comparison]] <- contrast
  }
  fit$gene_status[any_empty] <- "group_without_counts"
  fit
}

# ---- Multiple testing -------------------------------------------------------

# Stage-wise adjustment within one comparison (design section 9). stageR joins
# gene and PAC IDs with ":" and splits them again, so it receives simple keys
# and the results are mapped back by key.
stage_adjust <- function(gene_results, feature_results, alpha) {
  if (!nrow(feature_results)) return(numeric())
  if (anyDuplicated(gene_results$gene_id)) stop("stageR input has duplicate genes.")
  gene_index <- match(feature_results$gene_id, gene_results$gene_id)
  if (anyNA(gene_index)) stop("stageR input has PACs without a gene-level result.")
  if (any(table(feature_results$gene_id) < 2L)) {
    stop("stageR input has genes with fewer than 2 PACs.")
  }
  screen <- as.numeric(gene_results$pvalue)
  screen[!is.finite(screen)] <- NA_real_
  if (!any(is.finite(screen))) return(rep(NA_real_, nrow(feature_results)))
  gene_keys <- paste0("g", seq_len(nrow(gene_results)))
  names(screen) <- gene_keys
  feature_keys <- paste0("t", seq_len(nrow(feature_results)))
  confirmation <- as.numeric(feature_results$pvalue)
  missing <- !is.finite(confirmation)
  confirmation[missing] <- 1
  confirmation <- matrix(confirmation, ncol = 1L, dimnames = list(feature_keys, "contrast"))
  tx2gene <- data.frame(
    txID = feature_keys,
    geneID = gene_keys[gene_index],
    stringsAsFactors = FALSE
  )
  object <- stageRTx(
    pScreen = screen,
    pConfirmation = confirmation,
    pScreenAdjusted = FALSE,
    tx2gene = tx2gene
  )
  object <- stageWiseAdjustment(object, method = "dtu", alpha = alpha, allowNA = TRUE)
  adjusted <- suppressMessages(getAdjustedPValues(
    object,
    order = FALSE,
    onlySignificantGenes = FALSE
  ))
  if (!"transcript" %in% names(adjusted)) {
    stop("stageR did not return transcript-level adjusted p-values.")
  }
  ids <- if ("txID" %in% names(adjusted)) as.character(adjusted$txID) else rownames(adjusted)
  index <- match(feature_keys, ids)
  if (any(is.na(index) & is.finite(screen[gene_index]))) {
    stop("stageR returned no adjusted p-value for PACs of screened genes.")
  }
  result <- as.numeric(adjusted$transcript[index])
  result[missing] <- NA_real_
  result
}

# ---- Events -----------------------------------------------------------------

classify_event <- function(row, params) {
  number <- function(value) suppressWarnings(as.numeric(value))
  finite_number <- function(value) {
    value <- number(value)
    length(value) == 1L && is.finite(value)
  }
  at_least <- function(value, threshold) {
    value <- number(value)
    length(value) == 1L && is.finite(value) && value >= threshold
  }
  at_most <- function(value, threshold) {
    value <- number(value)
    length(value) == 1L && is.finite(value) && value <= threshold
  }
  truth <- function(value) {
    value <- tolower(trimws(as.character(value)))
    length(value) == 1L && !is.na(value) && value %in% c("true", "t", "1")
  }
  control_detected <- at_least(
    row$control_supporting_samples,
    params$event_min_supporting_samples
  )
  control_support_available <- finite_number(row$control_supporting_samples)
  treatment_detected <- at_least(
    row$treatment_supporting_samples,
    params$event_min_supporting_samples
  )
  treatment_support_available <- finite_number(row$treatment_supporting_samples)
  positive <- at_least(row$delta_pau, params$min_abs_delta_pau)
  negative <- at_most(row$delta_pau, -params$min_abs_delta_pau)
  significant <- at_most(row$gene_fdr, params$gene_fdr) &&
    at_most(row$pac_fdr, params$site_fdr)
  stable <- !truth(row$zero_boundary_unstable)
  confidence <- as.character(row$confidence)
  confident <- length(confidence) == 1L &&
    !is.na(confidence) &&
    confidence != "low" &&
    !truth(row$internal_priming_flag) &&
    !truth(row$exploratory_insufficient_replicates)
  gained_detection <- control_support_available && treatment_support_available &&
    !control_detected && treatment_detected &&
    at_most(row$fitted_control_pau, params$event_max_control_pau) &&
    at_least(row$fitted_treatment_pau, params$event_min_treatment_pau) && positive
  lost_detection <- control_support_available && treatment_support_available &&
    control_detected && !treatment_detected &&
    at_most(row$fitted_treatment_pau, params$event_max_control_pau) &&
    at_least(row$fitted_control_pau, params$event_min_treatment_pau) && negative
  # A PAC's absence from a group counts only when that group has enough reads
  # at the gene to show it (design section 9).
  control_covered <- at_least(row$control_gene_total, params$min_gene_total)
  treatment_covered <- at_least(row$treatment_gene_total, params$min_gene_total)
  if (gained_detection) {
    if (significant && stable && confident && control_covered) "gained" else "gained_candidate"
  } else if (lost_detection) {
    if (significant && stable && confident && treatment_covered) "lost" else "lost_candidate"
  } else if (control_detected && treatment_detected && significant && positive) {
    "increased_usage"
  } else if (control_detected && treatment_detected && significant && negative) {
    "decreased_usage"
  } else {
    "none"
  }
}

# Classify each row from typed column values rather than a character matrix.
classify_table_events <- function(table, params) {
  fields <- c(
    "control_supporting_samples", "treatment_supporting_samples", "delta_pau",
    "gene_fdr", "pac_fdr", "zero_boundary_unstable", "confidence",
    "internal_priming_flag", "exploratory_insufficient_replicates",
    "fitted_control_pau", "fitted_treatment_pau",
    "control_gene_total", "treatment_gene_total"
  )
  columns <- lapply(stats::setNames(fields, fields), function(field) table[[field]])
  vapply(
    seq_len(nrow(table)),
    function(index) classify_event(lapply(columns, `[[`, index), params),
    character(1)
  )
}

detection_candidates <- function(pacs, params) {
  gained <- pacs$fitted_control_pau <= params$event_max_control_pau &
    pacs$fitted_treatment_pau >= params$event_min_treatment_pau &
    pacs$treatment_supporting_samples >= params$event_min_supporting_samples &
    pacs$delta_pau >= params$min_abs_delta_pau
  lost <- pacs$fitted_treatment_pau <= params$event_max_control_pau &
    pacs$fitted_control_pau >= params$event_min_treatment_pau &
    pacs$control_supporting_samples >= params$event_min_supporting_samples &
    pacs$delta_pau <= -params$min_abs_delta_pau
  gained | lost
}

assign_events <- function(pacs, params) {
  row_groups <- split(seq_len(nrow(pacs)), pacs$gene_id)
  dominant_control <- vapply(row_groups, function(indices) {
    dominant_pac(pacs$feature_id[indices], pacs$fitted_control_pau[indices])
  }, character(1))
  dominant_treatment <- vapply(row_groups, function(indices) {
    dominant_pac(pacs$feature_id[indices], pacs$fitted_treatment_pau[indices])
  }, character(1))
  pacs$dominant_pac_control <- unname(dominant_control[pacs$gene_id])
  pacs$dominant_pac_treatment <- unname(dominant_treatment[pacs$gene_id])
  pacs$control_detected_complexity <- ave(
    pacs$fitted_control_pau >= params$event_min_treatment_pau,
    pacs$gene_id,
    FUN = function(values) sum(values, na.rm = TRUE)
  )
  pacs$treatment_detected_complexity <- ave(
    pacs$fitted_treatment_pau >= params$event_min_treatment_pau,
    pacs$gene_id,
    FUN = function(values) sum(values, na.rm = TRUE)
  )
  events <- classify_table_events(pacs, params)
  exploratory <- as.logical(pacs$exploratory_insufficient_replicates)
  events[exploratory & events == "gained"] <- "gained_candidate"
  events[exploratory & events == "lost"] <- "lost_candidate"
  # The descriptive labels below describe fitted usage, so they are given
  # only in genes that pass the comparison's gene-level screen.
  screened <- !is.na(pacs$gene_fdr) & pacs$gene_fdr <= params$gene_fdr
  dominant_switch <- !is.na(pacs$dominant_pac_control) &
    !is.na(pacs$dominant_pac_treatment) &
    pacs$dominant_pac_control != pacs$dominant_pac_treatment
  events[events == "none" & screened & dominant_switch &
    pacs$feature_id == pacs$dominant_pac_treatment] <- "dominant_switch"
  complexity_delta <- pacs$treatment_detected_complexity - pacs$control_detected_complexity
  # A gene without fitted usage in either group (a condition with no counts)
  # has no usage pattern to compare, so it gets no descriptive event.
  unfitted <- ave(
    !is.finite(pacs$fitted_control_pau) | !is.finite(pacs$fitted_treatment_pau),
    pacs$gene_id,
    FUN = any
  )
  complexity_delta[as.logical(unfitted)] <- 0
  events[events == "none" & screened & complexity_delta > 0] <- "complexity_gain"
  events[events == "none" & screened & complexity_delta < 0] <- "complexity_loss"
  pacs$event_type <- events
  pacs
}

# ---- Per-comparison tables --------------------------------------------------

comparison_outputs <- function(fit, comparison_row, counts, layout, atlas, empty, params, gene_names) {
  comparison <- comparison_row$comparison
  treatment <- comparison_row$treatment
  control <- layout$control
  contrast <- fit$contrasts[[comparison]]
  genes <- contrast$genes
  genes$gene_name <- unname(gene_names[genes$gene_id])
  genes$gene_fdr <- bh(genes$pvalue)
  features <- contrast$features
  features$pac_fdr <- stage_adjust(genes, features, params$site_fdr)
  gene_index <- match(features$gene_id, genes$gene_id)
  control_pau <- rowMeans(fit$proportions[, layout$groups[[control]], drop = FALSE])
  treatment_pau <- rowMeans(fit$proportions[, layout$groups[[treatment]], drop = FALSE])
  masked <- empty[features$gene_id, control] | empty[features$gene_id, treatment]
  control_pau[masked] <- NA_real_
  treatment_pau[masked] <- NA_real_
  precision <- unname(fit$precision[features$gene_id])
  annotation_index <- match(features$feature_id, atlas$pac_id)
  if (anyNA(annotation_index)) stop("The atlas lacks annotation for tested PACs in ", comparison, ".")
  annotation <- select_columns(
    atlas, c(LOCATION_COLUMNS, ATLAS_ANNOTATION_COLUMNS), "The atlas"
  )[annotation_index, , drop = FALSE]
  rownames(annotation) <- NULL
  pacs <- data.frame(
    feature_id = features$feature_id,
    gene_id = features$gene_id,
    gene_name = unname(gene_names[features$gene_id]),
    lr = features$lr,
    df = features$df,
    pvalue_pac = features$pvalue,
    pac_fdr = features$pac_fdr,
    pvalue_gene = genes$pvalue[gene_index],
    gene_fdr = genes$gene_fdr[gene_index],
    fitted_control_pau = unname(control_pau),
    fitted_treatment_pau = unname(treatment_pau),
    delta_pau = unname(treatment_pau - control_pau),
    control_supporting_samples = unname(supporting_samples(counts, layout$samples, control)),
    treatment_supporting_samples = unname(supporting_samples(counts, layout$samples, treatment)),
    control_gene_total = group_gene_totals(counts, layout$samples, control),
    treatment_gene_total = group_gene_totals(counts, layout$samples, treatment),
    raw_control_counts = unname(raw_count_text(counts, layout$samples, control)),
    raw_treatment_counts = unname(raw_count_text(counts, layout$samples, treatment)),
    observed_control_pau = unname(observed_pau_text(counts, layout$samples, control)),
    observed_treatment_pau = unname(observed_pau_text(counts, layout$samples, treatment)),
    model_status = features$model_status,
    stabilization_successes = features$stabilization_successes,
    stabilization_delta_pau_spread = features$stabilization_delta_pau_spread,
    zero_boundary_unstable = features$zero_boundary_unstable,
    zero_boundary_reason = features$zero_boundary_reason,
    precision = precision,
    alpha_control = unname(control_pau) * precision,
    alpha_treatment = unname(treatment_pau) * precision,
    annotation,
    pac_id = features$feature_id,
    condition = treatment,
    control_condition = control,
    exploratory_insufficient_replicates = layout$exploratory,
    delta_pau_ci_low = NA_real_,
    delta_pau_ci_high = NA_real_,
    bootstrap_successes = 0L,
    bootstrap_perturbed = 0L,
    bootstrap_status = "not_selected",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  pacs$effect_exceeds_threshold <- abs(pacs$delta_pau) >= params$min_abs_delta_pau
  pacs <- assign_events(pacs, params)
  genes$condition <- treatment
  genes$control_condition <- control
  genes$exploratory_insufficient_replicates <- layout$exploratory
  list(
    comparison = comparison,
    genes = genes,
    pacs = pacs,
    selected = bootstrap_gene_ids(
      pacs$gene_id,
      pacs$gene_fdr,
      detection_candidates(pacs, params),
      params$gene_fdr,
      params$dm_bootstrap_include_candidates
    )
  )
}

bootstrap_eligible_genes <- function(fit) {
  rows <- split(seq_len(nrow(fit$features)), factor(fit$features$gene_id, levels = fit$gene_ids))
  finite <- vapply(rows, function(index) all(is.finite(fit$proportions[index, , drop = FALSE])), logical(1))
  fit$gene_ids[finite & is.finite(fit$precision) & fit$precision > 0]
}

family_bootstrap_batches <- function(selections, counts, fit, layout, params, atlas_checksum, batch_size) {
  genes <- fit$gene_ids[fit$gene_ids %in% selections$gene_id]
  if (!length(genes)) return(list())
  chunks <- split(genes, ceiling(seq_along(genes) / batch_size))
  lapply(seq_along(chunks), function(index) {
    gene_ids <- chunks[[index]]
    rows <- counts$gene_id %in% gene_ids
    fitted <- data.frame(
      gene_id = counts$gene_id[rows],
      pac_id = counts$pac_id[rows],
      fit$proportions[rows, , drop = FALSE],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    rownames(fitted) <- NULL
    batch_counts <- counts[rows, c("gene_id", "pac_id", layout$sample_ids), drop = FALSE]
    rownames(batch_counts) <- NULL
    list(
      empty = FALSE,
      batch_id = sprintf("%s.batch-%03d", safe_file_component(layout$family), index),
      family = layout$family,
      control = layout$control,
      comparisons = layout$comparisons[, c("comparison", "treatment")],
      samples = layout$samples[, c("sample_id", "condition", layout$covariates), drop = FALSE],
      design = layout$design,
      params = params,
      atlas_checksum = atlas_checksum,
      counts = batch_counts,
      fitted = fitted,
      genes = data.frame(
        gene_id = gene_ids,
        precision = unname(fit$precision[gene_ids]),
        model_status = unname(fit$gene_status[gene_ids]),
        stringsAsFactors = FALSE
      ),
      selections = selections[selections$gene_id %in% gene_ids, , drop = FALSE]
    )
  })
}

empty_bootstrap_batch <- function(family) {
  list(
    empty = TRUE,
    batch_id = paste0(safe_file_component(family), ".batch-000"),
    family = family
  )
}

write_bootstrap_batches <- function(batches, output_dir, family) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!length(batches)) batches <- list(empty_bootstrap_batch(family))
  invisible(vapply(batches, function(batch) {
    path <- file.path(output_dir, paste0(batch$batch_id, ".rds"))
    saveRDS(batch, path, compress = FALSE)
    path
  }, character(1)))
}

write_empty_family_outputs <- function(layout, output_dir) {
  write_gzip_tsv(empty_table(OMNIBUS_COLUMNS), file.path(output_dir, paste0(layout$family, ".gene_omnibus.tsv.gz")))
  for (comparison in layout$comparisons$comparison) {
    write_gzip_tsv(empty_table(GENE_COLUMNS), file.path(output_dir, paste0(comparison, ".genes.tsv.gz")))
    write_gzip_tsv(empty_table(PAC_COLUMNS), file.path(output_dir, paste0(comparison, ".pacs.tsv.gz")))
  }
}

# ---- Fit mode ---------------------------------------------------------------

# Gene names by gene ID, from the atlas rows assigned to a single gene. The
# annotation step gives each gene one name, falling back to its ID.
atlas_gene_names <- function(atlas) {
  single <- !is.na(atlas$gene_id) & nzchar(atlas$gene_id) & !grepl(",", atlas$gene_id, fixed = TRUE)
  names <- unique(atlas[single, c("gene_id", "gene_name"), drop = FALSE])
  repeated <- names$gene_id[duplicated(names$gene_id)]
  if (length(repeated)) stop("Gene ", repeated[[1]], " has more than one gene_name in the atlas.")
  stats::setNames(as.character(names$gene_name), names$gene_id)
}

fit_family <- function(
  family,
  sample_rows,
  counts,
  atlas,
  params,
  output_dir,
  atlas_checksum,
  model_workers,
  bootstrap_batch_size
) {
  layout <- family_layout(family, sample_rows, params)
  if (is.null(layout)) return(NULL)
  empty_result <- list(
    layout = layout,
    precision = empty_table(GENE_PRECISION_COLUMNS),
    fitted_pau = empty_table(FITTED_PAU_COLUMNS),
    batches = list()
  )
  gene_names <- atlas_gene_names(atlas)
  filtering <- family_filter(counts, layout$sample_ids, params, family)
  write_gzip_tsv(
    select_columns(filtering$reasons, FILTERING_COLUMNS, "The filtering table"),
    file.path(output_dir, paste0(family, ".statistical_filtering.tsv.gz"))
  )
  filtered <- filtering$counts
  if (!nrow(filtered)) {
    warning("No testable genes in comparison family ", family)
    write_empty_family_outputs(layout, output_dir)
    return(empty_result)
  }
  filtered <- order_family_counts(filtered, atlas)
  features <- data.frame(
    gene_id = filtered$gene_id,
    feature_id = filtered$pac_id,
    stringsAsFactors = FALSE
  )
  bpparam <- drimseq_bpparam(
    model_workers,
    stable_seed(params$random_seed, atlas_checksum, family, "family_model")
  )

  # 1. Unmodified fit. dmPrecision estimates the common precision from a
  # random subset of genes, so the draw is seeded for reproducibility.
  set.seed(stable_seed(params$random_seed, atlas_checksum, family, "family_precision"))
  precision <- estimate_family_precision(filtered, layout, bpparam)
  fit <- fit_and_test(precision, layout, bpparam, features)
  fit <- initialize_status_columns(fit)

  # 2. Zero-count stabilization, once for the whole family.
  sets <- boundary_gene_sets(fit, filtered, layout)
  if (length(sets$perturbable)) {
    message(
      "Comparison family ", family, ": stabilizing ", length(sets$perturbable),
      " of ", length(fit$gene_ids), " gene(s) with zero-count groups."
    )
    repeats <- run_zero_sensitivity(
      filtered, sets$perturbable, layout, precision@common_precision,
      params, atlas_checksum, bpparam
    )
    summary <- summarize_zero_sensitivity(repeats, sets$perturbable, layout, params)
    fit <- write_back_boundary(fit, summary, sets)
  } else if (length(sets$unavailable)) {
    fit <- write_back_boundary(fit, NULL, sets)
  }
  empty <- empty_condition_groups(filtered, layout)
  fit <- mask_groups_without_counts(fit, empty, layout)

  # 3. Tables.
  omnibus <- fit$omnibus
  omnibus$gene_name <- unname(gene_names[omnibus$gene_id])
  omnibus$gene_fdr <- bh(omnibus$pvalue)
  omnibus$family <- family
  omnibus$exploratory_insufficient_replicates <- layout$exploratory
  write_gzip_tsv(
    select_columns(omnibus, OMNIBUS_COLUMNS, "The omnibus table"),
    file.path(output_dir, paste0(family, ".gene_omnibus.tsv.gz"))
  )

  outputs <- lapply(seq_len(nrow(layout$comparisons)), function(index) {
    comparison_outputs(
      fit, layout$comparisons[index, , drop = FALSE], filtered, layout, atlas, empty, params,
      gene_names
    )
  })
  settings <- bootstrap_settings(params)
  eligible <- bootstrap_eligible_genes(fit)
  selection_rows <- lapply(outputs, function(output) {
    genes <- output$selected[output$selected %in% eligible]
    if (!length(genes)) return(NULL)
    data.frame(gene_id = genes, comparison = output$comparison, stringsAsFactors = FALSE)
  })
  selections <- do.call(rbind, selection_rows)
  if (is.null(selections) || settings$replicates == 0L) {
    selections <- data.frame(gene_id = character(), comparison = character(), stringsAsFactors = FALSE)
  }
  fitted_rows <- list()
  for (output in outputs) {
    pacs <- output$pacs
    if (settings$replicates == 0L) {
      pacs$bootstrap_status <- "disabled"
    } else {
      pacs$bootstrap_status[pacs$gene_id %in% output$selected] <- "fit_unavailable"
      pending <- pacs$gene_id %in% selections$gene_id[selections$comparison == output$comparison]
      pacs$bootstrap_status[pending] <- "pending"
    }
    write_gzip_tsv(
      select_columns(output$genes, GENE_COLUMNS, "The genes table"),
      file.path(output_dir, paste0(output$comparison, ".genes.tsv.gz"))
    )
    write_gzip_tsv(
      select_columns(pacs, PAC_COLUMNS, "The PAC table"),
      file.path(output_dir, paste0(output$comparison, ".pacs.tsv.gz"))
    )
    fitted_rows[[output$comparison]] <- select_columns(pacs, FITTED_PAU_COLUMNS, "The fitted PAU table")
  }
  message(
    "Comparison family ", family, ": bootstrapping ", length(unique(selections$gene_id)),
    " gene(s) with ", settings$replicates, " replicate(s) each."
  )
  precision_table <- data.frame(
    gene_id = fit$gene_ids,
    gene_name = unname(gene_names[fit$gene_ids]),
    precision = unname(fit$precision),
    family = family,
    model_status = unname(fit$gene_status),
    stringsAsFactors = FALSE
  )
  list(
    layout = layout,
    precision = precision_table,
    fitted_pau = do.call(rbind, unname(fitted_rows)),
    batches = family_bootstrap_batches(
      selections, filtered, fit, layout, params, atlas_checksum, bootstrap_batch_size
    )
  )
}

# Motif usage (design section 9, "Motif Preference"). Each comparison is one
# limma fit on the family's main-model design (covariates, then condition with
# the family control as reference). Its rows are the score table's keys: a
# primary hexamer with its class, or a class alone. A row is tested in a
# comparison only when every sample of that comparison has at least
# motif_preference_min_genes informative genes for it.
fit_motif_preferences <- function(
  scores, layout, params, output_dir, suffix, keys = MOTIF_PREFERENCE_KEYS
) {
  scores <- scores[scores$sample_id %in% layout$sample_ids, , drop = FALSE]
  if (!nrow(scores)) return(invisible(NULL))
  rows <- unique(scores[, keys, drop = FALSE])
  # Classes first, then hexamers within a class.
  rows <- rows[do.call(order, c(unname(as.list(rows[rev(keys)])), method = "radix")), , drop = FALSE]
  rownames(rows) <- NULL
  key_text <- function(table) do.call(paste, c(unname(as.list(table[keys])), sep = "\r"))
  row_keys <- key_text(rows)
  score_cells <- cbind(match(key_text(scores), row_keys), match(scores$sample_id, layout$sample_ids))
  by_row <- function(column) {
    values <- matrix(NA_real_, length(row_keys), length(layout$sample_ids))
    colnames(values) <- layout$sample_ids
    values[score_cells] <- as.numeric(scores[[column]])
    values
  }
  usage <- by_row("motif_usage")
  transformed <- by_row("transformed_motif_usage")
  informative <- by_row("informative_genes")
  columns <- c(keys, setdiff(MOTIF_PREFERENCE_COLUMNS, MOTIF_PREFERENCE_KEYS))
  control_ids <- layout$groups[[layout$control]]
  for (index in seq_len(nrow(layout$comparisons))) {
    treatment <- layout$comparisons$treatment[[index]]
    treatment_ids <- layout$groups[[treatment]]
    comparison_ids <- c(control_ids, treatment_ids)
    minimum_genes <- apply(informative[, comparison_ids, drop = FALSE], 1, min)
    testable <- which(!is.na(minimum_genes) & minimum_genes >= params$motif_preference_min_genes)
    if (!length(testable)) next
    fit <- limma::eBayes(limma::lmFit(transformed[testable, , drop = FALSE], layout$design))
    table <- limma::topTable(
      fit, coef = layout$comparisons$coefficient[[index]], number = Inf, sort.by = "none"
    )
    control_mean <- rowMeans(usage[testable, control_ids, drop = FALSE])
    treatment_mean <- rowMeans(usage[testable, treatment_ids, drop = FALSE])
    output <- data.frame(
      rows[testable, keys, drop = FALSE],
      condition = treatment,
      control_condition = layout$control,
      control_mean = control_mean,
      treatment_mean = treatment_mean,
      delta_motif_usage = treatment_mean - control_mean,
      transformed_coefficient = table$logFC,
      p_value = table$P.Value,
      fdr = bh(table$P.Value),
      informative_genes = minimum_genes[testable],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    rownames(output) <- NULL
    write_gzip_tsv(
      select_columns(output, columns, "The motif preference table"),
      file.path(output_dir, paste0(layout$comparisons$comparison[[index]], ".", suffix, ".tsv.gz"))
    )
  }
  invisible(NULL)
}

run_fit_mode <- function(arguments) {
  require_args(arguments, c(
    "family", "samples", "counts", "atlas", "params", "motif_scores",
    "motif_sensitivity", "motif_class_scores", "motif_class_sensitivity", "output_dir"
  ))
  model_workers <- parse_model_workers(arguments$model_workers)
  batch_size <- parse_batch_size(arguments$bootstrap_batch_size)
  dir.create(arguments$output_dir, recursive = TRUE, showWarnings = FALSE)
  preliminary_dir <- file.path(arguments$output_dir, "preliminary")
  batch_dir <- file.path(arguments$output_dir, "bootstrap-batches")
  dir.create(preliminary_dir, recursive = TRUE, showWarnings = FALSE)
  params <- yaml::read_yaml(arguments$params)
  samples <- read_tsv(arguments$samples)
  counts <- read_tsv(arguments$counts)
  atlas <- read_tsv(arguments$atlas)

  available_families <- unique(samples$control_condition[
    samples$condition != samples$control_condition
  ])
  if (!arguments$family %in% available_families) {
    stop(
      "Requested comparison family ", arguments$family,
      " is absent from the normalized sample sheet."
    )
  }
  atlas_checksum <- unname(tools::md5sum(arguments$atlas))
  fitted <- fit_family(
    arguments$family,
    samples,
    counts,
    atlas,
    params,
    preliminary_dir,
    atlas_checksum,
    model_workers,
    batch_size
  )
  if (is.null(fitted)) stop("Comparison family ", arguments$family, " has no treatments.")
  write_gzip_tsv(
    select_columns(fitted$precision, GENE_PRECISION_COLUMNS, "The precision table"),
    file.path(preliminary_dir, "gene_precision.tsv.gz")
  )
  write_gzip_tsv(
    select_columns(fitted$fitted_pau, FITTED_PAU_COLUMNS, "The fitted PAU table"),
    file.path(preliminary_dir, "fitted_pau.tsv.gz")
  )
  write_bootstrap_batches(fitted$batches, batch_dir, arguments$family)
  fit_motif_preferences(
    read_tsv(arguments$motif_scores), fitted$layout, params, preliminary_dir, "preference"
  )
  fit_motif_preferences(
    read_tsv(arguments$motif_sensitivity), fitted$layout, params, preliminary_dir,
    "preference_known_rescue_sensitivity"
  )
  fit_motif_preferences(
    read_tsv(arguments$motif_class_scores), fitted$layout, params, preliminary_dir,
    "preference_class", keys = "primary_motif_class"
  )
  fit_motif_preferences(
    read_tsv(arguments$motif_class_sensitivity), fitted$layout, params, preliminary_dir,
    "preference_class_known_rescue_sensitivity", keys = "primary_motif_class"
  )
}

# ---- Versions mode ----------------------------------------------------------

# Appends R and statistics-package versions to a software-versions table,
# creating it with a header when it does not exist yet.
run_versions_mode <- function(arguments, versions) {
  require_args(arguments, "output")
  rows <- data.frame(
    software = c("R", names(versions)),
    version = c(paste(R.version$major, R.version$minor, sep = "."), unname(versions)),
    stringsAsFactors = FALSE
  )
  exists <- file.exists(arguments$output)
  utils::write.table(
    rows, arguments$output, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = !exists, append = exists
  )
}

# ---- Bootstrap mode ---------------------------------------------------------

validate_bootstrap_batch <- function(batch) {
  required <- c(
    "batch_id", "family", "control", "comparisons", "samples", "design", "params",
    "atlas_checksum", "counts", "fitted", "genes", "selections"
  )
  missing <- required[!required %in% names(batch)]
  if (length(missing)) {
    stop("Bootstrap batch is missing required fields: ", paste(missing, collapse = ", "), ".")
  }
  for (field in c("batch_id", "family", "control", "atlas_checksum")) {
    value <- batch[[field]]
    if (length(value) != 1L || is.na(value) || !nzchar(as.character(value))) {
      stop("Bootstrap batch has an invalid ", field, ".")
    }
  }
  label <- paste("Bootstrap batch", batch$batch_id)
  samples <- batch$samples
  if (!is.data.frame(samples) || anyDuplicated(samples$sample_id)) {
    stop(label, " has invalid samples.")
  }
  if (!is.matrix(batch$design) || nrow(batch$design) != nrow(samples)) {
    stop(label, " has an invalid design matrix.")
  }
  sample_ids <- as.character(samples$sample_id)
  for (table_name in c("counts", "fitted")) {
    table <- batch[[table_name]]
    if (!is.data.frame(table) || !identical(names(table), c("gene_id", "pac_id", sample_ids))) {
      stop(label, " has ", table_name, " columns that do not match its samples.")
    }
  }
  if (!identical(batch$counts$gene_id, batch$fitted$gene_id) ||
    !identical(batch$counts$pac_id, batch$fitted$pac_id)) {
    stop(label, " has counts and fitted proportions for different PACs.")
  }
  genes <- batch$genes
  if (!is.data.frame(genes) || !identical(as.character(genes$gene_id), unique(batch$counts$gene_id))) {
    stop(label, " lists genes that do not match its counts.")
  }
  if (any(table(batch$counts$gene_id) < 2L)) stop(label, " has genes with fewer than 2 PACs.")
  conditions <- levels(samples$condition)
  if (is.null(conditions) || !batch$control %in% conditions ||
    !all(batch$comparisons$treatment %in% conditions)) {
    stop(label, " has conditions that do not match its samples.")
  }
  selections <- batch$selections
  if (!is.data.frame(selections) ||
    !all(selections$gene_id %in% genes$gene_id) ||
    !all(selections$comparison %in% batch$comparisons$comparison) ||
    anyDuplicated(paste(selections$gene_id, selections$comparison, sep = "\r")) ||
    !all(genes$gene_id %in% selections$gene_id)) {
    stop(label, " has invalid gene selections.")
  }
  invisible(TRUE)
}

bootstrap_state <- function(batch) {
  samples <- batch$samples
  sample_ids <- as.character(samples$sample_id)
  genes <- as.character(batch$genes$gene_id)
  rows <- split(seq_len(nrow(batch$counts)), factor(batch$counts$gene_id, levels = genes))
  counts <- count_matrix(batch$counts, sample_ids)
  fitted <- count_matrix(batch$fitted, sample_ids)
  precision <- stats::setNames(as.numeric(batch$genes$precision), genes)
  totals <- lapply(genes, function(gene_id) {
    vapply(sample_ids, function(sample_id) {
      bootstrap_total(counts[rows[[gene_id]], sample_id], gene_id, sample_id)
    }, integer(1))
  })
  names(totals) <- genes
  shapes <- list()
  available <- stats::setNames(rep(FALSE, length(genes)), genes)
  for (gene_id in genes) {
    expected <- fitted[rows[[gene_id]], , drop = FALSE]
    usable <- is.finite(precision[[gene_id]]) && precision[[gene_id]] > 0 &&
      all(is.finite(expected)) && all(expected >= 0) &&
      all(abs(colSums(expected) - 1) <= 1e-6)
    available[[gene_id]] <- usable
    if (usable) shapes[[gene_id]] <- pmax(expected, 1e-10) * precision[[gene_id]]
  }
  covariates <- setdiff(names(samples), c("sample_id", "condition"))
  zero_levels <- c(
    list(condition = samples$condition),
    stats::setNames(lapply(covariates, function(covariate) samples[[covariate]]), covariates)
  )
  groups <- lapply(levels(samples$condition), function(condition) which(samples$condition == condition))
  names(groups) <- levels(samples$condition)
  list(
    batch = batch,
    genes = genes,
    rows = rows,
    feature_ids = as.character(batch$counts$pac_id),
    sample_ids = sample_ids,
    totals = totals,
    shapes = shapes,
    precision = precision,
    available = available,
    zero_levels = zero_levels,
    groups = groups,
    dm_samples = samples,
    design = batch$design,
    seed_parts = list(batch$params$random_seed, batch$atlas_checksum, batch$family)
  )
}

# One Dirichlet-multinomial draw per sample at the observed gene total.
simulate_gene_counts <- function(shapes, totals) {
  result <- matrix(0, nrow(shapes), ncol(shapes))
  for (column in seq_len(ncol(shapes))) {
    draw <- stats::rgamma(nrow(shapes), shape = shapes[, column], rate = 1)
    draw_total <- sum(draw)
    if (!is.finite(draw_total) || draw_total <= 0) return(NULL)
    if (totals[[column]] > 0L) {
      result[, column] <- stats::rmultinom(1L, totals[[column]], draw / draw_total)[, 1]
    }
  }
  result
}

has_zero_group <- function(values, zero_levels) {
  for (levels in zero_levels) {
    labels <- as.character(levels)
    for (level in unique(labels)) {
      if (any(rowSums(values[, labels == level, drop = FALSE]) == 0)) return(TRUE)
    }
  }
  FALSE
}

# Refit simulated genes with precision held at the family estimate. dmFit fits
# genes independently, so fitting them together changes no gene's result.
fit_fixed_precision_batch <- function(values, state) {
  gene_ids <- names(values)
  if (!length(gene_ids)) return(list())
  fit_genes <- function(ids) {
    table <- do.call(rbind, lapply(ids, function(gene_id) {
      matrix_values <- values[[gene_id]]
      colnames(matrix_values) <- state$sample_ids
      data.frame(
        gene_id = gene_id,
        pac_id = state$feature_ids[state$rows[[gene_id]]],
        matrix_values,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }))
    data <- dm_data(table, state$sample_ids, state$dm_samples)
    expression <- vapply(ids, function(gene_id) mean(colSums(values[[gene_id]])), numeric(1))
    object <- fixed_precision_object(
      data, state$design, state$precision, stats::median(state$precision[ids]), expression
    )
    fit <- DRIMSeq::dmFit(object, design = state$design, bb_model = FALSE, verbose = 0)
    proportions_table <- DRIMSeq::proportions(fit)
    result <- lapply(ids, function(gene_id) {
      rows <- proportions_table$gene_id == gene_id
      proportions <- as.matrix(proportions_table[rows, state$sample_ids, drop = FALSE])
      index <- match(state$feature_ids[state$rows[[gene_id]]], proportions_table$feature_id[rows])
      proportions[index, , drop = FALSE]
    })
    names(result) <- ids
    result
  }
  batched <- tryCatch(fit_genes(gene_ids), error = function(error) error)
  if (!inherits(batched, "error")) return(batched)
  # A DRIMSeq error in the joint call is retried gene by gene, so one gene
  # cannot fail the whole replicate.
  result <- lapply(gene_ids, function(gene_id) {
    single <- tryCatch(fit_genes(gene_id), error = function(error) NULL)
    if (is.null(single)) NULL else single[[gene_id]]
  })
  names(result) <- gene_ids
  result
}

finite_fit <- function(proportions) !is.null(proportions) && all(is.finite(proportions))

bootstrap_replicate <- function(replicate_number, state) {
  seed <- function(gene_id, purpose) {
    do.call(stable_seed, c(state$seed_parts, list(gene_id, purpose, replicate_number)))
  }
  genes <- state$genes[state$available]
  values <- list()
  for (gene_id in genes) {
    set.seed(seed(gene_id, "bootstrap"))
    simulated <- simulate_gene_counts(state$shapes[[gene_id]], state$totals[[gene_id]])
    if (!is.null(simulated)) values[[gene_id]] <- simulated
  }
  # Zero groups are perturbed before fitting, exactly as in the family fit.
  perturbed <- names(values)[vapply(values, has_zero_group, logical(1), zero_levels = state$zero_levels)]
  for (gene_id in perturbed) {
    values[[gene_id]] <- perturb_zero_cells(values[[gene_id]], seed(gene_id, "bootstrap_zero"))
  }
  fits <- fit_fixed_precision_batch(values, state)
  retry <- names(values)[vapply(names(values), function(gene_id) {
    !finite_fit(fits[[gene_id]]) && !gene_id %in% perturbed && any(values[[gene_id]] == 0)
  }, logical(1))]
  if (length(retry)) {
    for (gene_id in retry) {
      values[[gene_id]] <- perturb_zero_cells(values[[gene_id]], seed(gene_id, "bootstrap_zero"))
    }
    refits <- fit_fixed_precision_batch(values[retry], state)
    fits[retry] <- refits[retry]
    perturbed <- c(perturbed, retry)
  }
  result <- lapply(state$genes, function(gene_id) {
    proportions <- fits[[gene_id]]
    if (!finite_fit(proportions)) return(list(success = FALSE, perturbed = FALSE))
    means <- vapply(state$groups, function(columns) {
      rowMeans(proportions[, columns, drop = FALSE])
    }, numeric(nrow(proportions)))
    if (is.null(dim(means))) means <- matrix(means, nrow = nrow(proportions))
    colnames(means) <- names(state$groups)
    list(success = TRUE, perturbed = gene_id %in% perturbed, means = means)
  })
  names(result) <- state$genes
  result
}

bootstrap_interval_rows <- function(state, gene_id, comparison, treatment, status, replicate_results, settings) {
  feature_ids <- state$feature_ids[state$rows[[gene_id]]]
  rows <- data.frame(
    comparison = comparison,
    gene_id = gene_id,
    feature_id = feature_ids,
    delta_pau_ci_low = NA_real_,
    delta_pau_ci_high = NA_real_,
    bootstrap_successes = 0L,
    bootstrap_perturbed = 0L,
    bootstrap_status = status,
    stringsAsFactors = FALSE
  )
  if (!identical(status, "ok")) return(rows)
  control <- state$batch$control
  draws <- lapply(replicate_results, `[[`, gene_id)
  successful <- Filter(function(draw) isTRUE(draw$success), draws)
  rows$bootstrap_successes <- length(successful)
  rows$bootstrap_perturbed <- sum(vapply(successful, function(draw) isTRUE(draw$perturbed), logical(1)))
  if (length(successful) < settings$min_successes) {
    rows$bootstrap_status <- "insufficient_successes"
    return(rows)
  }
  deltas <- vapply(successful, function(draw) {
    draw$means[, treatment] - draw$means[, control]
  }, numeric(length(feature_ids)))
  if (is.null(dim(deltas))) deltas <- matrix(deltas, nrow = length(feature_ids))
  rows$delta_pau_ci_low <- apply(deltas, 1, stats::quantile, probs = 0.025, names = FALSE, type = 7)
  rows$delta_pau_ci_high <- apply(deltas, 1, stats::quantile, probs = 0.975, names = FALSE, type = 7)
  rows
}

run_bootstrap_batch <- function(batch_path, output_path, bootstrap_workers) {
  batch <- readRDS(batch_path)
  if (!is.list(batch)) {
    stop("Bootstrap batch payload must be a list.")
  }
  if (isTRUE(batch$empty)) {
    write_gzip_tsv(empty_table(INTERVAL_COLUMNS), output_path)
    return(invisible(NULL))
  }
  validate_bootstrap_batch(batch)
  settings <- bootstrap_settings(batch$params)
  state <- bootstrap_state(batch)
  replicate_results <- list()
  if (settings$replicates > 0L && any(state$available)) {
    message(
      "Bootstrap batch ", batch$batch_id, ": ", sum(state$available), " gene(s), ",
      settings$replicates, " replicate(s)."
    )
    replicate_results <- bootstrap_apply(
      seq_len(settings$replicates),
      bootstrap_workers,
      function(replicate_number) bootstrap_replicate(replicate_number, state)
    )
    broken <- vapply(replicate_results, function(result) {
      is.null(result) || inherits(result, "try-error") || !is.list(result)
    }, logical(1))
    if (any(broken)) {
      detail <- replicate_results[[which(broken)[[1]]]]
      stop(
        "A bootstrap worker failed in ", batch$batch_id, ": ",
        if (inherits(detail, "try-error")) as.character(detail)[[1]] else "no result returned"
      )
    }
  }
  treatments <- stats::setNames(batch$comparisons$treatment, batch$comparisons$comparison)
  rows <- lapply(seq_len(nrow(batch$selections)), function(index) {
    gene_id <- batch$selections$gene_id[[index]]
    comparison <- batch$selections$comparison[[index]]
    status <- if (settings$replicates == 0L) {
      "disabled"
    } else if (!state$available[[gene_id]]) {
      "fit_unavailable"
    } else {
      "ok"
    }
    bootstrap_interval_rows(
      state, gene_id, comparison, treatments[[comparison]], status, replicate_results, settings
    )
  })
  output <- if (length(rows)) do.call(rbind, rows) else empty_table(INTERVAL_COLUMNS)
  write_gzip_tsv(select_columns(output, INTERVAL_COLUMNS, "The interval table"), output_path)
}

# ---- Finalize mode ----------------------------------------------------------

read_bootstrap_intervals <- function(paths) {
  if (!length(paths)) {
    stop("No bootstrap interval files were provided for finalization.")
  }
  values <- lapply(paths, read_tsv)
  invalid <- vapply(
    values,
    function(value) !all(INTERVAL_COLUMNS %in% names(value)),
    logical(1)
  )
  if (any(invalid)) {
    stop("A bootstrap interval file is missing required columns.")
  }
  intervals <- do.call(rbind, lapply(values, function(value) value[, INTERVAL_COLUMNS, drop = FALSE]))
  if (!nrow(intervals)) return(intervals)
  key_columns <- c("comparison", "gene_id", "feature_id")
  invalid_keys <- vapply(
    intervals[, key_columns, drop = FALSE],
    function(value) any(is.na(value) | !nzchar(as.character(value))),
    logical(1)
  )
  if (any(invalid_keys)) {
    stop("Bootstrap interval rows require non-missing comparison, gene_id, and feature_id.")
  }
  for (column in c("delta_pau_ci_low", "delta_pau_ci_high")) {
    intervals[[column]] <- suppressWarnings(as.numeric(intervals[[column]]))
  }
  counts <- lapply(c("bootstrap_successes", "bootstrap_perturbed"), function(column) {
    suppressWarnings(as.numeric(intervals[[column]]))
  })
  if (any(vapply(counts, function(values) {
    any(!is.finite(values) | values < 0 | values != floor(values))
  }, logical(1)))) {
    stop("Bootstrap interval successes must be non-negative integers.")
  }
  if (any(counts[[2]] > counts[[1]])) {
    stop("Bootstrap interval rows report more perturbed replicates than successes.")
  }
  if (!all(intervals$bootstrap_status %in% BOOTSTRAP_STATUSES)) {
    stop("Bootstrap interval rows have an unknown bootstrap_status.")
  }
  lower <- intervals$delta_pau_ci_low
  upper <- intervals$delta_pau_ci_high
  if (any(!is.na(lower) & !is.na(upper) & lower > upper)) {
    stop("Bootstrap interval lower bounds cannot exceed upper bounds.")
  }
  interval_keys <- do.call(paste, c(intervals[, key_columns, drop = FALSE], sep = "\r"))
  if (anyDuplicated(interval_keys)) {
    stop("Bootstrap interval files contain duplicate comparison/gene/PAC rows.")
  }
  intervals$bootstrap_successes <- as.integer(counts[[1]])
  intervals$bootstrap_perturbed <- as.integer(counts[[2]])
  intervals
}

apply_bootstrap_intervals <- function(output, comparison, intervals) {
  matches <- intervals[intervals$comparison == comparison, , drop = FALSE]
  if (!nrow(matches)) return(output)
  interval_keys <- paste(matches$gene_id, matches$feature_id, sep = "\r")
  output_keys <- paste(output$gene_id, output$pac_id, sep = "\r")
  index <- match(interval_keys, output_keys)
  if (anyNA(index)) {
    stop("Bootstrap intervals for ", comparison, " refer to PACs absent from its table.")
  }
  if (any(output$bootstrap_status[index] != "pending")) {
    stop("Bootstrap intervals for ", comparison, " refer to PACs that were not selected.")
  }
  for (column in c(
    "delta_pau_ci_low", "delta_pau_ci_high", "bootstrap_successes",
    "bootstrap_perturbed", "bootstrap_status"
  )) {
    output[[column]][index] <- matches[[column]]
  }
  output
}

finalize_preliminary_outputs <- function(preliminary_dir, interval_paths, output_dir) {
  if (!dir.exists(preliminary_dir)) {
    stop("Preliminary statistics directory does not exist: ", preliminary_dir)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  intervals <- read_bootstrap_intervals(interval_paths)
  sources <- list.files(
    preliminary_dir,
    pattern = "\\.tsv\\.gz$",
    full.names = TRUE,
    recursive = FALSE
  )
  if (!length(sources)) {
    warning("No preliminary statistics tables found in ", preliminary_dir)
    return(invisible(NULL))
  }
  comparisons <- sub("\\.pacs\\.tsv\\.gz$", "", basename(sources[grepl("\\.pacs\\.tsv\\.gz$", sources)]))
  unknown <- setdiff(unique(intervals$comparison), comparisons)
  if (length(unknown)) {
    stop("Bootstrap intervals name unknown comparisons: ", paste(unknown, collapse = ", "), ".")
  }
  for (source in sources) {
    destination <- file.path(output_dir, basename(source))
    if (!grepl("\\.pacs\\.tsv\\.gz$", source)) {
      if (!file.copy(source, destination, overwrite = TRUE)) {
        stop("Could not copy preliminary statistics file: ", source)
      }
      next
    }
    comparison <- sub("\\.pacs\\.tsv\\.gz$", "", basename(source))
    output <- read_tsv(source)
    output <- apply_bootstrap_intervals(output, comparison, intervals)
    if (any(output$bootstrap_status == "pending", na.rm = TRUE)) {
      stop("Selected PACs in ", comparison, " received no bootstrap intervals.")
    }
    write_gzip_tsv(output, destination)
    events <- output[output$event_type != "none", , drop = FALSE]
    event_path <- file.path(
      output_dir,
      paste0(comparison, ".events.tsv.gz")
    )
    write_gzip_tsv(events, event_path)
  }
  invisible(NULL)
}

run_bootstrap_mode <- function(arguments) {
  require_args(arguments, c("batch", "output"))
  run_bootstrap_batch(
    arguments$batch,
    arguments$output,
    parse_bootstrap_workers(arguments$bootstrap_workers)
  )
}

run_finalize_mode <- function(arguments) {
  require_args(arguments, c("preliminary_dir", "intervals", "output_dir"))
  interval_paths <- strsplit(arguments$intervals, ",", fixed = TRUE)[[1]]
  finalize_preliminary_outputs(
    arguments$preliminary_dir,
    interval_paths,
    arguments$output_dir
  )
}

main <- function(argv) {
  arguments <- parse_args(argv)
  versions <- load_statistics_packages()
  switch(
    arguments$mode,
    versions = run_versions_mode(arguments, versions),
    fit = run_fit_mode(arguments),
    bootstrap = run_bootstrap_mode(arguments),
    finalize = run_finalize_mode(arguments),
    stop("Unknown --mode: ", arguments$mode)
  )
}

if (sys.nframe() == 0L) main(commandArgs(trailingOnly = TRUE))
