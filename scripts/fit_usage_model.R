#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  required <- c("DRIMSeq", "stageR", "limma", "yaml")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". Use the PACusage Conda environment or Apptainer image."
    )
  }
  library(DRIMSeq)
  library(stageR)
  library(limma)
})

parse_args <- function(arguments) {
  result <- list()
  index <- 1
  while (index <= length(arguments)) {
    key <- sub("^--", "", arguments[[index]])
    if (index == length(arguments)) stop("Missing value for --", key)
    result[[gsub("-", "_", key)]] <- arguments[[index + 1]]
    index <- index + 2
  }
  required <- c("samples", "counts", "atlas", "params", "output_dir")
  missing <- required[!vapply(required, function(key) !is.null(result[[key]]), logical(1))]
  if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))
  result
}

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_gzip_tsv <- function(value, path) {
  connection <- gzfile(path, "wt")
  on.exit(close(connection), add = TRUE)
  write.table(value, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

bh <- function(values) {
  result <- rep(NA_real_, length(values))
  finite <- is.finite(values)
  result[finite] <- p.adjust(values[finite], method = "BH")
  result
}

family_filter <- function(counts, sample_ids, params) {
  sample_counts <- counts[, sample_ids, drop = FALSE]
  gene_total <- ave(rowSums(sample_counts), counts$gene_id, FUN = sum)
  site_total <- rowSums(sample_counts)
  supporting <- rowSums(sample_counts > 0)
  site_usage <- ifelse(gene_total > 0, site_total / gene_total, 0)
  site_ok <- site_total >= params$min_site_count &
    supporting >= params$min_test_supporting_samples &
    site_usage >= params$min_site_usage
  eligible_per_gene <- ave(site_ok, counts$gene_id, FUN = sum)
  gene_ok <- gene_total >= params$min_gene_total & eligible_per_gene >= 2
  counts[site_ok & gene_ok, c("gene_id", "pac_id", sample_ids), drop = FALSE]
}

observed_group_pau <- function(counts, sample_rows, group_name) {
  ids <- sample_rows$sample_id[sample_rows$condition == group_name]
  values <- rowSums(counts[, ids, drop = FALSE])
  totals <- ave(values, counts$gene_id, FUN = sum)
  ifelse(totals > 0, values / totals, NA_real_)
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

stage_adjust <- function(gene_results, feature_results, alpha) {
  screen <- gene_results$pvalue
  names(screen) <- gene_results$gene_id
  missing_screen <- !is.finite(screen)
  screen[missing_screen] <- 1
  confirmation <- matrix(feature_results$pvalue, ncol = 1)
  missing_confirmation <- !is.finite(confirmation[, 1])
  confirmation[missing_confirmation, 1] <- 1
  rownames(confirmation) <- feature_results$feature_id
  colnames(confirmation) <- "contrast"
  tx2gene <- data.frame(
    txID = feature_results$feature_id,
    geneID = feature_results$gene_id,
    stringsAsFactors = FALSE
  )
  object <- stageRTx(
    pScreen = screen,
    pConfirmation = confirmation,
    pScreenAdjusted = FALSE,
    tx2gene = tx2gene
  )
  object <- stageWiseAdjustment(object, method = "dtu", alpha = alpha, allowNA = TRUE)
  adjusted <- getAdjustedPValues(
    object,
    order = FALSE,
    onlySignificantGenes = FALSE
  )
  match_index <- match(feature_results$feature_id, rownames(adjusted))
  result <- as.numeric(adjusted[match_index, ncol(adjusted)])
  result[missing_confirmation] <- NA_real_
  result
}

safe_precision <- function(fit, genes) {
  value <- tryCatch(DRIMSeq::genewise_precision(fit), error = function(error) NULL)
  if (is.null(value)) {
    return(data.frame(gene_id = genes, precision = NA_real_))
  }
  if (is.data.frame(value)) {
    gene_column <- intersect(c("gene_id", "gene"), colnames(value))[1]
    precision_column <- intersect(
      c("genewise_precision", "precision", "common_precision"),
      colnames(value)
    )[1]
    if (!is.na(gene_column) && !is.na(precision_column)) {
      return(data.frame(
        gene_id = value[[gene_column]],
        precision = value[[precision_column]]
      ))
    }
  }
  data.frame(gene_id = genes, precision = as.numeric(value)[seq_along(genes)])
}

fitted_group_pau <- function(fit, counts, sample_rows, group_name) {
  fitted <- proportions(fit)
  sample_ids <- sample_rows$sample_id[sample_rows$condition == group_name]
  values <- rowMeans(fitted[, sample_ids, drop = FALSE])
  names(values) <- paste(fitted$gene_id, fitted$feature_id, sep = "\r")
  keys <- paste(counts$gene_id, counts$pac_id, sep = "\r")
  as.numeric(values[keys])
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

stable_seed <- function(base_seed, ...) {
  text <- paste(..., collapse = "|")
  values <- utf8ToInt(text)
  offset <- sum(values * seq_along(values)) %% 1000000000
  as.integer((as.numeric(base_seed) + offset) %% .Machine$integer.max)
}

stabilize_boundary_gene <- function(
  gene_id,
  counts,
  sample_rows,
  design,
  coefficient,
  control,
  treatment,
  params,
  atlas_checksum,
  comparison
) {
  gene_counts <- counts[counts$gene_id == gene_id, , drop = FALSE]
  dm_counts <- gene_counts
  colnames(dm_counts)[colnames(dm_counts) == "pac_id"] <- "feature_id"
  dm_samples <- sample_rows[, c(
    "sample_id", "condition", unlist(params$model_covariates)
  ), drop = FALSE]
  runs <- list()
  for (repeat_number in seq_len(params$dm_zero_sensitivity_repeats)) {
    set.seed(stable_seed(
      params$random_seed, atlas_checksum, comparison, gene_id, repeat_number
    ))
    run <- tryCatch({
      value <- dmDSdata(counts = dm_counts, samples = dm_samples)
      value <- dmPrecision(value, design = design, verbose = 0)
      value <- dmFit(value, design = design, add_uniform = TRUE, verbose = 0)
      tested <- dmTest(value, coef = coefficient, verbose = 0)
      feature_results <- results(tested, level = "feature")
      fitted <- proportions(value)
      control_ids <- sample_rows$sample_id[sample_rows$condition == control]
      treatment_ids <- sample_rows$sample_id[sample_rows$condition == treatment]
      data.frame(
        feature_id = fitted$feature_id,
        p_value = feature_results$pvalue[
          match(fitted$feature_id, feature_results$feature_id)
        ],
        control_pau = rowMeans(fitted[, control_ids, drop = FALSE]),
        treatment_pau = rowMeans(fitted[, treatment_ids, drop = FALSE])
      )
    }, error = function(error) NULL)
    if (!is.null(run)) runs[[length(runs) + 1]] <- run
  }
  if (!length(runs)) {
    return(data.frame(
      feature_id = gene_counts$pac_id,
      stabilized_p_value = NA_real_,
      stabilized_control_pau = NA_real_,
      stabilized_treatment_pau = NA_real_,
      stabilization_successes = 0L,
      zero_boundary_unstable = TRUE
    ))
  }
  long <- do.call(rbind, runs)
  output <- lapply(gene_counts$pac_id, function(feature_id) {
    values <- long[long$feature_id == feature_id, , drop = FALSE]
    finite <- is.finite(values$p_value) &
      is.finite(values$control_pau) &
      is.finite(values$treatment_pau)
    values <- values[finite, , drop = FALSE]
    deltas <- values$treatment_pau - values$control_pau
    enough <- nrow(values) >= ceiling(
      params$dm_zero_sensitivity_repeats * params$dm_bootstrap_min_success_fraction
    )
    direction_stable <- !length(deltas) ||
      length(unique(sign(deltas[abs(deltas) > .Machine$double.eps]))) <= 1
    spread <- if (length(deltas)) diff(range(deltas)) else Inf
    data.frame(
      feature_id = feature_id,
      stabilized_p_value = if (enough) median(values$p_value) else NA_real_,
      stabilized_control_pau = if (enough) median(values$control_pau) else NA_real_,
      stabilized_treatment_pau = if (enough) median(values$treatment_pau) else NA_real_,
      stabilization_successes = nrow(values),
      zero_boundary_unstable = !enough || !direction_stable ||
        spread > params$dm_zero_max_delta_pau_spread
    )
  })
  do.call(rbind, output)
}

bootstrap_gene <- function(
  gene_id,
  counts,
  sample_rows,
  design,
  fitted_model,
  precision,
  control,
  treatment,
  params,
  atlas_checksum,
  comparison
) {
  gene_counts <- counts[counts$gene_id == gene_id, , drop = FALSE]
  fitted <- proportions(fitted_model)
  fitted <- fitted[fitted$gene_id == gene_id, , drop = FALSE]
  fitted <- fitted[match(gene_counts$pac_id, fitted$feature_id), , drop = FALSE]
  sample_ids <- sample_rows$sample_id
  control_ids <- sample_rows$sample_id[sample_rows$condition == control]
  treatment_ids <- sample_rows$sample_id[sample_rows$condition == treatment]
  runs <- list()
  for (repeat_number in seq_len(params$dm_bootstrap_replicates)) {
    set.seed(stable_seed(
      params$random_seed,
      atlas_checksum,
      comparison,
      gene_id,
      "bootstrap",
      repeat_number
    ))
    simulated <- gene_counts
    for (sample_id in sample_ids) {
      total <- sum(gene_counts[[sample_id]])
      expected <- pmax(as.numeric(fitted[[sample_id]]), 1e-10)
      shapes <- expected * max(as.numeric(precision), 1e-6)
      draw <- rgamma(length(shapes), shape = shapes, rate = 1)
      draw <- draw / sum(draw)
      simulated[[sample_id]] <- as.integer(rmultinom(1, total, draw)[, 1])
    }
    run <- tryCatch({
      dm_counts <- simulated
      colnames(dm_counts)[colnames(dm_counts) == "pac_id"] <- "feature_id"
      dm_samples <- sample_rows[, c(
        "sample_id", "condition", unlist(params$model_covariates)
      ), drop = FALSE]
      value <- dmDSdata(counts = dm_counts, samples = dm_samples)
      value <- dmPrecision(value, design = design, verbose = 0)
      value <- dmFit(value, design = design, verbose = 0)
      fitted_run <- proportions(value)
      data.frame(
        feature_id = fitted_run$feature_id,
        delta_pau =
          rowMeans(fitted_run[, treatment_ids, drop = FALSE]) -
          rowMeans(fitted_run[, control_ids, drop = FALSE])
      )
    }, error = function(error) NULL)
    if (!is.null(run) && all(is.finite(run$delta_pau))) {
      runs[[length(runs) + 1]] <- run
    }
  }
  minimum_successes <- ceiling(
    params$dm_bootstrap_replicates * params$dm_bootstrap_min_success_fraction
  )
  if (length(runs) < minimum_successes) {
    return(data.frame(
      feature_id = gene_counts$pac_id,
      delta_pau_ci_low = NA_real_,
      delta_pau_ci_high = NA_real_,
      bootstrap_successes = length(runs)
    ))
  }
  long <- do.call(rbind, runs)
  output <- lapply(gene_counts$pac_id, function(feature_id) {
    values <- long$delta_pau[long$feature_id == feature_id]
    data.frame(
      feature_id = feature_id,
      delta_pau_ci_low = unname(quantile(values, 0.025)),
      delta_pau_ci_high = unname(quantile(values, 0.975)),
      bootstrap_successes = length(values)
    )
  })
  do.call(rbind, output)
}

classify_event <- function(row, params) {
  number <- function(value) suppressWarnings(as.numeric(value))
  truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
  control_detected <- number(row$control_supporting_samples) >=
    params$event_min_supporting_samples
  treatment_detected <- number(row$treatment_supporting_samples) >=
    params$event_min_supporting_samples
  positive <- number(row$delta_pau) >= params$min_abs_delta_pau
  negative <- number(row$delta_pau) <= -params$min_abs_delta_pau
  significant <- is.finite(number(row$gene_fdr)) &&
    number(row$gene_fdr) <= params$gene_fdr &&
    is.finite(number(row$pac_fdr)) &&
    number(row$pac_fdr) <= params$site_fdr
  stable <- !truth(row$zero_boundary_unstable)
  confident <- as.character(row$confidence) != "low" &&
    !truth(row$internal_priming_flag) &&
    !truth(row$exploratory_insufficient_replicates)
  gained_detection <- !control_detected && treatment_detected &&
    number(row$fitted_control_pau) <= params$event_max_control_pau &&
    number(row$fitted_treatment_pau) >= params$event_min_treatment_pau && positive
  lost_detection <- control_detected && !treatment_detected &&
    number(row$fitted_treatment_pau) <= params$event_max_control_pau &&
    number(row$fitted_control_pau) >= params$event_min_treatment_pau && negative
  if (gained_detection) {
    if (significant && stable && confident) "gained" else "gained_candidate"
  } else if (lost_detection) {
    if (significant && stable && confident) "lost" else "lost_candidate"
  } else if (control_detected && treatment_detected && significant && positive) {
    "increased_usage"
  } else if (control_detected && treatment_detected && significant && negative) {
    "decreased_usage"
  } else {
    "none"
  }
}

fit_family <- function(
  family,
  sample_rows,
  counts,
  atlas,
  params,
  output_dir,
  atlas_checksum
) {
  control <- family
  treatments <- sort(unique(sample_rows$condition[
    sample_rows$control_condition == family &
      sample_rows$condition != family
  ]))
  if (!length(treatments)) return(NULL)
  family_samples <- sample_rows[
    sample_rows$condition %in% c(control, treatments),
    ,
    drop = FALSE
  ]
  family_samples$condition <- factor(
    family_samples$condition,
    levels = c(control, treatments)
  )
  covariates <- unlist(params$model_covariates)
  formula_text <- paste(
    "~",
    paste(c(covariates, "condition"), collapse = " + ")
  )
  design <- model.matrix(as.formula(formula_text), data = family_samples)
  if (qr(design)$rank < ncol(design)) {
    stop(
      "Comparison family ", family,
      " has a rank-deficient design. Check condition/covariate confounding."
    )
  }
  filtered <- family_filter(counts, family_samples$sample_id, params)
  if (!nrow(filtered)) {
    warning("No testable genes in comparison family ", family)
    return(NULL)
  }
  dm_counts <- filtered
  colnames(dm_counts)[colnames(dm_counts) == "pac_id"] <- "feature_id"
  dm_samples <- family_samples[, c("sample_id", "condition", covariates), drop = FALSE]
  precision_data <- dmDSdata(counts = dm_counts, samples = dm_samples)
  precision_data <- dmPrecision(precision_data, design = design, verbose = 0)
  data <- dmFit(precision_data, design = design, verbose = 0)
  condition_coefficients <- grep("^condition", colnames(design))
  omnibus <- dmTest(data, coef = condition_coefficients, verbose = 0)
  omnibus_results <- results(omnibus, level = "gene")
  omnibus_results$gene_fdr <- bh(omnibus_results$pvalue)
  write_gzip_tsv(
    omnibus_results,
    file.path(output_dir, paste0(family, ".gene_omnibus.tsv.gz"))
  )
  precision <- safe_precision(data, unique(filtered$gene_id))
  precision$family <- family

  fitted_list <- list()
  for (treatment in treatments) {
    comparison <- paste0(treatment, "_vs_", control)
    coefficient <- match(paste0("condition", treatment), colnames(design))
    if (is.na(coefficient)) {
      stop("No model coefficient found for ", comparison)
    }
    tested <- dmTest(data, coef = coefficient, verbose = 0)
    genes <- results(tested, level = "gene")
    features <- results(tested, level = "feature")
    control_pau <- fitted_group_pau(data, filtered, family_samples, control)
    treatment_pau <- fitted_group_pau(data, filtered, family_samples, treatment)
    model_status <- rep("drimseq", nrow(filtered))
    zero_boundary_unstable <- rep(FALSE, nrow(filtered))
    stabilization_successes <- rep(0L, nrow(filtered))
    feature_keys <- paste(features$gene_id, features$feature_id, sep = "\r")
    count_keys <- paste(filtered$gene_id, filtered$pac_id, sep = "\r")
    feature_match <- match(count_keys, feature_keys)
    boundary_genes <- unique(filtered$gene_id[
      !is.finite(control_pau) |
        !is.finite(treatment_pau) |
        is.na(feature_match) |
        !is.finite(features$pvalue[feature_match])
    ])
    for (gene_id in boundary_genes) {
      stabilized <- stabilize_boundary_gene(
        gene_id,
        filtered,
        family_samples,
        design,
        coefficient,
        control,
        treatment,
        params,
        atlas_checksum,
        comparison
      )
      count_index <- which(filtered$gene_id == gene_id)
      stabilization_match <- match(
        filtered$pac_id[count_index], stabilized$feature_id
      )
      feature_index <- feature_match[count_index]
      usable <- !stabilized$zero_boundary_unstable[stabilization_match]
      replace_indices <- count_index[usable]
      replace_stabilized <- stabilization_match[usable]
      replace_features <- feature_index[usable]
      if (length(replace_indices)) {
        control_pau[replace_indices] <-
          stabilized$stabilized_control_pau[replace_stabilized]
        treatment_pau[replace_indices] <-
          stabilized$stabilized_treatment_pau[replace_stabilized]
        features$pvalue[replace_features] <-
          stabilized$stabilized_p_value[replace_stabilized]
      }
      model_status[count_index] <- "drimseq_add_uniform"
      zero_boundary_unstable[count_index] <-
        stabilized$zero_boundary_unstable[stabilization_match]
      stabilization_successes[count_index] <-
        stabilized$stabilization_successes[stabilization_match]
      unstable_feature_indices <- feature_index[
        stabilized$zero_boundary_unstable[stabilization_match]
      ]
      unstable_feature_indices <- unstable_feature_indices[
        !is.na(unstable_feature_indices)
      ]
      features$pvalue[unstable_feature_indices] <- NA_real_
    }
    genes$gene_fdr <- bh(genes$pvalue)
    features$pac_fdr <- tryCatch(
      stage_adjust(genes, features, params$site_fdr),
      error = function(error) {
        warning("stageR adjustment failed for ", comparison, ": ", conditionMessage(error))
        rep(NA_real_, nrow(features))
      }
    )
    output <- merge(
      features,
      genes[, c("gene_id", "pvalue", "gene_fdr")],
      by = "gene_id",
      suffixes = c("_pac", "_gene"),
      all.x = TRUE
    )
    effect <- data.frame(
      gene_id = filtered$gene_id,
      feature_id = filtered$pac_id,
      fitted_control_pau = control_pau,
      fitted_treatment_pau = treatment_pau,
      delta_pau = treatment_pau - control_pau,
      control_supporting_samples = supporting_samples(filtered, family_samples, control),
      treatment_supporting_samples = supporting_samples(filtered, family_samples, treatment),
      control_gene_total = group_gene_totals(filtered, family_samples, control),
      treatment_gene_total = group_gene_totals(filtered, family_samples, treatment),
      raw_control_counts = raw_count_text(filtered, family_samples, control),
      raw_treatment_counts = raw_count_text(filtered, family_samples, treatment),
      observed_control_pau = observed_pau_text(filtered, family_samples, control),
      observed_treatment_pau = observed_pau_text(filtered, family_samples, treatment),
      model_status = model_status,
      stabilization_successes = stabilization_successes,
      zero_boundary_unstable = zero_boundary_unstable
    )
    effect <- merge(effect, precision, by = "gene_id", all.x = TRUE)
    effect$alpha_control <- effect$fitted_control_pau * effect$precision
    effect$alpha_treatment <- effect$fitted_treatment_pau * effect$precision
    output <- merge(output, effect, by = c("gene_id", "feature_id"), all.x = TRUE)
    output <- merge(
      output,
      atlas[, c(
        "pac_id", "assignment_class", "known_pac", "known_rescue_only",
        "confidence", "internal_priming_flag", "primary_pas_motif",
        "primary_pas_motif_rna", "primary_motif_class"
      )],
      by.x = "feature_id",
      by.y = "pac_id",
      all.x = TRUE
    )
    output$pac_id <- output$feature_id
    output$site_class <- output$assignment_class
    output$condition <- treatment
    output$control_condition <- control
    condition_counts <- table(family_samples$condition)
    output$exploratory_insufficient_replicates <- any(
      condition_counts < params$min_replicates_per_condition
    )
    output$delta_pau_ci_low <- NA_real_
    output$delta_pau_ci_high <- NA_real_
    output$bootstrap_successes <- 0L
    output$effect_exceeds_threshold <- abs(output$delta_pau) >= params$min_abs_delta_pau
    detection_candidate <- (
      output$fitted_control_pau <= params$event_max_control_pau &
        output$fitted_treatment_pau >= params$event_min_treatment_pau &
        output$treatment_supporting_samples >= params$event_min_supporting_samples &
        output$delta_pau >= params$min_abs_delta_pau
    ) | (
      output$fitted_treatment_pau <= params$event_max_control_pau &
        output$fitted_control_pau >= params$event_min_treatment_pau &
        output$control_supporting_samples >= params$event_min_supporting_samples &
        output$delta_pau <= -params$min_abs_delta_pau
    )
    bootstrap_genes <- unique(output$gene_id[
      (!is.na(output$gene_fdr) & output$gene_fdr <= params$gene_fdr) |
        detection_candidate
    ])
    for (gene_id in bootstrap_genes) {
      precision_value <- precision$precision[match(gene_id, precision$gene_id)]
      intervals <- bootstrap_gene(
        gene_id,
        filtered,
        family_samples,
        design,
        data,
        precision_value,
        control,
        treatment,
        params,
        atlas_checksum,
        comparison
      )
      output_indices <- which(output$gene_id == gene_id)
      interval_match <- match(output$feature_id[output_indices], intervals$feature_id)
      output$delta_pau_ci_low[output_indices] <-
        intervals$delta_pau_ci_low[interval_match]
      output$delta_pau_ci_high[output_indices] <-
        intervals$delta_pau_ci_high[interval_match]
      output$bootstrap_successes[output_indices] <-
        intervals$bootstrap_successes[interval_match]
    }
    row_groups <- split(seq_len(nrow(output)), output$gene_id)
    dominant_control <- vapply(
      row_groups,
      function(indices) {
        output$feature_id[indices[which.max(output$fitted_control_pau[indices])]]
      },
      character(1)
    )
    dominant_treatment <- vapply(
      row_groups,
      function(indices) {
        output$feature_id[indices[which.max(output$fitted_treatment_pau[indices])]]
      },
      character(1)
    )
    output$dominant_pac_control <- unname(dominant_control[output$gene_id])
    output$dominant_pac_treatment <- unname(dominant_treatment[output$gene_id])
    output$control_detected_complexity <- ave(
      output$fitted_control_pau >= params$event_min_treatment_pau,
      output$gene_id,
      FUN = sum
    )
    output$treatment_detected_complexity <- ave(
      output$fitted_treatment_pau >= params$event_min_treatment_pau,
      output$gene_id,
      FUN = sum
    )
    output$event_type <- apply(
      output,
      1,
      function(row) classify_event(as.list(row), params)
    )
    output$event_type[
      output$exploratory_insufficient_replicates &
        output$event_type == "gained"
    ] <- "gained_candidate"
    output$event_type[
      output$exploratory_insufficient_replicates &
        output$event_type == "lost"
    ] <- "lost_candidate"
    dominant_switch <- output$dominant_pac_control != output$dominant_pac_treatment
    output$event_type[
      output$event_type == "none" &
        dominant_switch &
        output$feature_id == output$dominant_pac_treatment
    ] <- "dominant_switch"
    complexity_delta <- output$treatment_detected_complexity -
      output$control_detected_complexity
    output$event_type[
      output$event_type == "none" & complexity_delta > 0
    ] <- "complexity_gain"
    output$event_type[
      output$event_type == "none" & complexity_delta < 0
    ] <- "complexity_loss"
    gene_path <- file.path(output_dir, paste0(comparison, ".genes.tsv.gz"))
    pac_path <- file.path(output_dir, paste0(comparison, ".pacs.tsv.gz"))
    event_path <- file.path(output_dir, paste0(comparison, ".events.tsv.gz"))
    write_gzip_tsv(genes, gene_path)
    write_gzip_tsv(output, pac_path)
    events <- output[output$event_type != "none", , drop = FALSE]
    write_gzip_tsv(events, event_path)
    fitted_list[[comparison]] <- output[, c(
      "gene_id", "feature_id", "condition", "control_condition",
      "fitted_control_pau", "fitted_treatment_pau", "delta_pau",
      "precision", "alpha_control", "alpha_treatment", "model_status"
    ), drop = FALSE]
  }
  list(fitted = fitted_list, precision = precision)
}

fit_motif_preferences <- function(scores, samples, params, output_dir, suffix = "preference") {
  if (is.null(scores) || !nrow(scores)) return(invisible(NULL))
  if (!"primary_pas_motif_rna" %in% names(scores)) {
    if ("primary_pas_motif" %in% names(scores)) {
      scores$primary_pas_motif_rna <- chartr("T", "U", scores$primary_pas_motif)
      scores$primary_pas_motif_rna[
        is.na(scores$primary_pas_motif_rna) |
          scores$primary_pas_motif_rna == ""
      ] <- "none"
    } else {
      scores$primary_pas_motif_rna <- "unresolved"
    }
  }
  merged <- merge(scores, samples, by = "sample_id")
  comparison_map <- unique(
    samples[
      samples$condition != samples$control_condition,
      c("condition", "control_condition"),
      drop = FALSE
    ]
  )
  for (family in unique(comparison_map$control_condition)) {
    treatments <- sort(comparison_map$condition[
      comparison_map$control_condition == family
    ])
    family_rows <- merged[
      merged$condition %in% c(family, treatments),
      ,
      drop = FALSE
    ]
    for (treatment in treatments) {
      subset_rows <- family_rows[
        family_rows$condition %in% c(family, treatment) &
          family_rows$informative_genes >= params$motif_preference_min_genes,
        ,
        drop = FALSE
      ]
      output_rows <- list()
      motifs <- unique(
        subset_rows[, c("primary_pas_motif_rna", "primary_motif_class"), drop = FALSE]
      )
      for (motif_index in seq_len(nrow(motifs))) {
        motif <- motifs$primary_pas_motif_rna[[motif_index]]
        motif_class <- motifs$primary_motif_class[[motif_index]]
        values <- subset_rows[
          subset_rows$primary_pas_motif_rna == motif &
            subset_rows$primary_motif_class == motif_class,
          ,
          drop = FALSE
        ]
        if (length(unique(values$condition)) < 2) next
        values$condition <- relevel(factor(values$condition), ref = family)
        design <- model.matrix(~ condition, values)
        fit <- eBayes(lmFit(matrix(values$transformed_motif_usage, nrow = 1), design))
        coefficient <- paste0("condition", treatment)
        table <- topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
        output_rows[[length(output_rows) + 1]] <- data.frame(
          primary_pas_motif_rna = motif,
          primary_motif_class = motif_class,
          condition = treatment,
          control_condition = family,
          control_mean = mean(values$motif_usage[values$condition == family]),
          treatment_mean = mean(values$motif_usage[values$condition == treatment]),
          delta_motif_usage =
            mean(values$motif_usage[values$condition == treatment]) -
            mean(values$motif_usage[values$condition == family]),
          transformed_coefficient = table$logFC[[1]],
          p_value = table$P.Value[[1]],
          informative_genes = min(values$informative_genes)
        )
      }
      if (length(output_rows)) {
        output <- do.call(rbind, output_rows)
        output$fdr <- bh(output$p_value)
        name <- paste0(treatment, "_vs_", family, ".", suffix, ".tsv.gz")
        write_gzip_tsv(output, file.path(output_dir, name))
      }
    }
  }
}

arguments <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(arguments$output_dir, recursive = TRUE, showWarnings = FALSE)
params <- yaml::read_yaml(arguments$params)
samples <- read_tsv(arguments$samples)
counts <- read_tsv(arguments$counts)
atlas <- read_tsv(arguments$atlas)

families <- unique(samples$control_condition[
  samples$condition != samples$control_condition
])
all_precision <- list()
all_fitted <- list()
atlas_checksum <- unname(tools::md5sum(arguments$atlas))
for (family in families) {
  fitted <- fit_family(
    family,
    samples,
    counts,
    atlas,
    params,
    arguments$output_dir,
    atlas_checksum
  )
  if (!is.null(fitted)) {
    all_precision[[family]] <- fitted$precision
    all_fitted <- c(all_fitted, fitted$fitted)
  }
}
if (length(all_precision)) {
  write_gzip_tsv(
    do.call(rbind, all_precision),
    file.path(arguments$output_dir, "gene_precision.tsv.gz")
  )
}
if (length(all_fitted)) {
  fitted_rows <- do.call(rbind, all_fitted)
  write_gzip_tsv(
    fitted_rows,
    file.path(arguments$output_dir, "fitted_pau.tsv.gz")
  )
}

if (!is.null(arguments$motif_scores) && file.exists(arguments$motif_scores)) {
  fit_motif_preferences(
    read_tsv(arguments$motif_scores),
    samples,
    params,
    arguments$output_dir,
    "preference"
  )
}
if (!is.null(arguments$motif_sensitivity) && file.exists(arguments$motif_sensitivity)) {
  fit_motif_preferences(
    read_tsv(arguments$motif_sensitivity),
    samples,
    params,
    arguments$output_dir,
    "preference_known_rescue_sensitivity"
  )
}
