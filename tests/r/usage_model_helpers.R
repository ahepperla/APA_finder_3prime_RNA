# Shared helpers for the R statistics tests. Base R only, so the tests add no
# dependencies beyond the pipeline's own R packages.

# ---- Test harness ---------------------------------------------------------

test_state <- new.env()
test_state$passed <- character()
test_state$failed <- character()

# Run one named case. Errors become failures so later cases still run.
test_case <- function(id, name, code) {
  message_text <- tryCatch(
    {
      force(code)
      NULL
    },
    error = function(error) conditionMessage(error)
  )
  label <- paste0(id, " ", name)
  if (is.null(message_text)) {
    test_state$passed <- c(test_state$passed, label)
    cat("PASS ", label, "\n", sep = "")
  } else {
    test_state$failed <- c(test_state$failed, label)
    cat("FAIL ", label, ": ", message_text, "\n", sep = "")
  }
  invisible(is.null(message_text))
}

check <- function(condition, ...) {
  if (!isTRUE(condition)) stop(paste0(...), call. = FALSE)
  invisible(TRUE)
}

finish_tests <- function() {
  cat(sprintf(
    "\n%d passed, %d failed\n",
    length(test_state$passed),
    length(test_state$failed)
  ))
  if (length(test_state$failed)) quit(save = "no", status = 1L)
  invisible(TRUE)
}

# Select exactly one row, so an assertion can never pass on an empty or
# ambiguous match.
one_row <- function(table, ...) {
  keys <- list(...)
  selected <- rep(TRUE, nrow(table))
  for (column in names(keys)) {
    check(column %in% names(table), "Column ", column, " is missing.")
    selected <- selected & !is.na(table[[column]]) &
      as.character(table[[column]]) == as.character(keys[[column]])
  }
  check(
    sum(selected) == 1L,
    "Expected exactly one row for ",
    paste(names(keys), unlist(keys), sep = "=", collapse = ", "),
    "; found ", sum(selected), "."
  )
  table[selected, , drop = FALSE]
}

as_number <- function(value) suppressWarnings(as.numeric(as.character(value)))

as_flag <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")

# ---- Loading the statistics script ---------------------------------------

# The script runs its command line only when Rscript runs it directly, so
# sourcing it just defines its functions.
load_usage_model <- function(path) {
  environment <- new.env(parent = globalenv())
  sys.source(path, envir = environment)
  environment$load_statistics_packages()
  environment
}

# ---- DRIMSeq call counters ------------------------------------------------

# The tracer runs inside DRIMSeq's methods and finds this environment through
# the global environment. Counts are only reliable in serial runs, because
# increments made in forked workers are lost.
assign(".pacusage_drimseq_calls", new.env(), envir = globalenv())

reset_drimseq_calls <- function() {
  assign("dmPrecision", 0L, envir = .pacusage_drimseq_calls)
  assign("dmFit", 0L, envir = .pacusage_drimseq_calls)
  invisible(NULL)
}

drimseq_calls <- function() {
  c(
    dmPrecision = get("dmPrecision", envir = .pacusage_drimseq_calls),
    dmFit = get("dmFit", envir = .pacusage_drimseq_calls)
  )
}

install_drimseq_counters <- function() {
  reset_drimseq_calls()
  suppressMessages({
    trace(
      "dmPrecision",
      signature = "dmDSdata",
      where = asNamespace("DRIMSeq"),
      print = FALSE,
      tracer = quote(
        .pacusage_drimseq_calls$dmPrecision <- .pacusage_drimseq_calls$dmPrecision + 1L
      )
    )
    trace(
      "dmFit",
      signature = "dmDSprecision",
      where = asNamespace("DRIMSeq"),
      print = FALSE,
      tracer = quote(.pacusage_drimseq_calls$dmFit <- .pacusage_drimseq_calls$dmFit + 1L)
    )
  })
  # Self-check on a two-gene object, so a broken tracer cannot make a cost
  # assertion pass vacuously.
  samples <- data.frame(sample_id = paste0("s", 1:4), condition = factor(c("a", "a", "b", "b")))
  counts <- data.frame(
    gene_id = rep(c("g1", "g2"), each = 2L),
    feature_id = paste0("f", 1:4),
    s1 = c(10, 20, 30, 40), s2 = c(12, 18, 33, 37),
    s3 = c(20, 10, 40, 30), s4 = c(22, 8, 42, 28)
  )
  design <- stats::model.matrix(~ condition, data = samples)
  suppressMessages({
    precision <- DRIMSeq::dmPrecision(
      DRIMSeq::dmDSdata(counts = counts, samples = samples), design = design, verbose = 0
    )
    DRIMSeq::dmFit(precision, design = design, verbose = 0)
  })
  if (!identical(unname(drimseq_calls()), c(1L, 1L))) {
    stop("DRIMSeq call counters are not working.")
  }
  reset_drimseq_calls()
}

# ---- Synthetic comparison family -----------------------------------------

rdirichlet_one <- function(alpha) {
  draw <- stats::rgamma(length(alpha), shape = alpha)
  draw / sum(draw)
}

statistics_params <- function(replicates = 40L) {
  list(
    model_covariates = list(),
    random_seed = 1729L,
    min_gene_total = 20L,
    min_site_count = 5L,
    min_site_usage = 0.01,
    min_site_usage_samples = NULL,
    min_site_usage_gene_reads = 10L,
    min_test_supporting_samples = 2L,
    min_replicates_per_condition = 2L,
    insufficient_replicates_policy = "error",
    gene_fdr = 0.05,
    site_fdr = 0.05,
    min_abs_delta_pau = 0.10,
    apa_pattern_min_change = 0.10,
    dm_bootstrap_replicates = as.integer(replicates),
    dm_bootstrap_min_success_fraction = 0.80,
    dm_bootstrap_include_candidates = TRUE,
    dm_zero_sensitivity_repeats = 5L,
    dm_zero_max_delta_pau_spread = 0.02,
    event_min_treatment_pau = 0.05,
    active_pac_min_pau = 0.05,
    event_max_control_pau = 0.01,
    event_min_supporting_samples = 2L,
    potential_internal_priming_withheld_calls = TRUE,
    motif_preference_min_genes = 50L
  )
}

# Replicate columns listed per condition in sample order, one vector per
# sample, each holding the counts of every PAC in genomic order.
literal_gene <- function(...) {
  columns <- list(...)
  do.call(cbind, columns)
}

# One Dirichlet-multinomial draw per sample: a negative-binomial gene total,
# sample proportions from a Dirichlet around the condition's proportions,
# then multinomial PAC counts.
simulate_dm_counts <- function(proportions_by_condition, precision, depth, sample_condition) {
  counts <- vapply(seq_along(sample_condition), function(index) {
    base <- proportions_by_condition[[sample_condition[[index]]]]
    total <- stats::rnbinom(1L, mu = depth, size = 20)
    as.numeric(stats::rmultinom(1L, total, rdirichlet_one(base * precision)))
  }, numeric(length(proportions_by_condition[[1]])))
  matrix(as.integer(counts), ncol = length(sample_condition))
}

# Build the statistics inputs from genes given as list(gene_id, contig,
# coordinates, counts), where counts has one column per sample.
assemble_usage_dataset <- function(genes, sample_ids, sample_condition, control, params) {
  # Genes are named after their IDs, the way the annotation step names genes
  # the annotation leaves unnamed.
  gene_name_of <- function(gene) if (is.null(gene$gene_name)) toupper(gene$gene_id) else gene$gene_name
  # The identity block of the pipeline's PAC tables: exact PACs on the plus
  # strand, with their BED interval and 1-based locus.
  identity <- function(gene) {
    data.frame(
      pac_id = sprintf("PACv1.test.%s.+.%d", gene$contig, gene$coordinates),
      gene_id = gene$gene_id,
      gene_name = gene_name_of(gene),
      chrom = gene$contig,
      start = gene$coordinates,
      end = gene$coordinates + 1L,
      strand = "+",
      locus = sprintf("%s:%d-%d", gene$contig, gene$coordinates + 1L, gene$coordinates + 1L),
      stringsAsFactors = FALSE
    )
  }
  counts <- do.call(rbind, lapply(genes, function(gene) {
    values <- as.data.frame(gene$counts)
    names(values) <- sample_ids
    data.frame(identity(gene), values, check.names = FALSE, stringsAsFactors = FALSE)
  }))
  atlas <- do.call(rbind, lapply(genes, function(gene) {
    data.frame(
      identity(gene),
      coordinate = gene$coordinates,
      candidate_status = "primary",
      known_pac = FALSE,
      known_rescue_only = FALSE,
      gene_region = "last_exon",
      # Every test gene has one last exon holding all of its PACs.
      last_exon_locus = sprintf(
        "%s:%d-%d", gene$contig, min(gene$coordinates) - 99L, max(gene$coordinates)
      ),
      ambiguous_gene_assignment = FALSE,
      upstream_sequence = "",
      primary_pas_motif = "AATAAA",
      primary_pas_motif_rna = "AAUAAA",
      primary_motif_class = "canonical",
      internal_priming_flag = FALSE,
      confidence = "high",
      stringsAsFactors = FALSE
    )
  }))
  samples <- data.frame(
    sample_id = sample_ids,
    alignment = paste0(sample_ids, ".bam"),
    condition = sample_condition,
    control = ifelse(sample_condition == control, "", control),
    control_condition = control,
    replicate = stats::ave(seq_along(sample_ids), sample_condition, FUN = seq_along),
    batch = "",
    donor = "",
    layout = "SE",
    strandedness = "forward",
    library_profile = "exact_boundary",
    evidence_source = "read_3p",
    stringsAsFactors = FALSE
  )
  list(counts = counts, atlas = atlas, samples = samples, params = params)
}

# Three conditions (C, T1, T2) with three replicates each. The family is C.
# Designed genes use literal counts; the rest are seeded Dirichlet-multinomial
# draws. Exactly two genes have a PAC with zero counts across a condition.
build_usage_dataset <- function(replicates = 40L) {
  previous_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(previous_kind)), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(20260927)

  conditions <- c("C", "T1", "T2")
  sample_ids <- as.vector(t(outer(conditions, 1:3, paste, sep = "_")))
  sample_condition <- rep(conditions, each = 3L)

  genes <- list()
  add_gene <- function(gene_id, counts, contig = "chrT", start = NULL) {
    index <- length(genes) + 1L
    if (is.null(start)) start <- 1000L * index
    coordinates <- start + 50L * (seq_len(nrow(counts)) - 1L)
    genes[[index]] <<- list(
      gene_id = gene_id,
      contig = contig,
      coordinates = coordinates,
      counts = counts
    )
  }
  simulated <- function(proportions_by_condition, precision, depth) {
    simulate_dm_counts(proportions_by_condition, precision, depth, sample_condition)
  }

  add_gene("gene_gain", literal_gene(
    c(195, 105, 0), c(180, 120, 0), c(210, 90, 0),
    c(135, 60, 105), c(130, 65, 105), c(140, 55, 105),
    c(200, 100, 0), c(185, 115, 0), c(190, 110, 0)
  ))
  add_gene("gene_unrelated_zero", literal_gene(
    c(150, 90, 60), c(140, 95, 65), c(160, 85, 55),
    c(75, 165, 60), c(80, 160, 60), c(70, 170, 60),
    c(180, 120, 0), c(175, 125, 0), c(185, 115, 0)
  ))
  add_gene("gene:colon", literal_gene(
    c(150, 90, 60), c(145, 95, 60), c(155, 85, 60),
    c(75, 90, 135), c(80, 85, 135), c(70, 95, 135),
    c(150, 90, 60), c(148, 92, 60), c(152, 88, 60)
  ), contig = "HLA-A*01:01:01:01", start = 100L)
  add_gene("shift_both_01", literal_gene(
    c(150, 90, 60), c(140, 95, 65), c(160, 85, 55),
    c(75, 90, 135), c(70, 95, 135), c(80, 85, 135),
    c(75, 90, 135), c(78, 88, 134), c(72, 92, 136)
  ))
  add_gene("shift_both_02", literal_gene(
    c(148, 92, 60), c(152, 88, 60), c(150, 90, 60),
    c(72, 93, 135), c(78, 87, 135), c(75, 90, 135),
    c(74, 91, 135), c(76, 89, 135), c(75, 90, 135)
  ))
  for (index in 1:6) {
    add_gene(sprintf("shift_T1_%02d", index), simulated(
      list(C = c(0.5, 0.3, 0.2), T1 = c(0.25, 0.3, 0.45), T2 = c(0.5, 0.3, 0.2)),
      precision = 100, depth = 300
    ))
  }
  three_pac_bases <- list(c(0.5, 0.3, 0.2), c(0.4, 0.35, 0.25), c(0.6, 0.25, 0.15))
  two_pac_bases <- list(c(0.6, 0.4), c(0.7, 0.3))
  depths <- c(100, 200, 300)
  precisions <- c(50, 100, 150)
  for (index in 1:28) {
    base <- if (index <= 18) {
      three_pac_bases[[(index - 1L) %% 3L + 1L]]
    } else {
      two_pac_bases[[(index - 1L) %% 2L + 1L]]
    }
    add_gene(sprintf("null_%02d", index), simulated(
      list(C = base, T1 = base, T2 = base),
      precision = precisions[[(index - 1L) %% 3L + 1L]],
      depth = depths[[(index %/% 3L) %% 3L + 1L]]
    ))
  }

  dataset <- assemble_usage_dataset(genes, sample_ids, sample_condition, "C", statistics_params(replicates))

  zero_group_genes <- unique(unlist(lapply(genes, function(gene) {
    has_zero <- any(vapply(conditions, function(condition) {
      any(rowSums(gene$counts[, sample_condition == condition, drop = FALSE]) == 0)
    }, logical(1)))
    if (has_zero) gene$gene_id else NULL
  })))
  if (!identical(sort(zero_group_genes), c("gene_gain", "gene_unrelated_zero"))) {
    stop("Synthetic data have unexpected zero groups: ", paste(zero_group_genes, collapse = ", "))
  }
  dataset
}

write_usage_inputs <- function(dataset, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  paths <- list(
    samples = file.path(directory, "normalized_samples.tsv"),
    counts = file.path(directory, "pac_counts.tsv.gz"),
    atlas = file.path(directory, "pacs.v1.metadata.tsv.gz"),
    params = file.path(directory, "resolved_params.yaml"),
    motif_scores = file.path(directory, "motif_scores.tsv"),
    motif_sensitivity = file.path(directory, "motif_scores_known_rescue_sensitivity.tsv"),
    motif_class_scores = file.path(directory, "motif_class_scores.tsv"),
    motif_class_sensitivity = file.path(directory, "motif_class_scores_known_rescue_sensitivity.tsv")
  )
  utils::write.table(dataset$samples, paths$samples, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "")
  connection <- gzfile(paths$counts, "wt")
  utils::write.table(dataset$counts, connection, sep = "\t", quote = FALSE, row.names = FALSE)
  close(connection)
  connection <- gzfile(paths$atlas, "wt")
  utils::write.table(dataset$atlas, connection, sep = "\t", quote = FALSE, row.names = FALSE)
  close(connection)
  yaml::write_yaml(dataset$params, paths$params)
  motif_scores <- if (is.null(dataset$motif_scores)) empty_motif_scores() else dataset$motif_scores
  for (path in c(paths$motif_scores, paths$motif_sensitivity)) {
    utils::write.table(motif_scores, path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  class_scores <- if (is.null(dataset$motif_class_scores)) {
    empty_motif_scores()[, setdiff(names(empty_motif_scores()), "primary_pas_motif_rna")]
  } else {
    dataset$motif_class_scores
  }
  for (path in c(paths$motif_class_scores, paths$motif_class_sensitivity)) {
    utils::write.table(class_scores, path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  paths
}

# The columns of pacusage motif-scores output, with no rows.
empty_motif_scores <- function() {
  data.frame(
    sample_id = character(), primary_pas_motif_rna = character(),
    primary_motif_class = character(), motif_usage = numeric(),
    transformed_motif_usage = numeric(), informative_genes = integer(),
    stringsAsFactors = FALSE
  )
}

# ---- Running the three statistics modes ----------------------------------

list_batches <- function(fit_directory) {
  sort(list.files(
    file.path(fit_directory, "bootstrap-batches"),
    pattern = "\\.rds$",
    full.names = TRUE
  ))
}

interval_path_for <- function(batch, directory) {
  file.path(directory, sub("\\.rds$", ".intervals.tsv.gz", basename(batch)))
}

# In-process run. With count = TRUE, DRIMSeq calls are counted separately for
# the fit and the bootstrap. Counting needs serial workers.
run_statistics_in_process <- function(
  model, paths, directory, family = "C",
  model_workers = 1L, batch_size = 1L, bootstrap_workers = 1L, count = FALSE
) {
  fit_directory <- file.path(directory, "fit")
  interval_directory <- file.path(directory, "intervals")
  final_directory <- file.path(directory, "final")
  dir.create(interval_directory, recursive = TRUE, showWarnings = FALSE)
  if (count) reset_drimseq_calls()
  model$run_fit_mode(list(
    mode = "fit",
    family = family,
    samples = paths$samples,
    counts = paths$counts,
    atlas = paths$atlas,
    params = paths$params,
    motif_scores = paths$motif_scores,
    motif_sensitivity = paths$motif_sensitivity,
    motif_class_scores = paths$motif_class_scores,
    motif_class_sensitivity = paths$motif_class_sensitivity,
    model_workers = as.character(model_workers),
    bootstrap_batch_size = as.character(batch_size),
    output_dir = fit_directory
  ))
  fit_calls <- if (count) drimseq_calls() else NULL
  batches <- list_batches(fit_directory)
  if (count) reset_drimseq_calls()
  interval_paths <- vapply(batches, function(batch) {
    output <- interval_path_for(batch, interval_directory)
    model$run_bootstrap_mode(list(
      mode = "bootstrap",
      batch = batch,
      bootstrap_workers = as.character(bootstrap_workers),
      output = output
    ))
    output
  }, character(1))
  bootstrap_calls <- if (count) drimseq_calls() else NULL
  model$run_finalize_mode(list(
    mode = "finalize",
    preliminary_dir = file.path(fit_directory, "preliminary"),
    intervals = paste(interval_paths, collapse = ","),
    output_dir = final_directory
  ))
  list(
    fit_directory = fit_directory,
    final_directory = final_directory,
    batches = batches,
    fit_calls = fit_calls,
    bootstrap_calls = bootstrap_calls
  )
}

run_rscript <- function(arguments) {
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    shQuote(arguments),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Rscript failed (", status, "): ", paste(utils::tail(output, 15L), collapse = "\n"))
  }
  invisible(output)
}

# Separate Rscript processes, exactly as Nextflow runs the modes.
run_statistics_subprocess <- function(
  script, paths, directory, family = "C",
  model_workers = 2L, batch_size = 500L, bootstrap_workers = 2L
) {
  fit_directory <- file.path(directory, "fit")
  interval_directory <- file.path(directory, "intervals")
  final_directory <- file.path(directory, "final")
  dir.create(interval_directory, recursive = TRUE, showWarnings = FALSE)
  run_rscript(c(
    script, "--mode", "fit", "--family", family,
    "--samples", paths$samples, "--counts", paths$counts,
    "--atlas", paths$atlas, "--params", paths$params,
    "--motif-scores", paths$motif_scores, "--motif-sensitivity", paths$motif_sensitivity,
    "--motif-class-scores", paths$motif_class_scores,
    "--motif-class-sensitivity", paths$motif_class_sensitivity,
    "--model-workers", as.character(model_workers),
    "--bootstrap-batch-size", as.character(batch_size),
    "--output-dir", fit_directory
  ))
  batches <- list_batches(fit_directory)
  interval_paths <- vapply(batches, function(batch) {
    output <- interval_path_for(batch, interval_directory)
    run_rscript(c(
      script, "--mode", "bootstrap", "--batch", batch,
      "--bootstrap-workers", as.character(bootstrap_workers),
      "--output", output
    ))
    output
  }, character(1))
  run_rscript(c(
    script, "--mode", "finalize",
    "--preliminary-dir", file.path(fit_directory, "preliminary"),
    "--intervals", paste(interval_paths, collapse = ","),
    "--output-dir", final_directory
  ))
  list(fit_directory = fit_directory, final_directory = final_directory, batches = batches)
}

read_result <- function(directory, name) {
  path <- file.path(directory, name)
  check(file.exists(path), "Missing output ", name, ".")
  utils::read.delim(path, colClasses = "character", check.names = FALSE, na.strings = "")
}

# Every table in a final directory, rows sorted by all columns, for exact
# comparison between runs.
final_tables <- function(directory) {
  names <- sort(list.files(directory, pattern = "\\.tsv(\\.gz)?$"))
  tables <- lapply(names, function(name) {
    table <- read_result(directory, name)
    if (!nrow(table)) return(table)
    ordering <- do.call(order, unname(as.list(table)))
    table <- table[ordering, , drop = FALSE]
    rownames(table) <- NULL
    table
  })
  names(tables) <- names
  tables
}
