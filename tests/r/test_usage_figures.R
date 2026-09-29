# Tests for scripts/plot_usage_figures.R.
#
# Usage: Rscript tests/r/test_usage_figures.R scripts/plot_usage_figures.R scripts/fit_usage_model.R
#
# The statistics script is sourced only to compare its column lists; its
# packages are not loaded.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: Rscript tests/r/test_usage_figures.R scripts/plot_usage_figures.R scripts/fit_usage_model.R")
}
script <- normalizePath(arguments[[1]], mustWork = TRUE)
statistics_script <- normalizePath(arguments[[2]], mustWork = TRUE)
test_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
source(file.path(dirname(normalizePath(test_file)), "usage_model_helpers.R"))

figures <- new.env(parent = globalenv())
sys.source(script, envir = figures)
figures$load_figure_packages()
params <- list(gene_fdr = 0.05, site_fdr = 0.05, min_abs_delta_pau = 0.1, min_gene_total = 20)
work <- tempfile("pacusage-figures-")
dir.create(work)

# ---- Synthetic tables ---------------------------------------------------------

# One PAC row, typed as read_pacs() returns it.
pac_row <- function(pac_id, gene_id, strand = "+", event_type = "none",
                    assignment_class = "terminal_exon", control = 0.5, treatment = 0.5,
                    pac_fdr = NA_real_, gene_fdr = 0.5, pvalue = 0.5,
                    control_total = 100, treatment_total = 100, condition = "T",
                    control_condition = "C") {
  data.frame(
    pac_id = pac_id, gene_id = gene_id, gene_name = toupper(gene_id), strand = strand,
    locus = paste0("chr1:", pac_id), condition = condition,
    control_condition = control_condition, event_type = event_type,
    fitted_control_pau = control, fitted_treatment_pau = treatment,
    delta_pau = treatment - control, pac_fdr = pac_fdr, gene_fdr = gene_fdr,
    pvalue_pac = pvalue, assignment_class = assignment_class,
    control_gene_total = control_total, treatment_gene_total = treatment_total,
    exploratory_insufficient_replicates = FALSE, stringsAsFactors = FALSE
  )
}

pac_table <- function(...) {
  rows <- list(...)
  if (!length(rows)) return(pac_row("x", "x")[0, , drop = FALSE])
  do.call(rbind, rows)
}

gene_row <- function(gene_id, dominant_switch = FALSE, complexity_change = "none",
                     gene_fdr = 0.5) {
  data.frame(
    gene_id = gene_id, gene_name = toupper(gene_id), condition = "T",
    control_condition = "C", dominant_switch = dominant_switch,
    complexity_change = complexity_change, gene_fdr = gene_fdr,
    exploratory_insufficient_replicates = FALSE, stringsAsFactors = FALSE
  )
}

comparison <- data.frame(
  stem = "T_vs_C", condition = "T", control_condition = "C", title = "T vs C",
  stringsAsFactors = FALSE
)

write_table <- function(table, path) {
  connection <- if (grepl("\\.gz$", path)) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(connection))
  utils::write.table(table, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

png_size <- function(path) {
  bytes <- readBin(path, "raw", 24L)
  check(identical(bytes[1:8], as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))),
    basename(path), " is not a PNG.")
  number <- function(offset) sum(as.integer(bytes[offset + 0:3]) * 256^(3:0))
  c(number(17L), number(21L))
}

layer_classes <- function(layers) {
  vapply(layers, function(layer) {
    if (inherits(layer, "Scale")) "scale" else class(layer$stat)[[1]]
  }, character(1))
}

# ---- Figure data --------------------------------------------------------------

test_case("F-01", "the distal PAC is the most 3' terminal-exon or downstream PAC", {
  pacs <- pac_table(
    pac_row("p100", "plus", "+", assignment_class = "terminal_exon"),
    pac_row("p200", "plus", "+", "increased_usage", "downstream", 0.2, 0.45, 0.001, 0.01),
    pac_row("p300", "plus", "+", "decreased_usage", "intronic"),
    pac_row("m500", "minus", "-", "none", "terminal_exon"),
    pac_row("m400", "minus", "-", "lost", "terminal_exon", 0.3, 0.0, 0.002, 0.02),
    pac_row("i1", "upstream_only", "+", assignment_class = "intronic"),
    pac_row("i2", "upstream_only", "+", assignment_class = "other_exon")
  )
  coordinates <- data.frame(
    pac_id = c("p100", "p200", "p300", "m500", "m400", "i1", "i2"),
    coordinate = c(100, 200, 300, 500, 400, 50, 60)
  )
  result <- figures$distal_usage_table(pacs, coordinates)
  table <- result$table
  check(identical(names(table), figures$DISTAL_COLUMNS), "columns: ", paste(names(table), collapse = ", "))
  check(identical(table$gene_id, c("plus", "minus")), "genes: ", paste(table$gene_id, collapse = ", "))
  check(identical(table$distal_pac_id, c("p200", "m400")), "distal PACs: ", paste(table$distal_pac_id, collapse = ", "))
  check(identical(table$direction, c("lengthened", "shortened")), "directions: ", paste(table$direction, collapse = ", "))
  check(isTRUE(all.equal(table$delta_distal_pau, c(0.25, -0.3))), "delta: ", paste(table$delta_distal_pau, collapse = ", "))
  check(identical(table$distal_pac_fdr, c(0.001, 0.002)), "distal pac_fdr was not the distal row's.")
  check(identical(table$distal_assignment_class, c("downstream", "terminal_exon")), "classes.")
  check(identical(table$tested_pacs, c(3L, 2L)), "tested_pacs: ", paste(table$tested_pacs, collapse = ", "))
  check(identical(result$without_distal, 1L), "without_distal: ", result$without_distal)
})

test_case("F-02", "the direction is the distal PAC's own call", {
  events <- c(
    "gained", "increased_usage", "lost", "decreased_usage", "gained_candidate",
    "lost_candidate", "dominant_switch", "complexity_gain", "complexity_loss", "none"
  )
  genes <- sprintf("g%02d", seq_along(events))
  pacs <- do.call(rbind, lapply(seq_along(events), function(index) {
    pac_row(paste0("pac", index), genes[[index]], event_type = events[[index]], gene_fdr = index / 100)
  }))
  coordinates <- data.frame(pac_id = pacs$pac_id, coordinate = seq_along(events) * 10)
  table <- figures$distal_usage_table(pacs, coordinates)$table
  directions <- stats::setNames(table$direction, table$distal_event_type)
  expected <- c(
    gained = "lengthened", increased_usage = "lengthened", lost = "shortened",
    decreased_usage = "shortened", gained_candidate = "lengthened_candidate",
    lost_candidate = "shortened_candidate", dominant_switch = "none",
    complexity_gain = "none", complexity_loss = "none", none = "none"
  )
  check(identical(directions, expected), "directions: ", paste(names(directions), directions, sep = "=", collapse = ", "))
})

test_case("F-03", "distal rows sort by gene FDR with missing FDRs last", {
  pacs <- pac_table(
    pac_row("a", "gene_b", gene_fdr = 0.2),
    pac_row("b", "gene_na", gene_fdr = NA_real_),
    pac_row("c", "gene_a", gene_fdr = 0.2),
    pac_row("d", "gene_c", gene_fdr = 0.01)
  )
  coordinates <- data.frame(pac_id = pacs$pac_id, coordinate = 1:4)
  table <- figures$distal_usage_table(pacs, coordinates)$table
  check(identical(table$gene_id, c("gene_c", "gene_a", "gene_b", "gene_na")), "order: ", paste(table$gene_id, collapse = ", "))
  empty <- figures$distal_usage_table(pac_table(), coordinates)
  check(nrow(empty$table) == 0L && identical(names(empty$table), figures$DISTAL_COLUMNS), "empty table.")
  check(identical(empty$without_distal, 0L), "empty without_distal: ", empty$without_distal)
  missing <- tryCatch(
    figures$distal_usage_table(pac_table(pac_row("absent", "g")), coordinates),
    error = function(error) conditionMessage(error)
  )
  check(grepl("atlas lacks coordinates for PACs: absent", missing, fixed = TRUE), "missing coordinate: ", missing)
})

test_case("F-04", "site classes count confirmed calls only, with every class present", {
  pacs <- pac_table(
    pac_row("t1", "g1", event_type = "increased_usage"),
    pac_row("t2", "g1", event_type = "lost"),
    pac_row("t3", "g2", event_type = "gained_candidate"),
    pac_row("i1", "g3", event_type = "gained", assignment_class = "intronic"),
    pac_row("i2", "g3", event_type = "dominant_switch", assignment_class = "intronic"),
    pac_row("o1", "g4", event_type = "complexity_gain", assignment_class = "other_exon")
  )
  counts <- figures$site_class_counts(pacs)
  check(nrow(counts) == 8L, "rows: ", nrow(counts))
  check(identical(as.character(counts$assignment_class), rep(names(figures$SITE_CLASSES), 2)), "classes.")
  check(identical(as.character(counts$direction), rep(c("down", "up"), each = 4)), "directions.")
  check(identical(counts$count, c(1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L)), "counts: ", paste(counts$count, collapse = ", "))
  check(identical(counts$tested, rep(c(3L, 1L, 2L, 0L), 2)), "tested: ", paste(counts$tested, collapse = ", "))
  unexpected <- tryCatch(
    figures$site_class_counts(pac_table(pac_row("x", "g", assignment_class = "intergenic"))),
    error = function(error) conditionMessage(error)
  )
  check(grepl("Unexpected assignment_class values: intergenic", unexpected, fixed = TRUE), unexpected)
})

test_case("F-05", "gene events come from the genes table, PAC events from the PAC rows", {
  comparisons <- data.frame(
    stem = c("A_vs_C", "B_vs_C"), condition = c("A", "B"), control_condition = "C",
    title = c("A vs C", "B vs C"), stringsAsFactors = FALSE
  )
  pacs <- pac_table(
    pac_row("p1", "g1", event_type = "gained"),
    pac_row("p2", "g1", event_type = "gained"),
    pac_row("p3", "g2", event_type = "increased_usage"),
    pac_row("p4", "g2", event_type = "decreased_usage"),
    pac_row("p5", "g3", event_type = "gained_candidate"),
    pac_row("p6", "g3", event_type = "lost_candidate"),
    # Labels whose genes the genes table does not flag are not counted.
    pac_row("p7", "g4", event_type = "dominant_switch"),
    pac_row("p8", "g5", event_type = "complexity_gain"),
    pac_row("p9", "g5", event_type = "complexity_gain")
  )
  genes <- rbind(
    gene_row("g1", dominant_switch = TRUE),
    gene_row("g2", complexity_change = "gain"),
    gene_row("g3", dominant_switch = TRUE, complexity_change = "loss"),
    gene_row("g4"),
    gene_row("g5")
  )
  counts <- figures$event_count_table(
    comparisons,
    list(A_vs_C = pacs, B_vs_C = pac_table()),
    list(A_vs_C = genes, B_vs_C = genes[0, , drop = FALSE])
  )
  events <- c(figures$PAC_EVENTS, figures$GENE_EVENTS)
  check(identical(levels(counts$event), events), "event levels.")
  check(identical(levels(counts$level), c("PAC events", "Gene events")), "level levels.")
  check(identical(levels(counts$comparison), c("A_vs_C", "B_vs_C")), "comparison levels.")
  first <- counts[counts$comparison == "A_vs_C", , drop = FALSE]
  observed <- stats::setNames(first$count, as.character(first$event))
  expected <- c(
    gained = 2L, increased_usage = 1L, decreased_usage = 1L, lost = 0L,
    gained_candidate = 1L, lost_candidate = 1L, dominant_switch = 2L,
    complexity_gain = 1L, complexity_loss = 1L
  )
  check(identical(observed, expected), "counts: ", paste(names(observed), observed, sep = "=", collapse = ", "))
  check(identical(as.character(first$level), rep(c("PAC events", "Gene events"), c(6, 3))), "levels.")
  second <- counts[counts$comparison == "B_vs_C", , drop = FALSE]
  check(nrow(second) == 9L && all(second$count == 0L), "the empty comparison is not all zeros.")
})

test_case("F-06", "volcano data caps zero p-values, counts missing ones, and labels one PAC per gene", {
  confirmed <- lapply(1:12, function(index) {
    pac_row(sprintf("c%02d", index), sprintf("gene%02d", index), event_type = "decreased_usage",
      control = 0.6, treatment = 0.4, pvalue = 10^-(index + 2))
  })
  pacs <- do.call(pac_table, c(confirmed, list(
    pac_row("zero", "gene00", event_type = "increased_usage", control = 0.2, treatment = 0.5, pvalue = 0),
    pac_row("second", "gene00", event_type = "increased_usage", control = 0.1, treatment = 0.3, pvalue = 1e-20),
    pac_row("uncalled", "other", pvalue = 0.5),
    pac_row("candidate", "other", event_type = "gained_candidate", control = 0, treatment = 0.3, pvalue = 1e-4),
    pac_row("unstable", "other", event_type = "lost_candidate", control = 0.3, treatment = 0, pvalue = NA_real_),
    pac_row("no_p", "other", pvalue = NA_real_),
    pac_row("no_delta", "other", control = NA_real_, pvalue = 0.3)
  )))
  data <- figures$volcano_data(pacs)
  points <- data$points
  check(data$missing == 3L, "missing: ", data$missing)
  check(data$missing_candidates == 1L, "missing candidates: ", data$missing_candidates)
  check(nrow(points) == nrow(pacs) - 3L, "points: ", nrow(points))
  cap <- max(1, -log10(points$pvalue[points$pvalue > 0])) * 1.05
  zero <- points[points$pac_id == "zero", , drop = FALSE]
  check(isTRUE(zero$capped) && zero$neg_log10_p == cap, "the zero p-value was not capped at ", cap)
  check(sum(points$capped) == 1L, "capped rows: ", sum(points$capped))
  check(identical(as.character(points$call[1]), "none"), "uncalled PACs are not drawn first.")
  check(identical(as.character(points$call[points$pac_id == "candidate"]), "up_candidate"), "candidate call.")
  # One label per gene, on its most significant confirmed PAC; candidates and
  # uncalled PACs are not labelled.
  expected_labels <- c("zero", sprintf("c%02d", 12:1))
  check(identical(data$labels$pac_id, expected_labels), "labels: ", paste(data$labels$pac_id, collapse = ", "))
})

test_case("F-18", "volcano labels stop at 20 genes in each direction", {
  up <- lapply(1:25, function(index) {
    pac_row(sprintf("u%02d", index), sprintf("up%02d", index), event_type = "increased_usage",
      control = 0.3, treatment = 0.6, pvalue = 10^-(2 * index))
  })
  down <- lapply(1:25, function(index) {
    pac_row(sprintf("d%02d", index), sprintf("down%02d", index), event_type = "lost",
      control = 0.4, treatment = 0, pvalue = 10^-(2 * index + 1))
  })
  labels <- figures$volcano_data(do.call(pac_table, c(up, down)))$labels
  check(sum(labels$call == "up") == 20L && sum(labels$call == "down") == 20L,
    "labels per direction: ", paste(table(labels$call), collapse = ", "))
  check(setequal(labels$pac_id, c(sprintf("u%02d", 6:25), sprintf("d%02d", 6:25))),
    "the most significant genes are not the ones labelled: ", paste(labels$pac_id, collapse = ", "))
  check(!is.unsorted(labels$pvalue), "labels are not in order of significance.")
})

test_case("F-19", "distal labels go on the 20 largest changes in each direction", {
  lengthened <- lapply(1:22, function(index) {
    pac_row(sprintf("l%02d", index), sprintf("long%02d", index), event_type = "increased_usage",
      control = 0.2, treatment = 0.2 + index / 100, pac_fdr = 0.01)
  })
  pacs <- do.call(pac_table, c(lengthened, list(
    pac_row("s1", "short1", event_type = "decreased_usage", control = 0.6, treatment = 0.1, pac_fdr = 0.01),
    pac_row("s2", "short2", event_type = "lost", control = 0.3, treatment = 0.0, pac_fdr = 0.01),
    pac_row("c1", "candidate", event_type = "gained_candidate", control = 0, treatment = 0.9),
    pac_row("n1", "unchanged", control = 0.1, treatment = 0.95)
  )))
  coordinates <- data.frame(pac_id = pacs$pac_id, coordinate = seq_len(nrow(pacs)))
  table <- figures$distal_usage_table(pacs, coordinates)$table
  labels <- figures$distal_labels(table)
  # By size of change: short1 (0.5), short2 (0.3), then long22 (0.22) down to
  # long03, the 20th lengthened gene; long01, long02, and the candidate are left.
  check(identical(labels$gene_id, c("short1", "short2", sprintf("long%02d", 22:3))),
    "labels: ", paste(labels$gene_id, collapse = ", "))
})

test_case("F-20", "the coverage figure counts the PACs it leaves out", {
  comparisons <- data.frame(stem = "T_vs_C", title = "T vs C", stringsAsFactors = FALSE)
  data <- figures$coverage_data(comparisons, list(T_vs_C = pac_table(pac_row("p1", "g1"))))
  subtitle <- function(dropped) {
    figures$plot_effect_vs_coverage(data, comparisons, params, dropped)$labels$subtitle
  }
  check(grepl("\n3 PACs without reads in both groups, not shown$", subtitle(3)), subtitle(3))
  check(grepl("\n1 PAC without reads in both groups", subtitle(1), fixed = TRUE), subtitle(1))
  check(!grepl("\n", subtitle(0), fixed = TRUE), subtitle(0))
})

test_case("F-07", "coverage rows need reads in both groups and a change, with no pseudocount", {
  comparisons <- data.frame(stem = c("A_vs_C", "B_vs_C"), stringsAsFactors = FALSE)
  pacs <- pac_table(
    pac_row("kept", "g1", event_type = "increased_usage", control_total = 10, treatment_total = 50),
    pac_row("kept_2", "g1", control_total = 400, treatment_total = 35),
    pac_row("no_control", "g2", control_total = 0, treatment_total = 100),
    pac_row("missing_total", "g3", control_total = NA_real_, treatment_total = 20),
    pac_row("no_delta", "g4", control = NA_real_)
  )
  data <- figures$coverage_data(comparisons, list(A_vs_C = pacs, B_vs_C = pac_table()))
  check(identical(data$pac_id, c("kept_2", "kept")), "rows: ", paste(data$pac_id, collapse = ", "))
  check(identical(data$depth, c(35, 10)), "depth: ", paste(data$depth, collapse = ", "))
  check(identical(levels(data$comparison), c("A_vs_C", "B_vs_C")), "comparison levels.")
  check(identical(as.character(data$call), c("none", "up")), "calls: ", paste(data$call, collapse = ", "))
})

test_case("F-08", "uncalled PACs switch from points to a density above 5,000 per panel", {
  uncalled <- function(count, panel = "A") {
    data.frame(
      delta_pau = seq_len(count) / max(count, 1), neg_log10_p = rep(1, count),
      call = factor(rep("none", count), levels = figures$CALLS),
      comparison = factor(rep(panel, count), levels = c("A", "B"))
    )
  }
  few <- layer_classes(figures$uncalled_layers(uncalled(5000)))
  check(identical(few, "StatIdentity"), "5,000 uncalled PACs: ", paste(few, collapse = ", "))
  many <- layer_classes(figures$uncalled_layers(uncalled(5001)))
  check(identical(many, c("StatBin2d", "scale")), "5,001 uncalled PACs: ", paste(many, collapse = ", "))
  mixed <- layer_classes(figures$uncalled_layers(rbind(uncalled(5001, "A"), uncalled(10, "B")), panel = "comparison"))
  check(identical(mixed, c("StatIdentity", "StatBin2d", "scale")), "per panel: ", paste(mixed, collapse = ", "))
  check(length(figures$uncalled_layers(uncalled(0))) == 0L, "no uncalled PACs still made a layer.")
})

test_case("F-17", "comparisons are each condition against its control, sorted as bytes", {
  samples <- data.frame(
    sample_id = paste0("s", 1:7),
    condition = c("C", "C", "b", "B", "a", "a", "B"),
    control_condition = c("C", "C", "A", "a", "B", "B", "a"),
    stringsAsFactors = FALSE
  )
  comparisons <- figures$comparisons_from_samples(samples)
  check(identical(comparisons$stem, c("B_vs_a", "a_vs_B", "b_vs_A")), "stems: ", paste(comparisons$stem, collapse = ", "))
  check(identical(comparisons$title, c("B vs a", "a vs B", "b vs A")), "titles.")
  roots <- tryCatch(
    figures$comparisons_from_samples(samples[1:2, , drop = FALSE]),
    error = function(error) conditionMessage(error)
  )
  check(grepl("no treatment-control comparison", roots, fixed = TRUE), roots)
})

# ---- Saving -------------------------------------------------------------------

test_case("F-09", "PDFs carry no dates or producer, and PNGs have the requested size", {
  stem <- file.path(work, "metadata")
  figures$save_figure(figures$placeholder_plot("Title", "Message"), stem, 3, 2)
  pdf <- readBin(paste0(stem, ".pdf"), "raw", file.size(paste0(stem, ".pdf")))
  for (field in c("/CreationDate", "/ModDate", "/Producer")) {
    check(!length(grepRaw(field, pdf, fixed = TRUE)), "the PDF records ", field)
  }
  check(length(grepRaw("/Title (metadata)", pdf, fixed = TRUE)) == 1L, "the PDF title is not the figure name.")
  check(identical(png_size(paste0(stem, ".png")), c(600, 400)), "PNG size: ", paste(png_size(paste0(stem, ".png")), collapse = "x"))
})

test_case("F-10", "separate processes write byte-identical PDFs and PNGs", {
  driver <- file.path(work, "render.R")
  writeLines(c(
    "arguments <- commandArgs(trailingOnly = TRUE)",
    "figures <- new.env(parent = globalenv())",
    "sys.source(arguments[[1]], envir = figures)",
    "figures$load_figure_packages()",
    "points <- data.frame(x = 1:20, y = (1:20)^2, label = paste0('gene', 1:20))",
    "plot <- ggplot2::ggplot(points, ggplot2::aes(x, y)) + ggplot2::geom_point() +",
    "  ggplot2::geom_text(ggplot2::aes(label = label), check_overlap = TRUE) +",
    "  ggplot2::labs(title = 'Rendering test') + figures$figure_theme()",
    "figures$save_figure(plot, arguments[[2]], 4, 3)"
  ), driver)
  stems <- file.path(work, c("first", "second"), "render")
  for (stem in stems) {
    dir.create(dirname(stem))
    run_rscript(c(driver, script, stem))
  }
  for (suffix in c(".pdf", ".png")) {
    files <- paste0(stems, suffix)
    bytes <- lapply(files, function(path) readBin(path, "raw", file.size(path)))
    check(identical(bytes[[1]], bytes[[2]]), suffix, " differs between processes.")
  }
})

test_case("F-11", "a percent sign in a comparison name is written literally", {
  stem <- file.path(work, "Drug_10%_vs_DMSO.volcano")
  figures$save_figure(figures$placeholder_plot("Drug 10% vs DMSO", "Message"), stem, 2, 2)
  check(file.exists(paste0(stem, ".pdf")) && file.exists(paste0(stem, ".png")), "files are missing.")
})

# ---- End to end ---------------------------------------------------------------

# Two treatments of one control, each tested against it.
write_inputs <- function(directory, filled) {
  statistics <- file.path(directory, "statistics")
  dir.create(statistics, recursive = TRUE)
  samples <- data.frame(
    sample_id = paste0("s", 1:6),
    condition = c("C", "C", "T1", "T1", "T2", "T2"),
    control_condition = c("C", "C", "C", "C", "C", "C"),
    stringsAsFactors = FALSE
  )
  write_table(samples, file.path(directory, "samples.tsv"))
  pacs <- pac_table(
    pac_row("g1_p1", "g1", "+", "decreased_usage", control = 0.7, treatment = 0.3, pac_fdr = 0.001, gene_fdr = 0.001, pvalue = 1e-6, condition = "T1"),
    pac_row("g1_p2", "g1", "+", "increased_usage", control = 0.3, treatment = 0.7, pac_fdr = 0.001, gene_fdr = 0.001, pvalue = 1e-6, condition = "T1"),
    pac_row("g2_p1", "g2", "-", "none", control = 0.5, treatment = 0.55, gene_fdr = 0.8, pvalue = 0.7, condition = "T1"),
    pac_row("g2_p2", "g2", "-", "none", assignment_class = "intronic", control = 0.5, treatment = 0.45, gene_fdr = 0.8, pvalue = 0.7, condition = "T1")
  )
  genes <- rbind(gene_row("g1", dominant_switch = TRUE, gene_fdr = 0.001), gene_row("g2", gene_fdr = 0.8))
  genes$condition <- "T1"
  model <- new.env(parent = globalenv())
  sys.source(statistics_script, envir = model)
  empty_pacs <- model$empty_table(model$PAC_COLUMNS)
  empty_genes <- model$empty_table(model$GENE_COLUMNS)
  for (stem in c("T1_vs_C", "T2_vs_C")) {
    use <- filled && stem == "T1_vs_C"
    write_table(if (use) pacs else empty_pacs, file.path(statistics, paste0(stem, ".pacs.tsv.gz")))
    write_table(if (use) genes else empty_genes, file.path(statistics, paste0(stem, ".genes.tsv.gz")))
  }
  write_table(
    data.frame(pac_id = pacs$pac_id, gene_id = pacs$gene_id, coordinate = c(100, 200, 900, 800),
      upstream_sequence = "ACGT"),
    file.path(directory, "atlas.tsv.gz")
  )
  yaml::write_yaml(c(params, list(outdir = "results")), file.path(directory, "params.yaml"))
  c(
    "--mode", "figures", "--samples", file.path(directory, "samples.tsv"),
    "--statistics-dir", statistics, "--atlas", file.path(directory, "atlas.tsv.gz"),
    "--params", file.path(directory, "params.yaml"), "--output-dir", file.path(directory, "figures")
  )
}

expected_files <- sort(c(
  as.vector(outer(c("T1_vs_C", "T2_vs_C"), c(
    ".volcano.pdf", ".volcano.png", ".distal_usage.pdf", ".distal_usage.png",
    ".distal_usage.tsv.gz", ".site_classes.pdf", ".site_classes.png"
  ), paste0)),
  "event_counts.pdf", "event_counts.png", "effect_vs_coverage.pdf", "effect_vs_coverage.png"
), method = "radix")

test_case("F-12", "the command line writes every figure and the distal-usage tables", {
  directory <- file.path(work, "filled")
  run_rscript(c(script, write_inputs(directory, filled = TRUE)))
  files <- sort(list.files(file.path(directory, "figures")), method = "radix")
  check(identical(files, expected_files), "files: ", paste(files, collapse = ", "))
  distal <- utils::read.delim(
    file.path(directory, "figures", "T1_vs_C.distal_usage.tsv.gz"),
    colClasses = "character", na.strings = character()
  )
  check(identical(names(distal), figures$DISTAL_COLUMNS), "columns: ", paste(names(distal), collapse = ", "))
  # g1 (+) ends at g1_p2, which gained usage; g2 (-) ends at g2_p1, since its
  # other PAC is intronic.
  check(identical(distal$gene_id, c("g1", "g2")), "genes: ", paste(distal$gene_id, collapse = ", "))
  check(identical(distal$distal_pac_id, c("g1_p2", "g2_p1")), "distal PACs.")
  check(identical(distal$direction, c("lengthened", "none")), "directions: ", paste(distal$direction, collapse = ", "))
  check(identical(distal$tested_pacs, c("2", "2")), "tested_pacs.")
  check(identical(distal$distal_pac_fdr, c("0.001", "")), "pac_fdr: ", paste(distal$distal_pac_fdr, collapse = ", "))
  empty <- utils::read.delim(file.path(directory, "figures", "T2_vs_C.distal_usage.tsv.gz"))
  check(nrow(empty) == 0L && identical(names(empty), figures$DISTAL_COLUMNS), "the empty comparison's table.")
  for (name in c("event_counts.png", "effect_vs_coverage.png", "T1_vs_C.volcano.png")) {
    check(file.size(file.path(directory, "figures", name)) > 1000, name, " is nearly empty.")
  }
  check(identical(png_size(file.path(directory, "figures", "event_counts.png")), c(1500, 550)), "event_counts size.")
  check(identical(png_size(file.path(directory, "figures", "effect_vs_coverage.png")), c(1500, 650)), "coverage size.")
})

test_case("F-13", "every comparison empty still writes the same files", {
  directory <- file.path(work, "empty")
  command <- write_inputs(directory, filled = FALSE)
  figures$run_figures_mode(figures$parse_args(command))
  files <- sort(list.files(file.path(directory, "figures")), method = "radix")
  check(identical(files, expected_files), "files: ", paste(files, collapse = ", "))
})

test_case("F-14", "versions mode appends ggplot2, or starts the table", {
  path <- file.path(work, "versions.tsv")
  writeLines(c("software\tversion", "R\t4.5.3", "DRIMSeq\t1.38.0"), path)
  version <- as.character(utils::packageVersion("ggplot2"))
  figures$run_versions_mode(list(output = path), version)
  table <- utils::read.delim(path, colClasses = "character")
  check(identical(table$software, c("R", "DRIMSeq", "ggplot2")), "rows: ", paste(table$software, collapse = ", "))
  check(identical(table$version[[3]], version), "version: ", table$version[[3]])
  fresh <- file.path(work, "fresh_versions.tsv")
  figures$run_versions_mode(list(output = fresh), version)
  check(identical(readLines(fresh), c("software\tversion", paste0("ggplot2\t", version))), "new table: ", paste(readLines(fresh), collapse = " | "))
})

test_case("F-15", "the figures read only columns the statistics tables have", {
  model <- new.env(parent = globalenv())
  sys.source(statistics_script, envir = model)
  missing_pacs <- setdiff(figures$PACS_REQUIRED, model$PAC_COLUMNS)
  missing_genes <- setdiff(figures$GENES_REQUIRED, model$GENE_COLUMNS)
  check(!length(missing_pacs), "not in PAC_COLUMNS: ", paste(missing_pacs, collapse = ", "))
  check(!length(missing_genes), "not in GENE_COLUMNS: ", paste(missing_genes, collapse = ", "))
})

test_case("F-16", "the figure script is ASCII only", {
  bytes <- readBin(script, "raw", file.size(script))
  check(all(as.integer(bytes) < 128L), sum(as.integer(bytes) >= 128L), " non-ASCII bytes.")
})

finish_tests()
