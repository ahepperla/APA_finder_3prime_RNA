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
                    gene_region = "last_exon", control = 0.5, treatment = 0.5,
                    pac_fdr = NA_real_, gene_fdr = 0.5, pvalue = 0.5,
                    control_total = 100, treatment_total = 100, condition = "T",
                    control_condition = "C") {
  data.frame(
    pac_id = pac_id, gene_id = gene_id, gene_name = toupper(gene_id), strand = strand,
    locus = paste0("chr1:", pac_id), condition = condition,
    control_condition = control_condition, event_type = event_type,
    fitted_control_pau = control, fitted_treatment_pau = treatment,
    delta_pau = treatment - control, pac_fdr = pac_fdr, gene_fdr = gene_fdr,
    pac_pvalue = pvalue, gene_region = gene_region,
    control_gene_total = control_total, treatment_gene_total = treatment_total,
    exploratory_insufficient_replicates = FALSE, stringsAsFactors = FALSE
  )
}

pac_table <- function(...) {
  rows <- list(...)
  if (!length(rows)) return(pac_row("x", "x")[0, , drop = FALSE])
  do.call(rbind, rows)
}

SHIFT_FIELDS <- c(
  "shift_direction", "shift_from_pac_id", "shift_from_gene_region", "shift_from_event_type",
  "shift_to_pac_id", "shift_to_gene_region", "shift_to_event_type"
)

# One genes row; shift holds its shift's fields in SHIFT_FIELDS order.
gene_row <- function(gene_id, dominant_switch = FALSE, active_pacs_change = "none",
                     gene_fdr = 0.5, apa_pattern = "none", shift = rep(NA_character_, 7)) {
  row <- data.frame(
    gene_id = gene_id, gene_name = toupper(gene_id), condition = "T",
    control_condition = "C", dominant_switch = dominant_switch,
    active_pacs_change = active_pacs_change, apa_pattern = apa_pattern, gene_fdr = gene_fdr,
    exploratory_insufficient_replicates = FALSE, stringsAsFactors = FALSE
  )
  row[SHIFT_FIELDS] <- as.list(shift)
  row
}

# A genes row whose usage moved from one PAC to another, with the regions and
# event types given.
shift_row <- function(gene_id, direction, from, to, from_region = "last_exon",
                      to_region = "last_exon", from_event = "decreased_usage",
                      to_event = "increased_usage") {
  gene_row(gene_id, gene_fdr = 0.01, shift = c(direction, from, from_region, from_event, to,
    to_region, to_event))
}

# A genes table for a PAC table's genes, with the given APA patterns by gene.
genes_of <- function(pacs, patterns = character()) {
  ids <- unique(pacs$gene_id)
  rows <- lapply(ids, function(id) {
    gene_row(id, apa_pattern = if (id %in% names(patterns)) patterns[[id]] else "none")
  })
  if (!length(rows)) return(gene_row("x")[0, , drop = FALSE])
  do.call(rbind, rows)
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

test_case("F-01", "the distal PAC is the most 3' last-exon or downstream PAC", {
  pacs <- pac_table(
    pac_row("p100", "plus", "+", gene_region = "last_exon"),
    pac_row("p200", "plus", "+", "increased_usage", "downstream_of_gene", 0.2, 0.45, 0.001, 0.01),
    pac_row("p300", "plus", "+", "decreased_usage", "intron"),
    pac_row("m500", "minus", "-", "none", "last_exon"),
    pac_row("m400", "minus", "-", "lost", "last_exon", 0.3, 0.0, 0.002, 0.02),
    pac_row("i1", "upstream_only", "+", gene_region = "intron"),
    pac_row("i2", "upstream_only", "+", gene_region = "internal_exon")
  )
  coordinates <- data.frame(
    pac_id = c("p100", "p200", "p300", "m500", "m400", "i1", "i2"),
    coordinate = c(100, 200, 300, 500, 400, 50, 60)
  )
  result <- figures$distal_usage_table(pacs, coordinates, genes_of(pacs, c(plus = "utr_lengthening")))
  table <- result$table
  check(identical(names(table), figures$DISTAL_COLUMNS), "columns: ", paste(names(table), collapse = ", "))
  check(identical(table$gene_id, c("plus", "minus")), "genes: ", paste(table$gene_id, collapse = ", "))
  check(identical(table$distal_pac_id, c("p200", "m400")), "distal PACs: ", paste(table$distal_pac_id, collapse = ", "))
  check(identical(table$direction, c("distal_up", "distal_down")), "directions: ", paste(table$direction, collapse = ", "))
  check(identical(table$apa_pattern, c("utr_lengthening", "none")), "patterns: ", paste(table$apa_pattern, collapse = ", "))
  check(isTRUE(all.equal(table$delta_distal_pau, c(0.25, -0.3))), "delta: ", paste(table$delta_distal_pau, collapse = ", "))
  check(identical(table$distal_pac_fdr, c(0.001, 0.002)), "distal pac_fdr was not the distal row's.")
  check(identical(table$distal_gene_region, c("downstream_of_gene", "last_exon")), "regions.")
  check(identical(table$tested_pacs, c(3L, 2L)), "tested_pacs: ", paste(table$tested_pacs, collapse = ", "))
  check(identical(result$without_distal, 1L), "without_distal: ", result$without_distal)
})

test_case("F-02", "the direction is the distal PAC's own call", {
  events <- c(
    "gained", "increased_usage", "lost", "decreased_usage", "gained_candidate",
    "lost_candidate", "none"
  )
  genes <- sprintf("g%02d", seq_along(events))
  pacs <- do.call(rbind, lapply(seq_along(events), function(index) {
    pac_row(paste0("pac", index), genes[[index]], event_type = events[[index]], gene_fdr = index / 100)
  }))
  coordinates <- data.frame(pac_id = pacs$pac_id, coordinate = seq_along(events) * 10)
  table <- figures$distal_usage_table(pacs, coordinates, genes_of(pacs))$table
  directions <- stats::setNames(table$direction, table$distal_event_type)
  expected <- c(
    gained = "distal_up", increased_usage = "distal_up", lost = "distal_down",
    decreased_usage = "distal_down", gained_candidate = "distal_up_candidate",
    lost_candidate = "distal_down_candidate", none = "none"
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
  table <- figures$distal_usage_table(pacs, coordinates, genes_of(pacs))$table
  check(identical(table$gene_id, c("gene_c", "gene_a", "gene_b", "gene_na")), "order: ", paste(table$gene_id, collapse = ", "))
  empty <- figures$distal_usage_table(pac_table(), coordinates, genes_of(pac_table()))
  check(nrow(empty$table) == 0L && identical(names(empty$table), figures$DISTAL_COLUMNS), "empty table.")
  check(identical(empty$without_distal, 0L), "empty without_distal: ", empty$without_distal)
  missing <- tryCatch(
    figures$distal_usage_table(pac_table(pac_row("absent", "g")), coordinates, gene_row("g")),
    error = function(error) conditionMessage(error)
  )
  check(grepl("atlas lacks coordinates for PACs: absent", missing, fixed = TRUE), "missing coordinate: ", missing)
})

test_case("F-04", "shifts count each gene once, by region pair and direction, with shared rows", {
  a <- rbind(
    shift_row("g1", "proximal", "g1_from", "g1_to"),
    shift_row("g2", "distal", "g2_from", "g2_to"),
    # The to side has no call of its own; the gene still counts once.
    shift_row("g3", "distal", "g3_from", "g3_to", to_event = "none"),
    shift_row("g4", "proximal", "g4_from", "g4_to", to_region = "intron"),
    gene_row("g5"),
    # An empty field, as read from a table, is no shift.
    gene_row("g6", shift = rep("", 7))
  )
  b <- rbind(
    shift_row("h1", "distal", "h1_from", "h1_to", from_region = "intron"),
    shift_row("h2", "proximal", "h2_from", "h2_to", from_region = "downstream_of_gene",
      to_region = "internal_exon")
  )
  # The pairs of both comparisons, in region order: from, then to.
  pairs <- figures$shift_pairs(list(A = a, B = b))
  expected_pairs <- c(
    "last_exon|last_exon", "last_exon|intron", "intron|last_exon",
    "downstream_of_gene|internal_exon"
  )
  check(identical(pairs, expected_pairs), "pairs: ", paste(pairs, collapse = ", "))
  counts <- figures$shift_region_counts(a, pairs)
  check(identical(as.character(counts$pair), rep(pairs, 2)), "pair rows.")
  check(identical(as.character(counts$direction), rep(c("proximal", "distal"), each = 4)), "directions.")
  check(identical(counts$count, c(1L, 1L, 0L, 0L, 2L, 0L, 0L, 0L)), "counts: ", paste(counts$count, collapse = ", "))
  none <- figures$shift_pairs(list(A = rbind(gene_row("x"), gene_row("y", shift = rep("", 7)))))
  check(identical(none, character()), "pairs without shifts: ", paste(none, collapse = ", "))
  check(nrow(figures$shift_region_counts(gene_row("x"), none)) == 0L, "counts without pairs.")
  unexpected <- tryCatch(
    figures$shift_pairs(list(A = shift_row("x", "distal", "p1", "p2", from_region = "intergenic"))),
    error = function(error) conditionMessage(error)
  )
  check(grepl("Unexpected shift gene_region values: intergenic", unexpected, fixed = TRUE), unexpected)
})

test_case("F-39", "the shifts figure draws every pair's genes either side of zero, and grows with the pairs", {
  pairs <- c("last_exon|last_exon", "last_exon|intron")
  genes <- rbind(
    shift_row("g1", "proximal", "a", "b"),
    shift_row("g2", "proximal", "c", "d"),
    shift_row("g3", "distal", "e", "f"),
    shift_row("g4", "proximal", "g", "h", to_region = "intron")
  )
  plot <- figures$plot_shifts_by_gene_region(figures$shift_region_counts(genes, pairs), comparison)
  data <- plot$data
  check(identical(levels(data$label), c("Last exon -> Intron", "Last exon -> Last exon")), "rows: ", paste(levels(data$label), collapse = ", "))
  check(identical(data$signed, c(-2L, -1L, 1L, 0L)), "signed counts: ", paste(data$signed, collapse = ", "))
  check(grepl("^4 genes with a confirmed shift", plot$labels$subtitle), "subtitle: ", plot$labels$subtitle)
  check(grepl("3 to a more proximal PAC, 1 to a more distal one$", plot$labels$subtitle), "subtitle: ", plot$labels$subtitle)
  # A comparison without shifts keeps the rows, with a note.
  quiet <- figures$plot_shifts_by_gene_region(figures$shift_region_counts(gene_row("x"), pairs), comparison)
  check(identical(levels(quiet$data$label), levels(data$label)), "the empty comparison lost its rows.")
  check(any(vapply(quiet$layers, function(layer) identical(layer$aes_params$label, "No gene has a confirmed shift"), logical(1))), "no note on the empty comparison.")
  # No shift in any comparison: no rows, and a placeholder.
  placeholder <- figures$plot_shifts_by_gene_region(figures$shift_region_counts(gene_row("x"), character()), comparison)
  note <- vapply(placeholder$layers, function(layer) as.character(layer$aes_params$label %||% ""), character(1))
  check(identical(unname(note), "No gene has a confirmed shift in any comparison"), "placeholder: ", paste(note, collapse = ", "))
  check(identical(figures$shift_figure_height(character()), 2.5), "height without pairs.")
  check(identical(figures$shift_figure_height(pairs), 2.75), "height for two pairs.")
})

test_case("F-05", "each gene counts once by its PAC calls; gene events come from the genes table", {
  comparisons <- data.frame(
    stem = c("A_vs_C", "B_vs_C"), condition = c("A", "B"), control_condition = "C",
    title = c("A vs C", "B vs C"), stringsAsFactors = FALSE
  )
  pacs <- pac_table(
    pac_row("p1", "g1", event_type = "gained"),
    pac_row("p2", "g1", event_type = "lost"),
    pac_row("p3", "g2", event_type = "gained"),
    pac_row("p4", "g2", event_type = "gained"),
    pac_row("p5", "g2", event_type = "decreased_usage"),
    pac_row("p6", "g3", event_type = "lost"),
    pac_row("p7", "g3", event_type = "increased_usage"),
    pac_row("p8", "g3", event_type = "gained_candidate"),
    pac_row("p9", "g4", event_type = "increased_usage"),
    pac_row("p10", "g4", event_type = "decreased_usage"),
    pac_row("p11", "g5", event_type = "decreased_usage"),
    pac_row("p12", "g5", event_type = "lost_candidate"),
    pac_row("p13", "g6", event_type = "gained_candidate"),
    pac_row("p14", "g6", event_type = "lost_candidate"),
    pac_row("p15", "g7")
  )
  genes <- rbind(
    gene_row("g1", dominant_switch = TRUE),
    gene_row("g2", active_pacs_change = "more"),
    gene_row("g3", dominant_switch = TRUE, active_pacs_change = "fewer"),
    gene_row("g4"),
    gene_row("g5")
  )
  counts <- figures$event_count_table(
    comparisons,
    list(A_vs_C = pacs, B_vs_C = pac_table()),
    list(A_vs_C = genes, B_vs_C = genes[0, , drop = FALSE])
  )
  events <- c(figures$GENE_CALLS, figures$GENE_EVENTS)
  check(identical(levels(counts$event), events), "event levels.")
  check(identical(levels(counts$level), c("Genes by PAC calls", "Gene events")), "level levels.")
  check(identical(levels(counts$comparison), c("A_vs_C", "B_vs_C")), "comparison levels.")
  first <- counts[counts$comparison == "A_vs_C", , drop = FALSE]
  observed <- stats::setNames(first$count, as.character(first$event))
  # g1 gains and loses; g2 gains twice and g3 loses, each once; g4 and g5
  # change usage, a candidate aside; g6 has candidates only; g7 nothing.
  expected <- c(
    gained_and_lost = 1L, gained = 1L, lost = 1L, changed_usage = 2L, candidates_only = 1L,
    dominant_switch = 2L, more_active_pacs = 1L, fewer_active_pacs = 1L
  )
  check(identical(observed, expected), "counts: ", paste(names(observed), observed, sep = "=", collapse = ", "))
  check(identical(as.character(first$level), rep(c("Genes by PAC calls", "Gene events"), c(5, 3))), "levels.")
  second <- counts[counts$comparison == "B_vs_C", , drop = FALSE]
  check(nrow(second) == 8L && all(second$count == 0L), "the empty comparison is not all zeros.")
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
  # Each geneNN loses usage at cNN to a PAC without a call; gene00 moves
  # usage to zero.
  genes <- do.call(rbind, c(
    lapply(1:12, function(index) {
      shift_row(sprintf("gene%02d", index), "distal", sprintf("c%02d", index),
        sprintf("u%02d", index), to_event = "none")
    }),
    list(shift_row("gene00", "distal", "elsewhere", "zero"), gene_row("other"))
  ))
  data <- figures$volcano_data(pacs, genes)
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
  # One label per gene, at the confirmed PAC its usage moved to, or else
  # moved from; candidates and uncalled PACs are not labelled.
  expected_labels <- c("zero", sprintf("c%02d", 12:1))
  check(identical(data$labels$pac_id, expected_labels), "labels: ", paste(data$labels$pac_id, collapse = ", "))
})

test_case("F-18", "volcano labels stop at 20 genes on each side", {
  up <- lapply(1:25, function(index) {
    pac_row(sprintf("u%02d", index), sprintf("up%02d", index), event_type = "increased_usage",
      control = 0.3, treatment = 0.6, pvalue = 10^-(2 * index))
  })
  down <- lapply(1:25, function(index) {
    pac_row(sprintf("d%02d", index), sprintf("down%02d", index), event_type = "lost",
      control = 0.4, treatment = 0, pvalue = 10^-(2 * index + 1))
  })
  # The upNN genes gain at uNN; the downNN genes lose at dNN, with no
  # confirmed gain, so they are labelled there.
  genes <- do.call(rbind, c(
    lapply(1:25, function(index) {
      shift_row(sprintf("up%02d", index), "distal", sprintf("x%02d", index), sprintf("u%02d", index),
        from_event = "none")
    }),
    lapply(1:25, function(index) {
      shift_row(sprintf("down%02d", index), "distal", sprintf("d%02d", index),
        sprintf("y%02d", index), from_event = "lost", to_event = "none")
    })
  ))
  labels <- figures$volcano_data(do.call(pac_table, c(up, down)), genes)$labels
  check(sum(labels$call == "up") == 20L && sum(labels$call == "down") == 20L,
    "labels per direction: ", paste(table(labels$call), collapse = ", "))
  check(setequal(labels$pac_id, c(sprintf("u%02d", 6:25), sprintf("d%02d", 6:25))),
    "the most significant genes are not the ones labelled: ", paste(labels$pac_id, collapse = ", "))
  check(!is.unsorted(labels$pvalue), "labels are not in order of significance.")
})

test_case("F-19", "distal labels go on the 20 largest changes in each direction", {
  rising <- lapply(1:22, function(index) {
    pac_row(sprintf("l%02d", index), sprintf("long%02d", index), event_type = "increased_usage",
      control = 0.2, treatment = 0.2 + index / 100, pac_fdr = 0.01)
  })
  pacs <- do.call(pac_table, c(rising, list(
    pac_row("s1", "short1", event_type = "decreased_usage", control = 0.6, treatment = 0.1, pac_fdr = 0.01),
    pac_row("s2", "short2", event_type = "lost", control = 0.3, treatment = 0.0, pac_fdr = 0.01),
    pac_row("c1", "candidate", event_type = "gained_candidate", control = 0, treatment = 0.9),
    pac_row("n1", "unchanged", control = 0.1, treatment = 0.95)
  )))
  coordinates <- data.frame(pac_id = pacs$pac_id, coordinate = seq_len(nrow(pacs)))
  table <- figures$distal_usage_table(pacs, coordinates, genes_of(pacs))$table
  labels <- figures$distal_labels(table)
  # By size of change: short1 (0.5), short2 (0.3), then long22 (0.22) down to
  # long03, the 20th rising gene; long01, long02, and the candidate are left.
  check(identical(labels$gene_id, c("short1", "short2", sprintf("long%02d", 22:3))),
    "labels: ", paste(labels$gene_id, collapse = ", "))
})

test_case("F-21", "pattern counts give each class its genes with all tested PACs, and apart by flagged support", {
  comparisons <- data.frame(stem = c("A_vs_C", "B_vs_C"), stringsAsFactors = FALSE)
  potential <- function(pattern) paste0(pattern, "_potential_internal_priming")
  genes <- rbind(
    gene_row("g1", apa_pattern = "intronic_gain;utr_shortening"),
    gene_row("g2", apa_pattern = "intronic_gain"),
    gene_row("g3", apa_pattern = "alternative_last_exon"),
    gene_row("g4", apa_pattern = "unclassified_change"),
    gene_row("g5", apa_pattern = "none"),
    gene_row("g6", apa_pattern = potential("intronic_gain")),
    gene_row("g7", apa_pattern = paste0("utr_shortening;", potential("intronic_loss")))
  )
  counts <- figures$pattern_count_table(
    comparisons, list(A_vs_C = genes, B_vs_C = genes[0, , drop = FALSE])
  )
  classes <- c(
    "intronic_gain", "intronic_loss", "alternative_last_exon", "utr_shortening",
    "utr_lengthening", "unclassified_change"
  )
  check(identical(levels(counts$pattern), classes), "pattern levels.")
  check(identical(levels(counts$support), c("all", "supported", "flagged")), "support levels.")
  first <- counts[counts$comparison == "A_vs_C", , drop = FALSE]
  observed <- stats::setNames(first$count, paste(first$support, first$pattern))
  # With all tested PACs, g6's flagged intronic gain joins g1's and g2's, and
  # g7 counts for both its UTR shortening and its flagged intronic loss.
  expected <- c(
    "all intronic_gain" = 3L, "all intronic_loss" = 1L, "all alternative_last_exon" = 1L,
    "all utr_shortening" = 2L, "all utr_lengthening" = 0L, "all unclassified_change" = 1L,
    "supported intronic_gain" = 2L, "supported intronic_loss" = 0L,
    "supported alternative_last_exon" = 1L, "supported utr_shortening" = 2L,
    "supported utr_lengthening" = 0L, "supported unclassified_change" = 1L,
    "flagged intronic_gain" = 1L, "flagged intronic_loss" = 1L,
    "flagged alternative_last_exon" = 0L, "flagged utr_shortening" = 0L,
    "flagged utr_lengthening" = 0L
  )
  check(identical(observed, expected), "counts: ", paste(names(observed), observed, sep = "=", collapse = ", "))
  # A gene never has both forms of a pattern, so each all-PACs count is the
  # sum of the other two.
  parts <- tapply(first$count[first$support != "all"], first$pattern[first$support != "all"], sum)
  check(identical(first$count[first$support == "all"], as.vector(parts)), "all is not the sum of the parts.")
  check(all(counts$count[counts$comparison == "B_vs_C"] == 0L), "the empty comparison.")
})

test_case("F-22", "a gene is colored by its first APA pattern, and marked when only flagged PACs support it", {
  potential <- "intronic_gain_potential_internal_priming"
  values <- c("intronic_gain;utr_shortening", "unclassified_change", "none", NA, "unknown", potential,
    paste0("utr_lengthening;", potential))
  colors <- figures$first_pattern(values)
  check(identical(as.character(colors),
    c("intronic_gain", "unclassified_change", "none", "none", "none", "intronic_gain", "utr_lengthening")),
    "first patterns: ", paste(colors, collapse = ", "))
  check(identical(levels(colors), names(figures$APA_CLASSES)), "levels.")
  flagged <- figures$first_pattern_flagged(values)
  check(identical(flagged, c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)), "flagged: ", paste(flagged, collapse = ", "))
})

test_case("F-20", "the coverage figure counts the PACs it leaves out", {
  comparisons <- data.frame(stem = "T_vs_C", title = "T vs C", stringsAsFactors = FALSE)
  data <- figures$coverage_data(comparisons, list(T_vs_C = pac_table(pac_row("p1", "g1"))))
  subtitle <- function(dropped) {
    figures$plot_effect_vs_coverage(data, comparisons, params, dropped)$labels$subtitle
  }
  check(grepl("\nNot shown: 3 PACs without reads in both groups$", subtitle(3)), subtitle(3))
  check(grepl("\nNot shown: 1 PAC without reads in both groups", subtitle(1), fixed = TRUE), subtitle(1))
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

# ---- Across comparisons -------------------------------------------------------

comparisons_of <- function(stems) {
  parts <- strsplit(stems, "_vs_", fixed = TRUE)
  data.frame(
    stem = stems, condition = vapply(parts, `[[`, character(1), 1L),
    control_condition = vapply(parts, `[[`, character(1), 2L),
    title = sub("_vs_", " vs ", stems, fixed = TRUE), stringsAsFactors = FALSE
  )
}

# A genes table from gene IDs, patterns, and gene FDRs.
pattern_genes <- function(ids, patterns, fdrs) {
  do.call(rbind, Map(function(id, pattern, fdr) {
    gene_row(id, apa_pattern = pattern, gene_fdr = fdr)
  }, ids, patterns, fdrs))
}

placeholder_message <- function(plot) {
  if (length(plot$layers) != 1L) return(NA_character_)
  label <- plot$layers[[1]]$aes_params$label
  if (is.null(label)) NA_character_ else label
}

test_case("F-23", "the pattern grid table orders genes by shared patterns, then best FDR", {
  comparisons <- comparisons_of(c("A_vs_C", "B_vs_C", "D_vs_A"))
  genes_tables <- list(
    A_vs_C = pattern_genes(
      c("g1", "g2", "g3", "g4", "g6"),
      c("utr_shortening", "none", "intronic_gain;utr_lengthening", "unclassified_change", "utr_lengthening"),
      c(0.01, 0.5, 0.2, 0.03, 0.01)
    ),
    B_vs_C = pattern_genes(
      c("g1", "g3", "g4", "g5", "g6"),
      c("utr_lengthening", "none", "unclassified_change", "alternative_last_exon", "utr_lengthening"),
      c(0.02, 0.001, 0.04, 0.3, 0.05)
    ),
    # g7's empty pattern counts as none, as the grid draws it.
    D_vs_A = pattern_genes(
      c("g1", "g4", "g5", "g7"), c("utr_shortening", "none", "intronic_loss", ""),
      c(0.5, 0.01, NA, NA)
    )
  )
  table <- figures$pattern_grid_table(comparisons, genes_tables)
  expected_columns <- c("gene_id", "gene_name", "patterned_comparisons", comparisons$stem)
  check(identical(names(table), expected_columns), "columns: ", paste(names(table), collapse = ", "))
  # g4 and g6 tie on two comparisons and a best FDR of 0.01, so the ID
  # decides; g2 and g7 have none, and g7's missing FDR puts it last.
  check(identical(table$gene_id, c("g1", "g4", "g6", "g5", "g3", "g2", "g7")), "order: ", paste(table$gene_id, collapse = ", "))
  check(identical(table$patterned_comparisons, c(3L, 2L, 2L, 2L, 1L, 0L, 0L)), "counts: ", paste(table$patterned_comparisons, collapse = ", "))
  check(identical(table$gene_name[[1]], "G1"), "names come from the genes tables.")
  row <- function(id) table[table$gene_id == id, , drop = FALSE]
  check(identical(unname(unlist(row("g5")[comparisons$stem])), c(NA, "alternative_last_exon", "intronic_loss")), "g5's cells.")
  check(identical(unname(unlist(row("g3")[comparisons$stem])), c("intronic_gain;utr_lengthening", "none", NA)), "g3's cells.")
  shown <- figures$grid_genes(table)
  check(identical(shown$gene_id, c("g1", "g4", "g6", "g5")), "grid genes: ", paste(shown$gene_id, collapse = ", "))
})

test_case("F-24", "the grid draws at most 50 genes, with untested and multi-pattern cells marked", {
  comparisons <- comparisons_of(c("A_vs_C", "B_vs_C"))
  ids <- sprintf("g%02d", 1:60)
  both <- pattern_genes(ids, rep("utr_shortening", 60), seq(0.001, 0.06, by = 0.001))
  first <- both
  first$apa_pattern[[1]] <- "intronic_gain;utr_shortening"
  second <- both[-2, , drop = FALSE]
  second <- rbind(second, pattern_genes("g61", "unclassified_change", 0.9))
  table <- figures$pattern_grid_table(comparisons, list(A_vs_C = first, B_vs_C = second))
  shown <- figures$grid_genes(table)
  check(nrow(shown) == 50L, "grid genes: ", nrow(shown))
  check(identical(shown$gene_id, ids[-2][1:50]), "the grid genes are not the table's first 50.")
  plot <- figures$plot_pattern_grid(table, comparisons)
  check(identical(plot$labels$subtitle, "59 genes with a pattern in two or more comparisons; the 50 shared by the most are shown"), "subtitle: ", plot$labels$subtitle)
  # The caption explains only the marks the grid draws: g01's +, and no *.
  caption <- strsplit(plot$labels$caption, "\n", fixed = TRUE)[[1]]
  check(identical(caption, c("+ Two or more patterns in the comparison; the cell shows the first.", "apa_patterns_by_comparison.tsv.gz has every tested gene's patterns.")), "caption: ", plot$labels$caption)
  check(identical(levels(plot$data$gene_id), rev(shown$gene_id)), "rows are not the grid genes, top first.")
  check(nrow(plot$data) == 100L, "cells: ", nrow(plot$data))
  marks <- plot$layers[[2]]$data
  check(identical(as.character(marks$gene_id), "g01") && identical(as.character(marks$fill), "intronic_gain"), "the + marks: ", paste(marks$gene_id, collapse = ", "))
  # A cell shows the first pattern, white for none, and grey where the gene
  # was not tested.
  three <- comparisons_of(c("A_vs_C", "B_vs_C", "D_vs_C"))
  small <- figures$pattern_grid_table(three, list(
    A_vs_C = first[1:3, ], B_vs_C = first[c(1, 3), ], D_vs_C = pattern_genes("g03", "none", 0.5)
  ))
  cells <- figures$plot_pattern_grid(small, three)$data
  fills <- stats::setNames(as.character(cells$fill), paste(cells$gene_id, cells$comparison))
  check(identical(unname(fills[c("g01 A_vs_C", "g03 B_vs_C", "g03 D_vs_C", "g01 D_vs_C")]),
    c("intronic_gain", "utr_shortening", "none", "not_tested")), "fills: ", paste(names(fills), fills, collapse = ", "))
  single <- figures$plot_pattern_grid(table, comparisons[1, , drop = FALSE])
  check(identical(placeholder_message(single), "Only one comparison"), "one comparison.")
  none <- figures$pattern_grid_table(comparisons, list(A_vs_C = first[1:2, ], B_vs_C = pattern_genes("g01", "none", 0.5)))
  check(identical(placeholder_message(figures$plot_pattern_grid(none, comparisons)), "No gene has a pattern in two or more comparisons"), "no shared genes.")
})

test_case("F-25", "pairs are related by a shared control or a chain, and capped at 15 panels", {
  pairs <- figures$concordance_pairs(comparisons_of(c("X_vs_Y", "Y_vs_Z", "W_vs_Z")))
  check(identical(pairs$pair, c("X_vs_Y|Y_vs_Z", "X_vs_Y|W_vs_Z", "Y_vs_Z|W_vs_Z")), "pairs: ", paste(pairs$pair, collapse = ", "))
  check(identical(pairs$relation, c("chained", "unrelated", "shared_control")), "relations: ", paste(pairs$relation, collapse = ", "))
  check(identical(pairs$title_a, c("X vs Y", "X vs Y", "Y vs Z")), "titles.")
  seven <- figures$concordance_pairs(comparisons_of(c("A_vs_C", "B_vs_C", "D_vs_A", "E_vs_F", "G_vs_H", "I_vs_J", "K_vs_L")))
  check(nrow(seven) == 21L, "pairs of seven: ", nrow(seven))
  shown <- figures$shown_pairs(seven)
  check(identical(shown$pair, c("A_vs_C|B_vs_C", "A_vs_C|D_vs_A")), "related pairs: ", paste(shown$pair, collapse = ", "))
  check(identical(shown$relation, c("shared_control", "chained")), "related relations.")
  six <- figures$concordance_pairs(comparisons_of(c("A_vs_C", "B_vs_C", "E_vs_F", "G_vs_H", "I_vs_J", "K_vs_L")))
  check(nrow(figures$shown_pairs(six)) == 15L, "15 pairs are all shown.")
  shared <- figures$concordance_pairs(comparisons_of(paste0(LETTERS[1:7], "_vs_Z")))
  capped <- figures$shown_pairs(shared)
  check(identical(capped$pair, shared$pair[1:15]), "the first 15 related pairs.")
  none <- figures$concordance_pairs(comparisons_of("A_vs_C"))
  check(nrow(none) == 0L && "relation" %in% names(none), "one comparison has no pairs.")
})

test_case("F-26", "concordance joins PACs tested in both, with calls and Pearson r", {
  comparisons <- comparisons_of(c("A_vs_C", "B_vs_C"))
  pairs <- figures$concordance_pairs(comparisons)
  a <- pac_table(
    pac_row("p1", "g1", event_type = "increased_usage", control = 0.4, treatment = 0.5),
    pac_row("p2", "g1", event_type = "lost_candidate", control = 0.6, treatment = 0.4),
    pac_row("p3", "g2", event_type = "lost", control = 0.3, treatment = 0.6),
    pac_row("p4", "g2", control = NA, treatment = 0.5),
    pac_row("p5", "g3")
  )
  b <- pac_table(
    pac_row("p3", "g2", control = 0.25, treatment = 0.5),
    pac_row("p1", "g1", event_type = "decreased_usage", control = 0.3, treatment = 0.5),
    pac_row("p2", "g1", event_type = "gained_candidate", control = 0.5, treatment = 0.4),
    pac_row("p4", "g2", control = 0.4, treatment = 0.5),
    pac_row("p6", "g4")
  )
  data <- figures$concordance_data(pairs, list(A_vs_C = a, B_vs_C = b))
  # Candidate calls are not confirmed calls; p4 has no change in A.
  check(identical(data$pac_id, c("p2", "p3", "p1")), "rows: ", paste(data$pac_id, collapse = ", "))
  check(identical(as.character(data$call), c("none", "one", "both")), "calls: ", paste(data$call, collapse = ", "))
  check(isTRUE(all.equal(data$delta_a, c(-0.2, 0.3, 0.1))) && isTRUE(all.equal(data$delta_b, c(-0.1, 0.25, 0.2))), "changes.")
  check(identical(levels(data$pair), pairs$pair), "pair levels.")
  summary <- figures$concordance_summary(pairs, data)
  check(identical(names(summary), c("comparison_a", "comparison_b", "relation", "shared_pacs", "called_in_both", "pearson_r")), "summary columns.")
  check(identical(summary$shared_pacs, 3L) && identical(summary$called_in_both, 1L), "summary counts.")
  expected <- stats::cor(c(-0.2, 0.3, 0.1), c(-0.1, 0.25, 0.2))
  check(isTRUE(all.equal(summary$pearson_r, expected)), "r: ", summary$pearson_r, " not ", expected)
  two <- figures$concordance_summary(pairs, data[data$pac_id != "p1", , drop = FALSE])
  check(is.na(two$pearson_r) && two$shared_pacs == 2L, "two PACs have no r.")
  flat <- data
  flat$delta_b <- 0.1
  check(is.na(figures$concordance_summary(pairs, flat)$pearson_r), "no spread has no r.")
  check(nrow(figures$concordance_summary(pairs[0, ], data[0, ])) == 0L, "no pairs.")
})

test_case("F-27", "the correlation matrix is symmetric with 1 on the diagonal", {
  comparisons <- comparisons_of(c("A_vs_C", "B_vs_C", "D_vs_A"))
  summary <- data.frame(
    comparison_a = c("A_vs_C", "A_vs_C", "B_vs_C"), comparison_b = c("B_vs_C", "D_vs_A", "D_vs_A"),
    relation = c("shared_control", "chained", "unrelated"), shared_pacs = c(10L, 10L, 2L),
    called_in_both = c(1L, 2L, 0L), pearson_r = c(0.4, -0.7, NA), stringsAsFactors = FALSE
  )
  grid <- figures$concordance_matrix(summary, comparisons)
  check(nrow(grid) == 9L, "cells: ", nrow(grid))
  value <- function(a, b) grid$r[grid$a == a & grid$b == b]
  check(all(vapply(comparisons$stem, function(stem) value(stem, stem), numeric(1)) == 1), "diagonal.")
  check(value("A_vs_C", "B_vs_C") == 0.4 && value("B_vs_C", "A_vs_C") == 0.4, "A and B.")
  check(value("D_vs_A", "A_vs_C") == -0.7 && value("A_vs_C", "D_vs_A") == -0.7, "A and D.")
  check(is.na(value("B_vs_C", "D_vs_A")) && is.na(value("D_vs_A", "B_vs_C")), "missing r.")
  check(identical(placeholder_message(figures$plot_concordance_matrix(summary[0, ], comparisons[1, , drop = FALSE])), "Only one comparison"), "one comparison.")
})

test_case("F-28", "the PAU PCA uses genes covered in every sample and fixes each sign", {
  samples <- data.frame(
    sample_id = c("S3", "S1", "S4", "S2"), condition = c("T", "C", "T", "C"),
    control_condition = c("C", "C", "C", "C"), stringsAsFactors = FALSE
  )
  sample_ids <- c("S1", "S2", "S3", "S4")
  rows <- function(gene, pac, values, totals) {
    data.frame(gene_id = gene, pac_id = pac, sample_id = sample_ids, gene_total = totals,
      pau = values, stringsAsFactors = FALSE)
  }
  pau <- rbind(
    rows("g1", "g1.a", c(0.2, 0.3, 0.7, 0.8), c(40, 30, 50, 20)),
    rows("g1", "g1.b", c(0.5, 0.4, 0.1, 0.2), c(40, 30, 50, 20)),
    rows("g1", "g1.c", c(0.3, 0.3, 0.2, 0.0), c(40, 30, 50, 20)),
    # g2 has 19 reads in S4, and g3 none in S2.
    rows("g2", "g2.a", c(0.9, 0.1, 0.5, 0.5), c(40, 30, 50, 19)),
    rows("g3", "g3.a", c(0.5, NA, 0.4, 0.6), c(40, 0, 50, 60)),
    # g4 lacks a row for one PAC in S3, which drops that PAC alone.
    rows("g4", "g4.a", c(0.6, 0.65, 0.2, 0.3), c(25, 25, 25, 25)),
    rows("g4", "g4.b", c(0.4, 0.35, 0.8, 0.7), c(25, 25, 25, 25))[-3, ]
  )
  pca <- figures$pau_pca(pau, samples, 20)
  check(identical(pca$sample_id, sample_ids) && identical(pca$condition, c("C", "C", "T", "T")), "samples.")
  check(identical(attr(pca, "genes"), 2L) && identical(attr(pca, "pacs"), 4L), "genes ", attr(pca, "genes"), ", PACs ", attr(pca, "pacs"))
  values <- rbind(
    c(0.2, 0.3, 0.7, 0.8), c(0.5, 0.4, 0.1, 0.2), c(0.3, 0.3, 0.2, 0.0), c(0.6, 0.65, 0.2, 0.3)
  )
  centered <- t(values - rowMeans(values))
  direct <- stats::prcomp(centered, center = FALSE)
  scores <- as.matrix(pca[, c("PC1", "PC2")])
  check(isTRUE(all.equal(abs(unname(scores)), abs(unname(direct$x[, 1:2])))), "scores differ from prcomp.")
  variance <- direct$sdev^2 / sum(direct$sdev^2)
  check(isTRUE(all.equal(pca$pc1_variance_fraction[[1]], variance[[1]])) && isTRUE(all.equal(pca$pc2_variance_fraction[[1]], variance[[2]])), "variance fractions.")
  # Each component's loadings, recovered from its scores, have their largest
  # entry positive.
  for (k in 1:2) {
    loading <- crossprod(centered, scores[, k])
    check(loading[[which.max(abs(loading))]] > 0, "PC", k, "'s largest loading is negative.")
  }
  # A two-PAC gene's loadings tie, and the first PAC's is made positive, so
  # the samples where g5.a is highest score positive, in either arrangement.
  tied <- rbind(
    rows("g5", "g5.a", c(0.9, 0.8, 0.3, 0.2), c(40, 30, 50, 20)),
    rows("g5", "g5.b", c(0.1, 0.2, 0.7, 0.8), c(40, 30, 50, 20))
  )
  for (reversed in c(FALSE, TRUE)) {
    arranged <- tied
    if (reversed) arranged$pau <- 1 - tied$pau
    tie <- figures$pau_pca(arranged, samples, 20)
    high <- if (reversed) c("S3", "S4") else c("S1", "S2")
    check(identical(tie$PC1 > 0, tie$sample_id %in% high), "tied loadings, reversed ", reversed, ": ", paste(tie$PC1, collapse = ", "))
  }
  check(is.null(figures$pau_pca(pau, samples[1, , drop = FALSE], 20)), "one sample.")
  check(is.null(figures$pau_pca(pau, samples, 1000)), "no covered gene.")
  check(identical(placeholder_message(figures$plot_pau_pca(NULL, params)), "No PAC has observed PAU in every sample"), "PCA placeholder.")
  plot <- figures$plot_pau_pca(pca, params)
  check(identical(plot$labels$subtitle, "Observed PAU at 4 PACs in 2 genes with at least 20 reads in every sample"), "subtitle: ", plot$labels$subtitle)
})

test_case("F-29", "the concordance figure says which pairs it leaves out", {
  seven <- comparisons_of(c("A_vs_C", "B_vs_C", "D_vs_A", "E_vs_F", "G_vs_H", "I_vs_J", "K_vs_L"))
  pairs <- figures$concordance_pairs(seven)
  data <- figures$concordance_data(pairs, stats::setNames(rep(list(pac_table(pac_row("p1", "g1"))), 7), seven$stem))
  summary <- figures$concordance_summary(pairs, data)
  plot <- figures$plot_concordance(data, figures$shown_pairs(pairs), summary, nrow(pairs), seven)
  check(grepl("\n2 of 21 pairs shown", plot$labels$subtitle, fixed = TRUE), "subtitle: ", plot$labels$subtitle)
  check(identical(levels(plot$data$pair), c("A_vs_C|B_vs_C", "A_vs_C|D_vs_A")), "panels.")
  # Both drawn pairs start with A_vs_C, so the matrix has one column.
  check(identical(levels(plot$data$column), "A_vs_C") && identical(levels(plot$data$row), c("B_vs_C", "D_vs_A")), "matrix: ", paste(levels(plot$data$column), collapse = ", "), " by ", paste(levels(plot$data$row), collapse = ", "))
  unrelated <- figures$concordance_pairs(comparisons_of(c("A_vs_B", "C_vs_D", "E_vs_F", "G_vs_H", "I_vs_J", "K_vs_L", "M_vs_N")))
  check(identical(placeholder_message(figures$plot_concordance(data[0, ], figures$shown_pairs(unrelated), summary, 21L, seven)), "No pair shares a control or a condition"), "no related pairs.")
  check(identical(placeholder_message(figures$plot_concordance(data[0, ], pairs[0, ], summary[0, ], 0L, seven[1, , drop = FALSE])), "Only one comparison"), "one comparison.")
})

test_case("F-30", "flagged patterns are marked on the grid and the distal plot", {
  potential <- "intronic_gain_potential_internal_priming"
  # The grid marks a flagged first pattern with *, and two patterns with +.
  comparisons <- comparisons_of(c("A_vs_C", "B_vs_C"))
  table <- figures$pattern_grid_table(comparisons, list(
    A_vs_C = pattern_genes(c("g1", "g2"), c(potential, "utr_shortening"), c(0.01, 0.02)),
    B_vs_C = pattern_genes(
      c("g1", "g2"), c(paste0(potential, ";utr_shortening"), "unclassified_change"), c(0.01, 0.02)
    )
  ))
  plot <- figures$plot_pattern_grid(table, comparisons)
  marks <- plot$layers[[2]]$data
  observed <- stats::setNames(marks$mark, paste(marks$gene_id, marks$comparison))
  check(identical(observed[order(names(observed))], c("g1 A_vs_C" = "*", "g1 B_vs_C" = "+*")),
    "grid marks: ", paste(names(observed), observed, collapse = ", "))
  # The distal plot draws a flagged pattern's gene as a diamond.
  pacs <- pac_table(pac_row("p1", "g1", event_type = "increased_usage", control = 0.2, treatment = 0.5, pac_fdr = 0.001))
  distal <- figures$distal_usage_table(pacs, data.frame(pac_id = "p1", coordinate = 100),
    genes_of(pacs, c(g1 = potential)))
  shapes <- figures$plot_distal_usage(distal, comparison, params)
  check(identical(as.character(shapes$data$call), "Confirmed *"), "distal mark: ", as.character(shapes$data$call))
  check(grepl("\n* The gene's APA pattern comes only from PACs flagged for possible internal priming.", shapes$labels$caption, fixed = TRUE), "distal caption: ", shapes$labels$caption)
  check(grepl("* The pattern comes only from PACs flagged for possible internal priming.", plot$labels$caption, fixed = TRUE), "grid caption: ", plot$labels$caption)
})

test_case("F-31", "counts carry thousands separators", {
  text <- figures$count_text(c(0, 999, 1000, 12345, 60610, 3e9, NA))
  check(identical(text, c("0", "999", "1,000", "12,345", "60,610", "3,000,000,000", NA)), "counts: ", paste(text, collapse = " "))
  check(identical(figures$plural(1, "PAC"), "1 PAC") && identical(figures$plural(12000, "gene"), "12,000 genes"), "plurals.")
})

test_case("F-32", "the density scale breaks at powers of ten, or finer steps in a narrow range", {
  breaks <- function(low, high) figures$density_breaks(c(low, high))
  check(identical(breaks(1, 2500), c(1, 10, 100, 1000)), "1 to 2,500: ", paste(breaks(1, 2500), collapse = ", "))
  check(identical(breaks(1, 60), c(1, 3, 10, 30)), "1 to 60: ", paste(breaks(1, 60), collapse = ", "))
  check(identical(breaks(1, 5), c(1, 2, 5)), "1 to 5: ", paste(breaks(1, 5), collapse = ", "))
  check(identical(breaks(2, 2), 2), "one count: ", paste(breaks(2, 2), collapse = ", "))
})

test_case("F-33", "the volcano subtitle and caption say what is left out and what marks mean", {
  pacs <- pac_table(
    pac_row("up", "g1", event_type = "increased_usage", control = 0.2, treatment = 0.5, pvalue = 1e-8),
    pac_row("down", "g2", event_type = "decreased_usage", control = 0.5, treatment = 0.2, pvalue = 0),
    pac_row("flat", "g2", pvalue = 0.4),
    pac_row("no_p", "g3", pvalue = NA_real_),
    pac_row("no_p_candidate", "g3", event_type = "gained_candidate", control = 0, treatment = 0.3, pvalue = NA_real_)
  )
  genes <- rbind(
    shift_row("g1", "distal", "elsewhere", "up", from_event = "none"),
    shift_row("g2", "proximal", "down", "flat", to_event = "none"),
    gene_row("g3")
  )
  plot <- figures$plot_volcano(pacs, genes, comparison, params, exploratory = FALSE)
  check(identical(plot$labels$subtitle, paste0(
    "5 PACs tested in 3 genes; confirmed calls: 1 up and 1 down, in 2 genes\n",
    "Not shown: 2 PACs without a p-value, including 1 candidate call"
  )), "subtitle: ", plot$labels$subtitle)
  check(grepl("^Dashed lines: a change of 0.1 either way", plot$labels$caption), "caption: ", plot$labels$caption)
  check(grepl("\nLabels: one per gene, at the PAC its usage moved to", plot$labels$caption, fixed = TRUE), "no label note: ", plot$labels$caption)
  check(grepl("\nEach shift shows on both sides", plot$labels$caption, fixed = TRUE), "no both-sides note: ", plot$labels$caption)
  check(grepl("\nTriangles: p-values of 0", plot$labels$caption, fixed = TRUE), "no triangle note: ", plot$labels$caption)
  # Without missing candidates or p-values of 0, neither note appears.
  plain <- figures$plot_volcano(pacs[pacs$pac_id != "no_p_candidate", , drop = FALSE], genes, comparison, params, FALSE)
  check(grepl("without a p-value$", plain$labels$subtitle), "subtitle: ", plain$labels$subtitle)
  finite <- pacs[!pacs$pac_id %in% c("down", "no_p_candidate"), , drop = FALSE]
  check(!grepl("Triangles", figures$plot_volcano(finite, genes, comparison, params, FALSE)$labels$caption, fixed = TRUE), "triangle note without p-values of 0.")
})

test_case("F-34", "PCA labels drop the condition's name only when every sample ID starts with it", {
  labels <- figures$sample_labels(c("AS_NT_DMSO_rep1", "C-2", "C.3"), c("AS_NT_DMSO", "C", "C"))
  check(identical(labels, c("rep1", "2", "3")), "every ID shortened: ", paste(labels, collapse = ", "))
  # One ID that does not start with its condition's name and a separator
  # keeps every label whole: S1, Cx_1 (no separator), C_ (nothing after it).
  for (odd in c("S1", "Cx_1", "C_", "C")) {
    ids <- c("C_1", "C_2", odd)
    labels <- figures$sample_labels(ids, c("C", "C", "C"))
    check(identical(labels, ids), odd, ": ", paste(labels, collapse = ", "))
  }
  check(identical(figures$sample_labels(character(), character()), character()), "no samples.")
})

test_case("F-35", "the concordance figure places each pair by its two comparisons, below the diagonal", {
  comparisons <- comparisons_of(c("A_vs_C", "B_vs_C", "D_vs_A"))
  pairs <- figures$concordance_pairs(comparisons)
  tables <- stats::setNames(rep(list(pac_table(
    pac_row("p1", "g1", control = 0.2, treatment = 0.4), pac_row("p2", "g1", control = 0.5, treatment = 0.3),
    pac_row("p3", "g2", control = 0.1, treatment = 0.6)
  )), 3), comparisons$stem)
  data <- figures$concordance_data(pairs, tables)
  summary <- figures$concordance_summary(pairs, data)
  plot <- figures$plot_concordance(data, pairs, summary, nrow(pairs), comparisons)
  # The first layer is each drawn panel's background: one per pair, at the
  # pair's first comparison's column and second's row; B_vs_C's own cell
  # (column B_vs_C, row B_vs_C) stays empty.
  cells <- plot$layers[[1]]$data
  placed <- paste(cells$column, cells$row)
  check(identical(placed, c("A_vs_C B_vs_C", "A_vs_C D_vs_A", "B_vs_C D_vs_A")), "cells: ", paste(placed, collapse = "; "))
  check(identical(levels(cells$column), c("A_vs_C", "B_vs_C")) && identical(levels(cells$row), c("B_vs_C", "D_vs_A")), "matrix levels.")
  check(identical(paste(plot$data$column, plot$data$row)[plot$data$pair == "A_vs_C|D_vs_A"], rep("A_vs_C D_vs_A", 3)), "points are not in their pair's panel.")
  size <- figures$concordance_size(pairs, comparisons)
  check(identical(size, c(width = 5.5, height = 5.25)), "size: ", paste(size, collapse = " x "))
  # The correlation matrix draws each pair once, in the same cells.
  matrix <- figures$plot_concordance_matrix(summary, comparisons)$data
  check(identical(paste(matrix$a, matrix$b), c("A_vs_C B_vs_C", "A_vs_C D_vs_A", "B_vs_C D_vs_A")), "matrix cells: ", paste(matrix$a, matrix$b, collapse = "; "))
})

test_case("F-38", "label placement has a fixed seed and no time limit", {
  # F-10 compares bytes, but its few labels settle long before ggrepel's
  # default half-second limit, so it would not notice the limit coming back.
  layer <- figures$repel_text(data.frame(x = 1, y = 1, label = "a"),
    aes(x = x, y = y, label = label), point_size = 1, nudge_y = 0)
  # ggrepel's default force and 10,000 iterations left crowded labels
  # overlapping at a real run's scale.
  settings <- layer$geom_params[c("max.time", "max.iter", "max.overlaps", "seed")]
  check(identical(settings, list(max.time = Inf, max.iter = 100000L, max.overlaps = Inf, seed = 1L)),
    "settings: ", paste(names(settings), unlist(settings), collapse = ", "))
  check(identical(layer$geom_params$force, 32), "force: ", layer$geom_params$force)
  check(inherits(layer$position, "PositionNudgeRepel"), "labels are not nudged off their points.")
})

test_case("F-37", "each gene is labelled once, at the PAC its usage moved to", {
  pacs <- pac_table(
    pac_row("switch_up", "switch", event_type = "increased_usage", control = 0.3, treatment = 0.7, pvalue = 1e-30),
    pac_row("switch_down", "switch", event_type = "decreased_usage", control = 0.7, treatment = 0.3, pvalue = 1e-30),
    pac_row("switch_down2", "switch", event_type = "lost", control = 0.2, treatment = 0, pvalue = 1e-10),
    # The label goes on the shift's PAC, not the gene's most significant one.
    pac_row("big_a", "big", event_type = "increased_usage", control = 0.2, treatment = 0.5, pvalue = 1e-5),
    pac_row("big_b", "big", event_type = "increased_usage", control = 0.2, treatment = 0.35, pvalue = 1e-30),
    pac_row("big_c", "big", event_type = "decreased_usage", control = 0.6, treatment = 0.15, pvalue = 1e-30),
    # Without a confirmed gain, the label goes where usage moved from.
    pac_row("loss_down", "loss", event_type = "lost", control = 0.3, treatment = 0, pvalue = 1e-8),
    pac_row("loss_up", "loss", control = 0.7, treatment = 1, pvalue = 0.2),
    pac_row("other_up", "other", event_type = "gained", control = 0, treatment = 0.4, pvalue = 1e-5)
  )
  genes <- rbind(
    shift_row("switch", "distal", "switch_down", "switch_up"),
    shift_row("big", "proximal", "big_c", "big_a"),
    shift_row("loss", "distal", "loss_down", "loss_up", from_event = "lost", to_event = "none"),
    shift_row("other", "distal", "elsewhere", "other_up", from_event = "none", to_event = "gained")
  )
  labels <- figures$volcano_data(pacs, genes)$labels
  observed <- paste(labels$pac_id, labels$call)
  check(identical(sort(observed), c("big_a up", "loss_down down", "other_up up", "switch_up up")), "labels: ", paste(observed, collapse = ", "))
})

test_case("F-36", "among volcano PACs with equal p-values, the larger changes are labelled", {
  zeros <- lapply(1:25, function(index) {
    pac_row(sprintf("z%02d", index), sprintf("zero%02d", index), event_type = "increased_usage",
      control = 0.1, treatment = 0.1 + (index + 10) / 100, pvalue = 0)
  })
  genes <- do.call(rbind, lapply(1:25, function(index) {
    shift_row(sprintf("zero%02d", index), "distal", "elsewhere", sprintf("z%02d", index))
  }))
  labels <- figures$volcano_data(do.call(pac_table, zeros), genes)$labels
  check(identical(labels$pac_id, sprintf("z%02d", 25:6)), "labels: ", paste(labels$pac_id, collapse = ", "))
})

# ---- Layout -------------------------------------------------------------------

# The thirteen conditions of the lead's design, four samples each, placed
# around a circle.
lead_pca <- function() {
  conditions <- c(
    "AS_CDK12_13_CR8", "AS_CDK12_13_DMSO", "AS_CDK12_13_PHA", "AS_CDK12_CR8", "AS_CDK12_DMSO",
    "AS_CDK12_PHA", "AS_CDK13_CR8", "AS_CDK13_DMSO", "AS_CDK13_PHA", "AS_NT_CR8", "AS_NT_DMSO",
    "AS_NT_PHA", "WT_NT_DMSO"
  )
  condition <- rep(conditions, each = 4)
  angle <- 2 * pi * (rep(seq_along(conditions), each = 4) + rep(0:3, 13) / 40) / 13
  pca <- data.frame(
    sample_id = paste0(condition, "_rep", 1:4), condition = condition,
    PC1 = 10 * cos(angle), PC2 = 6 * sin(angle), pc1_variance_fraction = 0.54,
    pc2_variance_fraction = 0.13, stringsAsFactors = FALSE
  )
  attr(pca, "pacs") <- 171270L
  attr(pca, "genes") <- 9470L
  pca
}

# Whether a plot draws anything outside its width and height. It is drawn in
# the middle of a larger page, at its own size, and the PNG is compared with
# one whose margins are painted over.
draws_outside <- function(plot, width, height) {
  render <- function(cover) {
    path <- tempfile(fileext = ".png")
    grDevices::png(path, width = (width + 2) * 200, height = (height + 2) * 200, res = 200)
    grid::grid.newpage()
    print(plot, newpage = FALSE,
      vp = grid::viewport(width = grid::unit(width, "in"), height = grid::unit(height, "in")))
    if (cover) {
      white <- grid::gpar(fill = "white", col = NA)
      inch <- grid::unit(1, "in")
      grid::grid.rect(x = 0, width = inch, just = "left", gp = white)
      grid::grid.rect(x = 1, width = inch, just = "right", gp = white)
      grid::grid.rect(y = 0, height = inch, just = "bottom", gp = white)
      grid::grid.rect(y = 1, height = inch, just = "top", gp = white)
    }
    grDevices::dev.off()
    readBin(path, "raw", file.size(path))
  }
  !identical(render(FALSE), render(TRUE))
}

test_case("F-40", "the PCA's legend sits beside the panel, in columns, and the figure widens for it", {
  pca <- lead_pca()
  plot <- figures$plot_pau_pca(pca, params)
  check(identical(plot$theme$legend.position, "right"), "legend position: ", plot$theme$legend.position)
  check(identical(plot$scales$get_scales("colour")$guide$params$ncol, 1L), "one column for 13 conditions.")
  size <- figures$pca_size(pca)
  check(identical(size, c(width = 7.125, height = 5.5)), "size: ", paste(size, collapse = " x "))
  # A bottom legend of these 13 names was wider than the figure.
  check(!draws_outside(plot, size[["width"]], size[["height"]]), "the PCA draws outside its page.")
  # 41 conditions take three columns, and the figure widens for each.
  many <- data.frame(sample_id = sprintf("S%02d", 1:41), condition = sprintf("C%02d", 1:41),
    PC1 = 1:41, PC2 = (1:41)^2, pc1_variance_fraction = 0.5, pc2_variance_fraction = 0.2)
  attr(many, "pacs") <- 10L
  attr(many, "genes") <- 5L
  check(identical(figures$pca_legend_columns(many$condition), 3L), "columns for 41 conditions.")
  check(identical(figures$pca_size(many), c(width = 7.375, height = 5.5)), "size for 41 conditions: ", paste(figures$pca_size(many), collapse = " x "))
  check(identical(figures$pca_size(NULL), c(width = 5.5, height = 5)), "placeholder size.")
})

test_case("F-41", "long comparison titles and condition names wrap", {
  check(identical(figures$wrap_name(c("AS_NT_DMSO", "MW_CDK12_13_PHA_LONG_TREATMENT_NAMEX")),
    c("AS_NT_DMSO", "MW_CDK12_13_PHA_LONG_\nTREATMENT_NAMEX")), "wrapped names.")
  # A name without a break point stays whole.
  check(identical(figures$wrap_name(strrep("A", 40)), strrep("A", 40)), "an unbreakable name.")
  comparisons <- comparisons_of(c("T_vs_C", "AS_CDK12_13_CR8_vs_AS_CDK12_13_DMSO"))
  titles <- figures$strip_titles(comparisons)
  check(identical(unname(titles), c("T vs C", "AS_CDK12_13_CR8\nvs AS_CDK12_13_DMSO")), "strip titles: ", paste(titles, collapse = " | "))
  # The count axes use the same titles, so a long one narrows the panels less.
  counts <- figures$event_count_table(comparisons, list(T_vs_C = pac_table(), AS_CDK12_13_CR8_vs_AS_CDK12_13_DMSO = pac_table()),
    list(T_vs_C = gene_row("x")[0, ], AS_CDK12_13_CR8_vs_AS_CDK12_13_DMSO = gene_row("x")[0, ]))
  axis <- ggplot2::ggplot_build(figures$plot_event_bars(counts, comparisons))$layout$panel_params[[1]]$y$get_labels()
  check("AS_CDK12_13_CR8\nvs AS_CDK12_13_DMSO" %in% axis, "event axis labels: ", paste(axis, collapse = " | "))
  check(all(lengths(strsplit(figures$PATTERN_SUPPORT_LABELS, "\n", fixed = TRUE)) == 2L), "pattern facet titles are not all two lines.")
})

test_case("F-42", "the last axis label at a panel's right edge stays on the page", {
  # 4,364 genes each way put the axis's last break, 6,000, on the panel's
  # right edge, so half the five-character label needs the right margin.
  regions <- names(figures$GENE_REGIONS)
  pairs <- paste(rep(regions, each = 4), rep(regions, 4), sep = "|")
  counts <- data.frame(
    pair = factor(rep(pairs, 2), levels = pairs),
    direction = factor(rep(c("proximal", "distal"), each = 16), levels = c("proximal", "distal")),
    count = c(4364L, 580L, 1870L, 260L, 290L, 50L, 160L, 10L, 910L, 90L, 380L, 40L, 160L, 20L,
      90L, 0L, 4364L, 360L, 890L, 90L, 630L, 110L, 70L, 20L, 2110L, 220L, 500L, 130L, 320L, 30L,
      30L, 10L)
  )
  plot <- figures$plot_shifts_by_gene_region(counts, comparison)
  check(!draws_outside(plot, 6.5, figures$shift_figure_height(pairs)), "the shifts figure draws outside its page.")
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
    "  figures$repel_text(points, ggplot2::aes(label = label), point_size = 1.5, nudge_y = 10) +",
    "  ggplot2::labs(title = 'Rendering test') + figures$figure_theme()",
    "figures$save_figure(plot, arguments[[2]], 4, 3)",
    "stack <- figures$figure_stack(list(plot, plot), heights = c(2, 1))",
    "figures$save_figure(stack, paste0(arguments[[2]], '_stack'), 4, 3)"
  ), driver)
  stems <- file.path(work, c("first", "second"), "render")
  for (stem in stems) {
    dir.create(dirname(stem))
    run_rscript(c(driver, script, stem))
  }
  for (suffix in c(".pdf", ".png", "_stack.pdf", "_stack.png")) {
    files <- paste0(stems, suffix)
    bytes <- lapply(files, function(path) readBin(path, "raw", file.size(path)))
    check(identical(bytes[[1]], bytes[[2]]), suffix, " differs between processes.")
  }
  stack_pdf <- readBin(paste0(stems[[1]], "_stack.pdf"), "raw", file.size(paste0(stems[[1]], "_stack.pdf")))
  check(!length(grepRaw("/CreationDate", stack_pdf, fixed = TRUE)), "the stacked PDF records a date.")
  check(identical(png_size(paste0(stems[[1]], "_stack.png")), c(800, 600)), "stacked PNG size.")
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
    pac_row("g2_p2", "g2", "-", "none", gene_region = "intron", control = 0.5, treatment = 0.45, gene_fdr = 0.8, pvalue = 0.7, condition = "T1")
  )
  genes <- rbind(
    gene_row("g1", dominant_switch = TRUE, gene_fdr = 0.001, apa_pattern = "utr_lengthening",
      shift = c("distal", "g1_p1", "last_exon", "decreased_usage", "g1_p2", "last_exon",
        "increased_usage")),
    gene_row("g2", gene_fdr = 0.8)
  )
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
  # Observed PAU as MERGE_COUNTS writes it: g1's usage moves in T1 and T2.
  shift <- c(0, 0, 0.3, 0.35, 0.1, 0.15)
  pau <- rbind(
    data.frame(pac_id = "g1_p1", sample_id = samples$sample_id, pau = 0.6 - shift),
    data.frame(pac_id = "g1_p2", sample_id = samples$sample_id, pau = 0.4 + shift),
    data.frame(pac_id = "g2_p1", sample_id = samples$sample_id, pau = c(0.5, 0.52, 0.49, 0.5, 0.51, 0.48)),
    data.frame(pac_id = "g2_p2", sample_id = samples$sample_id, pau = c(0.5, 0.48, 0.51, 0.5, 0.49, 0.52))
  )
  pau$gene_id <- sub("_p[0-9]$", "", pau$pac_id)
  pau$gene_name <- toupper(pau$gene_id)
  pau$gene_total <- 50
  pau$count <- pau$pau * pau$gene_total
  pau <- pau[, c("pac_id", "gene_id", "gene_name", "sample_id", "count", "gene_total", "pau")]
  write_table(if (filled) pau else pau[0, , drop = FALSE], file.path(directory, "observed_pau.tsv.gz"))
  yaml::write_yaml(c(params, list(outdir = "results")), file.path(directory, "params.yaml"))
  c(
    "--mode", "figures", "--samples", file.path(directory, "samples.tsv"),
    "--statistics-dir", statistics, "--atlas", file.path(directory, "atlas.tsv.gz"),
    "--pau", file.path(directory, "observed_pau.tsv.gz"),
    "--params", file.path(directory, "params.yaml"), "--output-dir", file.path(directory, "figures")
  )
}

expected_files <- sort(c(
  as.vector(outer(c("T1_vs_C", "T2_vs_C"), c(
    ".volcano.pdf", ".volcano.png", ".distal_usage.pdf", ".distal_usage.png",
    ".distal_usage.tsv.gz", ".shifts_by_gene_region.pdf", ".shifts_by_gene_region.png"
  ), paste0)),
  "event_counts.pdf", "event_counts.png", "effect_vs_coverage.pdf", "effect_vs_coverage.png",
  "apa_pattern_grid.pdf", "apa_pattern_grid.png", "apa_patterns_by_comparison.tsv.gz",
  "concordance.pdf", "concordance.png", "concordance.tsv.gz", "concordance_matrix.pdf",
  "concordance_matrix.png", "pau_pca.pdf", "pau_pca.png", "pau_pca.tsv"
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
  check(identical(distal$direction, c("distal_up", "none")), "directions: ", paste(distal$direction, collapse = ", "))
  check(identical(distal$apa_pattern, c("utr_lengthening", "none")), "patterns: ", paste(distal$apa_pattern, collapse = ", "))
  check(identical(distal$tested_pacs, c("2", "2")), "tested_pacs.")
  check(identical(distal$distal_pac_fdr, c("0.001", "")), "pac_fdr: ", paste(distal$distal_pac_fdr, collapse = ", "))
  empty <- utils::read.delim(file.path(directory, "figures", "T2_vs_C.distal_usage.tsv.gz"))
  check(nrow(empty) == 0L && identical(names(empty), figures$DISTAL_COLUMNS), "the empty comparison's table.")
  for (name in c("event_counts.png", "effect_vs_coverage.png", "T1_vs_C.volcano.png")) {
    check(file.size(file.path(directory, "figures", name)) > 1000, name, " is nearly empty.")
  }
  check(identical(png_size(file.path(directory, "figures", "event_counts.png")), c(1500, 1075)), "event_counts size.")
  # One region pair, last exon to last exon, in both comparisons' figures.
  for (stem in c("T1_vs_C", "T2_vs_C")) {
    size <- png_size(file.path(directory, "figures", paste0(stem, ".shifts_by_gene_region.png")))
    check(identical(size, c(1300, 500)), stem, " shifts size: ", paste(size, collapse = "x"))
  }
  check(identical(png_size(file.path(directory, "figures", "effect_vs_coverage.png")), c(1500, 800)), "coverage size.")
  sizes <- list(apa_pattern_grid = c(900, 525), concordance = c(1100, 750),
    concordance_matrix = c(950, 950), pau_pca = c(1225, 1100))
  for (name in names(sizes)) {
    size <- png_size(file.path(directory, "figures", paste0(name, ".png")))
    check(identical(size, sizes[[name]]), name, " size: ", paste(size, collapse = "x"))
  }
  read_output <- function(name) {
    utils::read.delim(file.path(directory, "figures", name), colClasses = "character",
      na.strings = character(), check.names = FALSE)
  }
  grid <- read_output("apa_patterns_by_comparison.tsv.gz")
  check(identical(names(grid), c("gene_id", "gene_name", "patterned_comparisons", "T1_vs_C", "T2_vs_C")), "grid columns.")
  check(identical(unname(as.matrix(grid)), rbind(c("g1", "G1", "1", "utr_lengthening", ""), c("g2", "G2", "0", "none", ""))), "grid rows.")
  concordance <- read_output("concordance.tsv.gz")
  check(identical(unname(unlist(concordance)), c("T1_vs_C", "T2_vs_C", "shared_control", "0", "0", "")), "concordance: ", paste(unlist(concordance), collapse = ", "))
  pca <- read_output("pau_pca.tsv")
  check(identical(names(pca), c("sample_id", "condition", "PC1", "PC2", "pc1_variance_fraction", "pc2_variance_fraction")), "PCA columns.")
  check(identical(pca$sample_id, paste0("s", 1:6)) && identical(pca$condition, c("C", "C", "T1", "T1", "T2", "T2")), "PCA samples.")
  # PC1 follows g1's shift. g1's two PACs tie for the largest loading, so
  # g1_p1's is positive, and the controls, where g1_p1 is highest, score
  # highest.
  check(identical(order(as.numeric(pca$PC1)), c(4L, 3L, 6L, 5L, 1L, 2L)), "PC1 order: ", paste(pca$PC1, collapse = ", "))
})

test_case("F-13", "every comparison empty still writes the same files", {
  directory <- file.path(work, "empty")
  command <- write_inputs(directory, filled = FALSE)
  figures$run_figures_mode(figures$parse_args(command))
  files <- sort(list.files(file.path(directory, "figures")), method = "radix")
  check(identical(files, expected_files), "files: ", paste(files, collapse = ", "))
  pca <- readLines(file.path(directory, "figures", "pau_pca.tsv"))
  check(identical(pca, "sample_id\tcondition\tPC1\tPC2\tpc1_variance_fraction\tpc2_variance_fraction"), "empty PCA table: ", paste(pca, collapse = " | "))
})

test_case("F-14", "versions mode appends ggplot2 and ggrepel, or starts the table", {
  path <- file.path(work, "versions.tsv")
  writeLines(c("software\tversion", "R\t4.5.3", "DRIMSeq\t1.38.0"), path)
  versions <- figures$load_figure_packages()
  expected <- vapply(c("ggplot2", "ggrepel"), function(package) {
    as.character(utils::packageVersion(package))
  }, character(1))
  check(identical(versions, expected), "loaded versions: ", paste(names(versions), versions, collapse = ", "))
  figures$run_versions_mode(list(output = path), versions)
  table <- utils::read.delim(path, colClasses = "character")
  check(identical(table$software, c("R", "DRIMSeq", "ggplot2", "ggrepel")), "rows: ", paste(table$software, collapse = ", "))
  check(identical(table$version[3:4], unname(expected)), "versions: ", paste(table$version, collapse = ", "))
  fresh <- file.path(work, "fresh_versions.tsv")
  figures$run_versions_mode(list(output = fresh), versions)
  check(identical(readLines(fresh), c("software\tversion", paste0(names(expected), "\t", expected))), "new table: ", paste(readLines(fresh), collapse = " | "))
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
