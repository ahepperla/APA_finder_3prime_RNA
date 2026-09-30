#!/usr/bin/env Rscript

# Treatment-control figures, drawn with ggplot2 from the tables that
# fit_usage_model.R writes. The figures show the calls in those tables; they
# never make calls of their own.
#
# Usage:
#   Rscript plot_usage_figures.R --mode figures --samples normalized_samples.tsv \
#     --statistics-dir DIR --atlas pacs.v1.metadata.tsv.gz --pau observed_pau.tsv.gz \
#     --params resolved_params.yaml --output-dir figures
#   Rscript plot_usage_figures.R --mode versions --output software_versions.tsv
#
# Every string in this file is ASCII: pdf() draws text in Latin-1 Helvetica,
# and the container runs under LC_ALL=C.

# ---- Tables -----------------------------------------------------------------

PACS_REQUIRED <- c(
  "pac_id", "gene_id", "gene_name", "strand", "locus", "condition",
  "control_condition", "event_type", "fitted_control_pau", "fitted_treatment_pau",
  "delta_pau", "pac_fdr", "gene_fdr", "pvalue_pac", "assignment_class",
  "control_gene_total", "treatment_gene_total", "exploratory_insufficient_replicates"
)
PACS_NUMERIC <- c(
  "fitted_control_pau", "fitted_treatment_pau", "delta_pau", "pac_fdr", "gene_fdr",
  "pvalue_pac", "control_gene_total", "treatment_gene_total"
)
GENES_REQUIRED <- c(
  "gene_id", "gene_name", "condition", "control_condition", "dominant_switch",
  "complexity_change", "apa_pattern", "gene_fdr", "exploratory_insufficient_replicates"
)
DISTAL_COLUMNS <- c(
  "gene_id", "gene_name", "condition", "control_condition", "direction", "apa_pattern",
  "distal_pac_id", "distal_locus", "distal_assignment_class",
  "fitted_control_distal_pau", "fitted_treatment_distal_pau", "delta_distal_pau",
  "distal_event_type", "gene_fdr", "distal_pac_fdr", "tested_pacs"
)

# Assignment classes of tested PACs, in display order. Intergenic PACs have no
# gene and are never tested.
SITE_CLASSES <- c(
  terminal_exon = "Terminal exon",
  other_exon = "Other exon",
  intronic = "Intron",
  downstream = "Downstream of gene end"
)
# A gene's distal PAC is its most 3' tested PAC in one of these classes.
DISTAL_CLASSES <- c("terminal_exon", "downstream")

# ---- Calls ------------------------------------------------------------------

UP_EVENTS <- c("gained", "increased_usage")
DOWN_EVENTS <- c("lost", "decreased_usage")
CALLS <- c("up", "down", "up_candidate", "down_candidate", "none")
DISTAL_DIRECTIONS <- c(
  up = "distal_up", down = "distal_down", up_candidate = "distal_up_candidate",
  down_candidate = "distal_down_candidate", none = "none"
)
# Gene-level APA patterns from the genes tables, in the order apa_pattern
# lists them; a gene is drawn with its first. A pattern that only calls on
# PACs flagged for possible internal priming support carries the suffix; it
# keeps its pattern's color and is marked instead, since a lighter tint of
# intronic gain would match intronic loss.
POTENTIAL_INTERNAL_PRIMING <- "_potential_internal_priming"
APA_CLASSES <- c(
  intronic_gain = "Intronic gain", intronic_loss = "Intronic loss",
  alternative_last_exon = "Alternative last exon", utr_shortening = "UTR shortening",
  utr_lengthening = "UTR lengthening", other = "Other change", none = "No pattern"
)
PAC_EVENTS <- c(
  "gained", "increased_usage", "decreased_usage", "lost", "gained_candidate",
  "lost_candidate"
)
GENE_EVENTS <- c("dominant_switch", "complexity_gain", "complexity_loss")
EVENT_LABELS <- c(
  gained = "Gained", increased_usage = "Increased usage",
  decreased_usage = "Decreased usage", lost = "Lost",
  gained_candidate = "Gained (candidate)", lost_candidate = "Lost (candidate)",
  dominant_switch = "Dominant switch", complexity_gain = "Complexity gain",
  complexity_loss = "Complexity loss"
)

# ---- Drawing ----------------------------------------------------------------

# Figure sizes are whole eighths of an inch, so their pixel sizes at 200 dpi
# are exact: the quartz device converts pixels to inches and truncates.
FIGURE_DPI <- 200
# Above this many uncalled PACs in a panel, they are drawn as a density.
UNCALLED_POINT_LIMIT <- 5000
# Gene labels per direction on the volcano and distal-usage plots.
LABEL_LIMIT <- 20
# Genes on the shared-genes grid, and panels on the concordance figure.
GRID_GENE_LIMIT <- 50
CONCORDANCE_PANEL_LIMIT <- 15
# Okabe-Ito colors for conditions on the PCA, with the pale yellow last;
# shapes change every eight conditions.
CONDITION_COLORS <- c(
  "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#F0E442"
)
CONDITION_SHAPES <- c(16, 17, 15, 18)
# Concordance points called in both comparisons of a pair, or in one.
CONCORDANCE_COLORS <- c(both = "#CC79A7", one = "#56B4E9")
COLOR_UP <- "#D55E00"
COLOR_DOWN <- "#0072B2"
COLOR_NONE <- "grey70"
EVENT_COLORS <- c(
  gained = "#D55E00", increased_usage = "#E69F00", decreased_usage = "#56B4E9",
  lost = "#0072B2", gained_candidate = "#EEBF99", lost_candidate = "#99C7E0",
  dominant_switch = "#CC79A7", complexity_gain = "#009E73", complexity_loss = "#80CEB9"
)
DIRECTION_COLORS <- c(up = COLOR_UP, down = COLOR_DOWN)
APA_COLORS <- c(
  intronic_gain = "#8C510A", intronic_loss = "#D8B365", alternative_last_exon = "#762A83",
  utr_shortening = "#0072B2", utr_lengthening = "#D55E00", other = "#555555", none = "grey75"
)
CALL_DIRECTION_LABELS <- c(up = "Increased or gained", down = "Decreased or lost")
MARK_SHAPES <- c(
  "Confirmed" = 16, "Candidate" = 1, "Confirmed, p = 0" = 17, "Candidate, p = 0" = 2
)
# The distal PAC's call on the distal-usage plot; diamonds mark a gene whose
# pattern only flagged PACs support.
DISTAL_SHAPES <- c(
  "Confirmed" = 16, "None or candidate" = 1, "Confirmed *" = 18, "None or candidate *" = 5
)

# ---- Command line -----------------------------------------------------------

parse_args <- function(arguments) {
  result <- list()
  index <- 1
  while (index <= length(arguments)) {
    key <- sub("^--", "", arguments[[index]])
    if (index == length(arguments)) stop("Missing value for --", key)
    result[[gsub("-", "_", key)]] <- arguments[[index + 1]]
    index <- index + 2
  }
  if (is.null(result$mode)) result$mode <- "figures"
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

# pdf(timestamp = FALSE) needs R 4.5, and the figures use the ggplot2 4 API.
load_figure_packages <- function() {
  remedy <- ". Use the PACusage Conda environment or Apptainer image."
  if (getRversion() < "4.5.0") {
    stop("The figures need R 4.5.0 or later, not ", getRversion(), remedy)
  }
  missing <- c("ggplot2", "yaml")[!vapply(
    c("ggplot2", "yaml"), requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), remedy)
  version <- utils::packageVersion("ggplot2")
  if (version < "4.0.0") stop("The figures need ggplot2 4.0.0 or later, not ", version, remedy)
  suppressPackageStartupMessages(library(ggplot2))
  invisible(as.character(version))
}

# ---- Reading ----------------------------------------------------------------

as_flag <- function(values) tolower(values) %in% c("true", "t", "1")

ascii_text <- function(values) iconv(values, from = "UTF-8", to = "ASCII", sub = "?")

# Every column is read as text, so identifiers keep their exact form, and only
# the columns the figures use are converted. An empty field becomes NA.
read_table <- function(path, required, numeric = character(), flags = character()) {
  if (!file.exists(path)) stop("Missing table: ", path)
  table <- utils::read.delim(
    path, colClasses = "character", na.strings = character(), quote = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  missing <- setdiff(required, names(table))
  if (length(missing)) {
    stop(basename(path), " lacks columns: ", paste(missing, collapse = ", "), ".")
  }
  for (column in numeric) table[[column]] <- suppressWarnings(as.numeric(table[[column]]))
  for (column in flags) table[[column]] <- as_flag(table[[column]])
  table
}

read_pacs <- function(path) {
  read_table(path, PACS_REQUIRED, PACS_NUMERIC, "exploratory_insufficient_replicates")
}

# Observed PAU without the count and name columns. An empty PAU, at a gene
# without reads in the sample, is NA.
read_pau <- function(path) {
  read_columns(path, c(
    gene_id = "character", pac_id = "character", sample_id = "character",
    gene_total = "numeric", pau = "numeric"
  ))
}

read_genes <- function(path) {
  read_table(
    path, GENES_REQUIRED, "gene_fdr",
    c("dominant_switch", "exploratory_insufficient_replicates")
  )
}

# Only the named columns, read with the given classes; the others are skipped.
read_columns <- function(path, classes) {
  if (!file.exists(path)) stop("Missing table: ", path)
  header <- names(utils::read.delim(path, nrows = 0, quote = "", check.names = FALSE))
  missing <- setdiff(names(classes), header)
  if (length(missing)) {
    stop(basename(path), " lacks columns: ", paste(missing, collapse = ", "), ".")
  }
  selected <- rep("NULL", length(header))
  selected[match(names(classes), header)] <- unname(classes)
  utils::read.delim(
    path, colClasses = selected, na.strings = character(), quote = "", check.names = FALSE
  )
}

# Only the PAC IDs and coordinates, skipping the atlas's sequence columns.
read_atlas_coordinates <- function(path) {
  read_columns(path, c(pac_id = "character", coordinate = "numeric"))
}

write_gzip_tsv <- function(value, path) {
  connection <- gzfile(path, "wt")
  on.exit(close(connection), add = TRUE)
  write.table(value, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

# ---- Figure data ------------------------------------------------------------

# Each condition against its direct control, keyed by the stem that
# fit_usage_model.R gives the comparison's tables.
comparisons_from_samples <- function(samples) {
  pairs <- unique(samples[
    samples$condition != samples$control_condition, c("condition", "control_condition"),
    drop = FALSE
  ])
  if (!nrow(pairs)) stop("The sample sheet has no treatment-control comparison.")
  pairs$stem <- paste0(pairs$condition, "_vs_", pairs$control_condition)
  pairs$title <- ascii_text(paste(pairs$condition, "vs", pairs$control_condition))
  pairs <- pairs[order(pairs$stem, method = "radix"), , drop = FALSE]
  rownames(pairs) <- NULL
  pairs[, c("stem", "condition", "control_condition", "title")]
}

pac_calls <- function(event_type) {
  call <- rep("none", length(event_type))
  call[event_type %in% UP_EVENTS] <- "up"
  call[event_type %in% DOWN_EVENTS] <- "down"
  call[event_type == "gained_candidate"] <- "up_candidate"
  call[event_type == "lost_candidate"] <- "down_candidate"
  factor(call, levels = CALLS)
}

# One row per tested gene with a distal PAC: its most 3' tested PAC in the
# terminal exon or downstream of the gene. The direction is that PAC's own
# call, so a change among upstream PACs alone leaves the gene at "none"; the
# gene's APA pattern comes from the genes table.
distal_usage_table <- function(pacs, coordinates, genes) {
  eligible <- pacs[pacs$assignment_class %in% DISTAL_CLASSES, , drop = FALSE]
  coordinate <- coordinates$coordinate[match(eligible$pac_id, coordinates$pac_id)]
  if (anyNA(coordinate)) {
    stop("The atlas lacks coordinates for PACs: ",
      paste(eligible$pac_id[is.na(coordinate)], collapse = ", "))
  }
  # In transcript orientation: the largest coordinate on the plus strand, the
  # smallest on the minus strand.
  position <- ifelse(eligible$strand == "-", -coordinate, coordinate)
  eligible <- eligible[order(eligible$gene_id, -position, eligible$pac_id, method = "radix"), ,
    drop = FALSE]
  distal <- eligible[!duplicated(eligible$gene_id), , drop = FALSE]
  table <- data.frame(
    gene_id = distal$gene_id,
    gene_name = distal$gene_name,
    condition = distal$condition,
    control_condition = distal$control_condition,
    direction = unname(DISTAL_DIRECTIONS[as.character(pac_calls(distal$event_type))]),
    apa_pattern = genes$apa_pattern[match(distal$gene_id, genes$gene_id)],
    distal_pac_id = distal$pac_id,
    distal_locus = distal$locus,
    distal_assignment_class = distal$assignment_class,
    fitted_control_distal_pau = distal$fitted_control_pau,
    fitted_treatment_distal_pau = distal$fitted_treatment_pau,
    delta_distal_pau = distal$delta_pau,
    distal_event_type = distal$event_type,
    gene_fdr = distal$gene_fdr,
    distal_pac_fdr = distal$pac_fdr,
    tested_pacs = tabulate(match(pacs$gene_id, distal$gene_id), nbins = nrow(distal)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  table$apa_pattern[is.na(table$apa_pattern)] <- "none"
  table <- table[order(is.na(table$gene_fdr), table$gene_fdr, table$gene_id, method = "radix"), ,
    drop = FALSE]
  rownames(table) <- NULL
  list(table = table, without_distal = length(unique(pacs$gene_id)) - nrow(table))
}

# Confirmed up and down calls by assignment class, with every class and both
# directions present. Candidates and the gene-level labels are not counted.
site_class_counts <- function(pacs) {
  unknown <- setdiff(pacs$assignment_class, names(SITE_CLASSES))
  if (length(unknown)) {
    stop("Unexpected assignment_class values: ",
      paste(sort(unknown, method = "radix"), collapse = ", "), ".")
  }
  classes <- factor(pacs$assignment_class, levels = names(SITE_CLASSES))
  call <- pac_calls(pacs$event_type)
  bins <- length(SITE_CLASSES)
  data.frame(
    assignment_class = factor(rep(names(SITE_CLASSES), 2), levels = names(SITE_CLASSES)),
    direction = factor(rep(c("down", "up"), each = bins), levels = c("down", "up")),
    count = c(tabulate(classes[call == "down"], bins), tabulate(classes[call == "up"], bins)),
    tested = rep(tabulate(classes, bins), 2)
  )
}

# PAC events count rows of .pacs. Gene events count genes flagged in .genes:
# the .pacs labels for them go only on PACs without a call of their own, so
# they would undercount.
event_count_table <- function(comparisons, pacs_tables, genes_tables) {
  rows <- lapply(comparisons$stem, function(stem) {
    pacs <- pacs_tables[[stem]]
    genes <- genes_tables[[stem]]
    genes <- genes[!duplicated(genes$gene_id), , drop = FALSE]
    counts <- c(
      vapply(PAC_EVENTS, function(event) sum(pacs$event_type == event), integer(1)),
      dominant_switch = sum(genes$dominant_switch),
      complexity_gain = sum(genes$complexity_change == "gain"),
      complexity_loss = sum(genes$complexity_change == "loss")
    )
    data.frame(comparison = stem, event = names(counts), count = unname(counts))
  })
  counts <- do.call(rbind, rows)
  data.frame(
    comparison = factor(counts$comparison, levels = comparisons$stem),
    level = factor(
      ifelse(counts$event %in% GENE_EVENTS, "Gene events", "PAC events"),
      levels = c("PAC events", "Gene events")
    ),
    event = factor(counts$event, levels = c(PAC_EVENTS, GENE_EVENTS)),
    count = counts$count
  )
}

# Genes per APA pattern in each comparison, with zeros, apart for patterns
# that only flagged PACs support; a gene with two patterns counts in both.
# "other" is not a pattern, so it counts with the supported ones.
pattern_count_table <- function(comparisons, genes_tables) {
  classes <- setdiff(names(APA_CLASSES), "none")
  labels <- c(classes, paste0(setdiff(classes, "other"), POTENTIAL_INTERNAL_PRIMING))
  rows <- lapply(comparisons$stem, function(stem) {
    genes <- genes_tables[[stem]]
    genes <- genes[!duplicated(genes$gene_id), , drop = FALSE]
    parts <- strsplit(genes$apa_pattern, ";", fixed = TRUE)
    counts <- vapply(labels, function(label) {
      sum(vapply(parts, function(part) label %in% part, logical(1)))
    }, integer(1))
    data.frame(comparison = stem, label = labels, count = unname(counts))
  })
  counts <- do.call(rbind, rows)
  potential <- endsWith(counts$label, POTENTIAL_INTERNAL_PRIMING)
  data.frame(
    comparison = factor(counts$comparison, levels = comparisons$stem),
    pattern = factor(pattern_base(counts$label), levels = classes),
    support = factor(ifelse(potential, "flagged", "supported"), levels = PATTERN_SUPPORT),
    count = counts$count
  )
}

# Facets of the pattern counts: patterns that unflagged PACs support, and
# those that only PACs flagged for possible internal priming do.
PATTERN_SUPPORT <- c("supported", "flagged")
PATTERN_SUPPORT_LABELS <- c(
  supported = "Supported by unflagged PACs",
  flagged = "Only flagged PACs: potential internal priming"
)

pattern_base <- function(apa_pattern) {
  sub(paste0(POTENTIAL_INTERNAL_PRIMING, "$"), "", apa_pattern)
}

# A gene's first APA pattern, which colors it, without the suffix.
first_pattern <- function(apa_pattern) {
  first <- pattern_base(sub(";.*$", "", apa_pattern))
  first[is.na(first) | !first %in% names(APA_CLASSES)] <- "none"
  factor(first, levels = names(APA_CLASSES))
}

# Whether that first pattern is one only flagged PACs support.
first_pattern_flagged <- function(apa_pattern) {
  first <- sub(";.*$", "", apa_pattern)
  !is.na(first) & endsWith(first, POTENTIAL_INTERNAL_PRIMING)
}

# PACs with a change and a p-value. A p-value that underflows to 0 is drawn
# just above the most significant finite one; this happens before binning,
# which would otherwise drop it. Labels go on each gene's most significant
# confirmed PAC, for at most LABEL_LIMIT genes in each direction.
volcano_data <- function(pacs) {
  call <- pac_calls(pacs$event_type)
  kept <- is.finite(pacs$delta_pau) & !is.na(pacs$pvalue_pac)
  points <- data.frame(
    pac_id = pacs$pac_id[kept],
    gene_id = pacs$gene_id[kept],
    gene_name = ascii_text(pacs$gene_name[kept]),
    delta_pau = pacs$delta_pau[kept],
    pvalue = pacs$pvalue_pac[kept],
    call = call[kept],
    stringsAsFactors = FALSE
  )
  points$capped <- points$pvalue == 0
  score <- -log10(points$pvalue)
  score[points$capped] <- max(1, score[!points$capped]) * 1.05
  points$neg_log10_p <- score
  # Uncalled PACs are drawn first, under the calls.
  points <- points[order(points$call != "none", points$pac_id, method = "radix"), , drop = FALSE]
  confirmed <- points[points$call %in% c("up", "down"), , drop = FALSE]
  confirmed <- confirmed[order(confirmed$pvalue, confirmed$pac_id, method = "radix"), ,
    drop = FALSE]
  labels <- limit_labels(confirmed[!duplicated(confirmed$gene_id), , drop = FALSE], "call")
  list(
    points = points,
    labels = labels,
    missing = sum(!kept),
    missing_candidates = sum(!kept & call %in% c("up_candidate", "down_candidate"))
  )
}

# The first LABEL_LIMIT rows of each direction, keeping the rows' order, which
# check_overlap also follows: an earlier label wins over a later one it covers.
limit_labels <- function(rows, direction) {
  rank <- stats::ave(seq_len(nrow(rows)), as.character(rows[[direction]]), FUN = seq_along)
  rows[rank <= LABEL_LIMIT, , drop = FALSE]
}

# Genes whose distal PAC has a confirmed call, largest change first, for at
# most LABEL_LIMIT genes in each direction. Candidates are not labelled.
distal_labels <- function(table) {
  called <- table[table$direction %in% c("distal_up", "distal_down") &
    is.finite(table$delta_distal_pau), , drop = FALSE]
  called <- called[order(-abs(called$delta_distal_pau), called$distal_pac_fdr, called$gene_id,
    method = "radix"), , drop = FALSE]
  limit_labels(called, "direction")
}

# Each PAC's change against the reads at its gene in the less-covered group. A
# group without reads has no fitted usage, so its rows have no change and drop
# out; no pseudocount is added.
coverage_data <- function(comparisons, pacs_tables) {
  rows <- lapply(comparisons$stem, function(stem) {
    pacs <- pacs_tables[[stem]]
    depth <- pmin(pacs$control_gene_total, pacs$treatment_gene_total)
    kept <- !is.na(depth) & depth > 0 & is.finite(pacs$delta_pau)
    data.frame(
      comparison = rep(stem, sum(kept)),
      pac_id = pacs$pac_id[kept],
      depth = depth[kept],
      delta_pau = pacs$delta_pau[kept],
      call = as.character(pac_calls(pacs$event_type[kept])),
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, rows)
  data$comparison <- factor(data$comparison, levels = comparisons$stem)
  data$call <- factor(data$call, levels = CALLS)
  data <- data[order(data$call != "none", data$comparison, data$pac_id, method = "radix"), ,
    drop = FALSE]
  rownames(data) <- NULL
  data
}

# ---- Across comparisons -----------------------------------------------------

# Every tested gene's pattern in each comparison, empty where the gene was not
# tested, ordered by the number of comparisons with a pattern, then by the
# gene's best FDR.
pattern_grid_table <- function(comparisons, genes_tables) {
  long <- do.call(rbind, c(
    list(data.frame(gene_id = character(), gene_name = character(), comparison = character(),
      apa_pattern = character(), gene_fdr = numeric(), stringsAsFactors = FALSE)),
    lapply(comparisons$stem, function(stem) {
      genes <- genes_tables[[stem]]
      genes <- genes[!duplicated(genes$gene_id), , drop = FALSE]
      data.frame(gene_id = genes$gene_id, gene_name = genes$gene_name,
        comparison = rep(stem, nrow(genes)), apa_pattern = genes$apa_pattern,
        gene_fdr = genes$gene_fdr, stringsAsFactors = FALSE)
    })
  ))
  ids <- unique(long$gene_id)
  table <- data.frame(
    gene_id = ids, gene_name = long$gene_name[match(ids, long$gene_id)],
    patterned_comparisons = integer(length(ids)), stringsAsFactors = FALSE, check.names = FALSE
  )
  for (stem in comparisons$stem) {
    rows <- long[long$comparison == stem, , drop = FALSE]
    table[[stem]] <- rows$apa_pattern[match(ids, rows$gene_id)]
    table$patterned_comparisons <- table$patterned_comparisons +
      as.integer(!is.na(table[[stem]]) & !table[[stem]] %in% c("", "none"))
  }
  finite <- long[is.finite(long$gene_fdr), , drop = FALSE]
  best <- unname(tapply(finite$gene_fdr, finite$gene_id, min)[ids])
  table <- table[order(-table$patterned_comparisons, is.na(best), best, table$gene_id,
    method = "radix"), , drop = FALSE]
  rownames(table) <- NULL
  table
}

# Pairs of comparisons with how they are related: sharing a control, or
# chained through a condition that is one's treatment and the other's control.
concordance_pairs <- function(comparisons) {
  pairs <- if (nrow(comparisons) < 2) {
    matrix(integer(), nrow = 2)
  } else {
    utils::combn(nrow(comparisons), 2)
  }
  a <- comparisons[pairs[1, ], , drop = FALSE]
  b <- comparisons[pairs[2, ], , drop = FALSE]
  relation <- as.character(ifelse(
    a$control_condition == b$control_condition, "shared_control",
    ifelse(a$condition == b$control_condition | b$condition == a$control_condition,
      "chained", "unrelated")
  ))
  data.frame(
    pair = paste(a$stem, b$stem, sep = "|"), comparison_a = a$stem, comparison_b = b$stem,
    title_a = a$title, title_b = b$title, relation = relation, stringsAsFactors = FALSE
  )
}

# The pairs drawn as panels: every pair up to the limit; beyond it, only the
# related ones, and at most the limit of those.
shown_pairs <- function(pairs) {
  if (nrow(pairs) <= CONCORDANCE_PANEL_LIMIT) return(pairs)
  related <- pairs[pairs$relation != "unrelated", , drop = FALSE]
  utils::head(related, CONCORDANCE_PANEL_LIMIT)
}

# Each PAC tested in both comparisons of a pair, with its change in each and
# whether it has a confirmed call in both, one, or neither.
concordance_data <- function(pairs, pacs_tables) {
  rows <- lapply(seq_len(nrow(pairs)), function(index) {
    a <- pacs_tables[[pairs$comparison_a[[index]]]]
    b <- pacs_tables[[pairs$comparison_b[[index]]]]
    shared <- intersect(a$pac_id, b$pac_id)
    x <- a$delta_pau[match(shared, a$pac_id)]
    y <- b$delta_pau[match(shared, b$pac_id)]
    called_a <- pac_calls(a$event_type[match(shared, a$pac_id)]) %in% c("up", "down")
    called_b <- pac_calls(b$event_type[match(shared, b$pac_id)]) %in% c("up", "down")
    kept <- is.finite(x) & is.finite(y)
    call <- ifelse(called_a & called_b, "both", ifelse(called_a | called_b, "one", "none"))
    data.frame(
      pair = rep(pairs$pair[[index]], sum(kept)), pac_id = shared[kept], delta_a = x[kept],
      delta_b = y[kept], call = call[kept], stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, c(list(data.frame(pair = character(), pac_id = character(),
    delta_a = numeric(), delta_b = numeric(), call = character(), stringsAsFactors = FALSE)),
    rows))
  data$pair <- factor(data$pair, levels = pairs$pair)
  data$call <- factor(data$call, levels = c("both", "one", "none"))
  # Drawing order: uncalled PACs first, and PACs called in both on top.
  data <- data[order(-as.integer(data$call), data$pair, data$pac_id, method = "radix"), ,
    drop = FALSE]
  rownames(data) <- NULL
  data
}

# One row per pair: shared PACs, PACs called in both, and the Pearson
# correlation of their changes (empty with fewer than 3 PACs or no spread).
concordance_summary <- function(pairs, data) {
  summary <- lapply(seq_len(nrow(pairs)), function(index) {
    rows <- data[data$pair == pairs$pair[[index]], , drop = FALSE]
    spread <- nrow(rows) >= 3L && stats::sd(rows$delta_a) > 0 && stats::sd(rows$delta_b) > 0
    data.frame(
      comparison_a = pairs$comparison_a[[index]], comparison_b = pairs$comparison_b[[index]],
      relation = pairs$relation[[index]], shared_pacs = nrow(rows),
      called_in_both = sum(rows$call == "both"),
      pearson_r = if (spread) stats::cor(rows$delta_a, rows$delta_b) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, c(list(data.frame(comparison_a = character(), comparison_b = character(),
    relation = character(), shared_pacs = integer(), called_in_both = integer(),
    pearson_r = numeric(), stringsAsFactors = FALSE)), summary))
}

# A PCA of observed PAU across samples, by the report's former rule: genes
# with at least minimum reads in every sample, all their PACs, and no
# zero-filling. Each component's sign makes its largest loading positive. A
# two-PAC gene's loadings tie in size, since its PAU sum to 1, so the first
# PAC by ID within a relative 1e-6 of the largest decides.
pau_pca <- function(pau, samples, minimum) {
  sample_ids <- sort(unique(samples$sample_id), method = "radix")
  totals <- tapply(pau$gene_total, pau$gene_id, min)
  covered <- names(totals)[!is.na(totals) & totals >= minimum]
  pau <- pau[pau$gene_id %in% covered & pau$sample_id %in% sample_ids, , drop = FALSE]
  pac_ids <- sort(unique(pau$pac_id), method = "radix")
  values <- matrix(NA_real_, nrow = length(pac_ids), ncol = length(sample_ids))
  values[cbind(match(pau$pac_id, pac_ids), match(pau$sample_id, sample_ids))] <- pau$pau
  complete <- stats::complete.cases(values)
  values <- values[complete, , drop = FALSE]
  genes <- length(unique(pau$gene_id[pau$pac_id %in% pac_ids[complete]]))
  if (length(sample_ids) < 2L || !nrow(values)) return(NULL)
  centered <- t(values) - matrix(rowMeans(values), nrow = length(sample_ids),
    ncol = nrow(values), byrow = TRUE)
  decomposition <- svd(centered)
  components <- min(2L, length(decomposition$d))
  scores <- matrix(0, nrow = length(sample_ids), ncol = 2)
  for (k in seq_len(components)) {
    loading <- decomposition$v[, k]
    largest <- which(abs(loading) >= max(abs(loading)) * (1 - 1e-6))[[1]]
    flip <- if (loading[[largest]] < 0) -1 else 1
    scores[, k] <- flip * decomposition$u[, k] * decomposition$d[[k]]
  }
  total <- sum(decomposition$d^2)
  variance <- if (total > 0) decomposition$d^2 / total else rep(0, length(decomposition$d))
  variance <- c(variance, 0, 0)[1:2]
  condition <- samples$condition[match(sample_ids, samples$sample_id)]
  result <- data.frame(
    sample_id = sample_ids, condition = condition, PC1 = scores[, 1], PC2 = scores[, 2],
    pc1_variance_fraction = variance[[1]], pc2_variance_fraction = variance[[2]],
    stringsAsFactors = FALSE
  )
  attr(result, "genes") <- genes
  attr(result, "pacs") <- nrow(values)
  result
}

# ---- Figures ----------------------------------------------------------------

figure_theme <- function() {
  theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.title.position = "plot",
      plot.subtitle = element_text(colour = "grey25"),
      legend.position = "bottom"
    )
}

plural <- function(count, word) paste(count, if (count == 1) word else paste0(word, "s"))

placeholder_plot <- function(title, message) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = message, size = 3.5, colour = "grey30") +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold")) +
    labs(title = title)
}

# Uncalled PACs as grey points, or as a grey density in any panel with more
# than UNCALLED_POINT_LIMIT of them. Capped points are triangles.
uncalled_layers <- function(data, panel = NULL) {
  uncalled <- data[data$call == "none", , drop = FALSE]
  if (!nrow(uncalled)) return(list())
  sizes <- if (is.null(panel)) {
    rep(nrow(uncalled), nrow(uncalled))
  } else {
    as.vector(table(uncalled[[panel]])[as.character(uncalled[[panel]])])
  }
  crowded <- sizes > UNCALLED_POINT_LIMIT
  capped <- if (is.null(uncalled$capped)) rep(FALSE, nrow(uncalled)) else uncalled$capped
  layers <- list()
  if (any(!crowded & !capped)) {
    layers <- c(layers, list(geom_point(
      data = uncalled[!crowded & !capped, , drop = FALSE], colour = COLOR_NONE, size = 0.8,
      shape = 16
    )))
  }
  if (any(!crowded & capped)) {
    layers <- c(layers, list(geom_point(
      data = uncalled[!crowded & capped, , drop = FALSE], colour = COLOR_NONE, size = 0.9,
      shape = 17
    )))
  }
  if (any(crowded)) {
    layers <- c(layers, list(
      geom_bin_2d(data = uncalled[crowded, , drop = FALSE], bins = c(80, 60)),
      scale_fill_gradient(low = "grey88", high = "grey35", name = "Uncalled PACs",
        guide = guide_colourbar(order = 3))
    ))
  }
  layers
}

# Colour by direction and shape by confirmed or candidate call, for the PACs
# with a call.
called_points <- function(data) {
  called <- data[data$call != "none", , drop = FALSE]
  candidate <- called$call %in% c("up_candidate", "down_candidate")
  capped <- if (is.null(called$capped)) rep(FALSE, nrow(called)) else called$capped
  mark <- ifelse(candidate, "Candidate", "Confirmed")
  mark[capped] <- paste0(mark[capped], ", p = 0")
  called$direction <- factor(
    ifelse(called$call %in% c("up", "up_candidate"), "up", "down"), levels = c("up", "down")
  )
  called$mark <- factor(mark, levels = names(MARK_SHAPES))
  called
}

# Legends in a fixed order: direction, then call type, then any density.
call_scales <- function() {
  list(
    scale_colour_manual(values = DIRECTION_COLORS, labels = CALL_DIRECTION_LABELS, name = NULL,
      guide = guide_legend(order = 1)),
    scale_shape_manual(values = MARK_SHAPES, name = NULL, guide = guide_legend(order = 2))
  )
}

plot_volcano <- function(pacs, comparison, params, exploratory) {
  if (!nrow(pacs)) return(placeholder_plot(comparison$title, "No tested PACs"))
  data <- volcano_data(pacs)
  calls <- pac_calls(pacs$event_type)
  subtitle <- sprintf(
    "%s tested in %s, %d with a confirmed call",
    plural(nrow(pacs), "PAC"), plural(length(unique(pacs$gene_id)), "gene"),
    sum(calls %in% c("up", "down"))
  )
  if (data$missing) {
    subtitle <- paste0(subtitle, sprintf(
      "\n%s without a p-value, not shown (%d of them candidate calls)",
      plural(data$missing, "PAC"), data$missing_candidates
    ))
  }
  if (exploratory) {
    subtitle <- paste0(subtitle, "\nExploratory comparison: gained and lost calls are candidates only")
  }
  if (!nrow(data$points)) return(placeholder_plot(comparison$title, "No tested PAC has a p-value"))
  largest <- max(abs(data$points$delta_pau))
  limit <- min(1, max(0.25, ceiling(round(largest * 10, 6)) / 10))
  threshold <- params$min_abs_delta_pau
  ggplot(data$points, aes(x = delta_pau, y = neg_log10_p)) +
    uncalled_layers(data$points) +
    geom_vline(xintercept = c(-threshold, threshold), linetype = "dashed", colour = "grey50",
      linewidth = 0.3) +
    geom_point(data = called_points(data$points), aes(colour = direction, shape = mark),
      size = 1.6) +
    geom_text(data = data$labels, aes(label = gene_name), size = 2.4, vjust = -0.7,
      check_overlap = TRUE) +
    call_scales() +
    coord_cartesian(xlim = c(-limit, limit)) +
    expand_limits(y = 0) +
    labs(
      title = comparison$title,
      subtitle = subtitle,
      x = ascii_text(sprintf(
        "Change in PAU (%s - %s)", comparison$condition, comparison$control_condition
      )),
      y = "-log10 PAC p-value"
    ) +
    figure_theme()
}

plot_distal_usage <- function(result, comparison, params) {
  table <- result$table
  if (!nrow(table)) {
    message <- "No tested gene has a PAC in its terminal exon or downstream"
    return(placeholder_plot(comparison$title, message))
  }
  fitted <- is.finite(table$fitted_control_distal_pau) &
    is.finite(table$fitted_treatment_distal_pau)
  shown <- table[fitted, , drop = FALSE]
  shown$pattern <- first_pattern(shown$apa_pattern)
  flagged <- first_pattern_flagged(shown$apa_pattern)
  call <- ifelse(shown$direction %in% c("distal_up", "distal_down"), "Confirmed", "None or candidate")
  shown$call <- factor(paste0(call, ifelse(flagged, " *", "")), levels = names(DISTAL_SHAPES))
  # Genes without a pattern first, under the others.
  shown <- shown[order(shown$pattern != "none", shown$gene_id, method = "radix"), ,
    drop = FALSE]
  subtitle <- plural(nrow(table), "gene")
  if (result$without_distal) {
    subtitle <- paste0(subtitle, sprintf(", and %d without a distal PAC", result$without_distal))
  }
  if (any(!fitted)) {
    subtitle <- paste0(subtitle, sprintf("; %d without fitted usage, not shown", sum(!fitted)))
  }
  subtitle <- paste0(
    subtitle, "\nDistal PAC: the most 3' tested PAC in a last exon or downstream"
  )
  if (any(flagged)) {
    subtitle <- paste0(subtitle, "\nDiamonds: pattern only from PACs flagged for internal priming")
  }
  labels <- distal_labels(shown)
  labels$gene_name <- ascii_text(labels$gene_name)
  counts <- sprintf(
    "Distal PAC up: %d\nDistal PAC down: %d",
    sum(table$direction == "distal_up"), sum(table$direction == "distal_down")
  )
  candidates <- sum(table$direction %in% c("distal_up_candidate", "distal_down_candidate"))
  if (candidates) counts <- paste0(counts, sprintf("\nCandidates: %d", candidates))
  threshold <- params$min_abs_delta_pau
  ggplot(shown, aes(x = fitted_control_distal_pau, y = fitted_treatment_distal_pau)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = c(-threshold, threshold), colour = "grey60",
      linetype = "dashed", linewidth = 0.3) +
    geom_point(aes(colour = pattern, shape = call), size = 1.6) +
    geom_text(data = labels, aes(label = gene_name), size = 2.4, vjust = -0.7,
      check_overlap = TRUE) +
    annotate("text", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 2.8, label = counts) +
    scale_colour_manual(
      values = APA_COLORS, labels = APA_CLASSES, name = "APA pattern",
      guide = guide_legend(order = 1, ncol = 2)
    ) +
    scale_shape_manual(
      values = DISTAL_SHAPES, name = "Distal PAC call", guide = guide_legend(order = 2, ncol = 2)
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = comparison$title,
      subtitle = subtitle,
      x = ascii_text(sprintf("Distal PAC usage in %s (fitted PAU)", comparison$control_condition)),
      y = ascii_text(sprintf("Distal PAC usage in %s (fitted PAU)", comparison$condition))
    ) +
    figure_theme() +
    # Two legends side by side would not fit the square figure; stacked, they
    # are kept tight so the panel stays large.
    theme(
      legend.box = "vertical", legend.spacing.y = grid::unit(2, "pt"),
      legend.key.spacing.y = grid::unit(0, "pt"), legend.margin = margin(0, 0, 0, 0)
    )
}

plot_site_classes <- function(counts, comparison) {
  tested <- counts$tested[counts$direction == "down"]
  if (!sum(tested)) return(placeholder_plot(comparison$title, "No tested PACs"))
  labels <- sprintf("%s (n = %d)", SITE_CLASSES, tested)
  # Terminal exon at the top.
  counts$label <- factor(labels[as.integer(counts$assignment_class)], levels = rev(labels))
  counts$signed <- ifelse(counts$direction == "down", -counts$count, counts$count)
  limit <- max(1, abs(counts$signed)) * 1.25
  shown <- counts[counts$count > 0, , drop = FALSE]
  whole_breaks <- function(limits) {
    breaks <- pretty(limits)
    breaks[breaks == round(breaks)]
  }
  plot <- ggplot(counts, aes(x = signed, y = label, fill = direction)) +
    geom_col(width = 0.6) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
    geom_text(data = shown, aes(label = count, hjust = ifelse(signed < 0, 1.3, -0.3)),
      size = 2.6) +
    scale_fill_manual(
      values = c(down = COLOR_DOWN, up = COLOR_UP),
      labels = c(down = "Lost or decreased", up = "Gained or increased"),
      name = NULL, drop = FALSE
    ) +
    scale_x_continuous(limits = c(-limit, limit), breaks = whole_breaks,
      labels = function(breaks) format(abs(breaks), trim = TRUE)) +
    labs(
      title = comparison$title,
      subtitle = "Confirmed PAC calls by site class: lost usage left of zero, gained right",
      x = "PACs", y = NULL
    ) +
    figure_theme()
  if (!nrow(shown)) {
    plot <- plot + annotate("text", x = 0, y = 2.5, label = "No confirmed PAC calls",
      size = 3, colour = "grey30")
  }
  plot
}

# A count axis: whole-number breaks, and at least 0 to 1, so an empty facet
# gets no fractional axis. The right-hand room holds the bar totals.
count_axis <- function() {
  scale_x_continuous(
    limits = function(range) c(0, max(range[[2]], 1)),
    breaks = function(limits) {
      values <- pretty(limits)
      # pretty() can miss a whole number by a rounding error.
      round(values[abs(values - round(values)) < 1e-9])
    },
    expand = expansion(mult = c(0, 0.2))
  )
}

# Two plots stacked in one figure: the PAC and gene events, and the genes per
# APA pattern, each with its own legend.
plot_event_counts <- function(counts, patterns, comparisons) {
  rows <- nrow(comparisons)
  figure_stack(
    list(plot_event_bars(counts, comparisons), plot_pattern_bars(patterns, comparisons)),
    heights = c(2 + 0.375 * rows, 1.5 + 0.375 * rows)
  )
}

figure_stack <- function(plots, heights) {
  structure(list(plots = plots, heights = heights), class = "figure_stack")
}

plot_pattern_bars <- function(patterns, comparisons) {
  titles <- stats::setNames(comparisons$title, comparisons$stem)
  totals <- stats::aggregate(count ~ comparison + support, data = patterns, FUN = sum)
  ggplot(patterns, aes(x = count, y = comparison, fill = pattern)) +
    geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
    geom_text(data = totals, aes(x = count, y = comparison, label = count),
      inherit.aes = FALSE, hjust = -0.3, size = 2.6) +
    facet_wrap(~support, scales = "free_x", labeller = as_labeller(PATTERN_SUPPORT_LABELS)) +
    scale_fill_manual(values = APA_COLORS, labels = APA_CLASSES, name = NULL, drop = FALSE) +
    scale_y_discrete(limits = rev(comparisons$stem), labels = titles) +
    count_axis() +
    guides(fill = guide_legend(nrow = 2)) +
    labs(
      title = "APA patterns per comparison",
      subtitle = "Genes per pattern, from the genes tables; a gene with two patterns counts in both",
      x = "Genes", y = NULL
    ) +
    figure_theme() +
    theme(panel.spacing.x = grid::unit(1, "lines"))
}

plot_event_bars <- function(counts, comparisons) {
  titles <- stats::setNames(comparisons$title, comparisons$stem)
  totals <- stats::aggregate(count ~ comparison + level, data = counts, FUN = sum)
  ggplot(counts, aes(x = count, y = comparison, fill = event)) +
    geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
    geom_text(data = totals, aes(x = count, y = comparison, label = count),
      inherit.aes = FALSE, hjust = -0.3, size = 2.6) +
    facet_wrap(~level, scales = "free_x") +
    scale_fill_manual(values = EVENT_COLORS, labels = EVENT_LABELS, name = NULL, drop = FALSE) +
    scale_y_discrete(limits = rev(comparisons$stem), labels = titles) +
    count_axis() +
    guides(fill = guide_legend(nrow = 3)) +
    labs(
      title = "Events per comparison",
      subtitle = "PAC events count PACs; gene events count genes that pass the gene-level screen",
      x = "Count", y = NULL
    ) +
    figure_theme() +
    theme(panel.spacing.x = grid::unit(1, "lines"))
}

plot_effect_vs_coverage <- function(data, comparisons, params, dropped) {
  title <- "Change in usage against gene coverage"
  if (!nrow(data)) return(placeholder_plot(title, "No tested PAC has reads in both groups"))
  titles <- stats::setNames(comparisons$title, comparisons$stem)
  subtitle <- "Dashed line: min_gene_total, the coverage a gained or lost call requires"
  if (dropped) {
    subtitle <- paste0(subtitle, sprintf(
      "\n%s without reads in both groups, not shown", plural(dropped, "PAC")
    ))
  }
  ggplot(data, aes(x = depth, y = delta_pau)) +
    uncalled_layers(data, panel = "comparison") +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
    geom_vline(xintercept = params$min_gene_total, linetype = "dashed", colour = "grey50",
      linewidth = 0.3) +
    geom_point(data = called_points(data), aes(colour = direction, shape = mark), size = 1.2) +
    call_scales() +
    scale_x_log10() +
    facet_wrap(~comparison, ncol = min(3, nrow(comparisons)), drop = FALSE,
      labeller = as_labeller(titles)) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(
      title = title, subtitle = subtitle,
      x = "Reads at the gene in the less-covered group", y = "Change in PAU"
    ) +
    figure_theme()
}

# Inches rounded up to a whole eighth.
eighths <- function(inches) ceiling(inches * 8 - 1e-9) / 8

# Comparison titles on two lines, for the grid's vertical column labels.
grid_titles <- function(comparisons) {
  stats::setNames(
    ascii_text(paste0(comparisons$condition, "\nvs ", comparisons$control_condition)),
    comparisons$stem
  )
}

# Room for the grid's vertical column labels, from their longest line.
grid_label_height <- function(comparisons) {
  lines <- unlist(strsplit(grid_titles(comparisons), "\n", fixed = TRUE))
  eighths(0.06 * max(nchar(lines)))
}

# The genes drawn on the grid: those with a pattern in at least two
# comparisons, up to GRID_GENE_LIMIT of them in the table's order.
grid_genes <- function(table) {
  utils::head(table[table$patterned_comparisons >= 2, , drop = FALSE], GRID_GENE_LIMIT)
}

plot_pattern_grid <- function(table, comparisons) {
  title <- "APA patterns shared between comparisons"
  if (nrow(comparisons) < 2) return(placeholder_plot(title, "Only one comparison"))
  shared <- sum(table$patterned_comparisons >= 2)
  if (!shared) {
    return(placeholder_plot(title, "No gene has a pattern in two or more comparisons"))
  }
  shown <- grid_genes(table)
  cells <- do.call(rbind, lapply(comparisons$stem, function(stem) {
    data.frame(gene_id = shown$gene_id, comparison = stem, value = shown[[stem]],
      stringsAsFactors = FALSE)
  }))
  fill <- as.character(first_pattern(cells$value))
  fill[is.na(cells$value)] <- "not_tested"
  cells$fill <- factor(fill, levels = c(names(APA_CLASSES), "not_tested"))
  cells$gene_id <- factor(cells$gene_id, levels = rev(shown$gene_id))
  cells$comparison <- factor(cells$comparison, levels = comparisons$stem)
  colors <- c(APA_COLORS[names(APA_COLORS) != "none"], none = "white", not_tested = "grey88")
  # The mark is white on a dark fill and black on a light one.
  dark <- colSums(grDevices::col2rgb(colors) * c(0.299, 0.587, 0.114)) < 128
  several <- !is.na(cells$value) & grepl(";", cells$value, fixed = TRUE)
  cells$mark <- paste0(ifelse(several, "+", ""), ifelse(first_pattern_flagged(cells$value), "*", ""))
  marked <- cells[cells$mark != "", , drop = FALSE]
  marked$ink <- ifelse(dark[as.character(marked$fill)], "white", "black")
  subtitle <- sprintf("%s with a pattern in 2 or more comparisons", plural(shared, "gene"))
  if (shared > GRID_GENE_LIMIT) {
    subtitle <- paste0(subtitle, sprintf("\nThe first %d are shown; the table lists all",
      GRID_GENE_LIMIT))
  }
  subtitle <- paste0(
    subtitle, "\n+ marks two or more patterns in one comparison",
    "\n* marks a pattern only from flagged PACs"
  )
  ggplot(cells, aes(x = comparison, y = gene_id)) +
    # show.legend = TRUE draws keys for patterns absent from the drawn genes.
    geom_tile(aes(fill = fill), colour = "grey60", linewidth = 0.2, show.legend = TRUE) +
    geom_text(data = marked, aes(label = mark, colour = ink), size = 3) +
    scale_fill_manual(values = colors, labels = c(APA_CLASSES, not_tested = "Not tested"),
      name = NULL, drop = FALSE) +
    scale_colour_identity() +
    scale_x_discrete(labels = grid_titles(comparisons), expand = c(0, 0)) +
    scale_y_discrete(labels = stats::setNames(ascii_text(shown$gene_name), shown$gene_id),
      expand = c(0, 0)) +
    guides(fill = guide_legend(nrow = 3)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    figure_theme() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 6.5)
    )
}

plot_concordance <- function(data, pairs, summary, all_pairs) {
  title <- "Concordance between comparisons"
  if (!all_pairs) return(placeholder_plot(title, "Only one comparison"))
  if (!nrow(pairs)) {
    return(placeholder_plot(title, "No pair shares a control or a condition"))
  }
  labels <- stats::setNames(paste0("x: ", pairs$title_a, "\ny: ", pairs$title_b), pairs$pair)
  shown <- data[data$pair %in% pairs$pair, , drop = FALSE]
  shown$pair <- factor(as.character(shown$pair), levels = pairs$pair)
  notes <- summary[match(pairs$pair, paste(summary$comparison_a, summary$comparison_b,
    sep = "|")), , drop = FALSE]
  notes$pair <- factor(pairs$pair, levels = pairs$pair)
  notes$label <- sprintf("r = %s\nn = %d",
    ifelse(is.na(notes$pearson_r), "NA", sprintf("%.2f", notes$pearson_r)), notes$shared_pacs)
  subtitle <- paste0(
    "Change in PAU at PACs tested in both comparisons. Comparisons that share a\n",
    "control correlate positively through its estimate alone, and chained ones negatively"
  )
  if (nrow(pairs) < all_pairs) {
    subtitle <- paste0(subtitle, sprintf(
      "\n%d of %d pairs shown, those sharing a control or a condition; the matrix has all",
      nrow(pairs), all_pairs
    ))
  }
  called <- shown[shown$call != "none", , drop = FALSE]
  called_layers <- if (nrow(called)) {
    list(
      geom_point(data = called, aes(colour = call), size = 1.2),
      scale_colour_manual(
        values = CONCORDANCE_COLORS,
        labels = c(both = "Confirmed call in both", one = "Confirmed call in one"),
        name = NULL, guide = guide_legend(order = 1)
      )
    )
  }
  ggplot(shown, aes(x = delta_a, y = delta_b)) +
    uncalled_layers(shown, panel = "pair") +
    geom_abline(slope = 1, intercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_abline(slope = -1, intercept = 0, colour = "grey60", linetype = "dashed",
      linewidth = 0.3) +
    called_layers +
    geom_text(data = notes, aes(x = -0.95, y = 0.95, label = label), inherit.aes = FALSE,
      hjust = 0, vjust = 1, size = 2.5) +
    facet_wrap(~pair, ncol = 3, drop = FALSE, labeller = as_labeller(labels)) +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
    labs(title = title, subtitle = subtitle, x = "Change in PAU, x comparison",
      y = "Change in PAU, y comparison") +
    figure_theme()
}

# Pearson r for every ordered pair of comparisons, with 1 on the diagonal.
concordance_matrix <- function(summary, comparisons) {
  stems <- comparisons$stem
  grid <- expand.grid(a = stems, b = stems, stringsAsFactors = FALSE)
  key <- paste(grid$a, grid$b, sep = "|")
  forward <- match(key, paste(summary$comparison_a, summary$comparison_b, sep = "|"))
  backward <- match(key, paste(summary$comparison_b, summary$comparison_a, sep = "|"))
  grid$r <- ifelse(grid$a == grid$b, 1,
    ifelse(!is.na(forward), summary$pearson_r[forward], summary$pearson_r[backward]))
  grid
}

plot_concordance_matrix <- function(summary, comparisons) {
  title <- "Correlation between comparisons"
  if (nrow(comparisons) < 2) return(placeholder_plot(title, "Only one comparison"))
  grid <- concordance_matrix(summary, comparisons)
  grid$label <- ifelse(is.na(grid$r), "NA", sprintf("%.2f", grid$r))
  grid$a <- factor(grid$a, levels = comparisons$stem)
  grid$b <- factor(grid$b, levels = rev(comparisons$stem))
  titles <- grid_titles(comparisons)
  ggplot(grid, aes(x = a, y = b, fill = r)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = label), size = 2.8) +
    scale_fill_gradient2(low = COLOR_DOWN, mid = "white", high = COLOR_UP, midpoint = 0,
      limits = c(-1, 1), name = "Pearson r", na.value = "grey88") +
    scale_x_discrete(labels = titles, expand = c(0, 0)) +
    scale_y_discrete(labels = titles, expand = c(0, 0)) +
    coord_equal() +
    labs(
      title = title,
      subtitle = paste0(
        "Pearson r of the change in PAU at shared PACs.\n",
        "Shared controls and chains correlate by design"
      ),
      x = NULL, y = NULL
    ) +
    figure_theme() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
}

plot_pau_pca <- function(pca, params) {
  title <- "PAU principal components"
  if (is.null(pca)) {
    return(placeholder_plot(title, "No PAC has observed PAU in every sample"))
  }
  conditions <- sort(unique(pca$condition), method = "radix")
  index <- seq_along(conditions) - 1L
  colors <- stats::setNames(CONDITION_COLORS[index %% length(CONDITION_COLORS) + 1L],
    conditions)
  shapes <- stats::setNames(
    CONDITION_SHAPES[(index %/% length(CONDITION_COLORS)) %% length(CONDITION_SHAPES) + 1L],
    conditions
  )
  pca$condition <- factor(pca$condition, levels = conditions)
  pca$label <- ascii_text(pca$sample_id)
  subtitle <- sprintf(
    "Observed PAU at %s in %s with at least %s reads in every sample",
    plural(attr(pca, "pacs"), "PAC"), plural(attr(pca, "genes"), "gene"),
    format(params$min_gene_total)
  )
  ggplot(pca, aes(x = PC1, y = PC2, colour = condition, shape = condition)) +
    geom_point(size = 2.2) +
    geom_text(aes(label = label), size = 2.3, vjust = -0.9, show.legend = FALSE) +
    scale_colour_manual(values = colors, labels = ascii_text(conditions), name = NULL) +
    scale_shape_manual(values = shapes, labels = ascii_text(conditions), name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    labs(
      title = title, subtitle = subtitle,
      x = sprintf("PC1 (%.1f%% of variance)", 100 * pca$pc1_variance_fraction[[1]]),
      y = sprintf("PC2 (%.1f%% of variance)", 100 * pca$pc2_variance_fraction[[1]])
    ) +
    figure_theme()
}

# Draws a ggplot, or a figure stack with each plot in its own row.
draw_figure <- function(plot) {
  if (!inherits(plot, "figure_stack")) return(print(plot))
  grid::grid.newpage()
  layout <- grid::grid.layout(length(plot$plots), 1, heights = grid::unit(plot$heights, "null"))
  grid::pushViewport(grid::viewport(layout = layout))
  for (index in seq_along(plot$plots)) {
    print(
      plot$plots[[index]], newpage = FALSE,
      vp = grid::viewport(layout.pos.row = index, layout.pos.col = 1)
    )
  }
  grid::popViewport()
  invisible(plot)
}

# Writes STEM.pdf and STEM.png. The PDF has no creation date, modification
# date, or producer, so reruns reproduce it byte for byte. The PNG uses the
# platform's default bitmap device: cairo on Linux, quartz on macOS.
save_figure <- function(plot, stem, width, height) {
  # Graphics devices read "%d" in a file name as a page-number format.
  escape <- function(path) gsub("%", "%%", path, fixed = TRUE)
  title <- gsub("[^A-Za-z0-9._-]", "_", basename(stem))
  grDevices::pdf(
    escape(paste0(stem, ".pdf")), width = width, height = height,
    timestamp = FALSE, producer = FALSE, title = title
  )
  tryCatch(draw_figure(plot), finally = grDevices::dev.off())
  grDevices::png(
    escape(paste0(stem, ".png")), width = round(width * FIGURE_DPI),
    height = round(height * FIGURE_DPI), units = "px", res = FIGURE_DPI
  )
  tryCatch(draw_figure(plot), finally = grDevices::dev.off())
  invisible(stem)
}

# ---- Modes ------------------------------------------------------------------

run_figures_mode <- function(arguments) {
  require_args(arguments, c("samples", "statistics_dir", "atlas", "pau", "params", "output_dir"))
  params <- yaml::read_yaml(arguments$params)
  samples <- read_table(arguments$samples, c("sample_id", "condition", "control_condition"))
  comparisons <- comparisons_from_samples(samples)
  coordinates <- read_atlas_coordinates(arguments$atlas)
  dir.create(arguments$output_dir, recursive = TRUE, showWarnings = FALSE)
  pacs_tables <- list()
  genes_tables <- list()
  for (index in seq_len(nrow(comparisons))) {
    comparison <- comparisons[index, , drop = FALSE]
    stem <- comparison$stem
    table_path <- function(suffix) file.path(arguments$statistics_dir, paste0(stem, suffix))
    pacs <- read_pacs(table_path(".pacs.tsv.gz"))
    genes <- read_genes(table_path(".genes.tsv.gz"))
    pacs_tables[[stem]] <- pacs
    genes_tables[[stem]] <- genes
    exploratory <- any(pacs$exploratory_insufficient_replicates) ||
      any(genes$exploratory_insufficient_replicates)
    output <- file.path(arguments$output_dir, stem)
    save_figure(plot_volcano(pacs, comparison, params, exploratory), paste0(output, ".volcano"),
      6.5, 5)
    distal <- distal_usage_table(pacs, coordinates, genes)
    write_gzip_tsv(distal$table[, DISTAL_COLUMNS], paste0(output, ".distal_usage.tsv.gz"))
    save_figure(plot_distal_usage(distal, comparison, params), paste0(output, ".distal_usage"),
      5.5, 5.5)
    save_figure(plot_site_classes(site_class_counts(pacs), comparison),
      paste0(output, ".site_classes"), 6.5, 3.5)
  }
  counts <- event_count_table(comparisons, pacs_tables, genes_tables)
  patterns <- pattern_count_table(comparisons, genes_tables)
  save_figure(plot_event_counts(counts, patterns, comparisons),
    file.path(arguments$output_dir, "event_counts"), 7.5, 3.5 + 0.75 * nrow(comparisons))
  coverage <- coverage_data(comparisons, pacs_tables)
  dropped <- sum(vapply(pacs_tables, nrow, integer(1))) - nrow(coverage)
  rows <- ceiling(nrow(comparisons) / min(3, nrow(comparisons)))
  save_figure(
    plot_effect_vs_coverage(coverage, comparisons, params, dropped),
    file.path(arguments$output_dir, "effect_vs_coverage"), 7.5, 0.75 + 2.5 * rows
  )
  run_across_comparisons(arguments, params, samples, comparisons, pacs_tables, genes_tables)
}

# The figures and tables that relate comparisons to each other, and the PCA
# of samples.
run_across_comparisons <- function(arguments, params, samples, comparisons, pacs_tables,
                                   genes_tables) {
  output <- function(name) file.path(arguments$output_dir, name)
  labels <- grid_label_height(comparisons)
  grid <- pattern_grid_table(comparisons, genes_tables)
  write_gzip_tsv(grid, output("apa_patterns_by_comparison.tsv.gz"))
  save_figure(plot_pattern_grid(grid, comparisons), output("apa_pattern_grid"),
    3 + 0.75 * nrow(comparisons), 2.375 + labels + 0.125 * nrow(grid_genes(grid)))
  pairs <- concordance_pairs(comparisons)
  pair_data <- concordance_data(pairs, pacs_tables)
  pair_summary <- concordance_summary(pairs, pair_data)
  write_gzip_tsv(pair_summary, output("concordance.tsv.gz"))
  shown <- shown_pairs(pairs)
  save_figure(plot_concordance(pair_data, shown, pair_summary, nrow(pairs)),
    output("concordance"), 7.5, 1 + 2.5 * max(1, ceiling(nrow(shown) / 3)))
  side <- max(4.5, 2.5 + 0.5 * nrow(comparisons)) + labels
  save_figure(plot_concordance_matrix(pair_summary, comparisons), output("concordance_matrix"),
    side, side)
  pca <- pau_pca(read_pau(arguments$pau), samples, params$min_gene_total)
  pca_table <- if (is.null(pca)) {
    data.frame(sample_id = character(), condition = character(), PC1 = numeric(),
      PC2 = numeric(), pc1_variance_fraction = numeric(), pc2_variance_fraction = numeric())
  } else {
    as.data.frame(pca)
  }
  utils::write.table(pca_table, output("pau_pca.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "")
  save_figure(plot_pau_pca(pca, params), output("pau_pca"), 5.5, 5)
}

# Appends the ggplot2 version to a software-versions table, creating it with a
# header when it does not exist yet.
run_versions_mode <- function(arguments, version) {
  require_args(arguments, "output")
  rows <- data.frame(software = "ggplot2", version = version, stringsAsFactors = FALSE)
  exists <- file.exists(arguments$output)
  utils::write.table(
    rows, arguments$output, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = !exists, append = exists
  )
}

main <- function(argv) {
  arguments <- parse_args(argv)
  version <- load_figure_packages()
  switch(
    arguments$mode,
    versions = run_versions_mode(arguments, version),
    figures = run_figures_mode(arguments),
    stop("Unknown --mode: ", arguments$mode)
  )
}

if (sys.nframe() == 0L) main(commandArgs(trailingOnly = TRUE))
