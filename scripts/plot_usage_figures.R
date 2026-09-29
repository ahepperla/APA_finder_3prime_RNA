#!/usr/bin/env Rscript

# Treatment-control figures, drawn with ggplot2 from the tables that
# fit_usage_model.R writes. The figures show the calls in those tables; they
# never make calls of their own.
#
# Usage:
#   Rscript plot_usage_figures.R --mode figures --samples normalized_samples.tsv \
#     --statistics-dir DIR --atlas pacs.v1.metadata.tsv.gz \
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
# lists them; a gene is drawn with its first.
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

read_genes <- function(path) {
  read_table(
    path, GENES_REQUIRED, "gene_fdr",
    c("dominant_switch", "exploratory_insufficient_replicates")
  )
}

# Only the PAC IDs and coordinates, skipping the atlas's sequence columns.
read_atlas_coordinates <- function(path) {
  header <- names(utils::read.delim(path, nrows = 0, quote = "", check.names = FALSE))
  missing <- setdiff(c("pac_id", "coordinate"), header)
  if (length(missing)) {
    stop(basename(path), " lacks columns: ", paste(missing, collapse = ", "), ".")
  }
  classes <- rep("NULL", length(header))
  classes[header == "pac_id"] <- "character"
  classes[header == "coordinate"] <- "numeric"
  utils::read.delim(
    path, colClasses = classes, na.strings = character(), quote = "", check.names = FALSE
  )
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

# Genes per APA pattern in each comparison, with zeros; a gene with two
# patterns counts in both.
pattern_count_table <- function(comparisons, genes_tables) {
  classes <- setdiff(names(APA_CLASSES), "none")
  rows <- lapply(comparisons$stem, function(stem) {
    genes <- genes_tables[[stem]]
    genes <- genes[!duplicated(genes$gene_id), , drop = FALSE]
    parts <- strsplit(genes$apa_pattern, ";", fixed = TRUE)
    counts <- vapply(classes, function(class) {
      sum(vapply(parts, function(part) class %in% part, logical(1)))
    }, integer(1))
    data.frame(comparison = stem, pattern = classes, count = unname(counts))
  })
  counts <- do.call(rbind, rows)
  data.frame(
    comparison = factor(counts$comparison, levels = comparisons$stem),
    pattern = factor(counts$pattern, levels = classes),
    count = counts$count
  )
}

# A gene's first APA pattern, which colors it.
first_pattern <- function(apa_pattern) {
  first <- sub(";.*$", "", apa_pattern)
  first[is.na(first) | !first %in% names(APA_CLASSES)] <- "none"
  factor(first, levels = names(APA_CLASSES))
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
  shown$call <- factor(
    ifelse(shown$direction %in% c("distal_up", "distal_down"), "Confirmed", "None or candidate"),
    levels = c("Confirmed", "None or candidate")
  )
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
      values = c("Confirmed" = 16, "None or candidate" = 1), name = "Distal PAC call",
      guide = guide_legend(order = 2, nrow = 1)
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
  totals <- stats::aggregate(count ~ comparison, data = patterns, FUN = sum)
  ggplot(patterns, aes(x = count, y = comparison, fill = pattern)) +
    geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
    geom_text(data = totals, aes(x = count, y = comparison, label = count),
      inherit.aes = FALSE, hjust = -0.3, size = 2.6) +
    scale_fill_manual(values = APA_COLORS, labels = APA_CLASSES, name = NULL, drop = FALSE) +
    scale_y_discrete(limits = rev(comparisons$stem), labels = titles) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(
      title = "APA patterns per comparison",
      subtitle = "Genes per pattern, from the genes tables; a gene with two patterns counts in both",
      x = "Genes", y = NULL
    ) +
    figure_theme()
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
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    guides(fill = guide_legend(nrow = 3)) +
    labs(
      title = "Events per comparison",
      subtitle = "PAC events count PACs; gene events count genes that pass the gene-level screen",
      x = "Count", y = NULL
    ) +
    figure_theme()
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
  require_args(arguments, c("samples", "statistics_dir", "atlas", "params", "output_dir"))
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
