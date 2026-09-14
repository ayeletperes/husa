# Conditional gene-pairing tables and panels shared by figure_dj_pairing, supp7 and
# supp_light_pairing: P(B|A) and P(A|B) per subject with the Spearman rank-order
# correlations of each axis, the label collapsing for slash-joined ASCs, and the boxplot /
# heatmap builders. Lifted from the manuscript tree unchanged.

normalize_group_label <- function(label) {
  parts <- trimws(unlist(strsplit(as.character(label), "/", fixed = TRUE)))
  parts <- parts[nzchar(parts)]
  if (!length(parts)) {
    return(NA_character_)
  }
  first_part <- parts[[1]]
  if (length(parts) == 1L) {
    return(first_part)
  }
  first_prefix <- sub("-.*$", "", first_part)
  remaining <- vapply(parts[-1], function(x) {
    if (grepl("^[A-Z]", x)) x else paste0(first_prefix, "-", x)
  }, character(1))
  c(first_part, remaining)[[1]]
}

short_axis_label <- function(x) {
  sub("[/].*$", "*", as.character(x))
}

read_gene_order <- function(gene_bed_path, gene_type, remove_prefix = FALSE) {
  gene_bed_dt <- data.table::fread(
    file = gene_bed_path,
    col.names = c("chrom", "start", "end", "gene")
  )

  d_gene_order <- gene_bed_dt[grepl(paste0("^", gene_type), gene), gene]
  if (remove_prefix) {
    d_gene_order <- gsub("^IG[HKL]", "", d_gene_order)
  }
  unique(d_gene_order)
}

order_labels_by_reference <- function(labels, reference_order) {
  label_dt <- data.table::data.table(
    label = as.character(labels),
    short_label = short_axis_label(labels),
    first_component = vapply(labels, normalize_group_label, character(1))
  )

  label_dt[, order_index := match(first_component, reference_order)]

  missing_n <- sum(is.na(label_dt$order_index))
  if (missing_n > 0L) {
    max_index <- suppressWarnings(max(label_dt$order_index, na.rm = TRUE))
    if (!is.finite(max_index)) max_index <- 0L
    label_dt[is.na(order_index), order_index := max_index + seq_len(.N)]
  }

  data.table::setorder(label_dt, order_index, short_label, label)
  label_dt
}

build_bias_input <- function(repertoire_dt, subject_col, gene_a_col, gene_b_col) {
  bias_dt <- data.table::copy(repertoire_dt[, .(
    subject = get(subject_col),
    gene_a = get(gene_a_col),
    gene_b = get(gene_b_col)
  )])

  bias_dt <- bias_dt[
    !is.na(subject) &
      !is.na(gene_a) &
      !is.na(gene_b) &
      nzchar(gene_a) &
      nzchar(gene_b)
  ]

  bias_dt
}

compute_bias_matrices <- function(
  repertoire_dt,
  subject_col,
  gene_a_col,
  gene_b_col,
  min_count_per_a = 10L,
  min_count_per_b = 10L,
  min_subjects_per_gene = 100L,
  top_n_subjects = NULL,
  combo_depth = NULL,
  require_complete_a = TRUE,
  eps = 1e-6
) {
  filtered_repertoire_dt <- build_bias_input(repertoire_dt, subject_col, gene_a_col, gene_b_col)

  filtered_repertoire_dt[, gene_a_count := .N, by = .(subject, gene_a)]
  filtered_repertoire_dt[, gene_b_count := .N, by = .(subject, gene_b)]

  gene_a_subject_dt <- filtered_repertoire_dt[
    ,
    .(subjects_above_min = data.table::uniqueN(subject[gene_a_count >= min_count_per_a])),
    by = gene_a
  ]

  kept_gene_a <- gene_a_subject_dt[subjects_above_min >= min_subjects_per_gene, gene_a]
  filtered_repertoire_dt <- filtered_repertoire_dt[gene_a %in% kept_gene_a]

  gene_b_subject_dt <- filtered_repertoire_dt[
    ,
    .(subjects_above_min = data.table::uniqueN(subject[gene_b_count >= min_count_per_b])),
    by = gene_b
  ]

  kept_gene_b <- gene_b_subject_dt[subjects_above_min >= min_subjects_per_gene, gene_b]
  filtered_repertoire_dt <- filtered_repertoire_dt[gene_b %in% kept_gene_b]

  if (min_count_per_a > 1L) {
    filtered_repertoire_dt <- filtered_repertoire_dt[gene_a_count >= min_count_per_a]
  }

  if (min_count_per_b > 1L) {
    filtered_repertoire_dt <- filtered_repertoire_dt[gene_b_count >= min_count_per_b]
  }

  if (!is.null(top_n_subjects)) {
    subject_coverage_dt <- filtered_repertoire_dt[
      ,
      .(has_all_gene_a = all(kept_gene_a %in% unique(gene_a))),
      by = subject
    ]

    candidate_subjects <- if (isTRUE(require_complete_a)) {
      subject_coverage_dt[has_all_gene_a == TRUE, subject]
    } else {
      subject_coverage_dt$subject
    }

    if (!length(candidate_subjects)) {
      candidate_subjects <- subject_coverage_dt$subject
    }

    subject_depth_rank_dt <- filtered_repertoire_dt[
      subject %in% candidate_subjects,
      .N,
      by = subject
    ][order(-N)]

    kept_subjects <- head(subject_depth_rank_dt$subject, top_n_subjects)
    filtered_repertoire_dt <- filtered_repertoire_dt[subject %in% kept_subjects]
  }

  pair_count_dt <- filtered_repertoire_dt[, .(pair_count = .N), by = .(subject, gene_a, gene_b)]

  if (!is.null(combo_depth)) {
    pair_count_dt <- pair_count_dt[pair_count > combo_depth]
  }

  subject_depth_dt <- pair_count_dt[, .(subject_depth = sum(pair_count)), by = subject]
  gene_a_total_dt <- pair_count_dt[, .(gene_a_total = sum(pair_count)), by = .(subject, gene_a)]
  gene_b_total_dt <- pair_count_dt[, .(gene_b_total = sum(pair_count)), by = .(subject, gene_b)]

  subject_probability_dt <- merge(pair_count_dt, subject_depth_dt, by = "subject")
  subject_probability_dt <- merge(subject_probability_dt, gene_a_total_dt, by = c("subject", "gene_a"))
  subject_probability_dt <- merge(subject_probability_dt, gene_b_total_dt, by = c("subject", "gene_b"))

  subject_probability_dt[, p_pair := pair_count / subject_depth]
  subject_probability_dt[, p_a := gene_a_total / subject_depth]
  subject_probability_dt[, p_b := gene_b_total / subject_depth]
  subject_probability_dt[, p_b_given_a := pair_count / pmax(gene_a_total, 1L)]
  subject_probability_dt[, p_a_given_b := pair_count / pmax(gene_b_total, 1L)]
  subject_probability_dt[, expected_pair := p_a * p_b]
  subject_probability_dt[, residual_pair := p_pair - expected_pair]
  subject_probability_dt[, log2_enrichment := log2((p_pair + eps) / (expected_pair + eps))]

  build_median_matrix <- function(value_col) {
    wide_dt <- data.table::dcast(
      subject_probability_dt[
        ,
        .(median_value = median(get(value_col), na.rm = TRUE)),
        by = .(gene_a, gene_b)
      ],
      gene_a ~ gene_b,
      value.var = "median_value",
      fill = NA_real_
    )

    row_ids <- wide_dt$gene_a
    wide_dt$gene_a <- NULL

    matrix_value <- as.matrix(wide_dt)
    rownames(matrix_value) <- row_ids
    matrix_value
  }

  p_b_given_a_matrix <- build_median_matrix("p_b_given_a")
  p_a_given_b_matrix <- build_median_matrix("p_a_given_b")
  observed_matrix <- build_median_matrix("p_pair")
  expected_matrix <- build_median_matrix("expected_pair")
  residual_matrix <- build_median_matrix("residual_pair")
  log2_enrichment_matrix <- build_median_matrix("log2_enrichment")

  spearman_b_given_a <- suppressWarnings(
    stats::cor(p_b_given_a_matrix, method = "spearman", use = "pairwise.complete.obs")
  )

  spearman_a_given_b <- suppressWarnings(
    stats::cor(t(p_a_given_b_matrix), method = "spearman", use = "pairwise.complete.obs")
  )

  gene_a_probability_check <- subject_probability_dt[
    ,
    .(sum_probability = sum(p_b_given_a)),
    by = .(subject, gene_a)
  ]

  gene_b_probability_check <- subject_probability_dt[
    ,
    .(sum_probability = sum(p_a_given_b)),
    by = .(subject, gene_b)
  ]

  list(
    subject_probability_dt = subject_probability_dt,
    median_matrices = list(
      observed_matrix = observed_matrix,
      expected_matrix = expected_matrix,
      residual_matrix = residual_matrix,
      log2_enrichment_matrix = log2_enrichment_matrix,
      p_b_given_a_matrix = p_b_given_a_matrix,
      p_a_given_b_matrix = p_a_given_b_matrix
    ),
    spearman_b_given_a = spearman_b_given_a,
    spearman_a_given_b = spearman_a_given_b,
    probability_checks = list(
      gene_a_probability_check = gene_a_probability_check,
      gene_b_probability_check = gene_b_probability_check,
      max_gene_a_deviation = max(abs(gene_a_probability_check$sum_probability - 1), na.rm = TRUE),
      max_gene_b_deviation = max(abs(gene_b_probability_check$sum_probability - 1), na.rm = TRUE)
    )
  )
}

build_correlation_long_dt <- function(correlation_matrix) {
  correlation_long_dt <- data.table::as.data.table(reshape2::melt(correlation_matrix))
  data.table::setnames(
    correlation_long_dt,
    c("Var1", "Var2", "value"),
    c("x_gene", "y_gene", "spearman")
  )

  correlation_long_dt[, x_short := short_axis_label(x_gene)]
  correlation_long_dt[, y_short := short_axis_label(y_gene)]
  correlation_long_dt[]
}

collapse_correlation_for_plot <- function(correlation_long_dt, x_levels, y_levels, set_diagonal_to_one = TRUE, remove_prefix = TRUE) {
  plot_dt <- data.table::copy(correlation_long_dt)

  plot_dt <- plot_dt[
    ,
    .(
      spearman = if (all(is.na(spearman))) {
        NA_real_
      } else {
        mean(spearman, na.rm = TRUE)
      }
    ),
    by = .(
      x_short = as.character(x_short),
      y_short = as.character(y_short)
    )
  ]

  if (isTRUE(set_diagonal_to_one)) {
    plot_dt[!is.na(x_short) & x_short == y_short, spearman := 1]
  }

  if (remove_prefix) {
    plot_dt[, x_short := gsub("^IG[HKL]", "", x_short)]
    plot_dt[, y_short := gsub("^IG[HKL]", "", y_short)]

    x_levels <- gsub("^IG[HKL]", "", x_levels)
    y_levels <- gsub("^IG[HKL]", "", y_levels)
  }

  plot_dt[, x_short := factor(x_short, levels = x_levels)]
  plot_dt[, y_short := factor(y_short, levels = y_levels)]

  plot_dt[!is.na(x_short) & !is.na(y_short)]
}

collapse_probability_for_plot <- function(subject_probability_dt, d_levels, j_levels) {
  plot_dt <- data.table::copy(subject_probability_dt)

  plot_dt[, gene_a_short := short_axis_label(gene_a)]
  plot_dt[, gene_b_short := short_axis_label(gene_b)]

  # Important: collapse counts first, then recompute conditional probabilities.
  # Otherwise two full labels that become the same short label are still treated as separate x positions.
  plot_dt <- plot_dt[
    ,
    .(pair_count = sum(pair_count, na.rm = TRUE)),
    by = .(subject, gene_a_short, gene_b_short)
  ]

  plot_dt[, subject_depth := sum(pair_count), by = subject]
  plot_dt[, gene_a_total := sum(pair_count), by = .(subject, gene_a_short)]
  plot_dt[, gene_b_total := sum(pair_count), by = .(subject, gene_b_short)]

  plot_dt[, p_pair := pair_count / subject_depth]
  plot_dt[, p_a := gene_a_total / subject_depth]
  plot_dt[, p_b := gene_b_total / subject_depth]
  plot_dt[, p_b_given_a := pair_count / pmax(gene_a_total, 1L)]
  plot_dt[, p_a_given_b := pair_count / pmax(gene_b_total, 1L)]

  plot_dt[, gene_a_short := factor(gene_a_short, levels = d_levels)]
  plot_dt[, gene_b_short := factor(gene_b_short, levels = j_levels)]

  plot_dt[!is.na(gene_a_short) & !is.na(gene_b_short)]
}

build_probability_order_dt <- function(subject_probability_dt, value_col, x_col, color_col) {
  order_dt <- subject_probability_dt[
    ,
    .(median_probability = median(get(value_col), na.rm = TRUE)),
    by = c(x_col, color_col)
  ]

  data.table::setorderv(order_dt, c(x_col, "median_probability"), c(1L, -1L))
  order_dt[, custom_order := seq_len(.N), by = x_col]
  order_dt
}


# ---- panels ----

make_named_palette <- function(levels) {
  levels <- as.character(levels)
  pal <- scales::hue_pal()(length(levels))
  stats::setNames(pal, levels)
}

build_probability_boxplot <- function(
  probability_long_dt,
  x_col,
  y_col,
  color_col,
  x_label,
  y_label,
  legend_title,
  color_values = NULL,
  base_size = 39
) {
  ggplot2::ggplot(
    probability_long_dt,
    ggplot2::aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      color = .data[[color_col]],
      group = interaction(.data[[x_col]], custom_order)
    )
  ) +
    ggplot2::geom_boxplot(
      #position = ggplot2::position_dodge2(width = 0.8, preserve = "single"),
      outlier.shape = NA,
      outliers = FALSE
    ) +
    ggplot2::scale_x_discrete(
      drop = FALSE,
      expand = ggplot2::expansion(add = 0)
    ) +
    {
      if (is.null(color_values)) {
        ggplot2::scale_color_discrete(drop = FALSE)
      } else {
        ggplot2::scale_color_manual(values = color_values, drop = FALSE)
      }
    } +
    ggpubr::theme_pubclean(base_size = base_size) +
    ggplot2::labs(y = y_label, x = x_label, color = legend_title)
}

build_correlation_heatmap <- function(
  correlation_long_dt,
  x_label,
  y_label,
  show_x_text = TRUE,
  show_y_text = TRUE,
  x_position_top = FALSE,
  base_size = 34
) {
  x_position <- if (x_position_top) "top" else "bottom"

  ggplot2::ggplot(
    correlation_long_dt,
    ggplot2::aes(x = x_short, y = y_short, fill = spearman)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_x_discrete(
      position = x_position,
      drop = FALSE,
      expand = ggplot2::expansion(add = 0)
    ) +
    ggplot2::scale_y_discrete(
      drop = FALSE,
      expand = ggplot2::expansion(add = 0)
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      oob = scales::squish,
      name = "Spearman"
    ) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::labs(x = x_label, y = y_label) +
    ggplot2::theme(
      axis.text.x = if (show_x_text) {
        ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      } else {
        ggplot2::element_blank()
      },
      axis.text.x.top = if (x_position_top && show_x_text) {
        ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0.5)
      } else {
        ggplot2::element_blank()
      },
      axis.text.y = if (show_y_text) {
        ggplot2::element_text(hjust = 1)
      } else {
        ggplot2::element_blank()
      },
      axis.title.x = ggplot2::element_text(),
      axis.title.y = ggplot2::element_text()
    )
}
