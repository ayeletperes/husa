#!/usr/bin/env Rscript
# Supplementary figure 7: IGHV/IGHD conditional usage, the V-D counterpart of the D-J pairing
# figure. A: P(D|V) per subject by IGHD ASC, coloured by IGHV ASC, over the IGHD rank-order
# Spearman heatmap. B: P(V|D), the mirror. The IGHD axis carries the same grouping and genomic
# order as the D-J figure; IGHV is seriated.

source("R/00_setup.R")
source("R/lib/dj_pairing.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(cowplot); library(seriation) })

rep_dt <- fread(need(file.path(OUT$repertoire, "gg_repertoire_data_IGH_genotype_corrected.csv.gz")))
res <- compute_bias_matrices(rep_dt, subject_col = "vdjbase_subject", gene_a_col = "d_gene_iuis", gene_b_col = "v_gene_iuis", top_n_subjects = uniqueN(rep_dt$vdjbase_subject))
d_seriated <- rownames(res$spearman_a_given_b)[get_order(seriate(dist(res$spearman_a_given_b)))]
v_seriated <- rownames(res$spearman_b_given_a)[get_order(seriate(dist(res$spearman_b_given_a)))]
d_short <- unique(order_labels_by_reference(d_seriated, read_gene_order(IN$gene_bed, "IGHD", remove_prefix = FALSE))$short_label)
v_short <- unique(short_axis_label(v_seriated))
prob <- collapse_probability_for_plot(res$subject_probability_dt, d_short, v_short)
ord_dv <- build_probability_order_dt(prob, "p_a_given_b", "gene_a_short", "gene_b_short")
ord_vd <- build_probability_order_dt(prob, "p_b_given_a", "gene_b_short", "gene_a_short")
prob <- merge(prob, ord_dv[, .(gene_a_short, gene_b_short, order_d_given_v = custom_order)], by = c("gene_a_short", "gene_b_short"))
prob <- merge(prob, ord_vd[, .(gene_b_short, gene_a_short, order_v_given_d = custom_order)], by = c("gene_b_short", "gene_a_short"))
corr_dv <- collapse_correlation_for_plot(build_correlation_long_dt(res$spearman_a_given_b), d_short, d_short)[, direction := "P(D|V)"]
corr_vd <- collapse_correlation_for_plot(build_correlation_long_dt(res$spearman_b_given_a), v_short, v_short)[, direction := "P(V|D)"]
fwrite(prob, file.path(OUT$source, "supp7_probabilities.csv.gz"))
fwrite(rbind(corr_dv, corr_vd), file.path(OUT$source, "supp7_correlations.csv.gz"))
fwrite(rbind(data.table(axis = "IGHD", short_label = d_short, plot_order = seq_along(d_short)),
             data.table(axis = "IGHV", short_label = v_short, plot_order = seq_along(v_short))), file.path(OUT$source, "supp7_axis_levels.csv"))
cat(sprintf("supp7: %d subjects, %d IGHD ASCs, %d IGHV ASCs\n", uniqueN(prob$subject), length(d_short), length(v_short)))

# ---- draw ----
d_short <- gsub("^IGH", "", d_short); v_short <- gsub("^IGH", "", v_short)
prob[, gene_a_short := factor(gsub("^IGH", "", as.character(gene_a_short)), levels = d_short)]
prob[, gene_b_short := factor(gsub("^IGH", "", as.character(gene_b_short)), levels = v_short)]
d_cols <- make_named_palette(d_short); v_cols <- make_named_palette(v_short)
no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# This figure's boxplot dodges the boxes per group, unlike the D-J figure's.
dodged_boxplot <- function(dt, x_col, y_col, color_col, x_label, y_label, legend_title, base_size = 39) {
  ggplot(dt, aes(x = get(x_col), y = get(y_col), color = get(color_col), group = interaction(get(x_col), custom_order))) +
    geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), outlier.shape = NA, outliers = FALSE) +
    ggpubr::theme_pubclean(base_size = base_size) + labs(y = y_label, x = x_label, color = legend_title)
}
# The boxplot ordering follows the merge-key sort of the tables it is drawn from.
dv_dt <- prob[, .(subject, gene_a_short, gene_b_short, probability_value = p_a_given_b, custom_order = order_d_given_v)]
setorder(dv_dt, gene_a_short, gene_b_short)
vd_dt <- prob[, .(subject, gene_a_short, gene_b_short, probability_value = p_b_given_a, custom_order = order_v_given_d)]
setorder(vd_dt, gene_b_short, gene_a_short)
p_d_given_v <- dodged_boxplot(dv_dt, "gene_a_short", "probability_value", "gene_b_short", "IGHD ASC", "P(D|V)", "IGHV ASC") +
  scale_x_discrete(drop = FALSE, expand = expansion(add = 0)) + scale_color_manual(values = v_cols, breaks = v_short, drop = FALSE) +
  no_x + theme(axis.title.x = element_text(), legend.position = "top")
p_v_given_d <- dodged_boxplot(vd_dt, "gene_b_short", "probability_value", "gene_a_short", "IGHV ASC", "P(V|D)", "IGHD ASC") +
  scale_x_discrete(drop = FALSE, expand = expansion(add = 0)) + scale_color_manual(values = d_cols, breaks = d_short, drop = FALSE) +
  no_x + theme(legend.position = "top")
tile_of <- function(d, levs, xlab) {
  d[, `:=`(x_short = factor(gsub("^IGH", "", as.character(x_short)), levels = levs), y_short = factor(gsub("^IGH", "", as.character(y_short)), levels = levs))]
  build_correlation_heatmap(d, xlab, NULL, show_x_text = TRUE, show_y_text = TRUE)
}
p_tile_dv <- tile_of(corr_dv, d_short, "IGHD ASC")
p_tile_vd <- tile_of(corr_vd, v_short, "IGHV ASC")

legend_patch <- wrap_elements(full = plot_grid(
  plot_grid(get_legend(p_d_given_v + theme(legend.position = "top", legend.direction = "horizontal") + guides(color = guide_legend(nrow = 3, byrow = TRUE))),
            get_legend(p_v_given_d + theme(legend.position = "top", legend.direction = "horizontal") + guides(color = guide_legend(nrow = 2, byrow = TRUE))),
            nrow = 1, rel_widths = c(1.4, 1)),
  get_legend(p_tile_vd + theme(legend.position = "top", legend.direction = "horizontal", legend.key.width = grid::unit(0.5, "inches"))),
  ncol = 1, rel_heights = c(1.4, 0.55)))
tile_theme <- theme(axis.text.y = element_text(color = NA), axis.ticks.y = element_blank(), axis.title.x = element_text(), axis.title.y = element_blank())
fig <- (legend_patch +
          (p_d_given_v + guides(color = "none") + xlab(NULL) + labs(tag = "A") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())) +
          (p_tile_dv + guides(fill = "none") + xlab("IGHD ASC") + tile_theme) +
          (p_v_given_d + guides(color = "none") + xlab(NULL) + labs(tag = "B") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())) +
          (p_tile_vd + guides(fill = "none") + xlab("IGHV ASC") + tile_theme) +
          plot_layout(design = "A\nB\nC\nD\nE", heights = c(1.1, 2, 1.3, 2, 1.3))) &
  theme(plot.tag = element_text(size = 58, face = "bold"), plot.tag.position = c(0.01, 0.98))
ggsave(file.path(OUT$figures, "supp7.pdf"), fig, width = 70, height = 32, device = grDevices::cairo_pdf, limitsize = FALSE)
