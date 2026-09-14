#!/usr/bin/env Rscript
# Light-chain pairing: P(J|V) and P(V|J) for IGK (A, B) and IGL (C, D), each over the
# Spearman heatmap of its rank orders. V is seriated then put back in genomic order; J is
# seriated only. Tables are computed only when missing from results/figures/source_data;
# the figure is always drawn from them.

source("R/00_setup.R")
source("R/lib/dj_pairing.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(cowplot) })

tables <- file.path(OUT$source, c("light_pairing_probabilities.csv.gz", "light_pairing_correlations.csv.gz", "light_pairing_axis_levels.csv"))
if (!all(file.exists(tables))) {
library(seriation)
locus_gene_order <- function(locus, segment) {
  b <- fread(need(IN$gene_bed), col.names = c("chrom", "start", "end", "gene"))
  unique(b[grepl(sprintf("^%s%s", locus, segment), gene), gsub("^IG[HKL]", "", gene)])
}

prob_all <- list(); corr_all <- list(); lev_all <- list()
for (locus in c("IGK", "IGL")) {
  rep_dt <- fread(need(file.path(OUT$repertoire, sprintf("gg_repertoire_data_%s_genotype_corrected.csv.gz", locus))))
  res <- compute_bias_matrices(rep_dt, subject_col = "vdjbase_subject", gene_a_col = "v_gene_iuis", gene_b_col = "j_gene_iuis", top_n_subjects = uniqueN(rep_dt$vdjbase_subject))
  v_seriated <- rownames(res$spearman_a_given_b)[get_order(seriate(dist(res$spearman_a_given_b)))]
  j_seriated <- rownames(res$spearman_b_given_a)[get_order(seriate(dist(res$spearman_b_given_a)))]
  v_short <- unique(order_labels_by_reference(v_seriated, locus_gene_order(locus, "V"))$short_label)
  j_short <- unique(short_axis_label(j_seriated))
  prob <- collapse_probability_for_plot(res$subject_probability_dt, v_short, j_short)
  ord_jv <- build_probability_order_dt(prob, "p_b_given_a", "gene_b_short", "gene_a_short")
  ord_vj <- build_probability_order_dt(prob, "p_a_given_b", "gene_a_short", "gene_b_short")
  prob <- merge(prob, ord_jv[, .(gene_b_short, gene_a_short, order_j_given_v = custom_order)], by = c("gene_b_short", "gene_a_short"))
  prob <- merge(prob, ord_vj[, .(gene_a_short, gene_b_short, order_v_given_j = custom_order)], by = c("gene_a_short", "gene_b_short"))
  prob[, locus := locus]
  prob_all[[locus]] <- prob
  corr_all[[locus]] <- rbind(
    collapse_correlation_for_plot(build_correlation_long_dt(res$spearman_b_given_a), j_short, j_short)[, `:=`(direction = "P(J|V)", locus = locus)],
    collapse_correlation_for_plot(build_correlation_long_dt(res$spearman_a_given_b), v_short, v_short)[, `:=`(direction = "P(V|J)", locus = locus)])
  lev_all[[locus]] <- rbind(data.table(locus = locus, axis = paste0(locus, "V"), short_label = v_short, plot_order = seq_along(v_short)),
                            data.table(locus = locus, axis = paste0(locus, "J"), short_label = j_short, plot_order = seq_along(j_short)))
  cat(sprintf("%s  %d subjects  %d V ASCs  %d J ASCs\n", locus, uniqueN(rep_dt$vdjbase_subject), length(v_short), length(j_short)))
}
fwrite(rbindlist(prob_all, use.names = TRUE), tables[1])
fwrite(rbindlist(corr_all, use.names = TRUE), tables[2])
fwrite(rbindlist(lev_all, use.names = TRUE), tables[3])
}

# ---- draw ----
prob <- fread(tables[1]); corr <- fread(tables[2]); lev <- fread(tables[3])
panels <- list()
for (lc in c("IGK", "IGL")) {
  v_short <- lev[locus == lc & axis == paste0(lc, "V")][order(plot_order), short_label]
  j_short <- lev[locus == lc & axis == paste0(lc, "J")][order(plot_order), short_label]
  p <- prob[locus == lc]
  p[, gene_a_short := factor(as.character(gene_a_short), levels = v_short)]
  p[, gene_b_short := factor(as.character(gene_b_short), levels = j_short)]
  panels[[paste0(lc, "_jv")]] <- build_probability_boxplot(p[, .(subject, gene_a_short, gene_b_short, probability_value = p_b_given_a, custom_order = order_j_given_v)],
                                                           "gene_b_short", "probability_value", "gene_a_short", paste0(lc, "J ASC"), "P(J|V)", paste0(lc, "V ASC"), make_named_palette(v_short))
  panels[[paste0(lc, "_vj")]] <- build_probability_boxplot(p[, .(subject, gene_a_short, gene_b_short, probability_value = p_a_given_b, custom_order = order_v_given_j)],
                                                           "gene_a_short", "probability_value", "gene_b_short", paste0(lc, "V ASC"), "P(V|J)", paste0(lc, "J ASC"), make_named_palette(j_short))
  for (dir in c("P(J|V)", "P(V|J)")) {
    lv <- if (dir == "P(J|V)") j_short else v_short
    d <- corr[locus == lc & direction == dir]
    d[, `:=`(x_short = factor(as.character(x_short), levels = lv), y_short = factor(as.character(y_short), levels = lv))]
    # The V matrices carry 22-23 labels and the x axis already names every row, so the y labels are dropped.
    panels[[paste0(lc, if (dir == "P(J|V)") "_tjv" else "_tvj")]] <- build_correlation_heatmap(d, paste0(lc, if (dir == "P(J|V)") "J" else "V", " ASC"), NULL, TRUE, dir == "P(J|V)")
  }
}
no_x <- theme(axis.text.x = element_blank())
stack <- function(box, tile) (box + xlab(NULL) + no_x) / (tile + guides(fill = "none")) + plot_layout(heights = c(2, 1))
fig <- wrap_elements(stack(panels$IGK_jv + guides(color = guide_legend(nrow = 2)), panels$IGK_tjv + theme(axis.text.x = element_text(angle = 0, vjust = 0.5)))) +
  wrap_elements(stack(panels$IGK_vj, panels$IGK_tvj)) +
  wrap_elements(stack(panels$IGL_jv + guides(color = guide_legend(nrow = 2)), panels$IGL_tjv + theme(axis.text.x = element_text(angle = 0, vjust = 0.5)))) +
  wrap_elements(stack(panels$IGL_vj, panels$IGL_tvj)) +
  plot_layout(design = "A\nB\nC\nD", heights = c(1, 1.35, 1, 1.35)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) & theme(plot.tag = element_text(size = 40, face = "bold"))
ggsave(file.path(OUT$figures, "supp_light_pairing.pdf"), fig, width = 40, height = 46, device = grDevices::cairo_pdf, limitsize = FALSE)
