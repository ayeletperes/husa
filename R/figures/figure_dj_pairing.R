#!/usr/bin/env Rscript
# IGHD/IGHJ pairing figure. A: P(J|D) per subject by IGHJ ASC, coloured by IGHD ASC, over
# the Spearman heatmap of the D rank orders. B: P(D|J), the mirror. C: the same conditional
# split by genotype at the strongest pairing QTL, with its marginals.

source("R/00_setup.R")
source("R/lib/dj_pairing.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(cowplot); library(seriation) })

# igh_37391 is the lowest p value in the enrichment MANOVA (J4 anchored); igh_38679 ties with it in
# perfect linkage disequilibrium.
variant <- "igh_37391"
qtl <- file.path(OUT$qtl, "source_data")

# ---- A and B: pairing bias ----
rep_dt <- fread(need(file.path(OUT$repertoire, "gg_repertoire_data_IGH_genotype_corrected.csv.gz")))
res <- compute_bias_matrices(rep_dt, subject_col = "vdjbase_subject", gene_a_col = "d_gene_iuis", gene_b_col = "j_gene_iuis",
                             top_n_subjects = uniqueN(rep_dt$vdjbase_subject))
# D is seriated then put back into genomic order; J is seriated only.
d_seriated <- rownames(res$spearman_a_given_b)[get_order(seriate(dist(res$spearman_a_given_b)))]
j_seriated <- rownames(res$spearman_b_given_a)[get_order(seriate(dist(res$spearman_b_given_a)))]
d_short <- unique(order_labels_by_reference(d_seriated, read_gene_order(IN$gene_bed, "IGHD", remove_prefix = FALSE))$short_label)
j_short <- unique(order_labels_by_reference(j_seriated, read_gene_order(IN$gene_bed, "IGHJ", remove_prefix = TRUE))$short_label)
# A slash-joined ASC ("IGHD2-2/IGHD2-15") collapses to "IGHD2-2*" on both the probabilities and the correlations.
prob <- collapse_probability_for_plot(res$subject_probability_dt, d_short, j_short)
ord_jd <- build_probability_order_dt(prob, "p_b_given_a", "gene_b_short", "gene_a_short")
ord_dj <- build_probability_order_dt(prob, "p_a_given_b", "gene_a_short", "gene_b_short")
prob <- merge(prob, ord_jd[, .(gene_b_short, gene_a_short, order_j_given_d = custom_order)], by = c("gene_b_short", "gene_a_short"))
prob <- merge(prob, ord_dj[, .(gene_a_short, gene_b_short, order_d_given_j = custom_order)], by = c("gene_a_short", "gene_b_short"))
corr_jd <- collapse_correlation_for_plot(build_correlation_long_dt(res$spearman_b_given_a), j_short, j_short)[, direction := "P(J|D)"]
corr_dj <- collapse_correlation_for_plot(build_correlation_long_dt(res$spearman_a_given_b), d_short, d_short)[, direction := "P(D|J)"]

# ---- C: the variant. Per-subject tables are joined from the complete QTL products; the
# cell fits are recomputed rather than read, so any filtered variant could be shown. ----
usage <- fread(file.path(qtl, "usage_associations_IGH.tsv.gz"), select = c("variant", "asc", "segment", "p_value", "significant"))
setnames(usage, c("p_value", "significant"), c("usage_p_value", "usage_significant"))
pairing <- fread(file.path(qtl, "pairing_associations.tsv.gz"), select = c("variant", "j_gene", "p_value", "significant"))
setnames(pairing, c("p_value", "significant"), c("pairing_by_j_p_value", "pairing_by_j_significant"))
enrichment <- fread(file.path(qtl, "dj_enrichment.tsv.gz"))
v <- variant  # not `variant`: inside `[` a bare name matching a column wins
dose <- fread(file.path(qtl, "dosage_long.tsv.gz"))[variant == v]
dose[, dosage := as.numeric(genotype)]
join_usage <- function(x, gene_col, seg) merge(x, usage[segment == seg, !"segment"], by.x = c("variant", gene_col), by.y = c("variant", "asc"), all.x = TRUE)
pj <- join_usage(merge(dose, unique(enrichment[, .(subject, j_gene, p_j)]), by = "subject", allow.cartesian = TRUE), "j_gene", "J")
pd <- join_usage(merge(dose, unique(enrichment[, .(subject, d_gene, p_d)]), by = "subject", allow.cartesian = TRUE), "d_gene", "D")
pjd <- merge(dose, enrichment[, .(subject, d_gene, j_gene, p_j_given_d, enrichment)], by = "subject", allow.cartesian = TRUE)
cells <- pjd[, {
  ok <- is.finite(enrichment) & is.finite(dosage)
  cf <- if (sum(ok) >= 10L && uniqueN(dosage[ok]) > 1L) summary(lm(enrichment[ok] ~ dosage[ok]))$coefficients
  .(cell_p_value = if (is.null(cf)) NA_real_ else cf[2L, "Pr(>|t|)"])
}, by = .(variant, d_gene, j_gene)]
cells <- merge(cells, pairing[variant == v], by = c("variant", "j_gene"))
# A cell is marked only inside a row the omnibus already called significant.
cells[, cell_significant := pairing_by_j_significant == TRUE & cell_p_value < 0.05]
pjd <- merge(pjd, cells, by = c("variant", "d_gene", "j_gene"))
stats <- fread(file.path(qtl, "pairing_associations.tsv.gz"))[variant == v][order(p_value)]
print(stats[, .(variant, j_gene, n, pillai, p_value, min_genotype_group, significant)], class = FALSE)

fwrite(prob, file.path(OUT$source, "dj_pairing_probabilities.csv.gz"))
fwrite(rbind(corr_jd, corr_dj), file.path(OUT$source, "dj_pairing_correlations.csv.gz"))
fwrite(rbind(data.table(axis = "IGHD", short_label = d_short, plot_order = seq_along(d_short)),
             data.table(axis = "IGHJ", short_label = j_short, plot_order = seq_along(j_short))), file.path(OUT$source, "dj_pairing_axis_levels.csv"))
fwrite(pj, file.path(OUT$source, "dj_snp_pj.csv.gz"))
fwrite(pd, file.path(OUT$source, "dj_snp_pd.csv.gz"))
fwrite(pjd, file.path(OUT$source, "dj_snp_pjd.csv.gz"))
fwrite(stats, file.path(OUT$source, "dj_snp_stats.tsv"), sep = "\t")

# ---- draw ----
geno_cols <- c(`0/0` = "#2a78d6", `0/1` = "#e34948", `1/1` = "#eda100")
d_short <- gsub("^IGH", "", d_short); j_short <- gsub("^IGH", "", j_short)
prob[, gene_a_short := factor(gsub("^IGH", "", as.character(gene_a_short)), levels = d_short)]
prob[, gene_b_short := factor(gsub("^IGH", "", as.character(gene_b_short)), levels = j_short)]
p_jd <- build_probability_boxplot(prob[, .(subject, gene_a_short, gene_b_short, probability_value = p_b_given_a, custom_order = order_j_given_d)],
                                  "gene_b_short", "probability_value", "gene_a_short", "IGHJ ASC", "P(J|D)", "IGHD ASC", make_named_palette(d_short))
p_dj <- build_probability_boxplot(prob[, .(subject, gene_a_short, gene_b_short, probability_value = p_a_given_b, custom_order = order_d_given_j)],
                                  "gene_a_short", "probability_value", "gene_b_short", "IGHD ASC", "P(D|J)", "IGHJ ASC", make_named_palette(j_short))
for (dt in list(corr_jd, corr_dj)) dt[, `:=`(x_short = gsub("^IGH", "", as.character(x_short)), y_short = gsub("^IGH", "", as.character(y_short)))]
corr_jd[, `:=`(x_short = factor(x_short, levels = j_short), y_short = factor(y_short, levels = j_short))]
corr_dj[, `:=`(x_short = factor(x_short, levels = d_short), y_short = factor(y_short, levels = d_short))]
tile_jd <- build_correlation_heatmap(corr_jd, "IGHJ ASC", NULL, TRUE, TRUE)
tile_dj <- build_correlation_heatmap(corr_dj, "IGHD ASC", NULL, TRUE, FALSE)

short <- function(x) ifelse(grepl("/", x), paste0(tstrsplit(x, "/")[[1L]], "*"), x)
pd[, d_gene_tag := factor(short(gsub("^IGH", "", d_gene)), levels = d_short)]
pjd[, d_gene_tag := factor(short(gsub("^IGH", "", d_gene)), levels = d_short)]
for (dt in list(pj, pd, pjd)) dt[, genotype := factor(c("0" = "0/0", "1" = "0/1", "2" = "1/1")[as.character(genotype)], levels = names(geno_cols))]
base_size <- 40  # panel C is drawn at the scale of A and B
base <- theme_bw(base_size = base_size) +
  theme(axis.title = element_text(size = base_size), axis.text = element_text(size = base_size), strip.text = element_text(size = base_size),
        legend.title = element_text(size = base_size), legend.text = element_text(size = base_size))
bare <- theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text = element_blank(), strip.background = element_blank())
geno_scale <- scale_colour_manual(values = geno_cols, drop = FALSE, name = paste0("Var: ", variant, " Genotype"))
# Stars sit above the tallest box (outliers excluded) of their group.
box_top <- function(x, f = 1.02) { q <- quantile(x, c(0.25, 0.75)); iqr <- q[2] - q[1]; max(x[x >= q[1] - 1.5 * iqr & x <= q[2] + 1.5 * iqr], na.rm = TRUE) * f }

p_geno <- ggplot(pjd[, .(n = uniqueN(subject)), by = genotype], aes(genotype, n, fill = genotype)) + geom_col(show.legend = FALSE) +
  scale_fill_manual(values = geno_cols, drop = FALSE) + labs(x = "Genotype", y = "# subjects") + base + theme(axis.ticks.x = element_blank())
# A star on a marginal panel means the variant is associated with that gene's own usage.
p_d <- ggplot(pd, aes(d_gene_tag, p_d, colour = genotype)) +
  geom_boxplot(outliers = FALSE, position = position_dodge(width = 0.9), width = 0.8, show.legend = FALSE) + geno_scale + facet_grid("1" ~ .) +
  scale_y_continuous(position = "right", limits = c(0, NA)) + labs(x = NULL, y = "P(D)") + base + bare +
  theme(strip.text = element_blank(), strip.background = element_blank(), strip.placement = "outside", strip.switch.pad.grid = unit(2.5, "cm"), axis.title.y.right = element_text(vjust = -14))
d_star <- unique(pd[usage_significant == TRUE, .(y = box_top(pd$p_d)), by = d_gene_tag])
if (nrow(d_star)) p_d <- p_d + geom_text(aes(d_gene_tag, y, label = "*"), data = d_star, colour = "#009E73", size = base_size / 2, vjust = 1, inherit.aes = FALSE)
p_j <- ggplot(pj, aes(genotype, p_j, colour = genotype)) + geom_boxplot(outliers = FALSE, show.legend = FALSE) + geno_scale +
  facet_grid(j_gene ~ ., scales = "free_y") + labs(x = NULL, y = "P(J)") + scale_y_continuous(limits = c(0, NA)) + base + bare
if (nrow(pj[usage_significant == TRUE])) {
  j_star <- pj[usage_significant == TRUE, .(y = box_top(p_j)), by = j_gene]
  p_j <- p_j + geom_text(aes(factor("1", levels = names(geno_cols)), y, label = "*"), data = j_star, colour = "#009E73", size = base_size / 2, vjust = 1, inherit.aes = FALSE)
}
# A green strip means that J gene's omnibus cleared the scan; the stars inside are the per-cell follow-up at nominal alpha.
strip_fill <- pjd[, .(fill = if (any(pairing_by_j_significant)) "#009E73" else "white"), by = j_gene][order(j_gene), fill]
star_cells <- unique(pjd[cell_significant == TRUE, .(j_gene, d_gene_tag, cell_p_value)])
if (nrow(star_cells)) {
  star_cells <- merge(star_cells, pjd[, .(y = box_top(p_j_given_d, 1.05)), by = j_gene], by = "j_gene")
  star_cells[, label := as.character(cut(cell_p_value, c(0, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "ns"), right = TRUE))]
}
p_cell <- ggplot(pjd, aes(d_gene_tag, p_j_given_d, colour = genotype)) +
  geom_boxplot(outliers = FALSE, position = position_dodge(width = 0.9), width = 0.8) + geno_scale + scale_y_continuous(position = "right", limits = c(0, NA)) +
  ggh4x::facet_grid2(j_gene ~ ., scales = "free_y", strip = ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = strip_fill))) +
  labs(x = NULL, y = "P(J|D)") + base +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.ticks.x = element_blank(), legend.position = "bottom",
        strip.placement = "outside", strip.switch.pad.grid = unit(2.5, "cm"), axis.title.y.right = element_text(vjust = -14))
if (nrow(star_cells)) p_cell <- p_cell + geom_text(aes(d_gene_tag, y, label = label), data = star_cells, colour = "#009E73", size = base_size / 2, vjust = 1, inherit.aes = FALSE)

keep_axis_space_y <- theme(axis.text.y = element_text(color = NA), axis.ticks.y = element_blank(), axis.title.x = element_text())
no_x <- theme(axis.text.x = element_blank())
# Design letters follow the order plots are added: A B C D E F, then I (spacer) and G H.
fig <- (p_jd + guides(color = guide_legend(nrow = 4)) +
          theme(axis.text.x = element_blank(), legend.position = "inside", legend.direction = "horizontal", legend.position.inside = c(0.87, 0.86)) + xlab(NULL) + no_x) +
  (tile_jd + theme(axis.text.x = element_text(angle = 0, vjust = 0.5), legend.position = "inside", legend.direction = "horizontal",
                   legend.position.inside = c(0.22, 3.45), legend.key.width = grid::unit(0.9, "inches")) + xlab("IGHJ ASC") + keep_axis_space_y) +
  (p_dj + guides(color = guide_legend(nrow = 1)) +
     theme(axis.text.x = element_blank(), legend.position = "inside", legend.direction = "horizontal", legend.position.inside = c(0.9, 0.93)) + xlab(NULL) + no_x) +
  (tile_dj + guides(fill = "none")) +
  p_j +
  (p_cell + theme(legend.position = "top", axis.text.x = element_blank())) +
  (p_geno + theme(axis.text.x = element_text(hjust = 0.5, vjust = 14), axis.title.x = element_blank())) +
  (p_d + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))) +
  plot_spacer() +
  plot_layout(design = "AA\nBB\nCC\nDD\nEF\nII\nGH", widths = c(1, 5), heights = c(2, 1, 2, 1.3, 8, -0.2, 2)) +
  plot_annotation(tag_levels = list(c("A", "", "B", "", "C", "", "", ""))) & theme(plot.tag = element_text(size = 58, face = "bold"))
ggsave(file.path(OUT$figures, "figure_dj_pairing.pdf"), fig, width = 42, height = 55, device = grDevices::cairo_pdf, limitsize = FALSE)
