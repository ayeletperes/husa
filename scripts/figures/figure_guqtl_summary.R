#!/usr/bin/env Rscript
# guQTL summary figure. A: a cartoon of the three steps of the scan (seeded synthetic values,
# not a result). B: significant gene-usage QTL SNP counts by location, per locus and segment,
# from the QTL variant summary; UTR is folded into intergenic.
# Tables are computed only when missing from figure_data; the figure is
# always drawn from them.

source("scripts/00_setup.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

feat_levels <- c("coding", "leader", "rss", "intergenic")
seg_levels <- c("V", "D", "DJ", "J")
loc_levels <- c("IGH", "IGK", "IGL")
tables <- file.path(OUT$source, paste0("figure_guqtl_summary_", c("counts", "chrom", "gene_model", "gene_snps", "assoc", "assoc_trend", "manhattan"), ".csv"))
names(tables) <- c("counts", "chrom", "gene_model", "gene_snps", "assoc", "assoc_trend", "manhattan")

if (!all(file.exists(tables))) {
  set.seed(42)
  vs <- fread(need(file.path(OUT$qtl, "reports", "variant_summary.tsv")))[analysis == "usage"]
  vs[feature == "utr", feature := "intergenic"]
  counts <- vs[feature %in% feat_levels, .(n_significant = sum(n_variants), n_kept = sum(n_variants_kept)), by = .(locus, segment, feature)]
  setorderv(counts, c("locus", "segment", "feature"))
  print(counts)

  # Schematic: a chromosome with widely spaced genes (V = leader | coding | rss, D = rss | coding | rss,
  # J = rss | coding), SNPs mostly in the intergenic gaps, an illustrative association and a kb-scale Manhattan.
  chrom <- data.table(xmin = 0, xmax = 132, ymin = 0.97, ymax = 1.03)
  gene_model <- rbindlist(list(
    data.table(feature = c("leader", "coding", "rss"), xmin = c(12, 14, 20), xmax = c(14, 20, 22)),
    data.table(feature = c("leader", "coding", "rss"), xmin = c(36, 38, 44), xmax = c(38, 44, 46)),
    data.table(feature = c("rss", "coding", "rss"), xmin = c(62, 64, 70), xmax = c(64, 70, 72)),
    data.table(feature = c("leader", "coding", "rss"), xmin = c(88, 90, 96), xmax = c(90, 96, 98)),
    data.table(feature = c("rss", "coding"), xmin = c(114, 116), xmax = c(116, 122))
  ))[, `:=`(ymin = 0.75, ymax = 1.25)]
  gene_snps <- data.table(x = c(5, 28, 30, 53, 56, 80, 83, 105, 128, 17, 41, 67, 93, 119, 13, 37, 89, 21, 63, 71, 115),
                          feature = c(rep("intergenic", 9), rep("coding", 5), rep("leader", 3), rep("rss", 4)))
  assoc <- rbindlist(lapply(0:2, function(g) data.table(genotype = g, usage = pmax(0, rnorm(18, 0.28 + 0.13 * g, 0.045)))))
  assoc_trend <- data.table(genotype = 0:2, usage = 0.28 + 0.13 * (0:2))
  manhattan <- rbind(data.table(pos = sort(runif(70, 0, 1500)), logp = abs(rnorm(70, 0, 0.7)), feature = "intergenic", significant = FALSE),
                     data.table(pos = c(340, 705, 780, 1020, 1215), logp = c(6.4, 9.1, 4.6, 7.2, 5.1),
                                feature = c("coding", "coding", "leader", "rss", "coding"), significant = TRUE))
  for (nm in names(tables)) fwrite(get(nm), tables[[nm]])
}

for (nm in names(tables)) assign(nm, fread(tables[[nm]]))
for (dt in list(gene_model, gene_snps, manhattan)) set(dt, j = "feature", value = factor(dt$feature, levels = feat_levels))
assoc[, genotype := factor(genotype, levels = 0:2)]
assoc_trend[, genotype := factor(genotype, levels = 0:2)]
counts[, `:=`(feature = factor(feature, levels = feat_levels), segment = factor(segment, levels = seg_levels), locus = factor(locus, levels = loc_levels))]

feat_cols <- c(coding = "#2a78d6", leader = "#eda100", rss = "#1baf7a", intergenic = "#8c6bb1")
geno_labels <- c(`0` = "0/0", `1` = "0/1", `2` = "1/1")
step_title <- function() theme(plot.title = element_text(size = 13, hjust = 0.5, face = "bold"))
set.seed(42)  # the jittered association points

g1 <- ggplot() +
  geom_rect(data = chrom, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = feat_cols[["intergenic"]], colour = NA) +
  geom_segment(data = gene_snps, aes(x = x, xend = x, y = 1.03, yend = 1.42), colour = "grey55") +
  geom_rect(data = gene_model, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = feature), colour = "grey25", linewidth = 0.3, show.legend = FALSE) +
  geom_point(data = gene_snps, aes(x = x, y = 1.45, colour = feature), size = 3.6, show.legend = FALSE) +
  scale_fill_manual(values = feat_cols) + scale_colour_manual(values = feat_cols) +
  coord_fixed(ratio = 12, xlim = c(-3, 135), ylim = c(0.55, 1.7), expand = FALSE) +
  ggtitle("1. Variants across the locus") + theme_void(base_size = 13) + step_title()
g2 <- ggplot(assoc, aes(genotype, usage)) +
  geom_line(data = assoc_trend, aes(genotype, usage, group = 1), colour = "grey40", linewidth = 0.7, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "grey15", width = 0.6) +
  geom_jitter(width = 0.14, size = 1, alpha = 0.6, colour = "grey15") +
  scale_x_discrete(labels = geno_labels) +
  labs(x = "Genotype", y = "Gene usage") + ggtitle("2. Test each SNP: usage ~ genotype") +
  ggpubr::theme_pubclean(base_size = 13) + theme(panel.grid.minor = element_blank()) + step_title()
g3 <- ggplot() +
  geom_hline(yintercept = 3.5, linetype = "dashed", colour = "grey55") +
  geom_point(data = manhattan[significant == FALSE], aes(pos, logp), colour = "grey70", size = 1) +
  geom_point(data = manhattan[significant == TRUE], aes(pos, logp, colour = feature), size = 3) +
  scale_colour_manual(values = feat_cols, drop = FALSE) +
  scale_x_continuous(labels = function(x) paste0(x, " kb")) +
  labs(x = "Position", y = expression(-log[10] * " p")) + ggtitle("3. Significant SNPs, by location") +
  ggpubr::theme_pubclean(base_size = 13) + theme(panel.grid.minor = element_blank(), legend.position = "none") + step_title()
arrow_plot <- function() ggplot() +
  annotate("segment", x = 0, xend = 1, y = 0, yend = 0, arrow = arrow(length = unit(0.18, "cm"), type = "closed"), linewidth = 0.8) +
  coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-1, 1)) + theme_void()
row_a <- g1 + arrow_plot() + g2 + arrow_plot() + g3 + plot_layout(nrow = 1, widths = c(2.4, 0.15, 1.0, 0.15, 1.15))

p_bar <- ggplot(counts, aes(segment, n_significant, fill = feature)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.1), width = 0.85) +
  geom_text(aes(label = n_significant), position = position_dodge2(width = 0.85, preserve = "single"), vjust = -0.35, size = 2.8, colour = "grey25") +
  facet_grid(. ~ locus, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = feat_cols, name = "SNP location") +
  scale_y_log10(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Segment", y = expression(log[10] * "(# significant gene-usage SNPs)")) +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

fig <- wrap_elements(row_a) / p_bar + plot_layout(heights = c(1.15, 1.5)) +
  plot_annotation(tag_levels = list(c("A", "B"))) & theme(plot.tag = element_text(size = 17, face = "bold"))
ggsave(file.path(OUT$figures, "figure_guqtl_summary.pdf"), fig, width = 14, height = 10, device = grDevices::cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT$figures, "figure_guqtl_summary.png"), fig, width = 14, height = 10, dpi = 150, limitsize = FALSE)
