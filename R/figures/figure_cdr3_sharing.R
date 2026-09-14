#!/usr/bin/env Rscript
# IGK public-CDR3 sharing and its genetic control. A-C reuse the supp6 tables (run supp6.R
# first). D-F show the IGKV1D-13 coding-allele example from the QTL stage: the usage
# Manhattan, gene usage by lead-SNP genotype, and lead-CDR3 usage against gene usage (zeros
# re-added so the through-origin fit anchors over the full range).

source("R/00_setup.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
set.seed(42)

v1d13 <- file.path(OUT$qtl, "igkv1d13")
manhattan <- fread(need(file.path(v1d13, "panelA_manhattan.tsv")))
usage <- fread(need(file.path(v1d13, "panelB_usage_by_genotype.tsv")))
cdr3 <- fread(need(file.path(v1d13, "panelC_cdr3_vs_usage.tsv")))
ccdf <- fread(need(file.path(OUT$source, "supp6_igk_public_cdr3_ccdf.csv")))[locus == "IGK"]
rgs <- fread(need(file.path(OUT$source, "supp6_igk_rgs_vs_cdr3_overlap.csv.gz")))[locus == "IGK"]
alle <- fread(need(file.path(OUT$source, "supp6_igk_cdr3_allele_sharing.csv")))[locus == "IGK"]

# Panel C of the QTL stage kept only lead-CDR3 users; every other depth-passing subject lands at zero.
cdr3_full <- merge(usage[low_depth == FALSE, .(subject, usage, dosage, genotype)], cdr3[, .(subject, cdr3_usage)], by = "subject", all.x = TRUE)
cdr3_full[is.na(cdr3_usage), cdr3_usage := 0]
# Within-gene frequency of the lead CDR3 among carriers normalises out gene abundance.
within_gene <- merge(usage[low_depth == FALSE, .(subject, n_v, dosage, genotype)], cdr3[, .(subject, n_cdr3)], by = "subject", all.x = TRUE)
within_gene[is.na(n_cdr3), n_cdr3 := 0L]
within_gene <- within_gene[n_v > 0][, within_freq := n_cdr3 / n_v]
print(within_gene[, .(n = .N, median_within_freq = round(median(within_freq), 4)), by = genotype][order(genotype)])
fwrite(manhattan, file.path(OUT$source, "figure_cdr3_sharing_manhattan.csv"))
fwrite(usage, file.path(OUT$source, "figure_cdr3_sharing_usage_by_genotype.csv"))
fwrite(cdr3_full, file.path(OUT$source, "figure_cdr3_sharing_cdr3_vs_usage.csv"))
fwrite(within_gene, file.path(OUT$source, "figure_cdr3_sharing_within_gene.csv"))

# ---- draw ----
geno_cols <- c(`0/0` = "#2a78d6", `0/1` = "#e34948", `1/1` = "#eda100")
diploid <- function(g) factor(as.character(g), levels = c("0", "1", "2"), labels = c("0/0", "0/1", "1/1"))
base_size <- 18

pA <- ggplot(ccdf, aes(n_samples, prop_ge_x)) + geom_step(direction = "hv", linewidth = 0.7) + scale_y_log10() +
  labs(x = "# individuals", y = expression(log[10] ~ "CCDF (shared CDR3s)")) + ggpubr::theme_pubclean(base_size = base_size)
bins <- unique(rgs$overlap_bin_quant)
bins <- bins[order(as.numeric(sub(",.*$", "", sub("^[\\[(]", "", bins))))]
rgs[, overlap_bin_quant := factor(overlap_bin_quant, levels = bins)]
pB <- ggplot(rgs, aes(overlap_bin_quant, overlap_coefficient_cdr3)) + geom_boxplot(outlier.shape = NA, outliers = FALSE) +
  labs(x = "Individual RGS overlap coefficient", y = "CDR3 AA overlap coefficient") + ggpubr::theme_pubclean(base_size = base_size) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
pC <- ggplot(alle, aes(share_number, cdr3_aa_median_frequency, colour = label)) + geom_point(size = 1.3, alpha = 0.7, shape = 16) +
  scale_colour_brewer(palette = "Dark2") + labs(x = "# subjects sharing CDR3", y = "Median CDR3 frequency", colour = "V allele") +
  ggpubr::theme_pubclean(base_size = base_size) + theme(legend.position = "inside", legend.position.inside = c(0.5, 0.8))

thr <- -log10(0.05 / nrow(manhattan))
lead <- manhattan[is_lead == TRUE][1]
top <- max(manhattan$logp)
pD <- ggplot(manhattan, aes(pos / 1e6, logp, colour = logp)) +
  geom_hline(yintercept = thr, linetype = "dashed", colour = "grey55") + geom_point(size = 1.3, show.legend = FALSE) +
  annotate("segment", x = lead$pos / 1e6 - 0.1, xend = lead$pos / 1e6 - 0.01, y = lead$logp - top * 0.05, yend = lead$logp,
           arrow = arrow(length = unit(0.22, "cm"), type = "closed"), colour = "black", linewidth = 0.7) +
  annotate("text", x = lead$pos / 1e6 - 0.11, y = lead$logp - top * 0.05, label = lead$variant, size = 4, fontface = "italic", hjust = 1) +
  scale_colour_gradient(low = "gray90", high = "#08306b", name = expression(-log[10] * " p")) +
  labs(x = "chr2 position (Mb)", y = expression(-log[10] * " p")) + ggpubr::theme_pubclean(base_size = base_size) + theme(legend.position = "right")

ug <- usage[low_depth == FALSE][, gt := diploid(genotype)]
ecount <- ug[, .(n = .N), by = gt][, y := max(ug$usage) * 1.06]
pE <- ggplot(ug, aes(gt, usage, colour = gt)) + geom_boxplot(outlier.shape = NA, fill = NA, show.legend = FALSE) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = ecount, aes(gt, y, label = paste0("n=", n)), inherit.aes = FALSE, size = 5.2, colour = "grey25") +
  scale_colour_manual(values = geno_cols) + labs(x = "IGKV1D-13 lead-SNP genotype", y = "IGKV1D-13 gene usage") + ggpubr::theme_pubclean(base_size = base_size)
cu <- copy(cdr3_full)[, gt := diploid(genotype)]
pF <- ggplot(cu, aes(usage, cdr3_usage, colour = gt)) +
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, colour = "grey40", linewidth = 0.6, linetype = "dashed") +
  geom_point(size = 1.6, alpha = 0.8, show.legend = FALSE) + scale_colour_manual(values = geno_cols, name = "Genotype") +
  labs(x = "IGKV1D-13 gene usage", y = "Lead-CDR3 AA usage") + ggpubr::theme_pubclean(base_size = base_size) + theme(legend.position = "right")

fig <- pA + pB + pC + pD + pE + pF + plot_layout(design = "ABC\nDEF", heights = c(1, 1)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))
ggsave(file.path(OUT$figures, "figure_cdr3_sharing.pdf"), fig, width = 18, height = 11, device = grDevices::cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT$figures, "figure_cdr3_sharing.png"), fig, width = 18, height = 11, dpi = 150, limitsize = FALSE)
