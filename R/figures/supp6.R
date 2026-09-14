#!/usr/bin/env Rscript
# Supplementary figure 6: IGK CDR3 amino acid sharing. A: public-CDR3 CCDF (how many
# individuals share a CDR3). B: RGS-allele overlap versus CDR3 overlap per pair of
# individuals, on depth-matched subsamples. C: CDR3 sharing against the V allele it is used
# with. Panel B is a full pairwise sweep and takes tens of minutes; its tables are reused by
# the CDR3-sharing main figure.

source("R/00_setup.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(ggpubr) })

igk_sample_min <- 5000L  # individuals need this many IGK sequences; pairwise overlap is measured on subsamples of this size
n_bins <- 4L
n_draws <- 10L
set.seed(42)

rep_dt <- fread(need(file.path(OUT$repertoire, "gg_repertoire_data_IGK_genotype_corrected.csv.gz")))
husa <- fread(need(file.path(OUT$husa, "husa.tsv")))

dt <- rep_dt[vdjbase_subject %in% unique(rep_dt[valid_both == TRUE, vdjbase_subject])]
dt[, depth := .N, by = vdjbase_subject]
dt[, cdr3_aa := substr(junction_aa, 2, nchar(junction_aa) - 1)]  # trim the conserved first and last junction residues
dt <- dt[depth >= igk_sample_min]
# V alleles are named by their ASC group; panel C labels them by the IUIS name, falling back to the ASC name.
iuis_of <- with(unique(husa[, .(allele, husa)])[!duplicated(allele)], setNames(husa, allele))
dt[, v_call_iuis := iuis_of[v_call_new]]
dt[is.na(v_call_iuis) | !nzchar(v_call_iuis), v_call_iuis := v_call_new]
overlap_coefficient <- function(x, y) { denom <- min(length(x), length(y)); if (denom == 0L) NA_real_ else length(intersect(x, y)) / denom }

# ---- A ----
prev <- unique(dt[, .(subject = vdjbase_subject, cdr3_aa, locus)])[, .(n_samples = .N), by = .(cdr3_aa, locus)]
ccdf <- prev[, { fn <- ecdf(n_samples); grid <- seq_len(max(n_samples)); prop <- 1 - fn(grid - 1e-9)
  data.table(n_samples = grid, prop_ge_x = prop, count_ge_x = as.integer(round(.N * prop))) }, by = locus]

# ---- B ----
igk_samples <- unique(dt$vdjbase_subject)
junctions <- dt[, .(cdr3_aa = list(unique(cdr3_aa))), by = vdjbase_subject]
setkey(junctions, vdjbase_subject)
cdr3_overlap <- CJ(id1 = igk_samples, id2 = igk_samples)[id1 < id2][, {
  j1 <- junctions[.(id1), cdr3_aa][[1]]; j2 <- junctions[.(id2), cdr3_aa][[1]]
  if (!length(j1) || !length(j2)) .(overlap_coefficient = NA_real_, locus = "IGK") else {
    oc <- numeric(n_draws)
    for (i in seq_len(n_draws)) oc[i] <- overlap_coefficient(j1[sample.int(length(j1), min(length(j1), igk_sample_min))], j2[sample.int(length(j2), min(length(j2), igk_sample_min))])
    .(overlap_coefficient = mean(oc, na.rm = TRUE), locus = "IGK")
  }
}, by = .(id1, id2)]
airr_igk <- husa[sample_count_AIRRseq > 0, .(allele, gene_type, samples_AIRRseq)][, .(sample = trimws(unlist(strsplit(samples_AIRRseq, ",", fixed = TRUE)))), by = .(allele, gene_type)
][nzchar(sample) & grepl("IGK", gene_type) & sample %in% igk_samples]
rgs_samples <- unique(airr_igk$sample)
rgs_overlap <- CJ(id1 = rgs_samples, id2 = rgs_samples)[id1 < id2][, .(overlap_coefficient = overlap_coefficient(airr_igk[sample == id1, unique(allele)], airr_igk[sample == id2, unique(allele)]), locus = "IGK"), by = .(id1, id2)]
rgs_vs_cdr3 <- merge(rgs_overlap, cdr3_overlap, by = c("id1", "id2", "locus"), suffixes = c("_airr_seq", "_cdr3"))
# Quartile bins of the genomic overlap, so each box holds the same number of pairs.
rgs_vs_cdr3[, overlap_bin_quant := {
  qs <- unique(unname(quantile(overlap_coefficient_airr_seq, probs = seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE, type = 7)))
  if (length(qs) < n_bins + 1L) qs <- seq(min(overlap_coefficient_airr_seq, na.rm = TRUE), max(overlap_coefficient_airr_seq, na.rm = TRUE), length.out = n_bins + 1L)
  cut(overlap_coefficient_airr_seq, breaks = qs, include.lowest = TRUE, right = FALSE)
}, by = locus]
bin_levels <- levels(rgs_vs_cdr3$overlap_bin_quant)
box_quartiles <- rgs_vs_cdr3[, .(q1 = quantile(overlap_coefficient_cdr3, 0.25, na.rm = TRUE), median = median(overlap_coefficient_cdr3, na.rm = TRUE),
                                 q3 = quantile(overlap_coefficient_cdr3, 0.75, na.rm = TRUE), n_pairs = .N), by = .(locus, overlap_bin_quant)][order(overlap_bin_quant)]

# ---- C ----
n_subjects <- uniqueN(dt$vdjbase_subject)
sharing <- dt[, .(share_number = uniqueN(vdjbase_subject)), by = .(locus, cdr3_aa, v_call_iuis)]
sharing[, unique_v_allele := uniqueN(v_call_iuis), by = .(locus, cdr3_aa)]
allele_frequency <- dt[, .(number_of_subjects = uniqueN(vdjbase_subject)), by = .(locus, v_call_iuis)]
cdr3_frequency <- dt[, .(cdr3_aa_frequency = .N / unique(depth)), by = .(locus, cdr3_aa, vdjbase_subject)][, .(
  cdr3_aa_median_frequency = median(cdr3_aa_frequency, na.rm = TRUE), cdr3_aa_25_frequency = quantile(cdr3_aa_frequency, probs = 0.25, na.rm = TRUE),
  cdr3_aa_75_frequency = quantile(cdr3_aa_frequency, probs = 0.75, na.rm = TRUE), subjects_with_cdr3 = uniqueN(vdjbase_subject)), by = .(locus, cdr3_aa)]
merged <- merge(merge(sharing, allele_frequency, by = c("locus", "v_call_iuis"), all.x = TRUE), cdr3_frequency, by = c("locus", "cdr3_aa"), all.x = TRUE)
# Alleles carried by more than 60% of individuals say nothing about sharing; a CDR3 must be shared by
# at least five people and tied to a single V allele.
allele_cdr3 <- merged[number_of_subjects <= round(n_subjects * 0.6) & share_number >= 5 & unique_v_allele == 1]
allele_cdr3[, label := v_call_iuis]
allele_cdr3[grepl("IGKV1D-39|IGKV1-39", v_call_iuis), label := "IGKV1-39/IGKV1D-39*01_†"]  # indistinguishable in the repertoire
allele_cdr3[, label := paste0(label, " (n=", number_of_subjects, ")")]
setorder(allele_cdr3, label, -share_number, cdr3_aa)

fwrite(data.table(axis = "overlap_bin", level = bin_levels, plot_order = seq_along(bin_levels)), file.path(OUT$source, "supp6_axis_levels.csv"))
fwrite(ccdf, file.path(OUT$source, "supp6_igk_public_cdr3_ccdf.csv"))
fwrite(rgs_vs_cdr3, file.path(OUT$source, "supp6_igk_rgs_vs_cdr3_overlap.csv.gz"))
fwrite(box_quartiles, file.path(OUT$source, "supp6_igk_overlap_bin_quartiles.csv"))
fwrite(allele_cdr3, file.path(OUT$source, "supp6_igk_cdr3_allele_sharing.csv"))
cat(sprintf("supp6: %d IGK individuals, %d pairs, %d V-allele-specific CDR3s\n", length(igk_samples), nrow(rgs_vs_cdr3), nrow(allele_cdr3)))

# ---- draw ----
base_size <- 24
p_public <- ggplot(ccdf, aes(x = n_samples, y = log10(prop_ge_x))) + geom_step(direction = "hv") +
  labs(x = "Number of Individuals", y = expression(log[10] ~ "CCDF (shared CDR3s)")) + facet_wrap(~locus, scales = "free") + theme_pubclean(base_size = base_size)
p_overlap <- ggplot(rgs_vs_cdr3, aes(x = overlap_bin_quant, y = overlap_coefficient_cdr3)) + geom_boxplot(outliers = FALSE) +
  coord_cartesian(ylim = c(min(box_quartiles$q1) - 0.02, max(box_quartiles$q3) + 0.02)) +
  labs(x = "Individual RGS Overlap Coefficient", y = "CDR3 AA Overlap Coefficient") + facet_wrap(~locus, scales = "free") +
  theme_pubclean(base_size = base_size) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
set.seed(42)
p_allele <- ggplot(allele_cdr3, aes(x = share_number, y = cdr3_aa_median_frequency, color = label)) +
  geom_point(alpha = 0.8, shape = 16, size = 2, position = position_jitter(width = 0.2, height = 0, seed = 42)) +
  labs(x = "Number of Subjects Sharing CDR3 AA", y = "Median Frequency of CDR3 AA", color = "V Allele") +
  scale_color_brewer(palette = "Dark2") + facet_wrap(~locus, scales = "free") + guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_pubclean(base_size = base_size) + theme(legend.position = "inside", legend.position.inside = c(0.4, 0.78), legend.background = element_rect(fill = "transparent"))
fig <- p_public + p_overlap + p_allele + plot_layout(design = "ABC", widths = c(1, 1, 1.3)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 32, face = "bold"))
ggsave(file.path(OUT$figures, "supp6.pdf"), fig, width = 26, height = 10, device = grDevices::cairo_pdf, limitsize = FALSE)
