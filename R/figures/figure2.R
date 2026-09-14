#!/usr/bin/env Rscript
# Figure 2: allele discovery across the baseline reference and the three genomic cohorts.
# A: per-gene-type Venn of allele membership (baseline, DS1 = 1KGP, DS2 = HPRC, DS3 = in-house).
# B: alleles ranked by carrier count, labels on the top ranks, novel/known pie inset.
# C: per-locus upset over ancestry, with the pairwise sample-to-sample allele overlap.

source("R/00_setup.R")
source("R/lib/upset_v3.R")
suppressPackageStartupMessages({ library(ggplot2); library(ggpubr); library(patchwork); library(ComplexUpset); library(cowplot)
  library(ggVennDiagram); library(ggvenn); library(ggrepel) })
set.seed(42)  # ggrepel places labels by random search

split_csv <- function(x) { v <- trimws(unlist(strsplit(x[!is.na(x) & nzchar(x)], ",", fixed = TRUE))); v[nzchar(v)] }
count_csv <- function(x) uniqueN(split_csv(x))

husa <- fread(need(file.path(OUT$husa, "husa.tsv")))
# The plotted name prefers the HUSA name, then the VDJbase name, then the ASC allele.
husa[, display_allele := fifelse(!is.na(husa) & nzchar(husa), husa, fifelse(!is.na(vdjbase_allele) & nzchar(vdjbase_allele), vdjbase_allele, allele))]
allele_counts <- husa[, .(sample_count_genomic_watson = count_csv(samples_genomic_watson),
                          sample_count_genomic_hprc = count_csv(samples_genomic_hprc),
                          sample_count_genomic_1kpg = count_csv(samples_genomic_1kpg),
                          in_baseline_reference = any(in_baseline_reference)), by = .(gene_type, allele, display_allele)]
allele_counts[, sample_count_genomic := sample_count_genomic_watson + sample_count_genomic_hprc + sample_count_genomic_1kpg]

# ---- A: source membership ----
panel_a <- allele_counts[, .(allele, gene_type, Baseline = in_baseline_reference, DS3 = sample_count_genomic_watson > 0L,
                             DS1 = sample_count_genomic_1kpg > 0L, DS2 = sample_count_genomic_hprc > 0L)]

# ---- B: carrier-count ranks. The three most carried alleles are labelled, plus the most
# carried one absent from the baseline; a label naming several alleles collapses to the first
# plus a dagger, listed separately for the caption. ----
rank_labels <- function(summary) {
  rank_parts <- list(); dagger_parts <- list()
  for (g in unique(summary$gene_type)) {
    sub <- summary[gene_type == g][order(order)]
    ids <- 1:min(3L, nrow(sub))
    if (sub[in_baseline_reference == FALSE, .N] > 0L) ids <- unique(c(ids, sub[in_baseline_reference == FALSE, min(order)]))
    sub[order %in% ids, label := gsub("IG[HKL]", "", display_allele)]
    multi <- sub[!is.na(label) & grepl(",", label)]
    if (nrow(multi) > 0L) {
      dagger_parts[[g]] <- data.table(gene_type = g, marker = paste0(sapply(strsplit(multi$label, ",", fixed = TRUE), function(x) trimws(x[[1]])), "†", seq_len(nrow(multi))),
                                      original_alleles = multi$display_allele)
    }
    sub[!is.na(label) & grepl(",", label), label := paste0(sapply(strsplit(label, ",", fixed = TRUE), function(x) trimws(x[[1]])), "†", seq_len(.N))]
    sub[!is.na(label) & nchar(label) > 30, label := paste0(sapply(strsplit(label, "_", fixed = TRUE), function(x) x[[1]]), "_‡", seq_len(.N))]
    sub[, ranked := !is.na(label)]
    rank_parts[[g]] <- sub
  }
  list(ranks = rbindlist(rank_parts), daggers = if (length(dagger_parts)) rbindlist(dagger_parts) else data.table())
}
genomic_summary <- allele_counts[sample_count_genomic > 0L]
genomic_summary[, sample_count := sample_count_genomic]
genomic_summary[order(sample_count, decreasing = TRUE), order := seq_len(.N), by = gene_type]
ranked <- rank_labels(genomic_summary)
rank_data <- ranked$ranks
rank_pie <- rank_data[, .(count = .N), by = .(gene_type, in_baseline_reference)]

# ---- C: which ancestries carry which allele, and how much two individuals share ----
ancestry_of <- list(
  DS3 = fread(need(IN$watson_metadata))[, setNames(ancestry_population, vdjbase_name)],
  DS2 = fread(need(IN$hprc_metadata))[, setNames(Ancestry, sample)],
  DS1 = fread(need(IN$kgp_metadata))[, setNames(Ancestry, sample)])
cohort <- function(count_col, samples_col, dataset) {
  husa[get(count_col) > 0, .(sample = split_csv(get(samples_col)), source = "GGS", genomic_dataset = dataset,
                             in_baseline_reference = any(in_baseline_reference)), by = .(allele, gene_type)
  ][, ancestry_population := ancestry_of[[dataset]][sample]]
}
genomic <- rbindlist(list(cohort("sample_count_genomic_watson", "samples_genomic_watson", "DS3"),
                          cohort("sample_count_genomic_hprc", "samples_genomic_hprc", "DS2"),
                          cohort("sample_count_genomic_1kpg", "samples_genomic_1kpg", "DS1")), fill = TRUE, use.names = TRUE)
genomic[, chain := substr(gene_type, 1, 3)]
ancestry_levels <- sort(unique(genomic$ancestry_population))

upset_membership <- rbindlist(lapply(unique(genomic$chain), function(ch) {
  d <- genomic[chain == ch]
  d <- d[!duplicated(paste(allele, sample, source, ancestry_population, in_baseline_reference))]
  d <- d[, .(present = uniqueN(sample) > 0, novel = any(!in_baseline_reference)), by = .(allele, gene_type, ancestry_population)]
  w <- dcast(d, gene_type + allele + novel ~ ancestry_population, value.var = "present", fill = FALSE)
  w[, chain := ch]
  w[, c("chain", "gene_type", "allele", "novel", ancestry_levels), with = FALSE]
}), use.names = TRUE, fill = TRUE)
for (col in ancestry_levels) upset_membership[is.na(get(col)), (col) := FALSE]

# Overlap coefficient and Jaccard distance for every pair of individuals within a locus and
# ancestry, counted as a matrix product.
overlap_key <- unique(genomic[!is.na(ancestry_population) & !is.na(allele), .(chain, ancestry_population, sample, allele)])
genomic_overlap <- rbindlist(lapply(split(overlap_key, by = c("chain", "ancestry_population")), function(grp) {
  samples <- unique(grp$sample); alleles <- unique(grp$allele)
  if (length(samples) < 2L) return(NULL)
  carried <- matrix(0, nrow = length(samples), ncol = length(alleles))
  carried[cbind(match(grp$sample, samples), match(grp$allele, alleles))] <- 1
  shared <- tcrossprod(carried); sizes <- diag(shared)
  ij <- which(upper.tri(shared), arr.ind = TRUE); i <- ij[, 1L]; j <- ij[, 2L]; n_shared <- shared[ij]
  data.table(locus = grp$chain[[1L]], ancestry_population = grp$ancestry_population[[1L]], Sample1 = samples[i], Sample2 = samples[j],
             jaccard_distance = 1 - n_shared / (sizes[i] + sizes[j] - n_shared), overlap_coefficient = n_shared / pmin(sizes[i], sizes[j]))
}))

fwrite(panel_a, file.path(OUT$source, "figure2_A.csv"))
fwrite(rank_data, file.path(OUT$source, "figure2_B.csv"))
fwrite(ranked$daggers, file.path(OUT$source, "figure2_dagger_alleles.csv"))
fwrite(rank_pie, file.path(OUT$source, "figure2_B_pie.csv"))
fwrite(genomic, file.path(OUT$source, "figure2_C_genomic_alleles.csv"))
fwrite(upset_membership, file.path(OUT$source, "figure2_C_upset.csv"))
fwrite(genomic_overlap, file.path(OUT$source, "figure2_C_genomic_overlap.csv"))
cat(sprintf("figure2: %d alleles, %d with genomic carriers, %d ancestries, %d sample pairs\n",
            nrow(panel_a), nrow(rank_data), length(ancestry_levels), nrow(genomic_overlap)))

# ---- draw ----
venn_sets <- c("Baseline", "DS3", "DS1", "DS2")
upset_sets <- rev(ancestry_levels)  # the upset lays its sets out bottom-up
venn_colors <- setNames(c("#D64A4A", "#eba525", "#8c4ac2", "#58A65C"), venn_sets)

venn_plots <- setNames(lapply(unique(panel_a$gene_type), function(g) {
  vdata <- process_data(Venn(data_frame_to_list(panel_a[gene_type == g, venn_sets, with = FALSE])))
  set_labels <- venn_setlabel(vdata)
  name_to_id <- setNames(vdata$setData$id, vdata$setData$name)
  set_labels$nudge_y <- ifelse(set_labels$Y < 0, -0.5, 0) + set_labels$Y
  p <- ggplot() +
    geom_polygon(aes(X, Y, group = id), data = venn_regionedge(vdata), fill = "transparent") +
    geom_path(aes(X, Y, group = id, color = id), data = venn_setedge(vdata), show.legend = FALSE) +
    geom_text(aes(X, Y, label = count), data = venn_regionlabel(vdata), size = 10) +
    coord_equal(clip = "off") +
    scale_color_manual(values = c(setNames(venn_colors, name_to_id[names(venn_colors)]), venn_colors)) +
    theme_void() + guides(color = "none") +
    theme(axis.title.x = element_text(angle = 0, size = 18, face = "bold", vjust = -3)) + xlab(g)
  if (g %in% c("IGHD", "IGKV", "IGLV")) p <- p + geom_label(aes(X, nudge_y, label = name, color = name), data = set_labels, size = 6)
  p
}), unique(panel_a$gene_type))
row1 <- plot_grid(venn_plots[["IGHV"]], venn_plots[["IGHD"]], venn_plots[["IGHJ"]], NULL, venn_plots[["IGKV"]], venn_plots[["IGKJ"]],
                  venn_plots[["IGLV"]], venn_plots[["IGLJ"]], nrow = 1, rel_widths = c(1 / 8, 1 / 8, 1 / 8, 1 / 32, 1 / 8, 1 / 8, 1 / 8, 1 / 8))

rank_data[ranked == FALSE, label := NA_character_]
rank_plots <- setNames(lapply(unique(rank_data$gene_type), function(g) {
  sub <- rank_data[gene_type == g][order(order)]
  ggplot(sub, aes(x = order, y = sample_count, label = label, color = !in_baseline_reference)) +
    geom_line(color = "black") +
    geom_point(aes(shape = ranked), size = ifelse(grepl("V", g), 1, 3)) +
    geom_label_repel(box.padding = unit(0.8, "lines"), point.padding = unit(0.8, "lines"), max.overlaps = Inf,
                     segment.color = "gray60", size = 8, na.rm = TRUE, alpha = 1, force_pull = 0.5, seed = 42) +
    geom_rug(sides = "b") +
    theme_pubclean(base_size = 28) +
    scale_x_continuous(breaks = function(x) pretty(x)[pretty(x) %% 1 == 0], limits = c(1, max(sub$order))) +
    scale_y_continuous(breaks = function(x) pretty(x)[pretty(x) %% 1 == 0], limits = c(1, max(sub$sample_count))) +
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
    labs(title = g, x = "Rank of GGS Alleles", y = "Number of Individuals", color = "Novel") +
    scale_color_manual(values = c("TRUE" = "steelblue3", "FALSE" = "gray60"))
}), unique(rank_data$gene_type))
pie_plots <- setNames(lapply(unique(rank_pie$gene_type), function(g) {
  ggplot(rank_pie[gene_type == g], aes(x = "", y = count, fill = in_baseline_reference, label = count)) +
    geom_bar(stat = "identity", width = 1, color = "gray60") +
    geom_text(position = position_stack(vjust = 0.5), color = "black", size = 8) +
    coord_polar("y", start = 0) + theme_void() +
    theme(legend.position = "none", panel.background = element_rect(fill = "white", color = "transparent")) +
    scale_fill_manual(values = c("FALSE" = "steelblue3", "TRUE" = "gray60"))
}), unique(rank_pie$gene_type))
add_pie <- function(p, pie, scale = 1, x = 0.75, y = 0.55) ggdraw(p + guides(color = "none", shape = "none")) + draw_plot(pie, x, y, .3, .3, scale = scale)
no_axis_titles <- theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
y_axis_only <- theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 24), legend.position = "none")
x_axis_only <- theme(axis.title.x = element_text(size = 24), axis.title.y = element_blank(), legend.position = "none")
both_axes <- theme(axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 24), legend.position = "none")
row2 <- plot_grid(
  plot_grid(add_pie(rank_plots[["IGHV"]] + y_axis_only + ylab(""), pie_plots[["IGHV"]], scale = 1.4),
            add_pie(rank_plots[["IGHD"]] + y_axis_only, pie_plots[["IGHD"]], scale = 1.4),
            add_pie(rank_plots[["IGHJ"]] + both_axes + ylab(""), pie_plots[["IGHJ"]], scale = 1.4), ncol = 1, axis = "r"),
  plot_grid(add_pie(rank_plots[["IGKV"]] + no_axis_titles, pie_plots[["IGKV"]], scale = 0.95, x = 0.7),
            add_pie(rank_plots[["IGKJ"]] + x_axis_only, pie_plots[["IGKJ"]], scale = 0.95, x = 0.7), ncol = 1),
  plot_grid(add_pie(rank_plots[["IGLV"]] + no_axis_titles, pie_plots[["IGLV"]], scale = 1, x = 0.7),
            add_pie(rank_plots[["IGLJ"]] + x_axis_only, pie_plots[["IGLJ"]], scale = 1, x = 0.7), ncol = 1),
  nrow = 1, rel_widths = c(3 / 8, 2 / 8, 2 / 8))

genomic_overlap[, ancestry_population := factor(ancestry_population, levels = ancestry_levels)]
overlap_plots <- setNames(lapply(unique(genomic_overlap$locus), function(ch) {
  ggplot(genomic_overlap[locus == ch], aes(x = ancestry_population, y = overlap_coefficient)) +
    geom_boxplot(position = "dodge") +
    upset_default_themes(text = element_text(size = 28)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
    scale_x_discrete(position = "top") + scale_y_continuous(position = "right") +
    labs(x = "", y = "Overlap Coefficient") + guides(color = "none")
}), unique(genomic_overlap$locus))
ancestry_upset <- function(membership, chain, corner) {
  upset_v3(data = membership, intersect = upset_sets, sort_sets = FALSE,
           base_annotations = list("Intersection size" = intersection_size(counts = FALSE, mapping = aes(fill = novel)) +
                                     labs(y = if (chain == "IGH") "Number of Alleles" else "") +
                                     scale_fill_manual(values = c("FALSE" = "gray60", "TRUE" = "steelblue3")) + guides(fill = "none")),
           set_sizes = upset_set_size(geom = geom_bar(aes(x = group, fill = novel), width = 0.5), position = "right") +
             scale_fill_manual(values = c("FALSE" = "gray60", "TRUE" = "steelblue3")) +
             scale_y_continuous(breaks = function(x) pretty(x)[pretty(x) %% 1 == 0]) + guides(fill = "none") +
             theme(axis.text.x = element_text(angle = 45, hjust = 1)),
           name = chain, guides = "collect", themes = upset_default_themes(text = element_text(size = 28), legend.position = "top"),
           wrap = FALSE, ignore_tag = TRUE, plot_for_spacer = corner + labs(x = "", y = if (chain == "IGL") "Overlap Coefficient" else "")) &
    theme(legend.position = "none")
}
upset_plots <- setNames(lapply(unique(upset_membership$chain), function(ch) ancestry_upset(upset_membership[chain == ch, !"chain"], ch, overlap_plots[[ch]])),
                        unique(upset_membership$chain))
row3 <- plot_grid(upset_plots[["IGH"]], upset_plots[["IGK"]], upset_plots[["IGL"]], rel_widths = c(3 / 8, 2 / 8, 2 / 8), nrow = 1)

final <- plot_grid(row1, row2, row3, nrow = 3, labels = list("A", "B", "C"), label_fontface = "bold", label_size = 36, rel_heights = c(.15, .55, .3))
ggsave(file.path(OUT$figures, "figure2.pdf"), final, width = 30, height = 30, units = "in", device = grDevices::cairo_pdf)
