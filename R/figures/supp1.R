#!/usr/bin/env Rscript
# Supplementary figure 1: Baseline / HUSA / IMGT sequence-set Venn per segment. Membership
# is exact sequence equality after degapping and uppercasing; IMGT is restricted to F/ORF
# alleles of genes present in HUSA (compared in the IMGT name space, via the husa column).

source("R/00_setup.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(Biostrings); library(alakazam); library(ggVennDiagram); library(ggvenn) })

segments <- c("IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ")
degap <- function(x) toupper(gsub("[.-]", "", as.character(x)))
join_alleles <- function(x) paste(sort(unique(x)), collapse = ",")

husa <- fread(need(file.path(OUT$husa, "husa.tsv")))[chain %in% CHAINS & gene_type %in% segments]
husa[, seq_norm := degap(seq)]
husa <- husa[!is.na(seq_norm) & nzchar(seq_norm)]
husa[, husa_gene := getGene(husa, first = TRUE, strip_d = FALSE, omit_nl = FALSE)]
husa[, in_baseline := as.logical(in_baseline_reference)]
husa[is.na(in_baseline), in_baseline := FALSE]
husa_genes <- unique(trimws(unlist(strsplit(getGene(husa$husa, first = FALSE, strip_d = FALSE, omit_nl = FALSE), ",", fixed = TRUE))))

imgt <- rbindlist(lapply(segments, function(seg) {
  dss <- readDNAStringSet(need(file.path(IN$imgt_vdj, sprintf("imgt_human_%s.fasta", seg))))
  parts <- strsplit(names(dss), "|", fixed = TRUE)
  data.table(gene_type = seg,
             allele = vapply(parts, function(p) if (length(p) >= 2L) p[[2]] else NA_character_, character(1)),
             functional = gsub("[][()]", "", vapply(parts, function(p) if (length(p) >= 4L) p[[4]] else NA_character_, character(1))),
             seq_norm = degap(dss))
}))[!is.na(allele) & nzchar(allele) & !is.na(seq_norm) & nzchar(seq_norm)]
imgt[, gene := getGene(allele, first = TRUE, strip_d = FALSE, omit_nl = FALSE)]
imgt_full <- imgt[functional %in% c("F", "ORF")]
imgt_inscope <- imgt_full[gene %in% husa_genes]

membership <- Reduce(function(x, y) merge(x, y, by = c("gene_type", "seq_norm"), all = TRUE), list(
  husa[in_baseline == TRUE, .(baseline_alleles = join_alleles(allele)), by = .(gene_type, seq_norm)],
  husa[, .(husa_alleles = join_alleles(allele)), by = .(gene_type, seq_norm)],
  imgt_inscope[, .(imgt_alleles = join_alleles(allele)), by = .(gene_type, seq_norm)]))
membership[, `:=`(in_baseline = !is.na(baseline_alleles), in_husa = !is.na(husa_alleles), in_imgt = !is.na(imgt_alleles))]
for (col in c("baseline_alleles", "husa_alleles", "imgt_alleles")) membership[is.na(get(col)), (col) := ""]
setorder(membership, gene_type, -in_baseline, -in_husa, -in_imgt, seq_norm)

overlap_counts <- membership[, .(
  n_baseline_only = sum(in_baseline & !in_husa & !in_imgt), n_husa_only = sum(!in_baseline & in_husa & !in_imgt),
  n_imgt_only = sum(!in_baseline & !in_husa & in_imgt), n_baseline_husa = sum(in_baseline & in_husa & !in_imgt),
  n_baseline_imgt = sum(in_baseline & !in_husa & in_imgt), n_husa_imgt = sum(!in_baseline & in_husa & in_imgt),
  n_all_three = sum(in_baseline & in_husa & in_imgt), n_baseline_total = sum(in_baseline), n_husa_total = sum(in_husa),
  n_imgt_in_scope_total = sum(in_imgt), n_husa_novel = sum(in_husa & !in_baseline),
  n_husa_novel_in_imgt = sum(in_husa & !in_baseline & in_imgt), n_husa_novel_not_imgt = sum(in_husa & !in_baseline & !in_imgt)
), by = gene_type]
overlap_counts <- merge(data.table(gene_type = segments), overlap_counts, by = "gene_type", all.x = TRUE)
for (col in setdiff(names(overlap_counts), "gene_type")) overlap_counts[is.na(get(col)), (col) := 0L]
overlap_counts <- overlap_counts[order(factor(gene_type, levels = segments))]
print(overlap_counts)

# HUSA-only sequences that IMGT does know under a gene absent from HUSA.
husa_only <- merge(unique(husa[, .(allele, husa_gene, gene_type, seq_norm)]), membership[in_husa & !in_imgt, .(gene_type, seq_norm)], by = c("gene_type", "seq_norm"))
flagged <- merge(husa_only, imgt_full[, .(imgt_alleles = join_alleles(allele), imgt_genes = join_alleles(gene)), by = .(gene_type, seq_norm)], by = c("gene_type", "seq_norm"))
if (nrow(flagged)) flagged <- flagged[vapply(imgt_genes, function(gs) !any(trimws(unlist(strsplit(gs, ",", fixed = TRUE))) %in% husa_genes), logical(1))]
flagged <- flagged[, .(husa_allele = allele, husa_gene, gene_type, seq_norm, imgt_alleles, imgt_genes)][order(gene_type, husa_allele)]

fwrite(membership[, .(gene_type, seq_norm, in_baseline, in_husa, in_imgt, baseline_alleles, husa_alleles, imgt_alleles)], file.path(OUT$source, "supp1_sequence_membership.csv"))
fwrite(overlap_counts, file.path(OUT$source, "supp1_overlap_counts.csv"))
fwrite(flagged, file.path(OUT$source, "supp1_gene_attribution_notes.csv"))

venn_colors <- c(Baseline = "#D64A4A", HUSA = "#58A65C", IMGT = "#8c4ac2")
panels <- lapply(segments, function(seg) {
  sub <- membership[gene_type == seg, .(Baseline = in_baseline, HUSA = in_husa, IMGT = in_imgt)]
  base_theme <- theme(plot.tag = element_text(face = "bold", size = 13, hjust = 0), plot.tag.position = c(0.07, 0.07), plot.margin = margin(4, 6, 4, 6))
  if (!nrow(sub) || !any(unlist(sub))) return(ggplot() + theme_void() + labs(tag = seg) + base_theme)
  data <- process_data(Venn(data_frame_to_list(as.data.frame(sub))))
  labels <- venn_setlabel(data)
  name_to_id <- setNames(data$setData$id, data$setData$name)
  labels$nudge_x <- ifelse(labels$X < 0, -0.08, 0.08) + labels$X
  labels$nudge_y <- ifelse(labels$Y < 0, -0.22, 0.22) + labels$Y
  ggplot() +
    geom_polygon(aes(X, Y, group = id), data = venn_regionedge(data), fill = "transparent") +
    geom_path(aes(X, Y, group = id, color = id), data = venn_setedge(data), linewidth = 0.8, show.legend = FALSE) +
    geom_text(aes(X, Y, label = count), data = venn_regionlabel(data), size = 4) +
    geom_label(aes(nudge_x, nudge_y, label = name, color = name), data = labels, size = 3.4, show.legend = FALSE) +
    coord_equal(clip = "off") +
    scale_color_manual(values = c(setNames(venn_colors, name_to_id[names(venn_colors)]), venn_colors)) +
    theme_void() + guides(color = "none") + labs(tag = seg) + base_theme
})

# Three heavy-chain panels centred over the four light-chain ones on an 8-half-column grid.
design <- c(area(1, 2, 1, 3), area(1, 4, 1, 5), area(1, 6, 1, 7), area(2, 1, 2, 2), area(2, 3, 2, 4), area(2, 5, 2, 6), area(2, 7, 2, 8))
fig <- wrap_plots(panels, design = design, heights = c(1, 1), widths = rep(1, 8))
suppressWarnings(ggsave(file.path(OUT$figures, "supp1.pdf"), fig, width = 12, height = 6))
