#!/usr/bin/env Rscript
# Supplementary figure 3: IMGT versus HUSA leaders for IGKV and IGLV. The HUSA leader is
# l_part1 + l_part2 of expressed alleles; IMGT stores the combined leader, so the two are
# matched by exact equality after degapping. A: logo over modal-length alleles, B: unique
# alleles versus unique leaders by consensus class, C: HUSA-vs-IMGT upset.
# Tables are computed only when missing from figure_data; the figure is
# always drawn from them.

source("scripts/00_setup.R")
suppressPackageStartupMessages({ library(ggplot2); library(ggseqlogo); library(patchwork); library(ComplexUpset) })
set.seed(42)

v_segments <- c("IGKV", "IGLV")
consensus_levels <- c("Consensus", "Both", "Non-Consensus")
tables <- file.path(OUT$source, c("supp3_overlap.csv", "supp3_overlap_counts.csv", "supp3_seqlogo_input.csv", "supp3_consensus_counts.csv"))

if (!all(file.exists(tables))) {
  suppressPackageStartupMessages({ library(Biostrings); library(alakazam) })
  leader_file <- c(IGKV = "imgt_human_IGKL.fasta", IGLV = "imgt_human_IGLL.fasta")
  degap <- function(x) toupper(gsub("[.-]", "", as.character(x)))
  modal_len <- function(x) { x <- x[nzchar(x)]; if (!length(x)) return(NA_integer_); tl <- table(nchar(x)); as.integer(names(tl)[which.max(tl)]) }
  count_samples <- function(s) { v <- trimws(unlist(strsplit(s[!is.na(s) & nzchar(s)], ",", fixed = TRUE))); uniqueN(v[nzchar(v)]) }

  husa <- fread(need(file.path(OUT$husa, "husa.tsv")))
  husa[, present := as.logical(present)]
  husa <- husa[present == TRUE & gene_type %in% v_segments]
  husa[, `:=`(l1 = toupper(trimws(l_part1)), l2 = toupper(trimws(l_part2)))]
  husa <- husa[nzchar(l1) & nzchar(l2) & l1 != "NA" & l2 != "NA"]
  husa[, leader := paste0(l1, l2)]
  husa_genes <- unique(trimws(unlist(strsplit(getGene(husa$husa, first = FALSE, strip_d = FALSE, omit_nl = FALSE), ",", fixed = TRUE))))

  imgt <- rbindlist(lapply(v_segments, function(seg) {
    dss <- readDNAStringSet(need(file.path(IN$imgt_leader, leader_file[[seg]])))
    parts <- strsplit(names(dss), "|", fixed = TRUE)
    data.table(gene_type = seg,
               allele = vapply(parts, function(p) if (length(p) >= 2L) p[[2]] else NA_character_, character(1)),
               functional = gsub("[][()]", "", vapply(parts, function(p) if (length(p) >= 4L) p[[4]] else NA_character_, character(1))),
               leader_norm = degap(dss))
  }))[!is.na(allele) & nzchar(allele) & nzchar(leader_norm)]
  imgt[, gene := getGene(allele, first = TRUE, strip_d = FALSE, omit_nl = FALSE)]
  imgt <- imgt[functional %in% c("F", "ORF") & gene %in% husa_genes]

  # Consensus over the modal-length leader-1 and leader-2; an allele or a leader is Consensus /
  # Non-Consensus when all / none of its rows match both, Both when they disagree.
  husa[, l1_modal := modal_len(l1), by = gene_type]
  husa[, l2_modal := modal_len(l2), by = gene_type]
  cons_l1 <- husa[nchar(l1) == l1_modal, .(cons_l1 = consensusString(DNAStringSet(l1), ambiguityMap = "?")), by = gene_type]
  cons_l2 <- husa[nchar(l2) == l2_modal, .(cons_l2 = consensusString(DNAStringSet(l2), ambiguityMap = "?")), by = gene_type]
  husa[cons_l1, l1_regex := gsub("?", "[ACGT]", i.cons_l1, fixed = TRUE), on = "gene_type"]
  husa[cons_l2, l2_regex := gsub("?", "[ACGT]", i.cons_l2, fixed = TRUE), on = "gene_type"]
  husa[, l1_match := nchar(l1) == l1_modal & mapply(function(s, p) grepl(p, s), l1, l1_regex)]
  husa[, l2_match := nchar(l2) == l2_modal & mapply(function(s, p) grepl(p, s), l2, l2_regex)]
  husa[, leader_consensus := l1_match & l2_match]
  cons_label <- function(flags) if (all(flags)) "Consensus" else if (!any(flags)) "Non-Consensus" else "Both"
  husa[, cons_allele_label := cons_label(leader_consensus), by = .(gene_type, allele)]
  husa[, cons_leader_label := cons_label(leader_consensus), by = .(gene_type, leader)]

  # NNN stands in for the intron; the logo is a probability logo over allele rows.
  logo <- husa[nchar(l1) == l1_modal & nchar(l2) == l2_modal, .(gene_type, leader_logo_seq = paste0(l1, "NNN", l2))]
  seqlogo_input <- logo[, .(n_alleles = .N), by = .(gene_type, leader_logo_seq)]
  consensus_counts <- rbind(
    husa[, .(value = uniqueN(allele), variable = "Alleles"), by = .(gene_type, consensus_label = cons_allele_label)],
    husa[, .(value = uniqueN(leader), variable = "Leaders"), by = .(gene_type, consensus_label = cons_leader_label)]
  )[, .(gene_type, variable, value, consensus_label)]

  comp <- merge(husa[, .(present = TRUE, count_genomic = count_samples(samples_genomic), count_husa = .N), by = .(gene_type, leader)],
                imgt[, .(count_imgt = .N), by = .(gene_type, leader = leader_norm)], by = c("gene_type", "leader"), all = TRUE)
  for (col in c("count_husa", "count_imgt", "count_genomic")) comp[is.na(get(col)), (col) := 0]
  comp[is.na(present), present := FALSE]
  comp[, `:=`(HUSA = count_husa > 0, IMGT = count_imgt > 0)]
  overlap_counts <- comp[, .(n_husa_only = sum(HUSA & !IMGT), n_shared = sum(HUSA & IMGT), n_imgt_only = sum(!HUSA & IMGT),
                             n_husa_total = sum(HUSA), n_imgt_total = sum(IMGT)), by = gene_type]
  print(overlap_counts)

  setkey(comp, gene_type, leader)
  setkey(seqlogo_input, gene_type, leader_logo_seq)
  fwrite(comp[, .(gene_type, leader, HUSA, IMGT, present, count_genomic)], tables[1])
  fwrite(overlap_counts, tables[2])
  fwrite(seqlogo_input, tables[3])
  fwrite(consensus_counts, tables[4])
}

comp <- fread(tables[1]); seqlogo_input <- fread(tables[3]); consensus_counts <- fread(tables[4])
consensus_counts[, consensus_label := factor(consensus_label, levels = consensus_levels)]
consensus_counts[, variable := factor(variable, levels = c("Alleles", "Leaders"))]
consensus_colors <- c(Consensus = "#B8860B", Both = "#008B8B", `Non-Consensus` = "#A63D40")
upset_colors <- setNames(c("#0072B2", "#aa0415ff"), c("TRUE", "FALSE"))

logos <- lapply(v_segments, function(seg) {
  d <- seqlogo_input[gene_type == seg]
  ggplot() + geom_logo(rep(d$leader_logo_seq, d$n_alleles), method = "probability", seq_type = "dna") + theme_logo() +
    labs(x = NULL, y = seg) +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_text(size = 28, vjust = -1.5), plot.margin = margin(4, 6, 4, 6))
})
counts <- lapply(v_segments, function(seg) {
  ggplot(consensus_counts[gene_type == seg], aes(x = variable, y = value, fill = consensus_label)) +
    geom_col() +
    scale_fill_manual(values = consensus_colors, name = "L1-L2", drop = FALSE) +
    scale_y_continuous(breaks = function(x) pretty(x)[pretty(x) %% 1 == 0]) +
    labs(x = NULL, y = "Unique Count") +
    ggpubr::theme_pubclean(base_size = 28) + theme(legend.position = "none")
})
upsets <- lapply(v_segments, function(seg) {
  d <- as.data.frame(comp[gene_type == seg][, above_zero := count_genomic > 0])
  upset(d, c("IMGT", "HUSA"), name = "",
        annotations = list("# Unique\nSamples" = (
          ggplot(mapping = aes(x = intersection, y = count_genomic, color = present, shape = above_zero)) +
            geom_boxplot(na.rm = TRUE, outlier.shape = NA) +
            geom_point(position = position_jitterdodge(seed = 42), size = 1.3) +
            scale_shape_manual(values = c(17, 16)) + scale_color_manual(values = upset_colors) +
            guides(color = "none", shape = "none"))),
        base_annotations = list("Intersection size" = intersection_size(
          counts = FALSE, text_colors = c(on_background = "black", on_bar = "black"), bar_number_threshold = 0.65,
          mapping = aes(fill = present)) + labs(y = "# Unique\nLeaders") + guides(fill = "none") + scale_fill_manual(values = upset_colors)),
        set_sizes = (upset_set_size() + scale_y_reverse(breaks = function(x) pretty(x)[pretty(x) %% 1 == 0])),
        guides = "over", themes = upset_default_themes(text = element_text(size = 28)))
})

hh <- rep(1, length(v_segments))
fig <- wrap_elements(wrap_plots(logos, ncol = 1, heights = hh)) +
  wrap_elements(wrap_plots(counts, ncol = 1, heights = hh)) +
  wrap_elements(wrap_plots(upsets, ncol = 1, heights = hh)) +
  plot_layout(widths = c(1.6, 0.85, 1.4)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 28, face = "bold"))
suppressWarnings(ggsave(file.path(OUT$figures, "supp3.pdf"), fig, width = 30, height = 15, limitsize = FALSE))
