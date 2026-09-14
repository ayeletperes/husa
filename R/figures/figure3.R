#!/usr/bin/env Rscript
# Figure 3: the AIRR-seq allele ranking (A) and the IGHV RSS and leader characterisation:
# sequence logos (B, C), unique allele and unique motif counts by agreement with the IGHV
# consensus (D, E), and the HUSA-vs-IMGT overlap of unique motifs (F, G).

source("R/00_setup.R")
source("R/lib/rss_helpers.R")
suppressPackageStartupMessages({ library(ggplot2); library(ggpubr); library(ggseqlogo); library(ComplexUpset); library(cowplot); library(ggrepel)
  library(Biostrings); library(alakazam) })
set.seed(42)

rank_gene_types <- c("IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ")
consensus_levels <- c("Both", "Either", "None")
split_csv <- function(x) { v <- trimws(unlist(strsplit(x[!is.na(x) & nzchar(x)], ",", fixed = TRUE))); v[nzchar(v)] }
count_csv <- function(x) uniqueN(split_csv(x))
degap <- function(x) toupper(gsub("[.-]", "", as.character(x)))
modal_len <- function(x) { x <- x[nzchar(x)]; if (!length(x)) return(NA_integer_); tl <- table(nchar(x)); as.integer(names(tl)[which.max(tl)]) }

husa <- fread(need(file.path(OUT$husa, "husa.tsv")))

# ---- A: alleles ranked by the number of AIRR-seq individuals carrying them ----
husa[, display_allele := fifelse(!is.na(husa) & nzchar(husa), husa, fifelse(!is.na(vdjbase_allele) & nzchar(vdjbase_allele), vdjbase_allele, allele))]
airr <- husa[, .(sample_count_AIRRseq = count_csv(samples_AIRRseq), in_baseline_reference = any(in_baseline_reference), present = any(present)),
             by = .(gene_type, allele, display_allele)][sample_count_AIRRseq > 0L]
airr[, sample_count := sample_count_AIRRseq]
airr[order(sample_count, decreasing = TRUE), order := seq_len(.N), by = gene_type]
rank_parts <- list(); dagger_parts <- list()
for (g in rank_gene_types) {
  sub <- airr[gene_type == g][order(order)]
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
rank_data <- rbindlist(rank_parts)
rank_pie <- rank_data[, .(count = .N), by = .(gene_type, in_baseline_reference)]
dagger_alleles <- if (length(dagger_parts)) rbindlist(dagger_parts) else data.table()

# ---- IMGT RSS reference, restricted to F/ORF alleles of genes in the baseline reference ----
reference_genes <- unlist(lapply(list.files(need(IN$watson_reference), pattern = "IG[HKL][VDJ]\\.fasta$", full.names = TRUE),
                                 function(f) unique(getGene(names(tigger::readIgFasta(f)), strip_d = FALSE, omit_nl = FALSE))))
imgt_label_map <- c("V-HEPTAMER" = "Heptamer", "V-NONAMER" = "Nonamer", "V-SPACER" = "Spacer", "V-RS" = "rss_aligned",
                    "J-HEPTAMER" = "Heptamer", "J-NONAMER" = "Nonamer", "J-SPACER" = "Spacer", "J-RS" = "rss_aligned",
                    "3'D-HEPTAMER" = "Heptamer", "3'D-NONAMER" = "Nonamer", "3'D-SPACER" = "Spacer", "3'D-RS" = "rss_aligned",
                    "5'D-HEPTAMER" = "Heptamer", "5'D-NONAMER" = "Nonamer", "5'D-SPACER" = "Spacer", "5'D-RS" = "rss_aligned")
imgt_rss <- fread(need(IN$imgt_rss))
imgt_rss[, gene := getGene(allele, strip_d = FALSE, omit_nl = FALSE)]
imgt_rss <- imgt_rss[functional %in% c("F", "ORF") & gene %in% reference_genes]
imgt_rss[, `:=`(label_new = imgt_label_map[label], gene_type = substr(allele, 1, 4))]
imgt_rss[, rss_type := fifelse(grepl("3", label), "IGHD_3", fifelse(grepl("5", label), "IGHD_5", gene_type))]
imgt_reshaped <- dcast(imgt_rss, allele + rss_type ~ label_new, value.var = "sequence")[!is.na(Heptamer) & !is.na(Nonamer) & !is.na(Spacer)]
imgt_reshaped[, rss_aligned := paste0(Heptamer, Spacer, Nonamer)]
imgt_reshaped[, length := nchar(rss_aligned)]

# ---- B/D/F: the IGHV RSS ----
husa_rss <- fread(need(file.path(OUT$husa, "husa_rss_filter.tsv")))
husa_rss[, iuis_allele := fifelse(!is.na(husa) & husa != "", husa, vdjbase_allele)]
rss_data <- prepare_aligned_rss_table(husa_rss)
rss_data <- merge(rss_data, unique(husa_rss[, .(allele, husa, vdjbase_allele, iuis_allele, iglabel_allele, iglabel, samples_genomic)]), by = "allele", all.x = TRUE)
for (nm in c("husa", "vdjbase_allele", "iuis_allele", "iglabel_allele", "iglabel", "samples_genomic")) {
  if (nm %in% names(rss_data)) next
  if (paste0(nm, ".x") %in% names(rss_data)) setnames(rss_data, paste0(nm, ".x"), nm) else if (paste0(nm, ".y") %in% names(rss_data)) setnames(rss_data, paste0(nm, ".y"), nm)
}
rss_data[, `:=`(rss_type = fifelse(!is.na(gene_type_plot) & gene_type_plot != "", gene_type_plot, gene_type),
                gene_type = fifelse(!is.na(gene_type_plot) & gene_type_plot != "", gene_type_plot, gene_type),
                Heptamer = heptamer, Spacer = spacer, Nonamer = nonamer, present = sample_count_AIRRseq > 0)]
rss_data <- rss_data[sure_subject_count > 0]
rss_data[, rss_sequence := paste0(Heptamer, "NNN", Nonamer)]  # the spacer is dropped from the plotted motif
rss_data <- rss_data[nchar(Nonamer) == 9 & nchar(Heptamer) == 7 & rss_type == "IGHV"]

# "?" marks a position with no majority base and matches anything.
rss_consensus <- rss_data[, .(consensus_heptamer = consensusString(DNAStringSet(Heptamer), ambiguityMap = "?"),
                              consensus_nonamer = consensusString(DNAStringSet(Nonamer), ambiguityMap = "?")), by = .(rss_type, gene_type)]
rss_consensus[, `:=`(consensus_heptamer_regex = gsub("?", "[ATCG]", consensus_heptamer, fixed = TRUE),
                     consensus_nonamer_regex = gsub("?", "[ATCG]", consensus_nonamer, fixed = TRUE))]
rss_data[rss_consensus, `:=`(consensus_heptamer_regex = i.consensus_heptamer_regex, consensus_nonamer_regex = i.consensus_nonamer_regex), on = .(rss_type, gene_type)]
rss_data[, `:=`(consensus_heptamer = mapply(function(s, p) grepl(p, s), Heptamer, consensus_heptamer_regex),
                consensus_nonamer = mapply(function(s, p) grepl(p, s), Nonamer, consensus_nonamer_regex))]
rss_data[, hn_pair_label := fcase(consensus_heptamer & consensus_nonamer, "Both", !consensus_heptamer & !consensus_nonamer, "None", default = "Either")]
rss_data[, rss_label := hn_pair_label[1], by = .(gene_type, rss_aligned)]
rss_data[, allele_label := fcase(all(hn_pair_label == "Both"), "Both", all(hn_pair_label == "None"), "None", default = "Either"), by = .(gene_type, allele)]
rss_unique_counts <- rbind(rss_data[, .(variable = "unique_allele_count", value = uniqueN(allele)), by = .(gene_type, consensus_label = allele_label)],
                           rss_data[, .(variable = "unique_rss_count", value = uniqueN(rss_aligned)), by = .(gene_type, consensus_label = rss_label)])
setcolorder(rss_unique_counts, c("gene_type", "variable", "value", "consensus_label"))
rss_by_sequence <- rss_data[, .(count = .N, present = any(present), count_genomic = uniqueN(split_csv(samples_genomic))), by = .(rss_aligned, gene_type, rss_type)]
rss_overlap <- merge(rss_by_sequence[, .(gene_type, rss_type, present, count_genomic, count, sequence = rss_aligned, length = nchar(rss_aligned))],
                     imgt_reshaped[rss_type == "IGHV" & length %in% c(39, 28, 38, 27), .(count = .N, gene_type = "IGHV", length = nchar(rss_aligned)), by = .(rss_type, rss_aligned)],
                     by.x = c("gene_type", "rss_type", "sequence", "length"), by.y = c("gene_type", "rss_type", "rss_aligned", "length"), all = TRUE, suffixes = c("_husa", "_imgt"))
rss_overlap[is.na(rss_overlap)] <- 0
rss_overlap[, `:=`(HUSA = count_husa > 0, IMGT = count_imgt > 0, present = present == 1, above_zero = count_genomic > 0)]

# ---- C/E/G: the IGHV leader (l_part1 + l_part2, matched to IMGT after degapping) ----
leader_data <- husa[as.logical(present) == TRUE & gene_type == "IGHV"]
leader_data[, `:=`(l1 = toupper(trimws(l_part1)), l2 = toupper(trimws(l_part2)))]
leader_data <- leader_data[nzchar(l1) & nzchar(l2) & l1 != "NA" & l2 != "NA"]
leader_data[, leader := paste0(l1, l2)]
husa_genes <- unique(trimws(unlist(strsplit(getGene(leader_data$husa, first = FALSE, strip_d = FALSE, omit_nl = FALSE), ",", fixed = TRUE))))
imgt_fasta <- readDNAStringSet(need(file.path(IN$imgt_leader, "imgt_human_IGHL.fasta")))
imgt_parts <- strsplit(names(imgt_fasta), "|", fixed = TRUE)
imgt_leader <- data.table(gene_type = "IGHV",
                          allele = vapply(imgt_parts, function(p) if (length(p) >= 2L) p[[2]] else NA_character_, character(1)),
                          functional = gsub("[][()]", "", vapply(imgt_parts, function(p) if (length(p) >= 4L) p[[4]] else NA_character_, character(1))),
                          leader_norm = degap(imgt_fasta))[!is.na(allele) & nzchar(allele) & nzchar(leader_norm)]
imgt_leader[, gene := getGene(allele, first = TRUE, strip_d = FALSE, omit_nl = FALSE)]
imgt_leader <- imgt_leader[functional %in% c("F", "ORF") & gene %in% husa_genes]

# The logo needs equal-length input, so the modal leader-1 and leader-2 lengths define both the
# logo subset and the consensus the classification is measured against.
leader_data[, l1_modal := modal_len(l1), by = gene_type]
leader_data[, l2_modal := modal_len(l2), by = gene_type]
cons_l1 <- leader_data[nchar(l1) == l1_modal, .(cons_l1 = consensusString(DNAStringSet(l1), ambiguityMap = "?")), by = gene_type]
cons_l2 <- leader_data[nchar(l2) == l2_modal, .(cons_l2 = consensusString(DNAStringSet(l2), ambiguityMap = "?")), by = gene_type]
leader_data[cons_l1, l1_regex := gsub("?", "[ACGT]", i.cons_l1, fixed = TRUE), on = "gene_type"]
leader_data[cons_l2, l2_regex := gsub("?", "[ACGT]", i.cons_l2, fixed = TRUE), on = "gene_type"]
leader_data[, l1_match := nchar(l1) == l1_modal & mapply(function(s, p) grepl(p, s), l1, l1_regex)]
leader_data[, l2_match := nchar(l2) == l2_modal & mapply(function(s, p) grepl(p, s), l2, l2_regex)]
leader_data[, leader_pair_label := fcase(l1_match & l2_match, "Both", !l1_match & !l2_match, "None", default = "Either")]
leader_data[, leader_label := leader_pair_label[1], by = .(gene_type, leader)]
leader_data[, allele_label := fcase(all(leader_pair_label == "Both"), "Both", all(leader_pair_label == "None"), "None", default = "Either"), by = .(gene_type, allele)]
leader_seqlogo_input <- leader_data[nchar(l1) == l1_modal & nchar(l2) == l2_modal, .(n_alleles = .N), by = .(gene_type, leader_logo_seq = paste0(l1, "NNN", l2))]
leader_unique_counts <- rbind(leader_data[, .(variable = "Alleles", value = uniqueN(allele)), by = .(gene_type, consensus_label = allele_label)],
                              leader_data[, .(variable = "Leaders", value = uniqueN(leader)), by = .(gene_type, consensus_label = leader_label)])
setcolorder(leader_unique_counts, c("gene_type", "variable", "value", "consensus_label"))
leader_overlap <- merge(leader_data[, .(present = TRUE, count_genomic = uniqueN(split_csv(samples_genomic)), count_husa = .N), by = .(gene_type, leader)],
                        imgt_leader[, .(count_imgt = .N), by = .(gene_type, leader = leader_norm)], by = c("gene_type", "leader"), all = TRUE)
for (col in c("count_husa", "count_imgt", "count_genomic")) leader_overlap[is.na(get(col)), (col) := 0]
leader_overlap[is.na(present), present := FALSE]
leader_overlap[, `:=`(HUSA = count_husa > 0, IMGT = count_imgt > 0, above_zero = count_genomic > 0)]

fwrite(rank_data, file.path(OUT$source, "figure3_A_rank.csv"))
fwrite(rank_pie, file.path(OUT$source, "figure3_A_pie.csv"))
fwrite(dagger_alleles, file.path(OUT$source, "figure3_A_dagger_alleles.csv"))
fwrite(rss_data[, .(gene_type, allele, rss_aligned, rss_sequence, Heptamer, Spacer, Nonamer)], file.path(OUT$source, "figure3_B_rss_seqlogo_input.csv"))
fwrite(rss_unique_counts, file.path(OUT$source, "figure3_B_rss_unique_counts.csv"))
fwrite(rss_overlap, file.path(OUT$source, "figure3_B_rss_imgt_overlap.csv"))
fwrite(leader_seqlogo_input, file.path(OUT$source, "figure3_C_leader_seqlogo_input.csv"))
fwrite(leader_unique_counts, file.path(OUT$source, "figure3_C_leader_unique_counts.csv"))
fwrite(leader_overlap[, .(gene_type, leader, HUSA, IMGT, present, count_genomic, above_zero)], file.path(OUT$source, "figure3_C_leader_imgt_overlap.csv"))
cat(sprintf("figure3: %d ranked alleles (%d labelled), %d IGHV RSS rows, %d IGHV leader alleles\n",
            nrow(rank_data), sum(rank_data$ranked), nrow(rss_data), nrow(leader_data)))

# ---- draw ----
rss_unique_counts[, `:=`(consensus_label = factor(consensus_label, levels = consensus_levels), variable = factor(variable, levels = c("unique_allele_count", "unique_rss_count")))]
leader_unique_counts[, `:=`(consensus_label = factor(consensus_label, levels = consensus_levels), variable = factor(variable, levels = c("Alleles", "Leaders")))]
consensus_colors <- setNames(c("#B8860B", "#008B8B", "#A63D40"), consensus_levels)
support_colors <- setNames(c("#0072B2", "#aa0415ff"), c("TRUE", "FALSE"))

rank_data[ranked == FALSE, label := NA_character_]
rank_plots <- setNames(lapply(rank_gene_types, function(g) {
  sub <- rank_data[gene_type == g][order(order)]
  ggplot(sub, aes(x = order, y = sample_count, label = label)) +
    geom_line(color = "black") +
    geom_point(aes(color = !in_baseline_reference, shape = ranked), size = ifelse(grepl("V", g), 1, 3)) +
    geom_rug(aes(color = !in_baseline_reference), sides = "b") +
    theme_pubclean(base_size = 38) +
    geom_label_repel(aes(color = !in_baseline_reference), box.padding = 0.85, point.padding = 0.5, max.overlaps = Inf, segment.color = "grey60", size = 10, na.rm = TRUE) +
    scale_x_continuous(breaks = function(x) pretty(x)[pretty(x) %% 1 == 0], limits = c(1, max(sub$order))) +
    scale_y_continuous(breaks = function(x) pretty(x)[pretty(x) %% 1 == 0], limits = c(1, max(sub$sample_count))) +
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
    labs(x = "Rank of GGS Alleles", y = "Number of Individuals", color = "Novel") +
    scale_color_manual(values = c("TRUE" = "steelblue3", "FALSE" = "gray60"))
}), rank_gene_types)
pie_plots <- setNames(lapply(rank_gene_types, function(g) {
  ggplot(rank_pie[gene_type == g], aes(x = "", y = count, fill = in_baseline_reference, label = count)) +
    geom_bar(stat = "identity", width = 1, color = "gray60") +
    geom_text(position = position_stack(vjust = 0.5), color = "black", size = 10) +
    coord_polar("y", start = 0) + theme_void() +
    theme(legend.position = "none", panel.background = element_rect(fill = "white", color = "transparent")) +
    scale_fill_manual(values = c("FALSE" = "steelblue3", "TRUE" = "gray60"))
}), rank_gene_types)
add_pie <- function(p, pie) ggdraw(p + guides(color = "none", shape = "none")) + draw_plot(pie, .75, .55, .3, .3)
no_axis_titles <- theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
y_axis_only <- theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 38), legend.position = "none")
x_axis_only <- theme(axis.title.x = element_text(size = 38), axis.title.y = element_blank(), legend.position = "none")
both_axes <- theme(axis.title.x = element_text(size = 38), axis.title.y = element_text(size = 38), legend.position = "none")
panel_a <- plot_grid(
  plot_grid(add_pie(rank_plots[["IGHV"]] + y_axis_only + labs(title = "IGHV", y = ""), pie_plots[["IGHV"]]),
            add_pie(rank_plots[["IGHD"]] + y_axis_only + labs(title = "IGHD"), pie_plots[["IGHD"]]),
            add_pie(rank_plots[["IGHJ"]] + both_axes + labs(title = "IGHJ", y = ""), pie_plots[["IGHJ"]]), ncol = 1),
  plot_grid(add_pie(rank_plots[["IGKV"]] + no_axis_titles + labs(title = "IGKV"), pie_plots[["IGKV"]]),
            add_pie(rank_plots[["IGKJ"]] + x_axis_only + labs(title = "IGKJ"), pie_plots[["IGKJ"]]), ncol = 1),
  plot_grid(add_pie(rank_plots[["IGLV"]] + no_axis_titles + labs(title = "IGLV"), pie_plots[["IGLV"]]),
            add_pie(rank_plots[["IGLJ"]] + x_axis_only + labs(title = "IGLJ"), pie_plots[["IGLJ"]]), ncol = 1),
  nrow = 1, rel_widths = c(1, 1, 1))

logo_theme <- theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
                    axis.title.x = element_blank(), axis.ticks.length.x = unit(0, "lines"), plot.margin = margin(0, 0, 0, 0, "pt"))
rss_logo_seqs <- rss_data$rss_sequence
p_rss_logo <- ggplot() + geom_logo(rss_logo_seqs, method = "probability", seq_type = "dna") + theme_logo() + logo_theme +
  scale_x_continuous(breaks = seq_len(max(nchar(rss_logo_seqs))), expand = c(0, 0))
leader_logo_seqs <- rep(leader_seqlogo_input$leader_logo_seq, leader_seqlogo_input$n_alleles)
p_leader_logo <- ggplot() + geom_logo(leader_logo_seqs, method = "probability", seq_type = "dna") + theme_logo() + logo_theme +
  scale_x_continuous(breaks = seq_len(max(nchar(leader_logo_seqs))), expand = c(0, 0))

counts_theme <- theme(axis.text.x = element_text(angle = 0), legend.position = "inside", legend.position.inside = c(0.75, 0.8),
                      legend.direction = "vertical", legend.box.background = element_blank())
p_rss_counts <- ggplot(rss_unique_counts, aes(x = variable, y = value, fill = consensus_label)) + geom_col() +
  labs(x = "", y = "Unique Count") + scale_x_discrete(labels = c("Alleles", "RSS")) +
  scale_y_continuous(breaks = function(x) pretty(x)[pretty(x) %% 1 == 0]) +
  scale_fill_manual(values = consensus_colors, name = "Hep-Non\ncons", drop = FALSE) + theme_pubclean(base_size = 42) + counts_theme
p_leader_counts <- ggplot(leader_unique_counts, aes(x = variable, y = value, fill = consensus_label)) + geom_col() +
  labs(x = NULL, y = "Unique Count") + scale_y_continuous(breaks = function(x) pretty(x)[pretty(x) %% 1 == 0]) +
  scale_fill_manual(values = consensus_colors, name = "L1-L2\ncons", drop = FALSE) + theme_pubclean(base_size = 42) + counts_theme

imgt_upset <- function(d, annotation_label, intersection_label) {
  upset(as.data.frame(d), c("IMGT", "HUSA"), name = "",
        annotations = setNames(list(
          ggplot(mapping = aes(x = intersection, y = count_genomic, color = present, shape = above_zero)) +
            geom_boxplot(na.rm = TRUE, outlier.shape = NA) + geom_point(position = position_jitterdodge(), size = 4) +
            scale_shape_manual(values = c(17, 16)) + scale_color_manual(values = support_colors) + guides(color = "none", shape = "none")), annotation_label),
        base_annotations = list("Intersection size" = intersection_size(counts = FALSE, text = list(size = 18), text_colors = c(on_background = "black", on_bar = "black"),
                                                                        bar_number_threshold = 0.65, mapping = aes(fill = present)) +
                                  labs(y = intersection_label) + guides(fill = guide_legend(title = "Expressed")) + scale_fill_manual(values = support_colors)),
        set_sizes = (upset_set_size() + scale_y_reverse(breaks = function(x) pretty(x)[pretty(x) %% 1 == 0])),
        guides = "over", themes = upset_default_themes(text = element_text(size = 44)))
}
p_rss_upset <- imgt_upset(rss_overlap, "# Carriers", "# Unique\nRSSs")
p_leader_upset <- imgt_upset(leader_overlap, "# Unique\nSamples", "# Unique\nLeaders")

rows_bc <- plot_grid(p_rss_logo, p_rss_counts, p_rss_upset, p_leader_logo, p_leader_counts, p_leader_upset,
                     ncol = 3, nrow = 2, rel_widths = c(1.6, 0.85, 1.4), rel_heights = c(1, 1), labels = c("B", "D", "F", "C", "E", "G"), label_size = 48)
final <- plot_grid(panel_a, rows_bc, ncol = 1, rel_heights = c(2.0, 2.4), labels = c("A", ""), label_size = 48)
ggsave(file.path(OUT$figures, "figure3.pdf"), final, width = 45, height = 32, device = grDevices::cairo_pdf, limitsize = FALSE)
