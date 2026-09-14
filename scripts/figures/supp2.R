#!/usr/bin/env Rscript
# Supplementary figure 2: the RSS landscape of every segment except IGHV. Per segment:
# A: sequence logo plus the nine most frequent RSSs with their motif tag and allele count,
# B: unique alleles and unique RSSs by agreement with the segment consensus,
# C: HUSA versus IMGT overlap with the per-RSS genomic sample counts.
# Tables are computed only when missing from figure_data; the figure is
# always drawn from them.

source("scripts/00_setup.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(ggrepel); library(ggseqlogo); library(ComplexUpset); library(ggpubr) })
set.seed(42)

gene_types <- c("IGLV", "IGKV", "IGHD_3", "IGHJ", "IGLJ", "IGKJ", "IGHD_5")  # panel order down the figure
tables <- file.path(OUT$source, paste0("supp2_", c("rss_aligned.csv.gz", "seqlogo_input.csv", "rss_letters.csv", "rss_counts.csv", "rss_tags.csv",
                                                  "imgt_overlap.csv", "imgt_summary.csv", "unique_counts.csv", "panel_layout.csv")))
names(tables) <- c("aligned", "seqlogo", "letters", "counts", "tags", "overlap", "summary", "unique", "layout")

if (!all(file.exists(tables))) {
source("scripts/lib/rss_helpers.R")

# ---- IMGT reference RSSs, F/ORF alleles of genes in the baseline reference ----
reference_genes <- unlist(lapply(list.files(need(IN$watson_reference), pattern = "IG[HKL][VDJ]\\.fasta$", full.names = TRUE),
                                 function(f) unique(alakazam::getGene(names(tigger::readIgFasta(f)), strip_d = FALSE, omit_nl = FALSE))))
imgt_label_map <- c("V-HEPTAMER" = "Heptamer", "V-NONAMER" = "Nonamer", "V-SPACER" = "Spacer", "V-RS" = "rss_aligned",
                    "J-HEPTAMER" = "Heptamer", "J-NONAMER" = "Nonamer", "J-SPACER" = "Spacer", "J-RS" = "rss_aligned",
                    "3'D-HEPTAMER" = "Heptamer", "3'D-NONAMER" = "Nonamer", "3'D-SPACER" = "Spacer", "3'D-RS" = "rss_aligned",
                    "5'D-HEPTAMER" = "Heptamer", "5'D-NONAMER" = "Nonamer", "5'D-SPACER" = "Spacer", "5'D-RS" = "rss_aligned")
imgt_rss <- fread(need(IN$imgt_rss))
imgt_rss[, gene := alakazam::getGene(allele, strip_d = FALSE, omit_nl = FALSE)]
imgt_rss <- imgt_rss[functional %in% c("F", "ORF") & gene %in% reference_genes]
imgt_rss[, `:=`(label_new = imgt_label_map[label], gene_type = substr(allele, 1, 4))]
imgt_rss[, rss_type := fifelse(grepl("3", label), "IGHD_3", fifelse(grepl("5", label), "IGHD_5", gene_type))]
imgt_reshaped <- dcast(imgt_rss, allele + rss_type ~ label_new, value.var = "sequence")[!is.na(Heptamer) & !is.na(Nonamer) & !is.na(Spacer)]
imgt_reshaped[, rss_aligned := paste0(Heptamer, Spacer, Nonamer)]
imgt_reshaped[, length := nchar(rss_aligned)]

# ---- HUSA RSSs ----
husa_rss <- fread(need(file.path(OUT$husa, "husa_rss_filter.tsv")))
husa_rss[, iuis_allele := fifelse(!is.na(husa) & husa != "", husa, vdjbase_allele)]
rss_data <- prepare_aligned_rss_table(husa_rss)
rss_data <- merge(rss_data, unique(husa_rss[, .(allele, husa, vdjbase_allele, iuis_allele, iglabel_allele, iglabel, samples_genomic)]), by = "allele", all.x = TRUE)
for (nm in c("husa", "vdjbase_allele", "iuis_allele", "iglabel_allele", "iglabel", "samples_genomic")) {
  if (nm %in% names(rss_data)) next
  if (paste0(nm, ".x") %in% names(rss_data)) setnames(rss_data, paste0(nm, ".x"), nm) else if (paste0(nm, ".y") %in% names(rss_data)) setnames(rss_data, paste0(nm, ".y"), nm)
}
rss_data[, `:=`(iuis_allele = fifelse(!is.na(husa) & husa != "", husa, vdjbase_allele),
                rss_type = fifelse(!is.na(gene_type_plot) & gene_type_plot != "", gene_type_plot, gene_type),
                gene_type = fifelse(!is.na(gene_type_plot) & gene_type_plot != "", gene_type_plot, gene_type),
                samples = sure_subject, sure_sample_count = sure_subject_count,
                Heptamer = heptamer, Spacer = spacer, Nonamer = nonamer, specie = "human", present = sample_count_AIRRseq > 0)]
rss_data <- rss_data[sure_sample_count > 0]
rss_data[, rss_sequence := paste0(Heptamer, "NNN", Nonamer)]  # the spacer is dropped from the plotted motif
rss_data <- rss_data[nchar(Nonamer) == 9 & nchar(Heptamer) == 7]

# ---- letter grid of the nine most frequent RSSs per segment ----
rss_counts <- rss_data[, .(count = .N), by = .(rss_type, specie, rss_aligned, gene_type, Nonamer, Heptamer, Spacer, rss_sequence)]
rss_counts[, color := ifelse(uniqueN(specie) > 1, "both", "single"), by = .(Nonamer, Heptamer, gene_type)]
rss_human <- rss_counts[specie == "human"]
rss_human[, underline := ifelse(length(unique(gene_type)) > 1, "underline", "no underline"), by = .(Nonamer, Heptamer)]
rss_human <- rss_human[order(gene_type, -count, rss_aligned)]
rss_human[, ig_tag := substr(gene_type, 4, 4)]
rss_human <- rss_human[order(ig_tag, -count)]
# Motif names V1, V2, ... J1, ... D1, ... by descending frequency within V, J, D.
rss_names <- c(); counters <- c(V = 0L, J = 0L, D = 0L)
for (segment in c("V", "J", "D")) for (seq in rss_human[ig_tag == segment, unique(rss_sequence)]) {
  if (seq %in% names(rss_names)) next
  counters[segment] <- counters[segment] + 1L
  rss_names[seq] <- paste0(segment, counters[segment])
}
rss_human[, tag := rss_names[rss_sequence]]
rss_human <- rss_human[order(tag, -count)]
rss_human[, rss_tag := paste0(tag, "-", seq_len(.N)), by = tag]
rss_human[, nseq := .N > 9, by = .(gene_type)]
rss_top <- rss_human[, head(.SD, 9), by = gene_type][order(gene_type, count, rss_sequence)]
rss_top[, order := seq(1, by = 1, length.out = .N), by = .(gene_type)]
rss_top[, count_char := as.character(count)]
# Segments with more than nine RSSs get three blank rows, drawn as an ellipsis in the count column.
for (segment in gene_types) {
  tmp <- rss_top[gene_type == segment]
  if (!nrow(tmp) || !tmp$nseq[1]) next
  blank <- tmp[rep(1, 3)]
  blank[, `:=`(order = seq_len(.N), rss_sequence = gsub("[ATCGN]", " ", rss_sequence), count_char = ".", underline = "no underline",
               color = "single", rss_tag = paste0(".", seq_len(.N)))]
  rss_top[gene_type == segment, order := order + 3]
  rss_top <- rbind(rss_top, blank)
}
grid_keys <- c("gene_type", "specie", "rss_sequence", "order", "color", "underline", "rss_tag")
letters_dt <- rss_top[, .(letter = unlist(strsplit(rss_sequence, ""))), by = grid_keys]
letters_dt[, position := seq_len(.N), by = .(gene_type, rss_tag)]
counts_dt <- rss_top[, .(letter = count_char, position = 1L), by = grid_keys]
tags_dt <- rss_top[, .(letter = rss_tag, position = 1L), by = grid_keys]
tags_dt[, letter_underline := gsub("[.][0-9]+$", "", fifelse(underline == "underline", paste0("*", letter), letter))]  # shared motifs get a star
y_max <- letters_dt[, max(order)]
for (dt in list(letters_dt, counts_dt, tags_dt)) dt[, order_update := order + (y_max - max(order)), by = gene_type]  # bottom-align every panel

# ---- HUSA versus IMGT overlap ----
rss_aligned_counts <- rss_data[, {
  genomic <- trimws(unlist(strsplit(samples_genomic[!is.na(samples_genomic) & nzchar(samples_genomic)], ",", fixed = TRUE)))
  .(count = .N, present = any(present), count_genomic = uniqueN(genomic[nzchar(genomic)]))
}, by = .(rss_aligned, gene_type, rss_type)]
overlap_list <- lapply(setNames(gene_types, gene_types), function(g) {
  merged <- merge(rss_aligned_counts[rss_type == g, .(gene_type, rss_type, present, count_genomic, count, sequence = rss_aligned, length = nchar(rss_aligned))],
                  imgt_reshaped[rss_type == g & length %in% c(39, 28, 38, 27), .(count = .N, gene_type = g, length = nchar(rss_aligned)), by = .(rss_type, rss_aligned)],
                  by.x = c("gene_type", "rss_type", "sequence", "length"), by.y = c("gene_type", "rss_type", "rss_aligned", "length"), all = TRUE, suffixes = c("_husa", "_imgt"))
  merged[is.na(merged)] <- 0
  merged[, `:=`(HUSA = count_husa > 0, IMGT = count_imgt > 0, present = present == 1, above_zero = count_genomic > 0)][]
})
overlap <- rbindlist(overlap_list)
overlap_summary <- data.table(gene_type = gene_types,
                              total_imgt = vapply(overlap_list, function(x) sum(x$IMGT), integer(1)),
                              total_husa = vapply(overlap_list, function(x) sum(x$HUSA), integer(1)),
                              total_both = vapply(overlap_list, function(x) sum(x$IMGT & x$HUSA), integer(1)),
                              total_neither = vapply(overlap_list, function(x) sum(!x$IMGT & !x$HUSA), integer(1)))

# ---- agreement with the segment consensus heptamer and nonamer ----
consensus_dt <- rss_data[, .(consensus_heptamer = Biostrings::consensusString(Biostrings::DNAStringSet(Heptamer), ambiguityMap = "?"),
                             consensus_nonamer = Biostrings::consensusString(Biostrings::DNAStringSet(Nonamer), ambiguityMap = "?")), by = .(rss_type, gene_type)]
consensus_dt[, `:=`(consensus_heptamer_regex = gsub("?", "[ATCG]", consensus_heptamer, fixed = TRUE),
                    consensus_nonamer_regex = gsub("?", "[ATCG]", consensus_nonamer, fixed = TRUE))]
rss_data[consensus_dt, `:=`(consensus_heptamer_regex = i.consensus_heptamer_regex, consensus_nonamer_regex = i.consensus_nonamer_regex), on = .(rss_type, gene_type)]
rss_data[, `:=`(consensus_heptamer = mapply(grepl, consensus_heptamer_regex, Heptamer, USE.NAMES = FALSE),
                consensus_nonamer = mapply(grepl, consensus_nonamer_regex, Nonamer, USE.NAMES = FALSE))]
rss_data[, consensus_heptamer_nonamer := consensus_heptamer & consensus_nonamer]
cons_label <- function(x) ifelse(all(x), "Consensus", ifelse(all(!x), "Non-Consensus", "Both"))
rss_data[, consensus_heptamer_nonamer_alele_label := cons_label(consensus_heptamer_nonamer), by = .(gene_type, allele)]
rss_data[, consensus_heptamer_nonamer_rss_label := cons_label(consensus_heptamer_nonamer), by = .(gene_type, rss_aligned)]
unique_counts <- rbind(
  rss_data[, .(variable = "unique_allele_count", value = uniqueN(allele)), by = .(gene_type, consensus_heptamer_nonamer_label = consensus_heptamer_nonamer_alele_label)],
  rss_data[, .(variable = "unique_rss_count", value = uniqueN(rss_aligned)), by = .(gene_type, consensus_heptamer_nonamer_label = consensus_heptamer_nonamer_rss_label)])
setcolorder(unique_counts, c("gene_type", "variable", "value", "consensus_heptamer_nonamer_label"))

# Per-panel layout: the spacer length label and where it sits.
panel_layout <- rbindlist(lapply(gene_types, function(g) {
  segment_rss <- rss_data[gene_type == g & specie == "human"]
  short_spacer <- if (g %in% c("IGHD_3", "IGHD_5", "IGKV", "IGLJ")) 12L else 23L
  panel_letters <- letters_dt[gene_type == g]
  truncated <- nrow(rss_human[gene_type == g]) > 9
  data.table(gene_type = g,
             spacer_label = if (any(nchar(segment_rss$Spacer) < short_spacer)) paste0(short_spacer - 1L, "/\n", short_spacer) else as.character(short_spacer),
             spacer_label_y = if (truncated) max(panel_letters$order_update) / 2 else max(panel_letters$order_update) - uniqueN(panel_letters$rss_tag) / 2,
             y_max = y_max)
}))

fwrite(rss_data, tables[["aligned"]])
fwrite(rss_data[specie == "human", .(gene_type, rss_sequence)], tables[["seqlogo"]])
fwrite(letters_dt, tables[["letters"]]); fwrite(counts_dt, tables[["counts"]]); fwrite(tags_dt, tables[["tags"]])
fwrite(overlap, tables[["overlap"]]); fwrite(overlap_summary, tables[["summary"]])
fwrite(unique_counts, tables[["unique"]]); fwrite(panel_layout, tables[["layout"]])
cat(sprintf("supp2: %d RSS rows, %d segments, %d overlap rows\n", nrow(rss_data), length(gene_types), nrow(overlap)))
}

# ---- draw ----
seqlogo_dt <- fread(tables[["seqlogo"]]); letters_dt <- fread(tables[["letters"]]); counts_dt <- fread(tables[["counts"]]); tags_dt <- fread(tables[["tags"]])
overlap <- fread(tables[["overlap"]]); unique_counts <- fread(tables[["unique"]]); panel_layout <- fread(tables[["layout"]])
y_max <- panel_layout$y_max[1]
unique_counts[, consensus_heptamer_nonamer_label := factor(consensus_heptamer_nonamer_label, levels = c("Consensus", "Both", "Non-Consensus"))]
letter_colors <- c("single" = "#000000", "both" = "#000000", "underline" = "#000000")
support_colors <- setNames(c("#0072B2", "#aa0415ff"), c("TRUE", "FALSE"))
consensus_colors <- setNames(c("#B8860B", "#008B8B", "#A63D40"), c("Consensus", "Both", "Non-Consensus"))
bare_axes <- theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
                   axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.ticks.length.x = unit(0, "lines"), panel.spacing = unit(0, "cm"))

rss_panels <- lapply(gene_types, function(g) {
  title <- if (grepl("D", g)) { parts <- strsplit(g, "_")[[1]]; paste0(parts[2], "' ", parts[1]) } else g  # IGHD_5 reads as "5' IGHD"
  p_tag <- ggplot() +
    geom_text_repel(data = tags_dt[gene_type == g], mapping = aes(as.numeric(position), order_update, label = letter_underline, color = color, hjust = -.1),
                    size = 8, bg.r = .1, force = 0, parse = FALSE) +
    scale_color_manual(values = letter_colors) + theme_logo() + bare_axes + theme(plot.margin = margin(0, 0, 0.1, 0, "cm")) +
    labs(y = title) + theme(axis.title.y = element_text(size = 34, angle = 90, hjust = 1))
  p_sequences <- ggplot() +
    geom_text(data = letters_dt[gene_type == g], mapping = aes(x = as.numeric(position), y = order_update, label = gsub("N", " ", letter), color = color,
                                                              size = ifelse(letter == ".", 16, 8))) +
    annotate("text", x = 9, y = panel_layout[gene_type == g, spacer_label_y], label = panel_layout[gene_type == g, spacer_label], size = 16) +
    scale_x_continuous(breaks = 1:17, expand = c(0.03, 0)) + scale_color_manual(values = letter_colors) + scale_size_identity() +
    theme_logo() + bare_axes + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
  p_count <- ggplot() +
    geom_text_repel(data = counts_dt[gene_type == g], mapping = aes(as.numeric(position), order_update, label = letter, color = color, hjust = 0,
                                                                   size = ifelse(letter == ".", 16, 8)), bg.r = .1, force = 0) +
    scale_color_manual(values = letter_colors) + scale_size_identity() + theme_logo() + bare_axes + theme(plot.margin = margin(0, -0.5, 0, 0, "cm"))
  p_seqlogo <- ggplot() + geom_logo(seqlogo_dt[gene_type == g, rss_sequence], method = "probability", seq_type = "dna") + theme_logo() +
    scale_x_continuous(breaks = 1:19, expand = c(0, 0)) + labs(title = title) +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
          axis.ticks.length.x = unit(0, "lines"), plot.title = element_text(hjust = .5, size = 40), plot.margin = margin(0, 0, 0, 0, "pt"))
  # J segments have far fewer RSSs than the tallest panel, so their grids are pinned to the shared y range.
  if (grepl("J", g)) {
    pinned <- coord_cartesian(expand = TRUE, ylim = c(1, y_max))
    p_tag <- p_tag + pinned; p_sequences <- p_sequences + pinned; p_count <- p_count + pinned
  }
  wrap_elements((p_seqlogo + labs(title = "")) + p_tag + p_sequences + p_count +
                  plot_layout(design = "#A#\nBCD", widths = c(1 / 3, 1, 1 / 5), heights = c(1 / 2, 1)))
})
unique_panels <- lapply(gene_types, function(g) {
  p <- ggplot(unique_counts[gene_type == g], aes(x = variable, y = value, fill = consensus_heptamer_nonamer_label)) + geom_col() +
    labs(x = "", y = "Unique Count", fill = "") + scale_x_discrete(labels = c("Alleles", "RSS")) +
    scale_y_continuous(breaks = function(x) pretty(x)[pretty(x) %% 1 == 0]) +
    scale_fill_manual(values = consensus_colors, name = "Heptamer-Nonamer") + theme_pubclean(base_size = 28) + theme(axis.text.x = element_text(angle = 0))
  if (g == gene_types[1]) p + theme(legend.position = "inside", legend.position.inside = c(0.9, 0.8), legend.direction = "vertical", legend.box.background = element_blank())
  else p + theme(legend.position = "none")
})
imgt_panels <- lapply(gene_types, function(g) {
  upset(overlap[rss_type == g], c("IMGT", "HUSA"), name = "",
        annotations = list("# Unique\nSamples" = (
          ggplot(mapping = aes(x = intersection, y = count_genomic, color = present, shape = above_zero)) +
            geom_boxplot(na.rm = TRUE, outlier.shape = NA) + geom_point(position = position_jitterdodge(seed = 42), size = 4) +
            scale_shape_manual(values = c(17, 16)) + scale_color_manual(values = support_colors) + guides(color = "none", shape = "none"))),
        base_annotations = list("Intersection size" = intersection_size(counts = FALSE, text = list(size = 14), text_colors = c(on_background = "black", on_bar = "black"),
                                                                        bar_number_threshold = 0.65, mapping = aes(fill = present)) +
                                  labs(y = "# Unique\nRSSs") + guides(fill = if (g == gene_types[1]) guide_legend(title = "Expressed") else "none") +
                                  scale_fill_manual(values = support_colors)),
        set_sizes = (upset_set_size() + scale_y_reverse(breaks = function(x) pretty(x)[pretty(x) %% 1 == 0])),
        guides = "over", themes = upset_default_themes(text = element_text(size = 28)))
})

row_heights <- c(1, 1, 1, 0.95, 0.95, 0.95, 1)
final <- wrap_elements(wrap_plots(rss_panels, ncol = 1, heights = row_heights) & theme(plot.margin = margin(-5, 0, 0, 0, "lines"))) +
  wrap_elements(wrap_plots(unique_panels, ncol = 1, heights = row_heights) & theme(plot.margin = margin(2, 0, 0, 0, "lines"))) +
  plot_spacer() +
  wrap_elements(wrap_plots(imgt_panels, ncol = 1, heights = row_heights) & theme(plot.margin = margin(2, 0, 0, 0, "lines"))) +
  plot_layout(widths = c(1, 0.7, 0.1, 1)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 70))
ggsave(file.path(OUT$figures, "supp2.pdf"), final, width = 30, height = 40, units = "in", dpi = 800, limitsize = FALSE)
