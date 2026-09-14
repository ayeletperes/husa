#!/usr/bin/env Rscript
# Figure 4. A: IGHV gene structure schematic. B: IGHV RSS heatmap over the IGHV leader
# heatmap, sharing IUIS-group columns and per-group usage boxplots. C: pairwise coding
# distance against pairwise RSS / leader distance per V locus, with Mantel statistics.
#
# Seeds are the ones the published figure used: 123 before each optimize_order_enhanced()
# call (it searches dendrogram flips stochastically) and 1234 before the Mantel block.

source("R/00_setup.R")
source("R/lib/figure56_helpers.R")
suppressPackageStartupMessages(library(ggplot2))

rss_leader <- fread(need(file.path(OUT$rss_leader, "rss_leader_iuis_data.csv.gz")))
husa_tsv <- need(file.path(OUT$husa, "husa.tsv"))
repertoire_igh <- need(file.path(OUT$repertoire, "gg_repertoire_data_IGH_genotype_corrected.csv.gz"))
permutations <- 4000L

# ---- A: schematic geometry (not to genomic scale) ----
schematic_colors <- c(Leader = "#C5B3E6", FR = "#9FA8DA", CDR = "#BBDEFB", Heptamer = "#FFB74D", Spacer = "#FFE0B2", Nonamer = "#FFA726")
zoom_y <- 0.20; block_h <- 0.24
block <- function(xmin, xmax, fill, label) data.table(xmin = xmin, xmax = xmax, ymin = zoom_y - block_h / 2, ymax = zoom_y + block_h / 2,
                                                       fill = fill, label = label, xmid = (xmin + xmax) / 2, y = zoom_y)
blocks <- rbindlist(list(
  block(1.0, 2.0, "Leader", "Exon 1"), block(3.2, 4.2, "Leader", "Exon 2"),
  block(4.2, 4.9, "FR", "FR1"), block(4.9, 5.7, "CDR", "CDR1"), block(5.7, 6.6, "FR", "FR2"), block(6.6, 7.4, "CDR", "CDR2"),
  block(7.4, 8.3, "FR", "FR3"), block(8.3, 9.1, "CDR", "CDR3"), block(9.1, 9.8, "FR", "FR4"),
  block(9.8, 11.0, "Heptamer", "Heptamer\n(7 bp)"), block(11.0, 12.2, "Spacer", "Spacer\n(23 bp)"), block(12.2, 13.4, "Nonamer", "Nonamer\n(9 bp)")))
blocks[, `:=`(layer = "block", colour = unname(schematic_colors[fill]), size = 11, fontface = "plain")]
regions <- data.table(xmin = c(1.0, 4.2, 9.8), xmax = c(4.2, 9.8, 13.4), y = 0.46, label = c("Leader", "Coding region", "RSS"))
regions[, `:=`(layer = "region", xmid = (xmin + xmax) / 2, size = 12, fontface = "bold")]
marks <- data.table(layer = c("backbone", "end_label", "end_label", "title"), xmin = c(0.8, NA, NA, NA), xmax = c(13.9, NA, NA, NA),
                    xmid = c(NA, 0.55, 14.05, 7.2), y = c(zoom_y, zoom_y, zoom_y, 0.75),
                    label = c(NA, "5'", "3'", "Schematics of IGHV gene structure"), size = c(NA, 10, 10, 15), fontface = c(NA, "plain", "plain", "bold"))
schematic <- rbindlist(list(blocks, regions, marks), use.names = TRUE, fill = TRUE)

# ---- B: IGHV RSS and leader heatmaps ----
genes_order_v <- gsub("IG[HKL]", "", fread(need(IN$gene_bed))[grepl("IGHV", V4)]$V4)
rss_v <- rss_leader[gene_type == "IGHV" & present == TRUE]
nonamer_start <- 30L  # heptamer (7) + 23 nt spacer, so the nonamer gap goes in at column 30
leader2_start <- max(nchar(rss_v$l_part1_aligned))

# Column tree: cluster coding sequences allele by allele, then collapse each IUIS group to one tip.
coding <- unique(rss_v[, .(vdjbase_allele, iuis_group, seq)])
tree <- ape::as.phylo(hclust(stringdist::stringdistmatrix(setNames(coding$seq, coding$vdjbase_allele), method = "lv", useNames = "names"), method = "complete"))
collapsed_tree <- Reduce(function(tr, node) collapse_clade(tr, node), tapply(coding$vdjbase_allele, coding$iuis_group, identity), init = tree)
collapsed_tree$tip.label <- coding$iuis_group[match(collapsed_tree$tip.label, coding$vdjbase_allele)]
hc_cols <- ape::as.hclust.phylo(collapsed_tree)
dend_cols <- as.dendrogram(hc_cols)

seqs_rss <- unique(rss_v$rss_aligned)
hc_rows_rss <- hclust(stringdist::stringdistmatrix(seqs_rss, method = "hamming", useNames = "strings"), method = "complete")
mat_rss <- table(rss_v$rss_aligned, rss_v$iuis_group)[hc_rows_rss$labels[hc_rows_rss$order], hc_cols$labels[hc_cols$order]]
set.seed(123)
result_rss <- optimize_order_enhanced(mat = mat_rss, row_dend = as.dendrogram(hc_rows_rss), col_dend = dend_cols,
                                      kernel_size = 3, max_iter_family = 50, max_iter_global = 200, max_steps = 3, verbose = FALSE,
                                      consider_mixed_flips = FALSE, consider_minimal_subtree = FALSE, consider_subgroups = FALSE,
                                      optimize_subgroups = TRUE, use_cache = TRUE)

# Optimised subgroup blocks, genes in chromosome order inside each block (absent from the bed sort last).
cols_dt <- data.table(label = labels(result_rss$col_dend))
cols_dt[, subgroup := gsub("D", "", sub("-.*$", "", label))]
cols_dt[, chr_order := match(label, genes_order_v)]
cols_dt[, subgroup_rank := match(subgroup, rev(c("V5", "V2", "V7", "V1", "V4", "V6", "V3")))]
setorder(cols_dt, subgroup_rank, chr_order, na.last = TRUE)

# Genes whose alleles cannot be told apart in the repertoire share one usage boxplot.
group_map <- build_allele_group_map(husa_path = husa_tsv, locus = "IGHV", keep_labels = cols_dt$label)
usage <- build_usage_objects(repertoire_csv = repertoire_igh, locus = "IGHV", call_col = "v_call_new", group_map = group_map, keep_labels = cols_dt$label)
usage_map <- copy(usage$map)[, merged := gsub("IG[KLH]", "", merged)]
cols_dt[, split_group := usage_map[match(cols_dt$label, iuis_group)]$merged]

hc_rss <- as.hclust(result_rss$row_dend)
rss_subgroup <- rss_v[, .(subgroup = paste0(unique(iuis_subgroup), collapse = ",")), by = rss_aligned][order(match(rss_aligned, hc_rss$labels))]
consensus_rss <- generate_consensus(seqs_rss, vec_onehot)

leader_v <- rss_leader[gene_type == "IGHV" & present == TRUE & !is.na(leader) & nzchar(leader)]
leader_seqs <- unique(leader_v$leader)
dist_method <- if (uniqueN(nchar(leader_seqs)) == 1L) "hamming" else "lv"
hc_rows_leader <- hclust(stringdist::stringdistmatrix(leader_seqs, method = dist_method, useNames = "strings"), method = "complete")
mat_leader <- table(leader_v$leader, leader_v$iuis_group)
missing_cols <- setdiff(cols_dt$label, colnames(mat_leader))
if (length(missing_cols)) mat_leader <- cbind(mat_leader, matrix(0L, nrow = nrow(mat_leader), ncol = length(missing_cols), dimnames = list(rownames(mat_leader), missing_cols)))
mat_leader <- mat_leader[hc_rows_leader$labels[hc_rows_leader$order], cols_dt$label, drop = FALSE]
set.seed(123)
result_leader <- optimize_order_enhanced(mat = mat_leader, row_dend = as.dendrogram(hc_rows_leader), col_dend = dend_cols,
                                         kernel_size = 3, max_iter_family = 100, max_iter_global = 200, max_steps = 4, verbose = FALSE,
                                         consider_mixed_flips = TRUE, consider_minimal_subtree = FALSE, consider_subgroups = TRUE,
                                         optimize_subgroups = TRUE, use_cache = TRUE)
hc_leader <- as.hclust(result_leader$row_dend)
leader_subgroup <- leader_v[, .(subgroup = paste0(unique(iuis_subgroup), collapse = ",")), by = leader][order(match(leader, hc_leader$labels))]
consensus_leader <- generate_consensus(leader_seqs, vec_onehot)

dense_counts <- function(dt, row_col, rows, cols) {
  counts <- dt[, .N, by = c(row_col, "iuis_group")]
  setnames(counts, c(row_col, "iuis_group", "N"), c("row_seq", "column_label", "count"))
  dense <- CJ(row_seq = rows, column_label = cols, unique = TRUE)
  dense[counts, count := i.count, on = .(row_seq, column_label)]
  dense[is.na(count), count := 0L][]
}
matrix_dt <- rbind(data.table(layer = "rss", dense_counts(rss_v, "rss_aligned", hc_rss$labels, cols_dt$label)),
                   data.table(layer = "leader", dense_counts(leader_v, "leader", hc_leader$labels, cols_dt$label)))
rows_dt <- rbind(data.table(layer = "rss", position = seq_along(hc_rss$labels), row_seq = hc_rss$labels, annotation_subgroup = rss_subgroup$subgroup),
                 data.table(layer = "leader", position = seq_along(hc_leader$labels), row_seq = hc_leader$labels, annotation_subgroup = leader_subgroup$subgroup))
m_box <- usage$m_box
usage_dt <- data.table(subject_index = rep(seq_len(nrow(m_box)), times = ncol(m_box)), merged = rep(colnames(m_box), each = nrow(m_box)), rel_usage = as.numeric(m_box))

# ---- C: coding distance against RSS / leader distance. Coding is Levenshtein; RSS is Hamming
# when every aligned RSS shares one length, else Levenshtein; leader is Levenshtein. Mantel is
# the Spearman correlation of the upper-triangle distances with a permutation p. ----
cophy <- rss_leader[present == TRUE & !is.na(seq) & nzchar(seq)]
lv_dist <- function(seqs) as.matrix(stringdist::stringdistmatrix(seqs, seqs, method = "lv"))
mantel <- function(mat_a, mat_b, np, idx) {
  a <- mat_a[idx]; b <- mat_b[idx]
  if (uniqueN(a) < 2L || uniqueN(b) < 2L) return(c(r = NA_real_, p = NA_real_))
  ra <- rank(a, ties.method = "average"); rb_c <- rank(b, ties.method = "average"); rb_c <- rb_c - mean(rb_c); ra_c <- ra - mean(ra)
  denom <- sqrt(sum(ra_c^2) * sum(rb_c^2))
  if (!is.finite(denom) || denom == 0) return(c(r = NA_real_, p = NA_real_))
  r_obs <- sum(ra_c * rb_c) / denom
  n <- nrow(mat_a); rank_a <- matrix(NA_real_, n, n)
  rank_a[cbind(idx[, 1L], idx[, 2L])] <- ra; rank_a[cbind(idx[, 2L], idx[, 1L])] <- ra
  ra_mean <- mean(ra); ii <- idx[, 1L]; jj <- idx[, 2L]; perm_ge <- 0L
  for (k in seq_len(np)) {
    o <- sample.int(n)
    perm_ge <- perm_ge + (sum((rank_a[cbind(o[ii], o[jj])] - ra_mean) * rb_c) / denom >= r_obs)
  }
  c(r = r_obs, p = perm_ge / np)
}
run_element <- function(seg, element, value_col, keep_pairs) {
  dt <- cophy[gene_type == seg]
  if (!nrow(dt)) return(NULL)
  sub <- unique(dt[!is.na(get(value_col)) & nzchar(get(value_col)), .(vdjbase_allele, seq, val = get(value_col))])
  n <- nrow(sub)
  if (n < 4L) return(list(summary = data.table(gene_type = seg, element = element, n_alleles = n, coding_dist = NA_character_, element_dist = NA_character_,
                                               mantel_r = NA_real_, mantel_p = NA_real_), pairs = NULL))
  idx <- which(upper.tri(matrix(FALSE, n, n)), arr.ind = TRUE)
  cod <- lv_dist(sub$seq)
  method <- if (element == "RSS" && uniqueN(nchar(sub$val)) == 1L) "hamming" else "lv"
  ne <- as.matrix(stringdist::stringdistmatrix(sub$val, sub$val, method = method))
  res <- mantel(cod, ne, permutations, idx)
  list(summary = data.table(gene_type = seg, element = element, n_alleles = n, coding_dist = "lv", element_dist = method, mantel_r = unname(res["r"]), mantel_p = unname(res["p"])),
       pairs = if (keep_pairs) data.table(gene_type = seg, element = element, coding_dist = cod[idx], element_dist = ne[idx]))
}
set.seed(1234)
locus_levels <- c("IGHV", "IGKV", "IGLV"); element_levels <- c("RSS", "Leader")
v_elements <- c(RSS = "rss_aligned", Leader = "leader")
jd_segments <- intersect(c("IGHJ", "IGKJ", "IGLJ", "IGHD_5", "IGHD_3"), unique(cophy$gene_type))
v_results <- list()
for (seg in locus_levels[locus_levels %in% unique(cophy$gene_type)]) for (el in element_levels) v_results[[paste(seg, el)]] <- run_element(seg, el, v_elements[[el]], TRUE)
# J and D are not drawn but belong to the FDR family, which is corrected once over every test.
jd_results <- lapply(jd_segments, run_element, element = "RSS", value_col = "rss_aligned", keep_pairs = FALSE)
mantel_dt <- rbindlist(lapply(c(v_results, jd_results), `[[`, "summary"), use.names = TRUE)
mantel_dt[, mantel_q := p.adjust(mantel_p, method = "BH")]
pairs_dt <- rbindlist(Filter(Negate(is.null), lapply(v_results, `[[`, "pairs")), use.names = TRUE)
trend_dt <- pairs_dt[, { xr <- range(coding_dist, na.rm = TRUE); fit <- lm(element_dist ~ coding_dist)
  data.table(coding_dist = xr, element_dist = predict(fit, newdata = data.frame(coding_dist = xr))) }, by = .(gene_type, element)]
ranges <- pairs_dt[, .(x = 0.02 * max(coding_dist), y = 0.97 * max(element_dist)), by = .(gene_type, element)]
mantel_dt[ranges, `:=`(x = i.x, y = i.y), on = .(gene_type, element)]
mantel_dt[!is.na(x), label := sprintf("Mantel r = %.2f, %s", mantel_r, fifelse(is.na(mantel_q), "q = NA", fifelse(mantel_q < 0.001, "q < 0.001", sprintf("q = %.3f", mantel_q))))]
print(mantel_dt[, .(gene_type, element, n_alleles, mantel_r, mantel_p, mantel_q)])

fwrite(schematic, file.path(OUT$source, "figure4_schematic.csv"))
fwrite(matrix_dt, file.path(OUT$source, "figure4_heatmap_matrix.csv.gz"))
fwrite(rows_dt, file.path(OUT$source, "figure4_heatmap_rows.csv.gz"))
fwrite(data.table(position = seq_len(nrow(cols_dt)), column_label = cols_dt$label, column_subgroup = cols_dt$subgroup, split_group = cols_dt$split_group),
       file.path(OUT$source, "figure4_heatmap_columns.csv"))
fwrite(data.table(layer = c("rss", "leader"), consensus = c(consensus_rss$consensus, consensus_leader$consensus), boundary = c(nonamer_start, leader2_start), rss = c(TRUE, FALSE)),
       file.path(OUT$source, "figure4_consensus.csv"))
fwrite(usage_dt, file.path(OUT$source, "figure4_usage.csv.gz"))
fwrite(pairs_dt, file.path(OUT$source, "figure4_cophylogeny_pairs.csv.gz"))
fwrite(mantel_dt, file.path(OUT$source, "figure4_cophylogeny_mantel.csv"))
fwrite(trend_dt, file.path(OUT$source, "figure4_cophylogeny_trend.csv"))

# ---- draw ----
p_schematic <- ggplot() +
  geom_text(data = marks[layer == "title"], aes(x = xmid, y = y, label = label, size = size, fontface = fontface)) +
  geom_segment(data = marks[layer == "backbone"], aes(x = xmin, xend = xmax, y = y, yend = y), linewidth = 0.5) +
  geom_text(data = marks[layer == "end_label"], aes(x = xmid, y = y, label = label, size = size)) +
  geom_rect(data = blocks, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = colour), color = "black", linewidth = 0.3) +
  geom_text(data = blocks, aes(x = xmid, y = y, label = label, size = size)) +
  geom_segment(data = regions, aes(x = xmin, xend = xmax, y = y, yend = y), linewidth = 0.5) +
  geom_segment(data = regions, aes(x = xmin, xend = xmin, y = y, yend = y - 0.05), linewidth = 0.5) +
  geom_segment(data = regions, aes(x = xmax, xend = xmax, y = y, yend = y - 0.05), linewidth = 0.5) +
  geom_text(data = regions, aes(x = xmid, y = y + 0.08, label = label, size = size, fontface = fontface)) +
  scale_fill_identity() + scale_size_identity() +
  coord_cartesian(xlim = c(0, 14.3), ylim = c(0.02, 0.85), clip = "off") + theme_void() +
  theme(legend.position = "none", plot.margin = margin(5, 5, 5, 5))

count_matrix <- function(dt, rows) {
  wide <- dcast(dt, row_seq ~ column_label, value.var = "count")
  m <- as.matrix(wide[, -"row_seq"]); rownames(m) <- wide$row_seq
  m[rows, cols_dt$label, drop = FALSE]
}
subgroups <- unique(c(cols_dt$subgroup, rss_subgroup$subgroup, leader_subgroup$subgroup))
graphics <- setNames(lapply(subgroups, subgroups_annotation_custom, colors_rss_subgroup = colors15), subgroups)
column_ha <- HeatmapAnnotation(Subgroup = anno_customize(setNames(cols_dt$subgroup, cols_dt$label), graphics = graphics),
                               height = unit(20, "mm"), annotation_name_gp = gpar(fontsize = 22))
ha_usage <- build_collapsed_box_anno(m_box = m_box, align_to_vec = cols_dt$split_group, ylab = "Usage", show_axis_at = 1L,
                                     axis_fontsize = 30, ylab_fontsize = 36, ylab_x = -38, ylab_x_npc = FALSE)
column_split <- factor(cols_dt$split_group, levels = unique(cols_dt$split_group))
heatmap_of <- function(mat, hc, row_labels, consensus, boundary, rss, name, top, bottom, names_side) {
  Heatmap(mat, cluster_rows = hc, cluster_columns = FALSE, column_order = cols_dt$label, show_row_names = FALSE, show_column_names = TRUE, col = col_fun,
          column_title = "", column_title_gp = gpar(fontsize = 20, col = "black"), border = FALSE, show_heatmap_legend = FALSE,
          row_dend_side = "right", row_names_side = "left", row_names_max_width = unit(16, "cm"),
          row_names_gp = gpar(fontsize = 22, fontfamily = "mono"), column_names_gp = gpar(fontsize = 32, fontfamily = "mono"), column_names_side = names_side,
          row_dend_width = unit(4, "cm"), column_dend_height = unit(4, "cm"), rect_gp = gpar(col = "gray80", lwd = 2),
          left_annotation = do.call(rowAnnotation, c(setNames(list(row_anno_text(gt_render(create_row_labels(hc$labels, consensus, nucleotide_colors, boundary, rss = rss, fontsize_px = 32)))), name),
                                                     list(show_annotation_name = FALSE))),
          right_annotation = rowAnnotation(Subgroup = anno_customize(row_labels, graphics = graphics), width = unit(10, "mm"), show_annotation_name = FALSE),
          top_annotation = top, bottom_annotation = bottom, column_split = column_split)
}
ht_rss <- heatmap_of(count_matrix(matrix_dt[layer == "rss"], hc_rss$labels), hc_rss, rss_subgroup$subgroup, consensus_rss$consensus, nonamer_start, TRUE,
                     "labels_rss_IGHV", column_ha, NULL, "top")
ht_leader <- heatmap_of(count_matrix(matrix_dt[layer == "leader"], hc_leader$labels), hc_leader, leader_subgroup$subgroup, consensus_leader$consensus, leader2_start, FALSE,
                        "labels_leader_IGHV", ha_usage, column_ha, "bottom")
# Row-side titles take the column-title style rotated; Heatmap() exposes no argument for it.
title_param <- ht_rss@column_title_param
ht_rss@column_title <- ""
title_param$rot <- 270; title_param$just <- c(0.5, 0.3)
ht_rss@row_title <- "RSS"; ht_rss@row_title_param <- title_param
ht_leader@row_title <- "Leader"; ht_leader@row_title_param <- title_param
consensus_rss_label <- create_consensus_label(consensus_rss$consensus, nucleotide_colors, nonamer_start, fontsize_px = 32)
consensus_leader_label <- create_consensus_label(consensus_leader$consensus, nucleotide_colors, leader2_start, rss = FALSE, fontsize_px = 32)
lgd <- Legend(col_fun = col_fun, title = "Count", direction = "horizontal", labels_gp = gpar(fontsize = 28), title_gp = gpar(fontsize = 32, fontface = "bold"),
              legend_width = unit(12, "cm"), grid_height = unit(1, "cm"), grid_width = unit(1, "cm"))

pairs_dt <- pairs_dt[gene_type %in% locus_levels]; mantel_dt <- mantel_dt[gene_type %in% locus_levels]; trend_dt <- trend_dt[gene_type %in% locus_levels]
for (dt in list(pairs_dt, mantel_dt, trend_dt)) {
  set(dt, j = "gene_type", value = factor(dt$gene_type, levels = locus_levels))
  set(dt, j = "element", value = factor(dt$element, levels = element_levels))
}
p_grid <- ggplot(pairs_dt, aes(x = coding_dist, y = element_dist)) +
  geom_hex(bins = 30) +
  geom_line(data = trend_dt, aes(x = coding_dist, y = element_dist), inherit.aes = FALSE, color = "#222222", linewidth = 0.6) +
  geom_text(data = mantel_dt[!is.na(label)], aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 8, color = "#222222") +
  scale_fill_gradient(name = "Pair count", low = "#e8eef5", high = "#2b5d8a", trans = "log10", labels = scales::label_number()) +
  facet_grid(element ~ gene_type, scales = "free") +
  labs(x = "Pairwise coding-sequence distance", y = "Pairwise non-coding distance") +
  ggpubr::theme_pubclean(base_size = 34) +
  theme(legend.position = "inside", legend.direction = "horizontal", legend.position.inside = c(0.13, 1.13), legend.key.width = unit(2.2, "cm"), plot.title = element_blank())

width <- 42; height <- 50
row_heights <- c(schematic = 0.5, heatmap = 3.0, cophylogeny = 1.9)
g_schematic <- grid.grabExpr(print(p_schematic), width = width, height = height * row_heights[["schematic"]] / sum(row_heights))
grDevices::cairo_pdf(file.path(OUT$figures, "figure4.pdf"), width = width, height = height)
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 1, heights = unit(row_heights, "null"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(g_schematic)
grid.text("A", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 80, fontface = "bold", fontfamily = "sans"))
popViewport()
# ComplexHeatmap is drawn straight into the row: grabbing it collapses the column-split slices.
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1, y = unit(0.95, "npc"), height = unit(0.95, "npc"), just = c("center", "top")))
draw(ht_rss %v% ht_leader, newpage = FALSE)
draw(lgd, x = unit(0.05, "npc"), y = unit(0.98, "npc"), just = c("left", "top"))
decorate_annotation("labels_rss_IGHV", grid.draw(richtext_grob(consensus_rss_label, y = unit(-4, "mm"), x = unit(0, "npc"), hjust = 0)))
decorate_annotation("labels_leader_IGHV", grid.draw(richtext_grob(consensus_leader_label, y = unit(-4, "mm"), x = unit(0, "npc"), hjust = 0)))
grid.text("B", x = unit(0, "npc"), y = unit(0.97, "npc"), just = c("left", "top"), gp = gpar(fontsize = 80, fontface = "bold", fontfamily = "sans"))
popViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1, width = unit(0.98, "npc"), height = unit(0.96, "npc")))
grid.draw(ggplotGrob(p_grid))
grid.text("C", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 42, fontface = "bold", fontfamily = "sans"))
popViewport()
dev.off()
