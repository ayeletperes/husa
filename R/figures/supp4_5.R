#!/usr/bin/env Rscript
# Supplementary figures 4 and 5: RSS (and leader) count heatmaps, unique aligned sequence by
# IUIS group, rows and columns ordered by the blockiness optimiser (seed 123 before each call).
#   supp4  IGKV and IGLV: an RSS heatmap stacked over a leader heatmap.
#   supp5  IGHD 5'/3' and IGHJ / IGKJ / IGLJ: RSS heatmaps only.

source("R/00_setup.R")
source("R/lib/figure56_helpers.R")

rss_leader <- fread(need(file.path(OUT$rss_leader, "rss_leader_iuis_data.csv.gz")))
husa_tsv <- need(file.path(OUT$husa, "husa.tsv"))
gene_order <- fread(need(IN$gene_bed))
setnames(gene_order, c("ch", "s", "e", "g"))

# Per panel: the column where the nonamer gap is drawn, and the hand-fixed subgroup / gene order
# along the x axis (kept back-to-front where the original wrote it that way).
panels <- list(
  IGKV = list(figure = "supp4", bed_tag = "IGKV", nonamer_start = 19, leader = TRUE, subgroup_from_group = FALSE, strip_d = TRUE,
              order_by = "subgroup", order_levels = rev(c("V5", "V3", "V6", "V4", "V2", "V1")),
              usage_locus = "IGKV", chain = "IGK", call_col = "v_call_new", allow_cartesian = FALSE),
  IGLV = list(figure = "supp4", bed_tag = "IGLV", nonamer_start = 30, leader = TRUE, subgroup_from_group = FALSE, strip_d = TRUE,
              order_by = "subgroup", order_levels = rev(c("V9", "V7", "V5", "V4", "V3", "V8", "V2", "V10", "V6", "V1")),
              usage_locus = "IGLV", chain = "IGL", call_col = "v_call_new", allow_cartesian = FALSE),
  IGHD_5 = list(figure = "supp5", bed_tag = "IGHD", nonamer_start = 19, leader = FALSE, subgroup_from_group = TRUE, strip_d = FALSE, merge_asc = TRUE,
                order_by = "subgroup", order_levels = c("D4", "D5", "D7", "D1", "D2", "D3", "D6"),
                usage_locus = "IGHD", chain = "IGH", call_col = "d_call_new", allow_cartesian = FALSE),
  # The 3' panel sits under the 5' one, so it takes the 5' column order and its rows start on the D4 block.
  IGHD_3 = list(figure = "supp5", bed_tag = "IGHD", nonamer_start = 19, leader = FALSE, subgroup_from_group = TRUE, strip_d = FALSE, merge_asc = TRUE,
                columns_from = "IGHD_5", rotate_subgroup_first = "D4",
                usage_locus = "IGHD", chain = "IGH", call_col = "d_call_new", allow_cartesian = FALSE),
  IGHJ = list(figure = "supp5", bed_tag = "IGHJ", nonamer_start = 29, leader = FALSE, subgroup_from_group = FALSE, strip_d = TRUE,
              order_by = "label", order_levels = c("J2", "J5", "J1", "J4", "J3", "J6"),
              usage_locus = "IGHJ", chain = "IGH", call_col = "j_call_new", allow_cartesian = FALSE),
  IGKJ = list(figure = "supp5", bed_tag = "IGKJ", nonamer_start = 29, leader = FALSE, subgroup_from_group = FALSE, strip_d = TRUE,
              order_by = "label", order_levels = c("J2", "J3", "J5", "J1", "J4"),
              usage_locus = "IGKJ", chain = "IGK", call_col = "j_call_new", allow_cartesian = FALSE),
  IGLJ = list(figure = "supp5", bed_tag = "IGLJ", nonamer_start = 19, leader = FALSE, subgroup_from_group = FALSE, strip_d = TRUE,
              order_by = "label", order_levels = c("J1", "J2", "J3-3", "J3-1", "J3-2", "J3-4", "J6", "J7"),
              usage_locus = "IGLJ", chain = "IGL", call_col = "j_call_new", allow_cartesian = TRUE)
)

# IGHD heatmap columns are ASC groups, as on the D axis of the D-J pairing figure: an ASC bundling
# several IUIS D genes ("IGHD5-18/IGHD5-5") collapses to one starred column ("D5-18*").
ighd_asc_map <- local({
  asc <- unique(fread(need(file.path(OUT$repertoire, "gg_repertoire_data_IGH_genotype_corrected.csv.gz")), select = "d_gene_iuis")$d_gene_iuis)
  asc <- asc[!is.na(asc) & nzchar(asc)]
  pairs <- unique(rbindlist(lapply(asc, function(a) data.table(member = gsub("^IGH", "", strsplit(a, "/", fixed = TRUE)[[1]]),
                                                               group = gsub("^IGH", "", sub("[/].*$", "*", a))))), by = "member")
  setNames(pairs$group, pairs$member)
})

dend_tables <- function(hc) list(
  merges = data.table(step = seq_len(nrow(hc$merge)), left = hc$merge[, 1], right = hc$merge[, 2], height = hc$height),
  leaves = data.table(leaf_index = seq_along(hc$labels), label = hc$labels, order_slot = match(seq_along(hc$labels), hc$order)))
dense_counts <- function(dt, row_col, rows, cols) {
  counts <- dt[, .N, by = c(row_col, "iuis_group")]
  setnames(counts, c(row_col, "iuis_group", "N"), c("row_seq", "column_label", "count"))
  dense <- CJ(row_seq = rows, column_label = cols, unique = TRUE)
  dense[counts, count := i.count, on = .(row_seq, column_label)]
  dense[is.na(count), count := 0L][]
}
optimise <- function(mat, row_dend, col_dend, leader) {
  set.seed(123)
  optimize_order_enhanced(mat = mat, row_dend = row_dend, col_dend = col_dend, kernel_size = 3,
                          max_iter_family = if (leader) 100 else 50, max_iter_global = 200, max_steps = if (leader) 4 else 3, verbose = FALSE,
                          consider_mixed_flips = leader, consider_minimal_subtree = FALSE, consider_subgroups = leader,
                          optimize_subgroups = TRUE, use_cache = TRUE)
}

matrix_out <- list(); rows_out <- list(); cols_out <- list(); merges_out <- list(); leaves_out <- list(); consensus_out <- list(); usage_out <- list()
column_order_cache <- list()
for (panel in names(panels)) {
  cfg <- panels[[panel]]
  message("supp45: ", panel)
  relabel_asc <- function(x) { if (!isTRUE(cfg$merge_asc)) return(x); m <- unname(ighd_asc_map[x]); fifelse(is.na(m), x, m) }
  d <- rss_leader[gene_type == panel & present == TRUE]
  if (cfg$subgroup_from_group) d[, iuis_subgroup := sub("-.*$", "", iuis_group)]  # IGHD subgroups come through as bare numbers

  # Column dendrogram from the coding sequences, collapsed to one tip per IUIS group.
  coding <- unique(d[, .(vdjbase_allele, iuis_group, seq)])
  hc_cols <- hclust(stringdist::stringdistmatrix(setNames(coding$seq, coding$vdjbase_allele), method = "lv", useNames = "names"), method = "complete")
  collapsed_tree <- Reduce(function(tr, node) collapse_clade(tr, node), tapply(coding$vdjbase_allele, coding$iuis_group, identity), init = ape::as.phylo(hc_cols))
  collapsed_tree$tip.label <- coding$iuis_group[match(collapsed_tree$tip.label, coding$vdjbase_allele)]
  hc_cols <- ape::as.hclust.phylo(collapsed_tree)
  dend_cols <- as.dendrogram(hc_cols)

  seqs_rss <- unique(d$rss_aligned)
  hc_rows_rss <- hclust(stringdist::stringdistmatrix(seqs_rss, method = "hamming", useNames = "strings"), method = "complete")
  mat_rss <- table(d$rss_aligned, d$iuis_group)[hc_rows_rss$labels[hc_rows_rss$order], hc_cols$labels[hc_cols$order]]
  result_rss <- optimise(mat_rss, as.dendrogram(hc_rows_rss), dend_cols, leader = FALSE)

  if (is.null(cfg$columns_from)) {
    genes_order <- gsub("IG[HKL]", "", gene_order[grepl(cfg$bed_tag, g)]$g)
    strip <- function(x) if (cfg$strip_d) gsub("D", "", sub("-.*$", "", x)) else sub("-.*$", "", x)
    cols_dt <- data.table(label = labels(result_rss$col_dend), subgroup = strip(labels(result_rss$col_dend)), group_order = match(labels(result_rss$col_dend), genes_order))
    cols_dt[, subgroup_order := as.numeric(factor(subgroup, levels = unique(subgroup)))]
    cols_dt[order(subgroup_order, group_order), order := 1:.N]
    cols_dt <- cols_dt[order(match(if (cfg$order_by == "subgroup") subgroup else label, cfg$order_levels))]
    if (isTRUE(cfg$merge_asc)) {  # collapse IUIS-gene columns to ASC groups, keeping each group's first slot
      cols_dt[, label := relabel_asc(label)]
      cols_dt <- unique(cols_dt, by = "label")
      cols_dt[, order := seq_len(.N)]
    }
  } else {
    cols_dt <- copy(column_order_cache[[cfg$columns_from]])
  }
  column_order_cache[[panel]] <- cols_dt

  # Usage is fed the member genes: its union-find derives the merged groups, so each drawn column's
  # split_group is the merged group of its base gene (the label minus any "*").
  keep_indiv <- if (isTRUE(cfg$merge_asc)) unique(d$iuis_group) else cols_dt$label
  group_map <- build_allele_group_map(husa_path = husa_tsv, locus = cfg$usage_locus, keep_labels = keep_indiv, allow.cartesian = cfg$allow_cartesian)
  usage <- build_usage_objects(repertoire_csv = file.path(OUT$repertoire, sprintf("gg_repertoire_data_%s_genotype_corrected.csv.gz", cfg$chain)),
                               locus = cfg$usage_locus, call_col = cfg$call_col, group_map = group_map, keep_labels = keep_indiv)
  usage_map <- copy(usage$map)[, merged := gsub("IG[KLH]", "", merged)]
  cols_dt[, split_group := usage_map[match(sub("[*]$", "", cols_dt$label), iuis_group), merged]]
  m_box <- usage$m_box
  usage_out[[panel]] <- data.table(figure = cfg$figure, panel = panel, subject_index = rep(seq_len(nrow(m_box)), times = ncol(m_box)),
                                   merged = rep(colnames(m_box), each = nrow(m_box)), rel_usage = as.numeric(m_box))
  cols_out[[panel]] <- data.table(figure = cfg$figure, panel = panel, position = seq_len(nrow(cols_dt)), column_label = cols_dt$label,
                                  column_subgroup = cols_dt$subgroup, split_group = cols_dt$split_group)

  # The right-hand subgroup bar is positional; for the rotated IGHD 3' panel it is built against the unrotated leaf order.
  rss_subgroup <- d[, .(subgroup = paste0(unique(iuis_subgroup), collapse = ",")), by = rss_aligned][order(match(rss_aligned, as.hclust(result_rss$row_dend)$labels))]
  hc_rss <- if (is.null(cfg$rotate_subgroup_first)) as.hclust(result_rss$row_dend) else {
    dend <- result_rss$row_dend
    leaf_subgroup <- unique(d[, .(iuis_subgroup, rss_aligned)])[order(match(rss_aligned, labels(dend)))]
    first <- which(leaf_subgroup$iuis_subgroup == cfg$rotate_subgroup_first)
    as.hclust(dendextend::rotate(dend, order = labels(dend)[c(first, setdiff(seq_len(nrow(leaf_subgroup)), first))]))
  }
  consensus_rss <- generate_consensus(seqs_rss, vec_onehot)
  d_counts <- if (isTRUE(cfg$merge_asc)) copy(d)[, iuis_group := relabel_asc(iuis_group)] else d  # sums an ASC group's members into one cell

  add_layer <- function(layer, dt, row_col, hc, subgroup_dt, consensus, boundary) {
    key <- paste(panel, layer)
    matrix_out[[key]] <<- data.table(figure = cfg$figure, panel = panel, layer = layer, dense_counts(dt, row_col, hc$labels, cols_dt$label))
    rows_out[[key]] <<- data.table(figure = cfg$figure, panel = panel, layer = layer, position = seq_along(hc$labels), row_seq = hc$labels, annotation_subgroup = subgroup_dt$subgroup)
    dd <- dend_tables(hc)
    merges_out[[key]] <<- data.table(figure = cfg$figure, panel = panel, layer = layer, dd$merges)
    leaves_out[[key]] <<- data.table(figure = cfg$figure, panel = panel, layer = layer, dd$leaves)
    consensus_out[[key]] <<- data.table(figure = cfg$figure, panel = panel, layer = layer, consensus = consensus$consensus, nonamer_start = boundary, rss = layer == "rss")
  }
  add_layer("rss", d_counts, "rss_aligned", hc_rss, rss_subgroup, consensus_rss, cfg$nonamer_start)
  if (!isTRUE(cfg$leader)) next

  # Leader rows. Rows whose L-PART1 is entirely missing are dropped; non-ATG leaders are kept because
  # some genes (IGLV8-61) carry only those.
  leader_dt <- d[!is.na(leader) & nzchar(leader) & !is.na(l_part1) & nzchar(l_part1)]
  leader_seqs <- unique(leader_dt$leader)
  leader2_start <- max(nchar(d$l_part1_aligned))
  hc_rows_leader <- hclust(stringdist::stringdistmatrix(leader_seqs, method = if (uniqueN(nchar(leader_seqs)) == 1L) "hamming" else "lv", useNames = "strings"), method = "complete")
  mat_leader <- table(leader_dt$leader, leader_dt$iuis_group)
  missing_cols <- setdiff(cols_dt$label, colnames(mat_leader))
  if (length(missing_cols)) mat_leader <- cbind(mat_leader, matrix(0L, nrow = nrow(mat_leader), ncol = length(missing_cols), dimnames = list(rownames(mat_leader), missing_cols)))
  mat_leader <- mat_leader[hc_rows_leader$labels[hc_rows_leader$order], cols_dt$label, drop = FALSE]
  result_leader <- optimise(mat_leader, as.dendrogram(hc_rows_leader), dend_cols, leader = TRUE)
  hc_leader <- as.hclust(result_leader$row_dend)
  leader_subgroup <- leader_dt[, .(subgroup = paste0(unique(iuis_subgroup), collapse = ",")), by = leader][order(match(leader, hc_leader$labels))]
  add_layer("leader", leader_dt, "leader", hc_leader, leader_subgroup, generate_consensus(leader_seqs, vec_onehot), leader2_start)
}

matrix_dt <- rbindlist(matrix_out); rows_dt <- rbindlist(rows_out); cols_all <- rbindlist(cols_out)
merges_dt <- rbindlist(merges_out); leaves_dt <- rbindlist(leaves_out); consensus_dt <- rbindlist(consensus_out); usage_dt <- rbindlist(usage_out)
fwrite(matrix_dt, file.path(OUT$source, "supp45_matrix.csv.gz"))
fwrite(rows_dt, file.path(OUT$source, "supp45_rows.csv.gz"))
fwrite(cols_all, file.path(OUT$source, "supp45_columns.csv"))
fwrite(merges_dt, file.path(OUT$source, "supp45_dendrogram_merges.csv.gz"))
fwrite(leaves_dt, file.path(OUT$source, "supp45_dendrogram_leaves.csv.gz"))
fwrite(consensus_dt, file.path(OUT$source, "supp45_consensus.csv"))
fwrite(usage_dt, file.path(OUT$source, "supp45_usage.csv.gz"))
print(rows_dt[, .(rows = .N), by = .(figure, panel, layer)])

# ---- draw ----
count_matrix <- function(dt, rows, cols) {
  wide <- dcast(dt, row_seq ~ column_label, value.var = "count")
  m <- as.matrix(wide[, -"row_seq"]); rownames(m) <- wide$row_seq
  m[rows, cols, drop = FALSE]
}
restore_hclust <- function(merges, leaves) {
  merges <- merges[order(step)]
  structure(list(merge = matrix(c(merges$left, merges$right), ncol = 2L), height = merges$height, order = leaves[order(order_slot), leaf_index],
                 labels = leaves[order(leaf_index), label], method = "complete", dist.method = "hamming"), class = "hclust")
}
# Everything a panel layer needs that does not depend on where it sits on the page.
panel_parts <- function(p, lyr, palette, annotation_name_side = "right") {  # lyr, not layer: the column would shadow it inside `[`
  cols <- cols_all[panel == p][order(position)]
  rows <- rows_dt[panel == p & layer == lyr][order(position)]
  cons <- consensus_dt[panel == p & layer == lyr]
  subgroups <- unique(c(cols$column_subgroup, rows_dt[panel == p]$annotation_subgroup))
  graphics <- setNames(lapply(subgroups, subgroups_annotation_custom, colors_rss_subgroup = palette), subgroups)
  m_box <- as.matrix(dcast(usage_dt[panel == p], subject_index ~ merged, value.var = "rel_usage")[, -"subject_index"])
  list(mat = count_matrix(matrix_dt[panel == p & layer == lyr], rows$row_seq, cols$column_label),
       hc = restore_hclust(merges_dt[panel == p & layer == lyr], leaves_dt[panel == p & layer == lyr]),
       column_order = cols$column_label, column_split = factor(cols$split_group, levels = unique(cols$split_group)),
       column_ha = HeatmapAnnotation(Subgroup = anno_customize(setNames(cols$column_subgroup, cols$column_label), graphics = graphics),
                                     height = unit(10, "mm"), annotation_name_side = annotation_name_side),
       usage_ha = build_collapsed_box_anno(m_box = m_box, align_to_vec = cols$split_group, ylab = "Usage", show_axis_at = 1L),
       subgroup_anno = anno_customize(rows$annotation_subgroup, graphics = graphics, which = "row"),
       labels_anno = row_anno_text(gt_render(create_row_labels(rows$row_seq, cons$consensus, nucleotide_colors, cons$nonamer_start, rss = cons$rss))),
       consensus_label = create_consensus_label(cons$consensus, nucleotide_colors, cons$nonamer_start, rss = cons$rss))
}
base_heatmap <- function(parts, name, ...) {
  Heatmap(parts$mat, cluster_rows = parts$hc, cluster_columns = FALSE, column_order = parts$column_order, show_row_names = FALSE, show_column_names = TRUE,
          col = col_fun, border = FALSE, show_heatmap_legend = FALSE, row_names_side = "left", row_names_max_width = unit(16, "cm"),
          row_names_gp = gpar(fontsize = 18, fontfamily = "mono"), column_names_gp = gpar(fontsize = 18, fontfamily = "mono"),
          row_dend_width = unit(4, "cm"), column_dend_height = unit(4, "cm"), rect_gp = gpar(col = "gray80", lwd = 2), column_split = parts$column_split, ...)
}
lgd <- Legend(col_fun = col_fun, title = "Count", direction = "horizontal", labels_gp = gpar(fontsize = 18), title_gp = gpar(fontsize = 20, fontface = "bold"))
consensus_strip <- function(name, label) decorate_annotation(name, grid.draw(richtext_grob(label, y = unit(-4, "mm"), x = unit(0, "npc"), hjust = 0)))

## supp4: IGKV and IGLV, RSS over leader
pdf(file.path(OUT$figures, "supp4.pdf"), width = 35, height = 40)
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = unit(c(1, 1.3), "null"))))
for (i in seq_along(c("IGKV", "IGLV"))) {
  p <- c("IGKV", "IGLV")[i]
  rss <- panel_parts(p, "rss", colors15); leader <- panel_parts(p, "leader", colors15)
  ht_rss <- base_heatmap(rss, column_title = "RSS", column_title_gp = gpar(fontsize = 0, col = "black", rotate = 90), row_dend_side = "right", column_names_side = "top",
                         left_annotation = do.call(rowAnnotation, c(setNames(list(rss$labels_anno), paste0("labels_rss_", p)), list(show_annotation_name = FALSE))),
                         right_annotation = rowAnnotation(Subgroup = rss$subgroup_anno, width = unit(5, "mm"), show_annotation_name = FALSE),
                         top_annotation = rss$column_ha)
  ht_leader <- base_heatmap(leader, column_title = p, column_title_gp = gpar(fontsize = 18, col = "black", rotate = 90), row_dend_side = "right",
                            left_annotation = do.call(rowAnnotation, c(setNames(list(leader$labels_anno), paste0("labels_leader_", p)), list(show_annotation_name = FALSE))),
                            right_annotation = rowAnnotation(Subgroup = leader$subgroup_anno, width = unit(5, "mm"), show_annotation_name = FALSE),
                            bottom_annotation = leader$column_ha, top_annotation = leader$usage_ha)
  title_param <- ht_rss@column_title_param
  title_param$gp$fontsize <- 20
  ht_rss@column_title <- p; ht_rss@column_title_param <- title_param
  title_param$rot <- 270; title_param$just <- c(0.5, 0.3)
  ht_rss@row_title <- "RSS"; ht_rss@row_title_param <- title_param
  ht_leader@row_title <- "Leader"; ht_leader@row_title_param <- title_param
  pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1, y = unit(0.95, "npc"), height = unit(0.95, "npc"), just = c("center", "top")))
  draw(ht_rss %v% ht_leader, newpage = FALSE)
  if (i == 1) draw(lgd, x = unit(0.1, "npc"), y = unit(0.98, "npc"), just = c("left", "top"))
  consensus_strip(paste0("labels_rss_", p), rss$consensus_label)
  consensus_strip(paste0("labels_leader_", p), leader$consensus_label)
  grid.text(c("A", "B")[i], x = unit(0, "npc"), y = unit(0.97, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", cex = 3))
  popViewport()
}
dev.off()

## supp5: IGHD 5' over 3', then the three J segments
colors_d <- setNames(unname(colors15), paste0("D", 1:15))
colors_j <- setNames(unname(colors15), paste0("J", 1:15))
d5 <- panel_parts("IGHD_5", "rss", colors_d, "right"); d3 <- panel_parts("IGHD_3", "rss", colors_d, "left")
ht_d5 <- base_heatmap(d5, column_title = "IGHD", row_title = "5' RSS", row_title_side = "right", column_title_gp = gpar(fontsize = 20, col = "black", rotate = 90),
                      row_title_gp = gpar(fontsize = 20, col = "black", rotate = 270), row_dend_side = "right", column_names_side = "top",
                      left_annotation = do.call(rowAnnotation, c(list(labels_rss_IGHD_5 = d5$labels_anno), list(Subgroup = d5$subgroup_anno), list(show_annotation_name = FALSE))),
                      top_annotation = d5$column_ha)
ht_d3 <- base_heatmap(d3, row_title = "3' RSS", column_title_gp = gpar(fontsize = 0, col = "black", rotate = 90), row_title_gp = gpar(fontsize = 20, col = "black", rotate = 270),
                      row_dend_side = "left", column_names_side = "bottom",
                      right_annotation = do.call(rowAnnotation, c(list(Subgroup = d3$subgroup_anno), list(labels_rss_IGHD_3 = d3$labels_anno), list(show_annotation_name = FALSE))),
                      bottom_annotation = d3$column_ha, top_annotation = d3$usage_ha)
j_panels <- c("IGHJ", "IGKJ", "IGLJ")
j_parts <- lapply(setNames(j_panels, j_panels), panel_parts, lyr = "rss", palette = colors_j)
ht_j <- lapply(j_panels, function(p) {
  parts <- j_parts[[p]]
  base_heatmap(parts, column_title = p, column_title_gp = gpar(fontsize = 20, col = "black", rotate = 90), row_dend_side = "right", column_names_side = "bottom",
               left_annotation = do.call(rowAnnotation, c(setNames(list(parts$labels_anno), paste0("labels_rss_", p)), list(show_annotation_name = FALSE))),
               right_annotation = rowAnnotation(Subgroup = parts$subgroup_anno, width = unit(5, "mm"), show_annotation_name = FALSE),
               bottom_annotation = parts$column_ha, top_annotation = parts$usage_ha)
})
pdf(file.path(OUT$figures, "supp5.pdf"), width = 35, height = 30)
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 3, heights = unit(c(1.1, 0.15, 0.7), "null"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:3, y = unit(0.95, "npc"), height = unit(0.95, "npc"), just = c("center", "top")))
draw(ht_d5 %v% ht_d3, newpage = FALSE)
draw(lgd, x = unit(0.1, "npc"), y = unit(0.98, "npc"), just = c("left", "top"))
consensus_strip("labels_rss_IGHD_5", d5$consensus_label)
consensus_strip("labels_rss_IGHD_3", d3$consensus_label)
grid.text("A", x = unit(0, "npc"), y = unit(0.95, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", cex = 2))
popViewport()
for (i in seq_along(j_panels)) {
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = i))
  draw(ht_j[[i]], newpage = FALSE)
  consensus_strip(paste0("labels_rss_", j_panels[i]), j_parts[[i]]$consensus_label)
  grid.text(c("B", "C", "D")[i], x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", cex = 2))
  popViewport()
}
dev.off()
