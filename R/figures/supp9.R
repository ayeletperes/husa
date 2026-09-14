#!/usr/bin/env Rscript
# Supplementary figure 9: unique alleles per ASC, split by cohort support (AIRR-seq only,
# genomic only, both, neither). ASCs that share a gene are merged so each gene appears once.

source("R/00_setup.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(igraph) })

husa <- fread(need(file.path(OUT$husa, "husa.tsv")))

genes_of <- function(alleles) {
  alleles <- unique(unlist(strsplit(unlist(strsplit(alleles, ",")), "/", fixed = TRUE)))
  paste0(unique(sapply(strsplit(alleles, "*", fixed = TRUE), "[[", 1)), collapse = ",")
}
unify <- function(genes) paste0(unique(unlist(strsplit(genes, ","))), collapse = ",")

alleles <- unique(husa[, {
  n_airr <- uniqueN(unlist(strsplit(samples_AIRRseq, ",")))
  n_genomic <- uniqueN(unlist(strsplit(samples_genomic_watson, ","))) + uniqueN(unlist(strsplit(samples_genomic_hprc, ","))) +
    uniqueN(unlist(strsplit(samples_genomic_1kpg, ",")))
  .(genes = genes_of(husa),
    tag = if (n_airr > 0 && n_genomic == 0) "airr_seq_only" else if (n_airr == 0 && n_genomic > 0) "genomic_only"
          else if (n_airr == 0 && n_genomic == 0) "not_in_cohorts" else "both")
}, by = .(gene_type, asc, seq)])

# ASCs sharing a gene are joined through a bipartite ASC-gene graph, one merged ASC per component.
asc_genes <- alleles[, .(genes = unify(genes)), by = .(gene_type, asc)]
edges <- asc_genes[, .(gene = unique(trimws(unlist(strsplit(genes, ",", fixed = TRUE))))), by = .(gene_type, asc)][!is.na(gene) & gene != ""]
merge_map <- edges[, {
  comp <- components(graph_from_data_frame(data.table(from = paste0("ASC__", asc), to = paste0("GENE__", gene)), directed = FALSE))$membership
  asc_nodes <- names(comp)[startsWith(names(comp), "ASC__")]
  tmp <- data.table(asc = sub("^ASC__", "", asc_nodes), component = comp[asc_nodes])
  tmp[, merged_asc := paste(sort(unique(asc)), collapse = "/"), by = component]
  tmp[, .(asc, merged_asc)]
}, by = gene_type]
alleles <- merge(alleles, merge_map, by = c("gene_type", "asc"), all.x = TRUE)
alleles[is.na(merged_asc), merged_asc := asc]

counts <- alleles[, .(
  original_ascs = paste(sort(unique(asc)), collapse = ","), genes = unify(genes), unique_alleles = uniqueN(seq),
  airr_seq_only = sum(tag == "airr_seq_only"), genomic_only = sum(tag == "genomic_only"),
  both = sum(tag == "both"), not_in_cohorts = sum(tag == "not_in_cohorts")
), by = .(gene_type, asc = merged_asc)]
plot_data <- melt(counts, id.vars = c("gene_type", "asc", "original_ascs", "genes", "unique_alleles"),
                  measure.vars = c("airr_seq_only", "genomic_only", "both", "not_in_cohorts"),
                  variable.name = "tag", value.name = "count")
tag_levels <- c("Both", "Genomic only", "AIRR-seq only", "Not in cohorts")
plot_data[, tag_display := factor(c(airr_seq_only = "AIRR-seq only", genomic_only = "Genomic only", both = "Both",
                                    not_in_cohorts = "Not in cohorts")[as.character(tag)], levels = tag_levels)]
fwrite(plot_data, file.path(OUT$source, "supp9_allele_counts.csv.gz"))

fill_values <- c("Both" = "#1b9e77", "Genomic only" = "#d95f02", "AIRR-seq only" = "#7570b3", "Not in cohorts" = "#e7298a")
panels <- sapply(unique(plot_data$gene_type), function(gt) {
  d <- plot_data[gene_type == gt][order(unique_alleles)]
  d[, genes := factor(gsub("IG[HLK]", "", genes), levels = unique(gsub("IG[HLK]", "", genes)))]
  p <- ggplot(d, aes(x = genes, y = count, fill = tag_display)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~gene_type, scales = "free") +
    labs(x = "", y = "# of unique alleles", fill = "Category") +
    ggpubr::theme_pubclean(base_size = 24) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    scale_y_continuous(breaks = function(x) pretty(x)[pretty(x) %% 1 == 0]) +
    scale_fill_manual(values = fill_values)
  if (grepl("J", gt)) p <- p + ylab(NULL)
  if (grepl("V", gt)) p <- p + xlab(NULL)
  if (gt == "IGHJ") p <- p + xlab("ASCs")
  p
})

final <- panels[["IGHV"]] + panels[["IGLV"]] + panels[["IGKV"]] +
  panels[["IGHD"]] + panels[["IGHJ"]] + panels[["IGLJ"]] + panels[["IGKJ"]] +
  plot_layout(design = "AAAA\nBBBB\nCCCC\nDEFG", guides = "collect", widths = c(2, 1, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", plot.tag = element_text(size = 28, face = "bold", vjust = -3))
ggsave(file.path(OUT$figures, "supp9.pdf"), final, width = 16, height = 28)
