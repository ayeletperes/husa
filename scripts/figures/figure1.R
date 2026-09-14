#!/usr/bin/env Rscript
# Figure 1: subjects per locus and ancestry in the three genomic cohorts.
# The table is computed only when missing from figure_data; the figure is
# always drawn from the table, so the shipped source_data alone regenerates it.

source("scripts/00_setup.R")
suppressPackageStartupMessages(library(ggplot2))

tables <- file.path(OUT$source, "figure1_subject_counts.csv")
if (!all(file.exists(tables))) {
  per_locus <- function(subjects_rds, ancestry, id_col, ancestry_col, source) {
    subjects <- readRDS(subjects_rds)
    dt <- rbindlist(lapply(CHAINS, function(locus) data.table(subject = subjects[[locus]], locus = locus)))
    dt <- merge(dt, ancestry[, .(subject = get(id_col), ancestry = get(ancestry_col))], by = "subject", all.x = TRUE)
    dt[, .(count = uniqueN(subject), source = source), by = .(locus, ancestry)]
  }
  counts <- rbindlist(list(
    per_locus(need(IN$subject_filter), fread(need(IN$watson_metadata)), "vdjbase_name", "ancestry_population", "In-House"),
    per_locus(need(IN$kgp_subjects), fread(need(IN$kgp_metadata)), "sample", "Ancestry", "1KGP"),
    per_locus(need(IN$hprc_subjects), fread(need(IN$hprc_metadata)), "sample", "Ancestry", "HPRC")
  ))
  fwrite(counts, tables[1])
}

counts <- fread(tables[1])
p <- ggplot(counts, aes(x = locus, y = count, fill = ancestry)) +
  geom_bar(stat = "identity") +
  labs(x = "Locus", y = "Number of Individuals", fill = "Ancestry Population") +
  ggpubr::theme_pubclean(base_size = 24) +
  theme(legend.position = "top") +
  facet_wrap(~source)
ggsave(file.path(OUT$figures, "figure1.png"), p, width = 12, height = 8, dpi = 300)
