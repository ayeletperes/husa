suppressPackageStartupMessages(library(data.table))

CHAINS <- c("IGH", "IGK", "IGL")

IN <- list(
  ggs_watson            = "data/ggs/Watson_GGS_selected.csv",
  ggs_watson_unfiltered = "data/ggs/Watson_GGS_unfiltered.csv.gz",
  ggs_hprc              = "data/ggs/HPRC_GGS_selected.csv",
  ggs_1kpg              = "data/ggs/1KPG_GGS_selected.csv",
  genotype              = "data/genotype/genotype_inference.csv",
  snp_dosage            = "data/genotype/snp_dosage.matrix",
  repertoires           = "data/repertoires/ggs_repertoires_valid_samples.tsv.gz",
  subject_filter        = "data/reference/subject_filter_unfiltered.rds",
  hprc_subjects         = "data/reference/hprc_subject_filter.rds",
  kgp_subjects          = "data/reference/1kpg_subject_filter.rds",
  baseline_asc          = "data/reference/baseline_asc",
  iglabel_labels        = "data/reference/iglabel_labels",
  bed_dir               = "data/reference/immune_receptor_genomics_251106",
  imgt_rss              = "data/reference/imgt_rss_lv.tsv",
  imgt_leader           = "data/reference/imgt/leader",
  imgt_vdj              = "data/reference/imgt/vdj",
  watson_reference      = "data/reference/watson_reference_alleles",
  watson_metadata       = "data/metadata/watson_metadata.csv",
  hprc_metadata         = "data/metadata/hprc_metadata.tsv",
  kgp_metadata          = "data/metadata/1kpg_metadata.tsv"
)
IN$gene_bed <- file.path(IN$bed_dir, "gene.bed")

OUT <- list(
  husa       = "results/husa",
  repertoire = "results/repertoire",
  rss_leader = "results/rss_leader",
  qtl        = "results/qtl",
  figures    = "results/figures",
  source     = "results/figures/source_data"
)
for (d in OUT) dir.create(d, recursive = TRUE, showWarnings = FALSE)

need <- function(path) {
  if (!file.exists(path)) stop("missing input: ", path, " (run fetch_data.sh)", call. = FALSE)
  path
}
