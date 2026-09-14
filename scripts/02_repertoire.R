#!/usr/bin/env Rscript
# Genotype-corrected repertoires. Every rearrangement is named by its ASC allele and its
# IUIS gene group, and kept only where the subject's inferred genotype carries every gene
# it was assigned to.
#
#   data/repertoire/gg_repertoire_data_<CHAIN>_genotype_corrected.csv.gz
#   data/repertoire/genotypes_all_loci_collapsed.csv.gz       one row per subject x ASC gene
#   data/repertoire/genotype_inference_with_husa_labels.csv.gz genotype calls with ASC and IUIS names

source("scripts/00_setup.R")
suppressPackageStartupMessages({ library(alakazam); library(piglet) })

husa <- fread(need(file.path(OUT$husa, "husa.tsv")))
husa$iuis_allele <- husa$husa  # `$`, not `[`: a column is called husa and would shadow the table inside `[`

rep_all <- fread(need(IN$repertoires), select = c("vdjbase_subject", "subject", "locus", "productive", "v_call", "d_call",
                                                  "j_call", "junction", "junction_aa", "junction_length", "v_mut", "valid_both"))

# ASC subgroup -> IUIS group labels. D is relabelled per allele at IUIS gene resolution,
# V and J keep the ASC grouping (see ?ascIUISVocabulary).
voc <- ascIUISVocabulary(husa)

# Deletion alleles are gapped in the inferred sequence but not in the reference, so the
# sequence lookup is ungapped on both sides.
ungap <- function(x) gsub("-", "", x, fixed = TRUE)
alleles <- unique(husa[, .(asc_allele = allele, iuis_allele, seq)])
asc_of_seq <- setNames(alleles$asc_allele, ungap(alleles$seq))
iuis_of_seq <- setNames(alleles$iuis_allele, ungap(alleles$seq))

genotypes <- fread(need(IN$genotype))
genotypes[, allele_iuis := iuis_of_seq[ungap(seq)]]
genotypes[, asc_allele := asc_of_seq[ungap(seq)]]
unnamed <- genotypes[is.na(asc_allele) | !nzchar(asc_allele)]
cat(sprintf("%d of %d genotype calls have no reference name and are dropped (%d with reads)\n",
            nrow(unnamed), nrow(genotypes), nrow(unnamed[round_count >= 1])))
fwrite(genotypes, file.path(OUT$repertoire, "genotype_inference_with_husa_labels.csv.gz"))

name_of <- with(genotypes[!duplicated(paste(allele, asc_allele))], setNames(asc_allele, allele))
genotypes[, asc_gene := gsub("[*].*$", "", asc_allele)]
genotypes[, iuis_tag := fifelse(substr(asc_allele, 4, 4) == "D", voc$d_allele_to_iuis[asc_allele], voc$asc_to_iuis[asc_gene])]
collapsed <- genotypes[z_score >= 0, {
  o <- order(asc_allele)
  .(genotypes_alleles = paste0(asc_allele[o], collapse = ","),
    genotypes_iuis_alleles = paste0(allele_iuis[o], collapse = ";"),
    z_score = paste0(round(z_score[o], 2), collapse = ","),
    counts = paste0(round_count[o], collapse = ","),
    alleles = paste0(gsub(paste0(asc_gene, "[*]"), "", asc_allele[o]), collapse = ","),
    total_count = sum(round_count),
    relative_usage = sum(round_count) / depth,
    tag = substr(asc_gene, 4, 4))
}, by = .(vdjbase_subject, asc_gene, iuis_tag, depth)]
fwrite(collapsed, file.path(OUT$repertoire, "genotypes_all_loci_collapsed.csv.gz"))

# A comma-joined call maps part by part; parts with no name are dropped.
map_calls <- function(voc, calls) {
  if (is.na(calls) || !nzchar(calls)) return(NA_character_)
  m <- voc[unlist(strsplit(calls, ",", fixed = TRUE))]
  m <- m[!is.na(m) & nzchar(m)]
  if (!length(m)) NA_character_ else paste0(m, collapse = ",")
}

rep_all[, chain := substr(fifelse(!is.na(v_call) & nzchar(v_call), v_call, locus), 1, 3)]
rep_all <- rep_all[is.na(v_mut) | v_mut <= 0]

named <- !is.na(collapsed$asc_gene) & nzchar(collapsed$asc_gene)
allowed <- unique(paste(collapsed$vdjbase_subject[named], collapsed$asc_gene[named]))
has_genotype <- unique(collapsed$vdjbase_subject)

for (ch in CHAINS) {
  d <- rep_all[chain == ch]
  n_in <- nrow(d)
  d[, v_call_new := name_of[v_call], by = v_call]
  if (ch == "IGH") {
    d[, d_call_new := name_of[d_call], by = d_call]
    d[is.na(d_call_new), d_call_new := map_calls(name_of, d_call), by = d_call]
  } else {
    d[, d_call_new := ""]
  }
  d[, j_call_new := name_of[j_call], by = j_call]
  d[is.na(j_call_new), j_call_new := map_calls(name_of, j_call), by = j_call]
  d[, `:=`(v_gene = getGene(v_call_new, first = FALSE, omit_nl = FALSE, strip_d = FALSE, collapse = TRUE),
           d_gene = getGene(d_call_new, first = FALSE, omit_nl = FALSE, strip_d = FALSE, collapse = TRUE),
           j_gene = getGene(j_call_new, first = FALSE, omit_nl = FALSE, strip_d = FALSE, collapse = TRUE))]

  # Ambiguous J calls (several genes) cannot be checked against the genotype. Subjects with
  # no genotype keep all their rows.
  subject_order <- unique(d$vdjbase_subject)
  d <- d[!grepl(",", j_gene)]
  in_genotype <- function(gene) paste(d$vdjbase_subject, gene) %chin% allowed
  d <- d[!(vdjbase_subject %chin% has_genotype) |
           (in_genotype(v_gene) & in_genotype(j_gene) & (ch != "IGH" | in_genotype(d_gene)))]
  d[, .subject_order := match(vdjbase_subject, subject_order)]
  setorder(d, .subject_order)
  d[, .subject_order := NULL]

  d <- annotateRepertoireIUIS(d, voc, chain = ch)
  fwrite(d, file.path(OUT$repertoire, sprintf("gg_repertoire_data_%s_genotype_corrected.csv.gz", ch)))
  cat(sprintf("%s: %s of %s rearrangements kept, %d subjects\n", ch,
              format(nrow(d), big.mark = ","), format(n_in, big.mark = ","), uniqueN(d$vdjbase_subject)))
}
