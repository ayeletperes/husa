#!/usr/bin/env Rscript
# HUSA allele table. Merges the genomic germline sets of the three cohorts with the
# AIRR-seq genotypes and the baseline ASC reference, attaches the frozen IgLabel labels,
# and re-clusters every segment into allele similarity clusters with piglet.
#
#   data/husa/husa.tsv             one row per allele (2587 in the manuscript run)
#   data/husa/husa_rss_filter.tsv  the subset with complete, consistent RSS evidence

source("scripts/00_setup.R")

collapse_unique <- function(x, sep = ", ") {
  x <- unique(x[!is.na(x)]); x <- x[x != ""]
  if (!length(x)) "" else paste(x, collapse = sep)
}
tokens <- function(x, sep = ",") {
  v <- trimws(unlist(strsplit(paste(x[!is.na(x)], collapse = sep), sep, fixed = TRUE)))
  unique(v[nzchar(v)])
}
collapse_tokens <- function(x) paste(tokens(x), collapse = ", ")
count_tokens <- function(x) vapply(as.character(x), function(v) length(tokens(v)), integer(1))
short_iglabel <- function(x) {
  x <- trimws(x); x[is.na(x)] <- ""
  sub("\\*[0-9]+$", "", sub("^IG[HKL][VDJ]1-", "", x))
}
select_allele <- function(a) {
  a <- trimws(unique(unlist(strsplit(a, ",")))); a <- a[nzchar(a)]
  if (length(a) <= 1L) return(a)
  known <- grep("_", a, invert = TRUE, value = TRUE)
  if (length(known)) known[1] else a[which.min(nchar(a))]
}
select_gap_sequence <- function(s, alleles) {
  s <- unique(s)
  if (length(s) == 1L) return(s)
  if (any(!grepl("_", alleles))) return(s[min(which(!grepl("_", alleles))[1], length(s))])
  no_del <- grep("-", s, value = TRUE, invert = TRUE)
  if (length(no_del)) no_del[1] else s[1]
}

# One GGS cohort, deduplicated, with the per-subject evidence flags the Watson cohort
# carries (a subject is "sure" for an allele when its functionality and every RSS/leader
# part agree across contigs).
read_genomic <- function(path, dataset, valid = NULL) {
  dt <- fread(path)
  setnames(dt, c("selected_allele", "gapped_sequence", "sequence"),
           c("allele_selected", "seq_gapped", "seq"), skip_absent = TRUE)
  if (!"seq_gapped" %in% names(dt)) dt[, seq_gapped := seq]
  dt[is.na(seq_gapped) | seq_gapped == "", seq_gapped := seq]
  dt <- dt[!is.na(allele_selected) & allele_selected != "" & !is.na(seq) & seq != ""]
  dt[, chain := substr(gene_type, 1, 3)]
  if (dataset == "watson") {
    dt[, subject_id := vdjbase_subject]
    dt[, valid := vdjbase_subject %in% valid[[chain]], by = chain]
    dt <- dt[valid == TRUE]
  } else {
    dt[, subject_id := subject]
    dt <- dt[!is.na(subject_id) & subject_id != ""]
    dt[, vdjbase_subject := ""]
  }
  dt[, genomic_dataset := dataset]
  key_cols <- intersect(c("genomic_dataset", "subject_id", "vdjbase_subject", "subject", "allele_selected",
                          "gene_type", "v_heptamer", "v_nonamer", "j_heptamer", "j_nonamer",
                          "d_3_heptamer", "d_3_nonamer", "d_5_heptamer", "d_5_nonamer",
                          "seq", "seq_gapped", "functional", "spacer_3", "spacer_5", "l_part1", "l_part2"), names(dt))
  dt <- dt[!duplicated(dt[, do.call(paste, c(.SD, sep = "|")), .SDcols = key_cols])]
  if (dataset == "watson") {
    evidence <- intersect(c("v_nonamer", "v_heptamer", "j_nonamer", "j_heptamer", "d_3_heptamer", "d_3_nonamer",
                            "d_5_heptamer", "d_5_nonamer", "spacer_3", "spacer_5", "l_part1", "l_part2"), names(dt))
    dt[, has_both := uniqueN(functional) != 1L, by = .(vdjbase_subject, allele_selected)]
    dt[, evidence_consistent := all(vapply(.SD, uniqueN, integer(1)) == 1L),
       by = .(vdjbase_subject, allele_selected), .SDcols = evidence]
    dt[, sure_subject := fifelse(!has_both & evidence_consistent, vdjbase_subject, "")]
    dt[, non_sure_subject := fifelse(has_both | !evidence_consistent, vdjbase_subject, "")]
  } else {
    dt[, `:=`(has_both = FALSE, evidence_consistent = TRUE, sure_subject = "", non_sure_subject = "")]
  }
  dt[, seq_gapped_selected := select_gap_sequence(seq_gapped, allele_selected), by = seq]
  dt
}

# ---- inputs -----------------------------------------------------------------------

samples <- fread(need(IN$repertoires), select = c("vdjbase_subject", "locus", "valid_both"))[valid_both == TRUE]
valid_airr <- tapply(samples$vdjbase_subject, samples$locus, unique)
# Genomic evidence is gated on the contamination-only subject set, not on valid_both:
# gating on valid_both silently drops the alleles seen only in genomic data.
valid_genomic <- readRDS(need(IN$subject_filter))

baseline <- rbindlist(lapply(CHAINS, function(ch) {
  fread(need(file.path(IN$baseline_asc, ch, paste0(ch, "_crossref.csv"))))
}), fill = TRUE)
label_col <- if ("iglabel_label" %in% names(baseline)) "iglabel_label" else "name"
baseline <- baseline[, .(allele = asc_allele, gene_type, iglabel = short_iglabel(get(label_col)),
                         wasp, seq, seq_gapped, in_baseline_reference = TRUE)]
baseline[seq_gapped == "", seq_gapped := seq]

genomic <- rbindlist(list(
  read_genomic(need(IN$ggs_watson), "watson", valid_genomic),
  read_genomic(need(IN$ggs_hprc), "hprc"),
  read_genomic(need(IN$ggs_1kpg), "1kpg")
), fill = TRUE, use.names = TRUE)

# An HPRC or 1KGP novel allele whose name encodes an indel is kept only when its
# sequence is seen in more than one subject across the two cohorts.
suffix <- ifelse(grepl("_", genomic$allele_selected, fixed = TRUE), sub("^.*?_", "", genomic$allele_selected), "")
deletion <- nzchar(suffix) & grepl("-", suffix, fixed = TRUE)
insertion <- nzchar(suffix) & (
  grepl("^i\\d+[acgt]+", suffix, ignore.case = TRUE) |
    (grepl("^[0-9a-f]{4}$", suffix, ignore.case = TRUE) & !grepl("^[acgt]\\d+[acgt]$", suffix, ignore.case = TRUE)) |
    grepl("\\d+[acgt]{2,}\\d+", suffix, ignore.case = TRUE))
genomic[, structural_subject := fifelse(genomic_dataset == "watson", vdjbase_subject, subject_id)]
genomic[, structural_subject := fifelse(!is.na(structural_subject) & structural_subject != "",
                                        paste(genomic_dataset, structural_subject, sep = ":"), "")]
genomic[, structural_candidate := genomic_dataset %in% c("hprc", "1kpg") & (deletion | insertion)]
support <- genomic[structural_candidate == TRUE, .(n_structural_subjects = uniqueN(structural_subject[nzchar(structural_subject)])), by = seq]
genomic <- merge(genomic, support, by = "seq", all.x = TRUE)
genomic[is.na(n_structural_subjects), n_structural_subjects := 0L]
cat(sprintf("structural novel candidates: %d rows, %d dropped\n", sum(genomic$structural_candidate),
            sum(genomic$structural_candidate & genomic$n_structural_subjects <= 1L)))
genomic <- genomic[!structural_candidate | n_structural_subjects > 1L]

geno <- fread(need(IN$genotype))[z_score >= 0]
geno[, chain := substr(gene, 1, 3)]
geno[, valid := vdjbase_subject %in% valid_airr[[chain]], by = chain]
geno <- geno[valid == TRUE]
geno[, repertoire_allele := allele]
geno[, allele_selected := select_allele(repertoire_allele), by = seq]
setnames(geno, "seq_gapped", "seq_gapped_repertoire")
geno[, seq := gsub("[.-]", "", seq_gapped_repertoire)]
geno[, seq_gapped_selected := select_gap_sequence(seq_gapped_repertoire, allele_selected), by = seq]
geno <- geno[z_score >= 1 & round_count >= 1]

# ---- per-subject evidence, one row per subject x allele ---------------------------

merged <- merge(
  geno[, .(vdjbase_subject, gene, allele_selected, repertoire_allele, count, threshold, z_score, seq_gapped_selected, seq)],
  genomic, by = c("vdjbase_subject", "allele_selected", "seq_gapped_selected", "seq"), all = TRUE)
merged[, seq_gapped := NULL]
setnames(merged, c("allele_selected", "seq_gapped_selected"), c("allele", "seq_gapped"))
merged[is.na(gene_type), gene_type := substr(allele, 1, 4)]
merged[, in_AIRRseq := !is.na(z_score)]
merged[, in_genomic_watson := genomic_dataset == "watson" & !is.na(subject_id) & subject_id != ""]
merged[, in_genomic_hprc := genomic_dataset == "hprc" & !is.na(subject_id) & subject_id != ""]
merged[, in_genomic_1kpg := genomic_dataset == "1kpg" & !is.na(subject_id) & subject_id != ""]
merged[, in_genomic := in_genomic_watson | in_genomic_hprc | in_genomic_1kpg]

# Rows that reached here without a sequence take it from another row of the same allele,
# then from the baseline.
first_nonempty <- function(x) { x <- x[!is.na(x) & x != ""]; if (length(x)) x[1] else "" }
for (src in list(merged, baseline)) {
  lookup <- src[!is.na(seq) & seq != "", .(seq_fill = seq[1], seq_gapped_fill = first_nonempty(seq_gapped)), by = .(allele, gene_type)]
  merged <- merge(merged, lookup, by = c("allele", "gene_type"), all.x = TRUE)
  merged[(is.na(seq) | seq == "") & !is.na(seq_fill) & seq_fill != "", seq := seq_fill]
  merged[(is.na(seq_gapped) | seq_gapped == "") & !is.na(seq_gapped_fill) & seq_gapped_fill != "", seq_gapped := seq_gapped_fill]
  merged[(is.na(seq_gapped) | seq_gapped == "") & !is.na(seq) & seq != "", seq_gapped := seq]
  merged[, c("seq_fill", "seq_gapped_fill") := NULL]
}
if (any(is.na(merged$seq) | merged$seq == "")) stop("rows without a sequence: ", sum(is.na(merged$seq) | merged$seq == ""))

keep_cols <- intersect(c("vdjbase_subject", "subject_id", "genomic_dataset", "allele", "repertoire_allele", "vdjbase_allele",
                         "gene_type", "v_heptamer", "v_nonamer", "j_heptamer", "j_nonamer", "d_3_heptamer", "d_3_nonamer",
                         "d_5_heptamer", "d_5_nonamer", "spacer_3", "spacer_5", "seq", "seq_gapped", "functional", "notes",
                         "z_score", "in_genomic", "in_genomic_watson", "in_genomic_hprc", "in_genomic_1kpg", "in_AIRRseq",
                         "sure_subject", "non_sure_subject", "l_part1", "l_part2"), names(merged))
merged <- merged[, ..keep_cols]
chr_cols <- names(merged)[vapply(merged, is.character, logical(1))]
merged[, (chr_cols) := lapply(.SD, function(x) fifelse(is.na(x), "", x)), .SDcols = chr_cols]

# ---- one row per allele x sequence (plus its RSS/leader variant) ------------------

by_allele <- setdiff(names(merged), c("vdjbase_subject", "subject_id", "genomic_dataset", "notes", "z_score", "sure_subject",
                                      "non_sure_subject", "in_genomic", "in_genomic_watson", "in_genomic_hprc", "in_genomic_1kpg",
                                      "in_AIRRseq", "repertoire_allele", "vdjbase_allele"))
per_allele <- merged[, .(
  repertoire_allele = collapse_unique(repertoire_allele),
  vdjbase_allele = collapse_unique(vdjbase_allele),
  samples_genomic = paste(unique(vdjbase_subject[in_genomic_watson == TRUE & vdjbase_subject != ""]), collapse = ", "),
  samples_genomic_watson = paste(unique(vdjbase_subject[in_genomic_watson == TRUE & vdjbase_subject != ""]), collapse = ", "),
  samples_genomic_hprc = paste(unique(subject_id[in_genomic_hprc == TRUE & subject_id != ""]), collapse = ", "),
  samples_genomic_1kpg = paste(unique(subject_id[in_genomic_1kpg == TRUE & subject_id != ""]), collapse = ", "),
  samples_AIRRseq = paste(unique(vdjbase_subject[in_AIRRseq == TRUE & vdjbase_subject != ""]), collapse = ", "),
  sure_subject = paste(unique(sure_subject[sure_subject != "" & in_genomic_watson == TRUE]), collapse = ", "),
  non_sure_subject = paste(unique(non_sure_subject[non_sure_subject != "" & in_genomic_watson == TRUE]), collapse = ", "),
  notes = collapse_unique(notes, sep = "; "),
  z_score = collapse_unique(as.character(z_score[!is.na(z_score)]), sep = "; ")
), by = by_allele]

baseline_collapsed <- baseline[, .(wasp = collapse_unique(wasp)), by = .(allele, gene_type, iglabel, seq, seq_gapped, in_baseline_reference)]
combined <- merge(per_allele, baseline_collapsed, by = c("allele", "gene_type", "seq"), all = TRUE)
setnames(combined, c("seq_gapped.x", "seq_gapped.y"), c("seq_gapped_genomic", "seq_gapped"))
combined[seq_gapped == "" | is.na(seq_gapped) | seq_gapped == "NA", seq_gapped := seq_gapped_genomic]
combined[, seq_gapped_genomic := NULL]

by_combined <- setdiff(names(combined), c("repertoire_allele", "vdjbase_allele", "samples_genomic", "samples_genomic_watson",
                                          "samples_genomic_hprc", "samples_genomic_1kpg", "samples_AIRRseq", "sure_subject",
                                          "non_sure_subject", "notes", "z_score", "iglabel", "wasp"))
combined <- combined[, .(
  repertoire_allele = collapse_unique(repertoire_allele),
  vdjbase_allele = collapse_unique(vdjbase_allele),
  samples_genomic = collapse_tokens(samples_genomic),
  samples_genomic_watson = collapse_tokens(samples_genomic_watson),
  samples_genomic_hprc = collapse_tokens(samples_genomic_hprc),
  samples_genomic_1kpg = collapse_tokens(samples_genomic_1kpg),
  samples_AIRRseq = collapse_tokens(samples_AIRRseq),
  sure_subject = collapse_tokens(sure_subject),
  non_sure_subject = collapse_tokens(non_sure_subject),
  notes = collapse_unique(notes, sep = "; "),
  z_score = collapse_unique(as.character(z_score), sep = "; "),
  iglabel = collapse_unique(iglabel),
  wasp = collapse_unique(wasp)
), by = by_combined]
combined[, `:=`(sample_count_genomic = count_tokens(samples_genomic),
                sample_count_genomic_watson = count_tokens(samples_genomic_watson),
                sample_count_genomic_hprc = count_tokens(samples_genomic_hprc),
                sample_count_genomic_1kpg = count_tokens(samples_genomic_1kpg),
                sample_count_AIRRseq = count_tokens(samples_AIRRseq),
                sure_subject_count = count_tokens(sure_subject),
                non_sure_subject_count = count_tokens(non_sure_subject))]
combined[is.na(in_baseline_reference), in_baseline_reference := FALSE]
combined[, iglabel := short_iglabel(iglabel)]
combined[, `:=`(asc = sub("^.*-", "", sub("\\*.*$", "", allele)), chain = substr(allele, 1, 3))]
combined[, `:=`(spacer_3_size = nchar(spacer_3), spacer_5_size = nchar(spacer_5))]
combined[, rss_notes := fifelse(
  !(spacer_3_size %in% c(0, 12, 23)) & !(spacer_5_size %in% c(0, 12, 23)), "Both spacer_3 and spacer_5 have variation at the size",
  fifelse(!(spacer_3_size %in% c(0, 12, 23)), "spacer_3 has variation at the size",
          fifelse(!(spacer_5_size %in% c(0, 12, 23)), "spacer_5 has variation at the size", "")))]
combined[, present := sample_count_AIRRseq > 0]
combined[, multi_genomic_location := uniqueN(vdjbase_allele) > 1, by = seq]

# ---- one row per sequence: the seed set for IgLabel and the ASC clustering ---------

# The IUIS-style name is the VDJbase/WASP name; a novel allele inherits its base
# allele's names, and an AIRR-seq-only novel allele is named by its SNPs from the base.
seed <- copy(combined)
seed[, iuis_allele := apply(.SD, 1, function(x) {
  v <- unique(trimws(unlist(strsplit(paste(x, collapse = ","), ","))))
  paste(v[nzchar(v) & !is.na(v) & v != "NA"], collapse = ",")
}), .SDcols = c("vdjbase_allele", "wasp")]
seed[in_baseline_reference == FALSE & iuis_allele == "", iuis_allele := sapply(allele, function(x) {
  base <- gsub("_.*", "", x)
  bases <- trimws(unique(unlist(strsplit(seed$iuis_allele[seed$allele == base], ","))))
  paste(sapply(bases[nzchar(bases)], function(a) gsub(base, a, x, fixed = TRUE)), collapse = ",")
})]
seed[sample_count_genomic == 0 & sample_count_AIRRseq > 0 & in_baseline_reference == FALSE, iuis_allele := sapply(allele, function(x) {
  base <- gsub("_.*", "", x)
  bases <- trimws(unique(unlist(strsplit(seed$iuis_allele[seed$allele == base], ","))))
  paste(sapply(bases[nzchar(bases)], function(a) {
    seqs_base <- unique(seed$seq_gapped[grepl(a, seed$iuis_allele, fixed = TRUE)])
    seqs_novel <- unique(seed$seq_gapped[seed$allele == x])
    paste0(a, "_", paste0(piglet::allele_diff_strings(germs = c(seqs_base, seqs_novel), X = 0), collapse = "_"))
  }), collapse = ",")
})]
seed[, seq_gapped_selected := select_gap_sequence(seq_gapped, iuis_allele), by = seq]
seed[, allele_selected := select_allele(allele), by = seq]
seed <- seed[, .(
  allele = collapse_unique(allele_selected), asc = collapse_unique(asc),
  seq_gapped = collapse_unique(seq_gapped_selected), chain = collapse_unique(chain),
  gene_type = collapse_unique(gene_type),
  samples_AIRRseq = collapse_tokens(samples_AIRRseq), samples_genomic = collapse_tokens(samples_genomic),
  samples_genomic_watson = collapse_tokens(samples_genomic_watson),
  samples_genomic_hprc = collapse_tokens(samples_genomic_hprc),
  samples_genomic_1kpg = collapse_tokens(samples_genomic_1kpg),
  sure_subject = collapse_tokens(sure_subject), non_sure_subject = collapse_tokens(non_sure_subject),
  vdjbase_allele = collapse_unique(vdjbase_allele), wasp = collapse_unique(wasp),
  rss_notes = collapse_unique(rss_notes), notes = collapse_unique(notes), iglabel = collapse_unique(iglabel),
  novel = any(grepl("_", allele))
), by = .(seq)]
seed[, `:=`(sample_count_AIRRseq = count_tokens(samples_AIRRseq), sample_count_genomic = count_tokens(samples_genomic),
            sample_count_genomic_watson = count_tokens(samples_genomic_watson),
            sample_count_genomic_hprc = count_tokens(samples_genomic_hprc),
            sample_count_genomic_1kpg = count_tokens(samples_genomic_1kpg),
            sure_subject_count = count_tokens(sure_subject), non_sure_subject_count = count_tokens(non_sure_subject))]
seed[!grepl("_", allele), allele_seed := allele]
cat(sprintf("seed set: %d sequences (%d novel)\n", nrow(seed), sum(seed$novel)))

# ---- IgLabel labels (frozen review set) ---------------------------------------------

labels <- rbindlist(lapply(CHAINS, function(ch) {
  fread(need(file.path(IN$iglabel_labels, paste0(ch, ".labels.csv"))))[, .(allele = seq_id, iglabel = short_iglabel(iglabel_raw))]
}))
labels <- labels[nzchar(iglabel)][!duplicated(allele)]
label_of <- setNames(labels$iglabel, labels$allele)

labelled <- copy(seed)
labelled[(is.na(iglabel) | iglabel == "") & allele %in% names(label_of), iglabel := label_of[allele]]
labelled[, vdjbase_allele := setNames(combined$vdjbase_allele, combined$seq)[seq]]
fallback <- seed[!grepl("_", allele) & nzchar(iglabel)][!duplicated(seq)]
labelled[is.na(iglabel) | iglabel == "", iglabel := setNames(fallback$iglabel, fallback$seq)[seq]]
labelled[, iglabel := short_iglabel(iglabel)]
chr_cols <- names(labelled)[vapply(labelled, is.character, logical(1))]
labelled[, (chr_cols) := lapply(.SD, function(x) fifelse(is.na(x), "", x)), .SDcols = chr_cols]

labelled[, iuis_allele := apply(.SD, 1, function(x) {
  v <- unique(trimws(unlist(strsplit(paste(x, collapse = ","), ","))))
  paste(v[nzchar(v)], collapse = ",")
}), .SDcols = c("vdjbase_allele", "wasp")]
seed_iuis <- setNames(labelled$iuis_allele[!grepl("_", labelled$allele)], labelled$allele[!grepl("_", labelled$allele)])
labelled[grepl("_", allele) & iuis_allele == "", iuis_allele := sapply(allele, function(x) {
  base <- gsub("_.*", "", x)
  bases <- trimws(unique(unlist(strsplit(seed_iuis[base], ","))))
  paste(sapply(bases[nzchar(bases)], function(a) gsub(base, a, x, fixed = TRUE)), collapse = ",")
})]
labelled[sample_count_genomic == 0 & sample_count_AIRRseq > 0 & novel == TRUE, iuis_allele := sapply(allele, function(x) {
  base <- gsub("_.*", "", x)
  bases <- trimws(unique(unlist(strsplit(labelled$iuis_allele[labelled$allele == base], ","))))
  paste(sapply(bases[nzchar(bases)], function(a) {
    seqs_base <- unique(labelled$seq_gapped[grepl(a, labelled$iuis_allele, fixed = TRUE)])
    seqs_novel <- unique(labelled$seq_gapped[labelled$allele == x])
    paste0(a, "_", paste0(piglet::allele_diff_strings(germs = c(seqs_base, seqs_novel), X = 0), collapse = "_"))
  }), collapse = ",")
})]
labelled[, seq_gapped_selected := seq_gapped]
labelled[(duplicated(seq) & !duplicated(seq_gapped)) | (duplicated(seq, fromLast = TRUE) & !duplicated(seq_gapped, fromLast = TRUE)),
         seq_gapped_selected := select_gap_sequence(seq_gapped, iuis_allele), by = seq]

usofa <- labelled[, .(
  samples_AIRRseq = collapse_tokens(samples_AIRRseq), samples_genomic = collapse_tokens(samples_genomic),
  samples_genomic_watson = collapse_tokens(samples_genomic_watson),
  samples_genomic_hprc = collapse_tokens(samples_genomic_hprc),
  samples_genomic_1kpg = collapse_tokens(samples_genomic_1kpg),
  sure_subject = collapse_tokens(sure_subject), non_sure_subject = collapse_tokens(non_sure_subject),
  iglabel = collapse_unique(iglabel, sep = ","), rss_notes = collapse_unique(rss_notes, sep = ",")
), by = .(allele, asc, gene_type, chain, seq, seq_gapped = seq_gapped_selected, vdjbase_allele, wasp, notes, novel, allele_seed, iuis_allele)]
usofa[, `:=`(sample_count_AIRRseq = count_tokens(samples_AIRRseq), sample_count_genomic = count_tokens(samples_genomic),
             sample_count_genomic_watson = count_tokens(samples_genomic_watson),
             sample_count_genomic_hprc = count_tokens(samples_genomic_hprc),
             sample_count_genomic_1kpg = count_tokens(samples_genomic_1kpg),
             sure_subject_count = count_tokens(sure_subject), non_sure_subject_count = count_tokens(non_sure_subject))]
usofa[, iglabel_allele := paste0(gene_type, "1-", iglabel, "*00")]
if (any(!nzchar(trimws(usofa$iglabel)))) stop("alleles without an IgLabel label: ", sum(!nzchar(trimws(usofa$iglabel))))
if (any(grepl(",", usofa$iglabel, fixed = TRUE))) stop("alleles with more than one IgLabel label")

# ---- ASC clustering per segment ------------------------------------------------------

usofa[, iglabel_label := vapply(seq_len(.N), function(i) {
  x <- trimws(iglabel[i]); if (!nzchar(x)) x <- trimws(iglabel_allele[i]); if (!nzchar(x)) x <- trimws(allele[i])
  x <- sub("\\*[0-9]+$", "", sub("^IG[HKL][VDJ]1-", "", x))
  if (!nzchar(x)) x <- sub("\\*.*$", "", sub("^IG[HKL][VDJ]1-", "", trimws(allele[i])))
  x
}, character(1))]
usofa[, iglabel := iglabel_label]
usofa[, iglabel_allele := sprintf("%s1-%s*00", gene_type, iglabel_label)]
usofa[, seq_gapped_selected := seq_gapped]
usofa[is.na(seq_gapped_selected) | !nzchar(trimws(seq_gapped_selected)), seq_gapped_selected := seq]

wasp_of <- unique(combined[, .(seq, chain, gene_type, wasp)])[!duplicated(paste(seq, chain, gene_type))]
usofa <- merge(usofa, wasp_of, by = c("seq", "chain", "gene_type"), all.x = TRUE, suffixes = c("", ".combined"))
usofa[(is.na(wasp) | wasp == "") & !is.na(wasp.combined) & wasp.combined != "", wasp := wasp.combined]
usofa[, wasp.combined := NULL]

first_nonempty_trim <- function(x) { x <- unique(trimws(x[!is.na(x)])); x <- x[nzchar(x)]; if (length(x)) x[1] else "" }
cross <- rbindlist(lapply(intersect(c("IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ"), unique(usofa$gene_type)), function(seg) {
  seg_dt <- usofa[gene_type == seg]
  if (any(seg_dt[, uniqueN(seq), by = allele]$V1 > 1L)) stop("allele with two sequences in ", seg)
  seg_dt <- seg_dt[!duplicated(paste(allele, seq))]
  germline <- setNames(seg_dt$seq_gapped_selected, seg_dt$allele)
  asc <- piglet::inferAlleleClusters(
    germline_set = germline, locus = seg, clustering_method = "leiden",
    distance_method = if (substr(seg, 4, 4) == "V") "hamming" else "lv",
    family_threshold = max(nchar(gsub("[.]", "", germline))) * 0.25,
    trim_3prime_side = NULL, mask_5prime_side = 0, optimize_silhouette = TRUE,
    family_prefix = FALSE, ncores = 1L, quiet = TRUE)
  tab <- asc$alleleClusterTable
  allele_col <- intersect(c("imgt_allele", "iuis_allele"), colnames(tab))[1]
  info <- seg_dt[, .(iglabel = first_nonempty_trim(iglabel_label), iglabel_allele = first_nonempty_trim(iglabel_allele),
                     husa = collapse_unique(trimws(iuis_allele), sep = ","), wasp = collapse_unique(trimws(wasp), sep = ","),
                     seq = first_nonempty_trim(seq), seq_gapped = first_nonempty_trim(seq_gapped_selected)), by = allele]
  out <- data.table(allele = tab[[allele_col]], asc_allele_new = gsub("F", "", tab$new_allele), gene_type = seg, chain = substr(seg, 1, 3))
  cat(sprintf("%s: %d alleles -> %d ASC alleles\n", seg, length(germline), nrow(out)))
  merge(out, info, by = "allele", sort = FALSE)
}))

# Final names: HUSA family from the IUIS names, subgroup from the IgLabel label of the
# *01 member of each ASC subgroup, allele numbers renumbered within the subgroup.
gene_from_multi <- function(group, subgroup) {
  if (is.na(subgroup) || subgroup == "") return(subgroup)
  if (grepl("D-", subgroup, fixed = TRUE)) paste0(gsub(paste0(unique(group), "D?-"), "", subgroup), "D")
  else gsub(paste0(unique(group), "-"), "", subgroup)
}
cross[, husa_group := alakazam::getFamily(husa, strip_d = TRUE, omit_nl = FALSE, first = FALSE, collapse = TRUE)]
cross[, husa_subgroup := alakazam::getGene(husa, strip_d = FALSE, omit_nl = FALSE, first = FALSE, collapse = TRUE)]
cross[!(gene_type %in% c("IGHJ", "IGKJ", "IGLJ")), husa_subgroup := gene_from_multi(unique(husa_group), husa_subgroup), by = husa_subgroup]
cross[, asc_group := alakazam::getFamily(asc_allele_new, strip_d = FALSE, omit_nl = FALSE)]
cross[, asc_subgroup := alakazam::getGene(asc_allele_new, strip_d = FALSE, omit_nl = FALSE)]
cross[, asc_subgroup := gsub(paste0(unique(asc_group), "-"), "", asc_subgroup), by = asc_subgroup]
cross[, resolved_subgroup_name := iglabel[which(stringr::str_detect(asc_allele_new, "\\*01"))][1], by = .(asc_subgroup, gene_type)]
cross[is.na(resolved_subgroup_name) | resolved_subgroup_name == "", resolved_subgroup_name := iglabel[1], by = .(asc_subgroup, gene_type)]
cross[, final_group_subgroup_name := paste0(husa_group, "-", resolved_subgroup_name)]
cross[grepl("IGLJ[23],", final_group_subgroup_name), final_group_subgroup_name := paste0(gene_type, "2/3-", resolved_subgroup_name)]
cross[, allele_number := stringr::str_remove(stringr::str_extract(asc_allele_new, "\\*\\d+$"), "\\*")]
setorder(cross, chain, gene_type, final_group_subgroup_name, asc_allele_new, iglabel, husa, seq_gapped)
cross[, allele_number := {
  n <- suppressWarnings(as.integer(trimws(allele_number)))
  if (anyNA(n) || anyDuplicated(n[!is.na(n)]) > 0) sprintf("%02d", seq_len(.N)) else sprintf("%02d", n)
}, by = final_group_subgroup_name]
cross[, final_allele := paste0(final_group_subgroup_name, "*", allele_number)]

coalesce_chr <- function(...) { v <- list(...); out <- rep(NA_character_, length(v[[1]])); for (x in v) { i <- is.na(out) | out == ""; out[i] <- x[i] }; out }
assignments <- unique(cross[, .(seq, seq_gapped, chain, gene_type, final_allele, asc = resolved_subgroup_name)])
assignments[, seq_gapped_match := coalesce_chr(seq_gapped, seq)]
husa <- unique(usofa[, .(seq, chain, gene_type, seq_gapped_selected, iuis_allele, iglabel, iglabel_allele,
                         sample_count_AIRRseq, samples_AIRRseq, novel, allele_seed)])
husa[, seq_gapped_match := coalesce_chr(seq_gapped_selected, seq)]
husa <- merge(husa, assignments[, .(seq, seq_gapped_match, chain, gene_type, final_allele, asc)],
              by = c("seq", "seq_gapped_match", "chain", "gene_type"), all.x = TRUE)
annot <- copy(combined)[, intersect(c("allele", "asc", "sample_count_AIRRseq", "samples_AIRRseq", "iglabel", "novel", "allele_seed"), names(combined)) := NULL]
husa <- merge(husa, annot, by = c("seq", "chain", "gene_type"), all.x = TRUE)
setnames(husa, c("final_allele", "iuis_allele"), c("allele", "husa"))
husa[, seq_gapped := do.call(coalesce_chr, .SD), .SDcols = c("seq_gapped_selected", "seq_gapped", "seq")]
husa[, c("seq_gapped_selected", "seq_gapped_match") := NULL]
husa[, asc := coalesce_chr(asc, iglabel)]
lead <- c("allele", "asc", "gene_type", "chain", "functional", "novel", "present", "husa", "vdjbase_allele", "wasp", "iglabel_allele", "seq", "seq_gapped")
setcolorder(husa, c(lead, setdiff(names(husa), lead)))
husa <- husa[order(gene_type, allele, na.last = TRUE)]

fwrite(husa, file.path(OUT$husa, "husa.tsv"), sep = "\t", quote = FALSE)
rss_ok <- !grepl("NONAMER not found|RSS not found", husa$notes, ignore.case = TRUE) &
  husa$sample_count_genomic > 0 & husa$sample_count_AIRRseq > 0 & husa$sure_subject_count > 0
fwrite(husa[rss_ok], file.path(OUT$husa, "husa_rss_filter.tsv"), sep = "\t")
cat(sprintf("husa.tsv: %d alleles (%d with RSS evidence)\n", nrow(husa), sum(rss_ok)))
