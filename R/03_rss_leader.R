#!/usr/bin/env Rscript
# RSS and leader table behind every RSS/leader figure: one row per allele x distinct
# RSS/leader variant seen in the in-house genomic cohort, with short spacers gap-aligned
# to the canonical length and leaders repaired for known annotation artefacts.
#
#   results/rss_leader/rss_leader_iuis_data.csv.gz

source("R/00_setup.R")
source("R/lib/rss_helpers.R")
suppressPackageStartupMessages({ library(stringi); library(Biostrings) })

complement_reverse <- function(seqs) {
  seqs <- as.character(seqs); seqs[is.na(seqs) | !nzchar(seqs)] <- NA_character_
  out <- rep(NA_character_, length(seqs)); keep <- !is.na(seqs)
  if (any(keep)) out[keep] <- as.character(reverseComplement(DNAStringSet(seqs[keep])))
  out
}

gd <- fread(need(IN$ggs_watson))

# Leader-1 fixes: an ATG start shifted by one base, a trailing base, or two ATGs in the
# annotated exon; leader-2 fixes for three recurrent mis-annotations.
gd[gene_type == "IGKV" & grepl("G$", l_part1), l_part1 := substr(l_part1, 1, nchar(l_part1) - 1)]
gd[, l_part1_length := nchar(l_part1)]
gd[, starts_tg := grepl("^TG", l_part1)]
gd[, starts_catg := grepl("^CATG", l_part1)]
gd[, multi_atg := vapply(stri_locate_all(regex = "ATG", str = l_part1), nrow, integer(1)) > 1L]
gd[, l_part1_fixed := l_part1]
gd[l_part1_length == 46 & starts_tg & gene_type == "IGHV", l_part1_fixed := paste0("A", substr(l_part1, 1, 45))]
gd[l_part1_length == 45 & starts_tg & gene_type == "IGLV", l_part1_fixed := paste0("A", l_part1)]
gd[l_part1_length == 44 & starts_tg & gene_type == "IGLV", l_part1_fixed := paste0("A", l_part1)]
gd[l_part1_length == 47 & starts_catg & gene_type == "IGHV", l_part1_fixed := substr(l_part1, 2, 47)]
gd[l_part1_length == 49 & starts_catg & gene_type == "IGKV", l_part1_fixed := substr(l_part1, 2, nchar(l_part1))]
gd[l_part1_length == 47 & !grepl("A$", l_part1) & gene_type == "IGKV", l_part1_fixed := paste0(l_part1, "A")]
trim_multi_atg <- function(seq, l = 46) {
  loc <- stri_locate_all(regex = "ATG", str = seq)[[1]]
  for (i in seq_len(nrow(loc))) { sq <- substr(seq, loc[i, 1], nchar(seq)); if (nchar(sq) == l) return(sq) }
  seq
}
gd[l_part1_length > 46 & multi_atg & gene_type == "IGHV", l_part1_fixed := trim_multi_atg(l_part1), by = l_part1]
gd[l_part1_length > 48 & multi_atg & gene_type == "IGKV", l_part1_fixed := trim_multi_atg(l_part1, l = 48), by = l_part1]
gd[nchar(l_part1_fixed) > 46 & gene_type == "IGHV", l_part1_fixed := substr(l_part1_fixed, 1, 46)]
gd[nchar(l_part1_fixed) > 46 & gene_type == "IGLV", l_part1_fixed := substr(l_part1_fixed, 1, 46)]
gd[, l_part2_fixed := l_part2]
gd[gene_type == "IGLV" & l_part2 == "GGTCCTGGGC",    l_part2_fixed := "GGTCCTGGGCC"]
gd[gene_type == "IGLV" & l_part2 == "GGTCTCTCTCTCC", l_part2_fixed := "GGTCTCTCTCC"]
gd[gene_type == "IGKV" & l_part2 == "ACTACCAGATGT",  l_part2_fixed := "CTACCAGATGT"]
gd[, spacer_3_extended := as.character(spacer_3_extended)]
gd[, spacer_5_extended := as.character(spacer_5_extended)]
gd[gene_type %in% c("IGKV", "IGLV"), spacer_3_extended := spacer_3]
gd[gene_type %in% c("IGKJ", "IGLJ"), spacer_5_extended := spacer_5]

# One row per subject and allele; J and 5' D RSSs are reverse complemented so every RSS
# reads heptamer -> spacer -> nonamer.
mk <- function(dt, seg, hept, nona, spac, spacx, l1 = NA_character_, l2 = NA_character_) {
  data.table(gene_type = seg, iuis_subgroup = sub("-.*$", "", sub("\\*.*$", "", dt$vdjbase_allele)),
             iuis_group = sub("\\*.*$", "", dt$vdjbase_allele),
             heptamer = hept, nonamer = nona, spacer = spac, spacer_extended = spacx,
             l_part1 = l1, l_part2 = l2, vdjbase_subject = dt$vdjbase_subject, vdjbase_allele = dt$vdjbase_allele)
}
V <- gd[gene_type %in% c("IGHV", "IGKV", "IGLV")]; J <- gd[gene_type %in% c("IGHJ", "IGKJ", "IGLJ")]; D <- gd[gene_type == "IGHD"]
raw <- rbindlist(list(
  mk(V, V$gene_type, V$v_heptamer, V$v_nonamer, V$spacer_3, V$spacer_3_extended, V$l_part1_fixed, V$l_part2_fixed),
  mk(J, J$gene_type, complement_reverse(J$j_heptamer), complement_reverse(J$j_nonamer),
     complement_reverse(J$spacer_5), complement_reverse(J$spacer_5_extended)),
  mk(D, "IGHD_5", complement_reverse(D$d_5_heptamer), complement_reverse(D$d_5_nonamer),
     complement_reverse(D$spacer_5), complement_reverse(D$spacer_5_extended)),
  mk(D, "IGHD_3", D$d_3_heptamer, D$d_3_nonamer, D$spacer_3, D$spacer_3_extended)
), use.names = TRUE)
n_raw <- nrow(raw)
raw <- raw[nchar(heptamer) == 7 & nchar(nonamer) == 9]
cat(sprintf("RSS rows: %d of %d with a 7-mer heptamer and 9-mer nonamer\n", nrow(raw), n_raw))

rss <- raw[, .(n_samples = uniqueN(vdjbase_subject)),
           by = .(gene_type, vdjbase_allele, iuis_subgroup, iuis_group, heptamer, nonamer, spacer, spacer_extended, l_part1, l_part2)
][order(gene_type, vdjbase_allele, iuis_subgroup, iuis_group)]

# Spacers shorter than the canonical 12/23 are gap-aligned against the canonical spacers
# of their own subgroup; the two IGHV5 spacers have a fixed placement.
rss <- rss[, {
  spacer_len <- nchar(spacer)
  if (uniqueN(spacer_len) == 1L) {
    spacer_aligned <- spacer
  } else {
    is_canon <- spacer_len %in% c(12L, 23L)
    refs_by_sub <- split(spacer[is_canon], iuis_subgroup[is_canon])
    short_tbl <- unique(.SD[!is_canon, .(spacer, iuis_subgroup)])
    aligned_map <- setNames(short_tbl$spacer, short_tbl$spacer)
    for (i in seq_len(nrow(short_tbl))) {
      tgt <- short_tbl$spacer[i]; sub <- short_tbl$iuis_subgroup[i]
      if (!sub %in% names(refs_by_sub)) {
        if (gene_type == "IGHV" & grepl("V5", sub)) refs <- refs_by_sub[[grep("V[2]", names(refs_by_sub))]] else refs <- refs_by_sub
      } else refs <- refs_by_sub[[sub]]
      if (!is.null(refs) && length(refs) > 0 & gene_type == "IGHV" & grepl("V5", sub)) {
        if (tgt == "AGAGAAACCAGCACTGAGCCCG") aligned_map["AGAGAAACCAGCACTGAGCCCG"] <- "AGAGAAACCAGCACTGAGCCC-G"
        if (tgt == "AGAGAAACCAGCCCCGAGCCCG") aligned_map["AGAGAAACCAGCCCCGAGCCCG"] <- "AGAGAAACCAGCCCCGAGCCC-G"
      } else if (!is.null(refs) && length(refs) > 0) {
        al <- align_target_to_references(target = tgt, references = refs)
        aligned_map[tgt] <- if (!is.null(al$target)) al$target else tgt
      }
    }
    spacer_aligned <- spacer; sel <- spacer %in% names(aligned_map); spacer_aligned[sel] <- aligned_map[spacer[sel]]
  }
  .(vdjbase_allele, iuis_subgroup, iuis_group, heptamer, nonamer, spacer = spacer_aligned, n_samples,
    rss_aligned = paste0(heptamer, spacer_aligned, nonamer), spacer_extended, l_part1, l_part2,
    leader = ifelse(is.na(l_part1), NA_character_, paste0(l_part1, l_part2)))
}, by = gene_type]

# Leader-1 sequences off the canonical length (46 for IGHV, the longest seen otherwise)
# are padded or truncated to it against the closest canonical-length leader.
align_to_length <- function(target, references, target_length) {
  target_split <- strsplit(target, "")[[1]]; best_alignment <- NULL; best_mismatches <- Inf
  for (reference in references) {
    reference_split <- strsplit(reference, "")[[1]]
    if (length(reference_split) < target_length) reference_split <- c(reference_split, rep("-", target_length - length(reference_split)))
    else if (length(reference_split) > target_length) reference_split <- reference_split[1:target_length]
    n <- target_length
    padded_target <- if (length(target_split) >= n) target_split[1:n] else c(target_split, rep("-", n - length(target_split)))
    mismatches <- sum(padded_target != reference_split)
    if (mismatches < best_mismatches) { best_alignment <- padded_target; best_mismatches <- mismatches }
  }
  list(target = paste(best_alignment, collapse = ""))
}
rss <- rss[, {
  canon_len <- if (gene_type %in% c("IGHV")) 46 else max(nchar(l_part1))
  l1_aligned <- l_part1
  if (any(nchar(l_part1) != canon_len, na.rm = TRUE)) {
    refs <- l_part1[nchar(l_part1) == canon_len & !is.na(l_part1)]
    if (length(refs) > 0) {
      for (i in which(nchar(l_part1) != canon_len & !is.na(l_part1))) l1_aligned[i] <- align_to_length(l_part1[i], refs, target_length = canon_len)$target
    } else {
      l1_aligned <- vapply(l_part1, function(seq) {
        if (is.na(seq)) return(NA_character_)
        if (nchar(seq) >= canon_len) substr(seq, 1, canon_len) else paste0(seq, strrep("-", canon_len - nchar(seq)))
      }, character(1))
    }
  }
  .(vdjbase_allele, iuis_subgroup, iuis_group, heptamer, nonamer, spacer, n_samples, rss_aligned, spacer_extended,
    l_part1, l_part1_aligned = l1_aligned, l_part2, leader = ifelse(is.na(l1_aligned), NA_character_, paste0(l1_aligned, l_part2)))
}, by = gene_type]

# Flag the alleles expressed in the AIRR-seq cohort and attach the coding sequence.
husa <- fread(need(file.path(OUT$husa, "husa.tsv")))[, .(gene_type, seq, iuis_allele = husa, present)]
husa <- husa[, .(iuis_allele = unlist(strsplit(iuis_allele, ","))), by = .(gene_type, seq, present)]
husa <- husa[, .(present = any(present)), by = .(gene_type, iuis_allele, seq)]
rss <- merge(rss, husa[, .(iuis_allele, present, seq)], by.x = "vdjbase_allele", by.y = "iuis_allele", all.x = TRUE)
rss[, iuis_group := gsub("IG[KLH]", "", iuis_group)]
rss[, iuis_subgroup := gsub("D", "", gsub("IG[KLH]", "", iuis_subgroup))]
rss[, l_part1_artifact := !is.na(l_part1) & !grepl("^ATG", l_part1)]

fwrite(rss, file.path(OUT$rss_leader, "rss_leader_iuis_data.csv.gz"))
cat(sprintf("rows: %d | genes with >1 expressed leader: %d\n", nrow(rss),
            rss[!is.na(leader) & present == TRUE, uniqueN(leader), by = .(gene_type, iuis_group)][V1 > 1, .N]))
