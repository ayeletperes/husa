# RSS alignment helpers shared by the rss_leader stage and the RSS figures.

complement_reverse <- function(seqs) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seqs)))
}

compute_reference_distance <- function(sequences, method = "lv", trim_prim_3 = NULL, type = "both", quiet = FALSE) {
  if (!is.null(trim_prim_3)) sequences <- substr(sequences, 1, trim_prim_3)
  if (method == "hamming") {
    max_len <- max(nchar(sequences))
    sequences <- gsub("\\s", "N", format(sequences, width = max_len))
  }
  dist <- stringdist::stringdistmatrix(sequences, method = method, useNames = "names")
  if (type == "both") return(list(dist = dist, matrix = as.matrix(dist)))
  if (type == "matrix") return(as.matrix(dist))
  dist
}

# Place one gap in a short spacer where it best matches a canonical-length reference.
align_target_to_references <- function(target, references) {
  target_split <- strsplit(target, "", fixed = TRUE)[[1]]
  best_alignment <- NULL
  best_mismatches <- Inf
  best_reference <- NULL
  for (reference in references) {
    reference_split <- strsplit(reference, "", fixed = TRUE)[[1]]
    n <- max(length(target_split), length(reference_split))
    padded_target <- c(target_split, rep("-", n - length(target_split)))
    padded_reference <- c(reference_split, rep("-", n - length(reference_split)))
    best_alignment_for_ref <- padded_target
    best_mismatches_for_ref <- sum(padded_target != padded_reference)
    for (gap_pos in seq_len(n)) {
      temp_target <- c(padded_target[1:(gap_pos - 1)], "-", padded_target[gap_pos:n])[1:n]
      temp_mismatches <- sum(temp_target != padded_reference)
      if (temp_mismatches < best_mismatches_for_ref) {
        best_alignment_for_ref <- temp_target
        best_mismatches_for_ref <- temp_mismatches
      }
    }
    if (best_mismatches_for_ref < best_mismatches) {
      best_alignment <- best_alignment_for_ref
      best_mismatches <- best_mismatches_for_ref
      best_reference <- reference_split
    }
  }
  list(target = paste(best_alignment, collapse = ""), reference = paste(best_reference, collapse = ""))
}

# One segment's RSS table from the husa rows: heptamer, spacer, nonamer split out, and
# short spacers gap-aligned to the canonical length within their family.
extract_rss_table <- function(data, cols, spacer_size = c(23, 22), canonical_size = NULL) {
  if (is.null(canonical_size)) canonical_size <- intersect(c(23, 12), spacer_size)
  rss <- data[!grepl("NONAMER not found|RSS not found", notes, ignore.case = TRUE) &
                sample_count_genomic > 0 & sample_count_AIRRseq > 0]
  if (nrow(rss) == 0L) return(data.table::data.table())

  rss[, rss := do.call(paste, c(.SD, sep = ",")), .SDcols = cols]
  keep_cols <- intersect(c("allele", "vdjbase_allele", "iglabel", "gene_type", "chain", "rss", "sure_subject", "sure_subject_count",
                           "l_part1", "l_part2", "sample_count_genomic", "sample_count_AIRRseq"), names(rss))
  rss <- data.table::melt(rss[, ..keep_cols], id.vars = setdiff(keep_cols, "rss"), measure.vars = "rss",
                          variable.name = "rss_type", value.name = "rss_sequence")
  rss <- rss[rss_sequence != "," & nchar(rss_sequence) %in% (spacer_size + 7 + 9 + 2)]
  if (nrow(rss) == 0L) return(data.table::data.table())

  rss[, c("heptamer", "spacer", "nonamer") := data.table::tstrsplit(rss_sequence, ",", fixed = TRUE)]
  rss <- unique(rss[nchar(heptamer) == 7 & nchar(nonamer) == 9])
  rss[, rss := gsub(",", "", rss_sequence)]
  rss[, asc := allele]
  rss[, family := sub("-.*$", "", asc)]
  rss[, rss_aligned := rss]
  if (all(nchar(rss$spacer) %in% canonical_size)) return(rss)

  spacers <- setNames(unique(rss$spacer), unique(rss$spacer))
  spacer_dist <- compute_reference_distance(spacers, method = "lv", type = "both", quiet = TRUE)
  short_spacers <- names(spacers)[nchar(spacers) < canonical_size]
  if (length(short_spacers) == 0L) return(rss)

  # Each family's reference spacers are the canonical-length ones nearest to its short ones.
  short_by_family <- unique(rss[spacer %in% short_spacers, .(spacer, family)])
  references <- lapply(unique(short_by_family$family), function(fam) {
    candidates <- setdiff(colnames(spacer_dist$matrix), short_spacers)
    if (length(candidates) == 0L) return(character())
    hits <- unlist(lapply(short_by_family[family == fam, spacer], function(sp) {
      x <- spacer_dist$matrix[sp, candidates]
      names(x)[x <= (min(x) + 2)]
    }), use.names = FALSE)
    tab <- table(hits)
    names(tab[tab == max(tab)])
  })
  names(references) <- unique(short_by_family$family)

  aligned <- vapply(seq_len(nrow(short_by_family)), function(i) {
    refs <- references[[short_by_family$family[[i]]]]
    if (length(refs) == 0L) return(short_by_family$spacer[[i]])
    align_target_to_references(target = short_by_family$spacer[[i]], references = refs)$target
  }, character(1))
  names(aligned) <- short_by_family$spacer

  rss[, spacer_aligned := spacer]
  rss[spacer %in% names(aligned), spacer_aligned := aligned[spacer]]
  rss[, rss_aligned := paste0(heptamer, spacer_aligned, nonamer)]
  rss[]
}

# Every segment's aligned RSS table, 3' RSS for V, 5' RSS (reverse complemented) for J
# and both sides for D, with the spacer length class each segment is drawn at.
prepare_aligned_rss_table <- function(husa_dt) {
  dt <- data.table::copy(husa_dt)
  dt[, chain := toupper(chain)]
  dt[, gene_type := toupper(gene_type)]
  igh <- dt[chain == "IGH"]
  igk <- dt[chain == "IGK"]
  igl <- dt[chain == "IGL"]
  rc_5prime <- function(x, hept, nona) {
    x[, (hept) := complement_reverse(get(hept))]
    x[, spacer_5 := complement_reverse(spacer_5)]
    x[, (nona) := complement_reverse(get(nona))]
    x
  }
  expressed <- function(x, seg) x[gene_type == seg & sample_count_genomic > 0 & sample_count_AIRRseq > 0]

  ighv <- extract_rss_table(igh, c("v_heptamer", "spacer_3", "v_nonamer"), spacer_size = c(23, 22, 21), canonical_size = 23)
  if (nrow(ighv)) ighv[, `:=`(rss_side = "3prime", gene_type_plot = "IGHV", spacer_length_class = "23")]

  ighj <- expressed(igh, "IGHJ")
  if (nrow(ighj)) ighj <- rc_5prime(ighj, "j_heptamer", "j_nonamer")
  ighj <- extract_rss_table(ighj, c("j_heptamer", "spacer_5", "j_nonamer"), spacer_size = c(23, 22), canonical_size = 23)
  if (nrow(ighj)) ighj[, `:=`(rss_side = "5prime", gene_type_plot = "IGHJ", spacer_length_class = "23")]

  ighd5 <- expressed(igh, "IGHD")
  if (nrow(ighd5)) ighd5 <- rc_5prime(ighd5, "d_5_heptamer", "d_5_nonamer")
  ighd5 <- extract_rss_table(ighd5, c("d_5_heptamer", "spacer_5", "d_5_nonamer"), spacer_size = c(12, 11, 13), canonical_size = 12)
  if (nrow(ighd5)) ighd5[, `:=`(rss_side = "5prime", gene_type_plot = "IGHD_5", spacer_length_class = "12", gene_type = "IGHD_5")]

  ighd3 <- extract_rss_table(igh, c("d_3_heptamer", "spacer_3", "d_3_nonamer"), spacer_size = c(12, 11, 13), canonical_size = 12)
  if (nrow(ighd3)) ighd3[, `:=`(rss_side = "3prime", gene_type_plot = "IGHD_3", spacer_length_class = "12", gene_type = "IGHD_3")]

  igkv <- extract_rss_table(igk, c("v_heptamer", "spacer_3", "v_nonamer"), spacer_size = c(12, 11), canonical_size = 12)
  if (nrow(igkv)) igkv[, `:=`(rss_side = "3prime", gene_type_plot = "IGKV", spacer_length_class = "12")]

  igkj <- expressed(igk, "IGKJ")
  if (nrow(igkj)) igkj <- rc_5prime(igkj, "j_heptamer", "j_nonamer")
  igkj <- extract_rss_table(igkj, c("j_heptamer", "spacer_5", "j_nonamer"), spacer_size = c(23, 22, 21), canonical_size = 23)
  if (nrow(igkj)) igkj[, `:=`(rss_side = "5prime", gene_type_plot = "IGKJ", spacer_length_class = "23")]

  iglj <- expressed(igl, "IGLJ")
  if (nrow(iglj)) iglj <- rc_5prime(iglj, "j_heptamer", "j_nonamer")
  iglj <- extract_rss_table(iglj, c("j_heptamer", "spacer_5", "j_nonamer"), spacer_size = c(12, 11), canonical_size = 12)
  if (nrow(iglj)) iglj[, `:=`(rss_side = "5prime", gene_type_plot = "IGLJ", spacer_length_class = "12")]

  iglv <- extract_rss_table(igl, c("v_heptamer", "spacer_3", "v_nonamer"), spacer_size = c(23, 22, 21), canonical_size = 23)
  if (nrow(iglv)) iglv[, `:=`(rss_side = "3prime", gene_type_plot = "IGLV", spacer_length_class = "23")]

  out <- data.table::rbindlist(list(ighv, ighj, ighd5, ighd3, igkv, igkj, iglv, iglj), fill = TRUE)
  if (nrow(out) == 0L) return(out)
  out[, spacer_length := nchar(spacer)]
  out[]
}
