#!/usr/bin/env Rscript
# Gene-usage QTL (all loci) and D-J pairing QTL (IGH) with piglet, variant annotation by
# gene feature, and the IGKV1D-13 example behind the CDR3-sharing figure.
#
#   data/qtl/source_data/   usage_associations_<LOCUS>.tsv.gz, asc_usage.tsv.gz, dj_enrichment.tsv.gz,
#                              pairing_associations.tsv.gz, pairing_associations_by_d.tsv.gz, dosage_long.tsv.gz
#   data/qtl/reports/       usage_leads, usage_per_asc, pairing_leads, pairing_leads_by_d, pairing_per_partner,
#                              thresholds, variant_features, feature_summary_*, feature_enrichment,
#                              significant_variants, variant_summary  (.tsv)
#   data/qtl/igkv1d13/      panelA_manhattan, panelB_usage_by_genotype, panelC_cdr3_vs_usage  (.tsv)

source("scripts/00_setup.R")
suppressPackageStartupMessages(library(piglet))
for (d in c("source_data", "reports", "igkv1d13")) dir.create(file.path(OUT$qtl, d), showWarnings = FALSE)
w <- function(x, sub, name) fwrite(x, file.path(OUT$qtl, sub, name), sep = "\t")

min_subjects <- 60L; min_maf <- 0.05; alpha <- 0.05; r2_clump <- 0.8; min_genotype_group <- 5L; cis_window <- 5e4
contig_of <- c(IGH = "igh", IGK = "chr2", IGL = "chr22")
segments_of <- list(IGH = c(V = "v_gene_iuis", D = "d_gene_iuis", J = "j_gene_iuis"),
                    IGK = c(V = "v_gene_iuis", J = "j_gene_iuis"), IGL = c(V = "v_gene_iuis", J = "j_gene_iuis"))

# ---- inputs -----------------------------------------------------------------------

# Subjects whose personalised reference lost an IGHD or IGHJ gene: their reads for that
# gene land on the remaining genes and corrupt the partner profile, so they are dropped.
u <- fread(need(IN$ggs_watson_unfiltered), select = c("subject", "gene_type", "gene", "keep"))[gene_type %in% c("IGHD", "IGHJ")]
incomplete <- u[, .(any_keep = any(keep)), by = .(subject, gene_type, gene)][any_keep == FALSE, sort(unique(subject))]

g <- fread(need(IN$snp_dosage))
setnames(g, 1L, "snp"); g[, snp := as.character(snp)]
snp_contig <- sub("_.*$", "", g$snp)
snp_pos <- suppressWarnings(as.integer(sub("^[^_]*_([0-9]+).*$", "\\1", g$snp)))
dosage_for <- function(subjects, contig) {
  keep <- sort(intersect(subjects, names(g)))
  if (length(keep) < min_subjects) stop("only ", length(keep), " subjects overlap the repertoire and the genotypes")
  m <- as.matrix(g[, ..keep]); rownames(m) <- g$snp; storage.mode(m) <- "numeric"
  freq <- rowMeans(m, na.rm = TRUE) / 2
  maf <- pmin(freq, 1 - freq)
  ok <- snp_contig %chin% contig & is.finite(maf) & maf >= min_maf
  list(dosage = m[ok, , drop = FALSE],
       variants = data.table(variant = g$snp[ok], contig = snp_contig[ok], pos = snp_pos[ok], maf = maf[ok],
                             missing_rate = rowMeans(is.na(m))[ok]))
}

# Gene midpoints per locus, with the segment letter so IGHD3-9 and IGHV3-9 stay apart.
bed <- fread(need(IN$gene_bed), col.names = c("chr", "start", "end", "gene"))
positions_for <- function(locus) {
  b <- bed[chr == contig_of[[locus]] & startsWith(gene, locus)]
  b[, `:=`(segment = substr(gene, nchar(locus) + 1L, nchar(locus) + 1L),
           short = sub(paste0("^", locus, "[VDJ]"), "", gene), mid = (start + end) / 2)][]
}
# Every member midpoint of a slash-joined ASC label ("D5-12/5-18").
member_positions <- function(label, positions) {
  members <- sub("^[VDJ]", "", strsplit(label, "/", fixed = TRUE)[[1]])
  positions$mid[positions$segment == substr(label, 1, 1) & positions$short %chin% members]
}

# Productive, unmutated rearrangements of one locus.
read_repertoire <- function(locus) {
  d <- fread(need(file.path(OUT$repertoire, sprintf("gg_repertoire_data_%s_genotype_corrected.csv.gz", locus))),
             select = c("subject", unname(segments_of[[locus]]), "v_mut", "productive"))
  d[, subject := as.character(subject)]
  n_in <- nrow(d)
  d <- d[toupper(as.character(productive)) %chin% c("TRUE", "T")][v_mut <= 0]
  excluded <- intersect(incomplete, unique(d$subject))
  d <- d[!subject %chin% excluded]
  cat(sprintf("%s: %s of %s rearrangements productive and unmutated, %d subjects excluded for an incomplete reference\n",
              locus, format(nrow(d), big.mark = ","), format(n_in, big.mark = ","), length(excluded)))
  list(data = d, excluded = excluded)
}

# Written form of an association table: the exact-LD key and the columns that are
# constant or derivable are dropped, keys first.
share <- function(a) {
  a <- copy(a)
  drop <- intersect(c("ld_group", "contig", "chrom", "missing_rate", "well_powered", "n_response", "conditional"), names(a))
  if (length(drop)) a[, (drop) := NULL]
  setcolorder(a, intersect(c("variant", "pos", "maf", "j_gene", "d_gene", "asc", "segment", "n"), names(a)))
  a
}

# ---- gene usage, per locus ---------------------------------------------------------

usage_res <- list()
dosages <- list()
reps <- list()
for (locus in names(contig_of)) {
  rep <- reps[[locus]] <- read_repertoire(locus)
  pheno <- ascUsagePhenotype(rep$data, segments_of[[locus]])
  geno <- dosage_for(unique(pheno$subject), contig_of[[locus]])
  dosages[[locus]] <- geno
  res <- runGeneUsageQTL(rep$data, geno$dosage, geno$variants, segments = segments_of[[locus]],
                         positions = positions_for(locus), locus = locus, min_subjects = min_subjects,
                         alpha = alpha, r2_clump = r2_clump, min_genotype_group = min_genotype_group, cis_window = cis_window)
  if (nrow(res$leads)) {
    res$per_asc <- merge(res$per_asc, res$leads[, .(n_leads = .N, top_variant = variant[1L], top_beta = beta[1L]), by = asc],
                         by = "asc", all.x = TRUE)
    setnames(res$leads, "contig", "chrom")
  }
  res$thresholds[, n_excluded := length(rep$excluded)]
  setcolorder(res$thresholds, c("locus", "analysis", "n_subjects", "n_excluded"))
  usage_res[[locus]] <- res
  w(share(res$associations), "source_data", sprintf("usage_associations_%s.tsv.gz", locus))
}
asc_usage <- rbindlist(lapply(usage_res, `[[`, "phenotype"), use.names = TRUE)
asc_usage[, locus := rep(names(usage_res), vapply(usage_res, function(r) nrow(r$phenotype), integer(1)))]
setcolorder(asc_usage, c("subject", "asc", "count", "total", "segment", "locus", "n_asc"))
w(asc_usage, "source_data", "asc_usage.tsv.gz")
w(rbindlist(lapply(usage_res, `[[`, "leads"), use.names = TRUE, fill = TRUE), "reports", "usage_leads.tsv")
w(rbindlist(lapply(usage_res, `[[`, "per_asc"), use.names = TRUE, fill = TRUE), "reports", "usage_per_asc.tsv")
usage_thresholds <- rbindlist(lapply(usage_res, `[[`, "thresholds"), use.names = TRUE)
w(usage_thresholds, "reports", "usage_thresholds.tsv")

# ---- D-J pairing, IGH: both anchorings, each its own scan ----------------------------

igh <- reps$IGH
pairs <- pairingTable(igh$data, anchor = "j_gene_iuis", partner = "d_gene_iuis")
geno <- dosage_for(unique(pairs$subject), "igh")
positions <- positions_for("IGH")
d_position <- data.table(partner_gene = sort(unique(pairs$partner_gene)))
d_position[, partner_position := vapply(partner_gene, function(a) mean(member_positions(a, positions)), 1)]
ancestry <- fread(need(IN$watson_metadata), select = c("subject_id", "ancestry_population", "vdjbase_name"))
ancestry <- unique(rbind(ancestry[, .(subject = as.character(vdjbase_name), ancestry = ancestry_population)],
                         ancestry[, .(subject = as.character(subject_id), ancestry = ancestry_population)]))[nzchar(ancestry)]

jd <- runPairingQTL(igh$data, geno$dosage, geno$variants, anchor = "j_gene_iuis", partner = "d_gene_iuis",
                    conditional = "P(J|D)", locus = "IGH", min_subjects = min_subjects, alpha = alpha,
                    r2_clump = r2_clump, min_genotype_group = min_genotype_group)
dj <- runPairingQTL(igh$data, geno$dosage, geno$variants, anchor = "d_gene_iuis", partner = "j_gene_iuis",
                    conditional = "P(D|J)", locus = "IGH", min_subjects = min_subjects, alpha = alpha,
                    r2_clump = r2_clump, min_genotype_group = min_genotype_group)

# The J-anchored leads are adjudicated: rigid shift or reallocation, one subject, and
# whether either margin moved.
lead_detail <- pairingLeadCharacter(jd$leads, jd$pairs, geno$dosage, positions = d_position, ancestry = ancestry,
                                    min_genotype_group = min_genotype_group, min_subjects = min_subjects)
pairing_leads <- lead_detail$leads
setnames(pairing_leads, c("anchor_gene", "contig", "best_partner_gene", "n_partner_up", "n_partner_down", "n_partner_nominal",
                          "n_partner_separated", "anchor_usage_beta", "anchor_usage_p", "partner_usage_p_min",
                          "partner_usage_n_moved", "partner_usage_genes"),
         c("j_gene", "chrom", "best_d_gene", "n_d_up", "n_d_down", "n_d_nominal", "n_d_separated", "j_usage_beta",
           "j_usage_p", "d_usage_p_min", "d_usage_n_moved", "d_usage_genes"))
verdict_names <- c(pairing_plus_anchor_usage = "pairing_plus_j_usage", pairing_plus_partner_usage = "pairing_plus_d_usage")
pairing_leads[verdict %chin% names(verdict_names), verdict := verdict_names[verdict]]
pairing_leads[, conditional := NULL]
per_partner <- setnames(lead_detail$per_partner, c("anchor_gene", "partner_gene", "partner_position"), c("j_gene", "d_gene", "d_position"))
leads_d <- copy(dj$leads)[, conditional := NULL]
setnames(leads_d, c("anchor_gene", "contig"), c("d_gene", "chrom"))
setcolorder(leads_d, c("variant", "n", "n_response", "pillai", "f_stat", "p_value", "min_genotype_group", "well_powered", "d_gene", "chrom"))

enrichment <- copy(jd$pairs)
setnames(enrichment, c("anchor_gene", "partner_gene", "anchor_total", "partner_total", "p_anchor", "p_partner",
                       "p_partner_given_anchor", "p_anchor_given_partner"),
         c("j_gene", "d_gene", "j_total", "d_total", "p_j", "p_d", "p_d_given_j", "p_j_given_d"))
setcolorder(enrichment, c("subject", "d_gene", "j_gene", "count", "depth", "d_total", "j_total", "expected", "enrichment",
                          "p_d", "p_j", "p_j_given_d", "p_d_given_j"))
setorder(enrichment, subject, d_gene, j_gene)

thresholds <- rbind(usage_thresholds,
                    rbind(jd$thresholds, dj$thresholds)[, `:=`(n_subjects = ncol(geno$dosage), n_excluded = length(igh$excluded), grouped_by = c("J", "D"))],
                    fill = TRUE)
setorder(thresholds, analysis, locus)

w(enrichment, "source_data", "dj_enrichment.tsv.gz")
w(share(setnames(copy(jd$associations), "anchor_gene", "j_gene")), "source_data", "pairing_associations.tsv.gz")
w(share(setnames(copy(dj$associations), "anchor_gene", "d_gene")), "source_data", "pairing_associations_by_d.tsv.gz")
w(pairing_leads, "reports", "pairing_leads.tsv")
w(leads_d, "reports", "pairing_leads_by_d.tsv")
w(per_partner, "reports", "pairing_per_partner.tsv")
w(thresholds, "reports", "thresholds.tsv")

# Every filtered IGH variant by subject, so any variant's genotype groups can be looked up.
m <- geno$dosage
dose <- data.table(variant = rep(rownames(m), times = ncol(m)), subject = rep(colnames(m), each = nrow(m)),
                   dosage = as.numeric(m))[is.finite(dosage)]
if (any(dose$dosage != round(dose$dosage))) stop("genotype calls are not whole numbers; write dosage as its own column")
dose[, genotype := as.integer(round(dosage))][, dosage := NULL]
setorder(dose, variant, subject)
w(dose, "source_data", "dosage_long.tsv.gz")

# ---- variant annotation: which gene feature does each tested variant sit in ----------

feature_levels <- c("rss", "coding", "leader", "utr", "intergenic")
bed_of <- function(n) fread(need(file.path(IN$bed_dir, paste0(n, ".bed"))), col.names = c("chr", "start", "end", "gene"))
lp2 <- bed_of("l-part2")
# exon_1 coincides with the region for J segments, so it counts as leader only for genes
# that also have an L-PART2, which is the V segments.
features <- rbindlist(list(
  bed_of("heptamer")[, `:=`(feature = "rss", sub_feature = "heptamer")],
  bed_of("spacer")[, `:=`(feature = "rss", sub_feature = "spacer")],
  bed_of("nonamer")[, `:=`(feature = "rss", sub_feature = "nonamer")],
  bed_of("region")[, `:=`(feature = "coding", sub_feature = "region")],
  bed_of("exon_1")[gene %chin% unique(lp2$gene)][, `:=`(feature = "leader", sub_feature = "l_part1")],
  bed_of("intron")[, `:=`(feature = "leader", sub_feature = "leader_intron")],
  lp2[, `:=`(feature = "leader", sub_feature = "l_part2")],
  bed_of("utr")[, `:=`(feature = "utr", sub_feature = "utr")]
), use.names = TRUE)
features[, priority := match(feature, feature_levels)]

annotate <- function(variants, locus) {
  fl <- features[chr == contig_of[[locus]] & startsWith(gene, locus)]
  gl <- bed[chr == contig_of[[locus]] & startsWith(gene, locus)]
  v <- data.table(variant = variants, pos = suppressWarnings(as.integer(sub("^[^_]*_([0-9]+).*$", "\\1", variants))))[is.finite(pos)]
  setkey(fl, start, end)
  hits <- foverlaps(v[, .(variant, pos, start = pos, end = pos)], fl, by.x = c("start", "end"), type = "within", nomatch = NA)
  hits <- merge(hits, gl[, .(gene, gene_mid = (start + end) / 2)], by = "gene", all.x = TRUE)
  hits[, gene_gap := abs(pos - gene_mid)]
  setorder(hits, variant, priority, gene_gap)
  best <- hits[!is.na(feature), .SD[1L], by = variant, .SDcols = c("pos", "gene", "feature", "sub_feature")]
  best[, distance_to_gene := 0]
  loose <- v[!variant %chin% best$variant, .(variant, pos)]
  if (nrow(loose)) {
    loose[, gene := vapply(pos, function(p) gl$gene[which.min(pmin(abs(gl$start - p), abs(gl$end - p)))], character(1))]
    loose[, distance_to_gene := vapply(pos, function(p) min(pmin(abs(gl$start - p), abs(gl$end - p))), 1)]
    loose[, `:=`(feature = "intergenic", sub_feature = "intergenic")]
  }
  out <- rbind(best, loose, use.names = TRUE, fill = TRUE)
  out[, `:=`(locus = locus, feature = factor(feature, levels = feature_levels))][]
}

usage_leads <- rbindlist(lapply(usage_res, `[[`, "leads"), use.names = TRUE, fill = TRUE)
ann <- list(); by_feature <- list()
for (loc in names(usage_res)) {
  assoc <- usage_res[[loc]]$associations
  ann[[loc]] <- annotate(unique(assoc$variant), loc)
  sig_v <- unique(assoc[significant == TRUE, variant])
  lead_v <- unique(usage_leads[locus == loc, variant])
  by_feature[[loc]] <- ann[[loc]][, .(n_tested = .N, n_significant = sum(variant %chin% sig_v), n_lead = sum(variant %chin% lead_v)),
                                  by = .(locus, gene, feature, sub_feature)]
}
ann <- rbindlist(ann, use.names = TRUE)
long <- rbindlist(by_feature, use.names = TRUE)
by_gene <- long[, .(n_tested = sum(n_tested), n_significant = sum(n_significant), n_lead = sum(n_lead)), by = .(locus, gene, feature)]
wide <- dcast(by_gene, locus + gene ~ feature, value.var = c("n_tested", "n_significant"), fill = 0L, drop = FALSE)
# Enrichment of independent leads per feature class against the rest of the locus. Leads
# rather than significant variants: with locus-wide LD the raw count is a property of the
# haplotype block, not of the variant.
enrich <- rbindlist(lapply(unique(by_gene$locus), function(loc) {
  x <- by_gene[locus == loc, .(tested = sum(n_tested), lead = sum(n_lead)), by = feature]
  tt <- sum(x$tested); tl <- sum(x$lead)
  x[, `:=`(locus = loc, background_rate = tl / tt, lead_rate = lead / tested)]
  x[, odds_ratio := (lead / pmax(tested - lead, 0.5)) / ((tl - lead) / pmax((tt - tested) - (tl - lead), 0.5))]
  x[, p_value := mapply(function(a, b) fisher.test(matrix(c(a, b - a, tl - a, (tt - b) - (tl - a)), 2))$p.value, lead, tested)]
  x[]
}), use.names = TRUE)
w(ann, "reports", "variant_features.tsv")
w(long, "reports", "feature_summary_long.tsv")
w(wide, "reports", "feature_summary_by_gene.tsv")
w(enrich, "reports", "feature_enrichment.tsv")

# ---- significant variants by locus, segment and feature ------------------------------
# Distinct variants within each stratum: a variant significant for several ASCs counts
# once per segment. Variants whose smallest genotype class is below min_genotype_group
# are reported but not counted as kept.
per_variant <- rbindlist(c(
  lapply(names(usage_res), function(locus) {
    sig <- usage_res[[locus]]$associations[significant == TRUE, .(variant, segment, asc, p_value)]
    if (!nrow(sig)) return(NULL)
    mgc <- qtlSmallestGenotypeClass(dosages[[locus]]$dosage)
    x <- sig[, .(analysis = "usage", n_asc = uniqueN(asc), best_p = min(p_value)), by = .(variant, segment)]
    x[, `:=`(locus = locus, min_genotype_group = mgc[variant])]
  }),
  list(jd$associations[significant == TRUE][, .(analysis = "pairing", n_asc = uniqueN(anchor_gene), best_p = min(p_value),
                                                 locus = "IGH", segment = "DJ", min_genotype_group = min_genotype_group[1L]), by = variant])
), use.names = TRUE, fill = TRUE)
per_variant <- merge(per_variant, ann[, .(variant, locus, gene, feature, sub_feature)], by = c("variant", "locus"), all.x = TRUE)
min_group_at <- min_genotype_group  # the column would shadow the setting inside `[`
per_variant[, passes_min_group := is.finite(min_genotype_group) & min_genotype_group >= min_group_at]
setorder(per_variant, locus, analysis, segment, best_p)
variant_summary <- per_variant[, .(n_variants = uniqueN(variant), n_variants_kept = uniqueN(variant[passes_min_group]),
                                   median_min_group = as.numeric(median(min_genotype_group, na.rm = TRUE)), best_p = min(best_p)),
                               by = .(locus, analysis, segment, feature)]
setorder(variant_summary, locus, analysis, segment, feature)
w(per_variant, "reports", "significant_variants.tsv")
w(variant_summary, "reports", "variant_summary.tsv")
print(dcast(variant_summary, locus + analysis + segment ~ feature, value.var = "n_variants_kept", fill = 0L))

# ---- IGKV1D-13: usage against its lead variant, and the public CDR3 it carries -------

asc_name <- "V1-13/1D-13"
assoc <- share(usage_res$IGK$associations)[asc == asc_name]
lead <- assoc[order(p_value, pos)][1]  # tie-break on position so the lead is deterministic
cat(sprintf("IGKV1D-13 lead %s  p=%.3g  beta=%.3f\n", lead$variant, lead$p_value, lead$beta))
manhattan <- assoc[startsWith(variant, "chr2_"), .(variant, pos, p = p_value)]  # chr2_, not chr2: that also takes chr22
manhattan[, `:=`(logp = -log10(p), is_lead = variant == lead$variant)]
w(manhattan, "igkv1d13", "panelA_manhattan.tsv")

husa <- unique(fread(need(file.path(OUT$husa, "husa.tsv")), select = c("allele", "husa")))[!duplicated(allele)]
gene_of <- setNames(sub("\\*.*$", "", sub("_.*$", "", husa$husa)), husa$allele)
rp <- fread(need(file.path(OUT$repertoire, "gg_repertoire_data_IGK_genotype_corrected.csv.gz")),
            select = c("subject", "junction_aa", "v_call_new"))
rp[, subject := as.character(subject)]
rp[, gene := gene_of[v_call_new]]
depth <- rp[, .(total = .N), by = subject]
subj <- sort(intersect(depth$subject, names(g)))
dos <- data.table(subject = subj, dosage = as.numeric(unlist(g[snp == lead$variant, ..subj])))
vrows <- rp[gene == "IGKV1D-13"]
per <- merge(depth[subject %chin% subj], vrows[, .(n_v = .N), by = subject], by = "subject", all.x = TRUE)
per[is.na(n_v), n_v := 0L]
per <- merge(per, vrows[junction_aa == "CQQFNSYPFTF", .(n_cdr3 = .N), by = subject], by = "subject", all.x = TRUE)
per[is.na(n_cdr3), n_cdr3 := 0L]
per <- merge(per, dos, by = "subject")[!is.na(dosage)]
per[, `:=`(usage = n_v / total, cdr3_usage = n_cdr3 / total, low_depth = total < 5000L,
           genotype = factor(round(dosage), levels = 0:2))]
w(per[, .(subject, total, n_v, usage, dosage, genotype, low_depth)], "igkv1d13", "panelB_usage_by_genotype.tsv")
w(per[n_cdr3 > 0 & !low_depth, .(subject, total, n_v, usage, n_cdr3, cdr3_usage, dosage, genotype)], "igkv1d13", "panelC_cdr3_vs_usage.tsv")
cat(sprintf("IGKV1D-13: %d genotyped subjects, %d depth-passing carriers\n", nrow(per), per[n_v > 0 & !low_depth, .N]))
