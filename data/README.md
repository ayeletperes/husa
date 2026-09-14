# Input data

Populate this directory with `bash ../fetch_data.sh` (Zenodo archive) or
`bash ../fetch_data.sh --link-local` (symlinks into a sibling `husa_manuscript` checkout).
Data is git-ignored; only this README is tracked.

```
data/
  ggs/
    Watson_GGS_selected.csv        in-house cohort GGS, selected columns + extended spacers
    Watson_GGS_unfiltered.csv.gz   in-house GGS before filtering (its `keep` flag marks
                                   subjects whose reference lost an IGHD/IGHJ gene)
    HPRC_GGS_selected.csv          HPRC cohort GGS
    1KPG_GGS_selected.csv          1KGP cohort GGS
  genotype/
    genotype_inference.csv         PIgLET per-subject genotype (in-house AIRR-seq)
    snp_dosage.matrix              in-house SNP dosage matrix (variants x subjects, 0/1/2)
  repertoires/
    ggs_repertoires_valid_samples.tsv.gz   in-house AIRR-seq rearrangements, with valid_both
  reference/
    subject_filter_unfiltered.rds  in-house subjects passing the contamination gate, per chain
    hprc_subject_filter.rds        HPRC subjects passing genomic processing, per chain
    1kpg_subject_filter.rds        1KGP subjects passing genomic processing, per chain
    baseline_asc/<CHAIN>/<CHAIN>_crossref.csv   baseline ASC reference (2026-02-24)
    iglabel_labels/<CHAIN>.labels.csv           frozen IgLabel label set (seq_id, iglabel_raw)
    immune_receptor_genomics_251106/            gene.bed plus the feature BEDs
                                   (heptamer, spacer, nonamer, region, exon_1, intron, l-part2, utr)
    imgt_rss_lv.tsv                IMGT RSS reference (snapshot 2026-03-16)
    imgt/leader/imgt_human_IG?L.fasta           IMGT leader references
    imgt/vdj/imgt_human_IG??.fasta              IMGT V/D/J references
    watson_reference_alleles/Homo_sapiens_IG??.fasta   baseline reference FASTAs
  metadata/
    watson_metadata.csv            in-house subject metadata (vdjbase_name, ancestry_population)
    hprc_metadata.tsv              HPRC sample ancestry (sample, Ancestry)
    1kpg_metadata.tsv              1KGP sample ancestry (sample, Ancestry)
```

## Provenance (the manuscript run each input comes from)

| data/ file | source in husa_manuscript |
|---|---|
| ggs/Watson_GGS_selected.csv, Watson_GGS_unfiltered.csv.gz | results_v3/genomic_data_processing/2026-09-10 |
| ggs/HPRC_GGS_selected.csv | results_v3/hprc_genomic_data_processing/2026-09-10 |
| ggs/1KPG_GGS_selected.csv | results_v3/1kpg_genomic_data_processing/2026-09-10 |
| genotype/genotype_inference.csv | results_v3/piglet_inference/2026-09-09 |
| genotype/snp_dosage.matrix | external/watson_genomic_data/vcf.geno.updated.matrix |
| repertoires/ggs_repertoires_valid_samples.tsv.gz | results_v3/repertoires/2026-09-09 |
| reference/subject_filter_unfiltered.rds | results_v3/genomic_data_processing/2026-09-10/subject_filter_all_loci_unfiltered.rds |
| reference/hprc_subject_filter.rds, 1kpg_subject_filter.rds | results_v3/{hprc,1kpg}_genomic_data_processing/2026-09-10 |
| reference/baseline_asc/ | results/asc_baseline/2026-02-24 |
| reference/iglabel_labels/ | scratch_rerun/gldb_merged/human/<CHAIN>/husa_final_merged_2026-09-10/final/final.labels.csv |
| reference/immune_receptor_genomics_251106/ | external/immune_receptor_genomics/251106 |
| reference/imgt_rss_lv.tsv | results/reference_data/imgt_rss/2026-03-16 |
| reference/imgt/ | immunogenomic_db/var/imgt_cache/imgt_ref/human |
| reference/watson_reference_alleles/ | external/watson_reference_alleles |
| metadata/ | external/watson_metadata, external/hprc_genomic_data, external/1KPG_genomic_data |
