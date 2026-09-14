#!/usr/bin/env bash
# Populate data/ with the pipeline inputs (layout in data/README.md).
#   bash fetch_data.sh              # download + unpack the Zenodo archive
#   bash fetch_data.sh --link-local # symlink from a sibling husa_manuscript checkout (dev)
set -euo pipefail
cd "$(dirname "$0")"

ZENODO_URL="${ZENODO_URL:-}"                 # set to the Zenodo record archive URL for release
LOCAL_ROOT="${LOCAL_ROOT:-../husa_manuscript}"  # for --link-local
RUN=2026-09-10

mkdir -p data/ggs data/genotype data/repertoires data/reference/iglabel_labels data/reference/imgt data/metadata

if [[ "${1:-}" == "--link-local" ]]; then
  # Relative links so the same tree resolves on every machine that has the sibling checkout.
  up=../..
  ln -sfn "$up/$LOCAL_ROOT/results_v3/genomic_data_processing/$RUN/GGS_all_loci_selected_columns_extended_spacers.csv" data/ggs/Watson_GGS_selected.csv
  ln -sfn "$up/$LOCAL_ROOT/results_v3/genomic_data_processing/$RUN/GGS_all_loci_unfiltered.csv.gz"                    data/ggs/Watson_GGS_unfiltered.csv.gz
  ln -sfn "$up/$LOCAL_ROOT/results_v3/hprc_genomic_data_processing/$RUN/HPRC_all_loci_selected_columns_extended_spacers.csv" data/ggs/HPRC_GGS_selected.csv
  ln -sfn "$up/$LOCAL_ROOT/results_v3/1kpg_genomic_data_processing/$RUN/1KPG_all_loci_selected_columns_extended_spacers.csv" data/ggs/1KPG_GGS_selected.csv
  ln -sfn "$up/$LOCAL_ROOT/results_v3/piglet_inference/2026-09-09/genotype_inference.csv"        data/genotype/genotype_inference.csv
  ln -sfn "$up/$LOCAL_ROOT/external/watson_genomic_data/vcf.geno.updated.matrix"                  data/genotype/snp_dosage.matrix
  ln -sfn "$up/$LOCAL_ROOT/results_v3/repertoires/2026-09-09/ggs_repertoires_dt_valid_samples.tsv.gz" data/repertoires/ggs_repertoires_valid_samples.tsv.gz
  up=../..
  ln -sfn "$up/$LOCAL_ROOT/results_v3/genomic_data_processing/$RUN/subject_filter_all_loci_unfiltered.rds" data/reference/subject_filter_unfiltered.rds
  ln -sfn "$up/$LOCAL_ROOT/results_v3/hprc_genomic_data_processing/$RUN/hprc_subject_filter_all_loci.rds"   data/reference/hprc_subject_filter.rds
  ln -sfn "$up/$LOCAL_ROOT/results_v3/1kpg_genomic_data_processing/$RUN/1kpg_subject_filter_all_loci.rds"   data/reference/1kpg_subject_filter.rds
  ln -sfn "$up/$LOCAL_ROOT/results/asc_baseline/2026-02-24"                        data/reference/baseline_asc
  ln -sfn "$up/$LOCAL_ROOT/external/immune_receptor_genomics/251106"               data/reference/immune_receptor_genomics_251106
  ln -sfn "$up/$LOCAL_ROOT/results/reference_data/imgt_rss/2026-03-16/imgt_rss_lv.tsv" data/reference/imgt_rss_lv.tsv
  ln -sfn "$up/$LOCAL_ROOT/external/watson_reference_alleles"                     data/reference/watson_reference_alleles
  up=../../..
  ln -sfn "$up/$LOCAL_ROOT/immunogenomic_db/var/imgt_cache/imgt_ref/human/leader" data/reference/imgt/leader
  ln -sfn "$up/$LOCAL_ROOT/immunogenomic_db/var/imgt_cache/imgt_ref/human/vdj"    data/reference/imgt/vdj
  for ch in IGH IGK IGL; do
    ln -sfn "$up/$LOCAL_ROOT/scratch_rerun/gldb_merged/human/$ch/husa_final_merged_$RUN/final/final.labels.csv" "data/reference/iglabel_labels/$ch.labels.csv"
  done
  up=../..
  ln -sfn "$up/$LOCAL_ROOT/external/watson_metadata/P28_IGK_metadata_with_vdjbase_names.csv" data/metadata/watson_metadata.csv
  ln -sfn "$up/$LOCAL_ROOT/external/hprc_genomic_data/sample_metadata.tsv"               data/metadata/hprc_metadata.tsv
  ln -sfn "$up/$LOCAL_ROOT/external/1KPG_genomic_data/sample_metadata.tsv"               data/metadata/1kpg_metadata.tsv
  missing=$(find data -xtype l | sort)
  if [[ -n "$missing" ]]; then echo "dangling links:"; echo "$missing"; exit 1; fi
  echo "linked $(find data -type l | wc -l) inputs from $LOCAL_ROOT"
  exit 0
fi

[[ -n "$ZENODO_URL" ]] || { echo "ZENODO_URL not set. Set it, or use --link-local for a local tree." >&2; exit 1; }
curl -L "$ZENODO_URL" -o data/husa_pipeline_inputs.tar.gz
tar -xzf data/husa_pipeline_inputs.tar.gz -C data/
rm -f data/husa_pipeline_inputs.tar.gz
echo "data ready under data/"
