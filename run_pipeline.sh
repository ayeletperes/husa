#!/usr/bin/env bash
# Run the HUSA pipeline, GGS -> figures. Each stage is skipped when its output exists;
# delete the output (or pass --force) to redo it.
#   bash run_pipeline.sh                 # every stage in order
#   bash run_pipeline.sh qtl             # one stage
#   bash run_pipeline.sh figures supp6   # one figure script
set -euo pipefail
cd "$(dirname "$0")"
RSCRIPT="${RSCRIPT:-Rscript}"

force=""; args=()
for a in "$@"; do [[ "$a" == "--force" ]] && force=1 || args+=("$a"); done
want="${args[0]:-}"; only="${args[1]:-}"

declare -A SCRIPT=(
  [husa_table]=R/01_husa_table.R
  [repertoire]=R/02_repertoire.R
  [rss_leader]=R/03_rss_leader.R
  [qtl]=R/04_qtl.R
)
declare -A OUTPUT=(
  [husa_table]=results/husa/husa.tsv
  [repertoire]=results/repertoire/gg_repertoire_data_IGL_genotype_corrected.csv.gz
  [rss_leader]=results/rss_leader/rss_leader_iuis_data.csv.gz
  [qtl]=results/qtl/igkv1d13/panelC_cdr3_vs_usage.tsv
)

run() {  # run <label> <script> <output>
  if [[ -z "$force" && -e "$3" ]]; then echo "[skip] $1 ($3 exists)"; return; fi
  echo "[run ] $1"; "$RSCRIPT" "$2"; echo "[done] $1"
}

# Figure scripts, in dependency order (figure_cdr3_sharing reads supp6's tables).
FIGURES=(figure1 figure2 figure3 figure4 supp1 supp2 supp3 supp4_5 supp6 figure_cdr3_sharing
         figure_dj_pairing figure_guqtl_summary supp7 supp9 supp_light_pairing)
figure_output() {
  case "$1" in
    figure1) echo results/figures/figure1.png ;;
    supp4_5) echo results/figures/supp5.pdf ;;
    *) echo "results/figures/$1.pdf" ;;
  esac
}
run_figures() {
  for f in "${FIGURES[@]}"; do
    [[ -n "$only" && "$only" != "$f" ]] && continue
    run "figures/$f" "R/figures/$f.R" "$(figure_output "$f")"
  done
}

if [[ -z "$want" ]]; then
  for s in husa_table repertoire rss_leader qtl; do run "$s" "${SCRIPT[$s]}" "${OUTPUT[$s]}"; done
  run_figures
  echo "pipeline complete -> results/"
elif [[ "$want" == "figures" ]]; then
  run_figures
elif [[ -n "${SCRIPT[$want]:-}" ]]; then
  run "$want" "${SCRIPT[$want]}" "${OUTPUT[$want]}"
else
  echo "unknown stage: $want (husa_table repertoire rss_leader qtl figures)" >&2; exit 1
fi
