# HUSA pipeline

R scripts that reproduce the analysis behind the HUSA manuscript from the genomic germline
sets (GGS) to the figures: the HUSA allele table, the genotype-corrected repertoires, the
RSS and leader table, the gene-usage and D-J pairing QTL scans, and every manuscript figure.

The upstream steps (raw-read annotation, IgLabel review, PIgLET genotype training) are not
part of this pipeline. Their outputs are the inputs here, hosted on Zenodo (see **Data**).

## Requirements

R (>= 4.2) with `data.table`, `alakazam`, `tigger`, `Biostrings`, `stringdist`, `stringi`,
`igraph`, `ape`, `dendextend`, `seriation`, `ggplot2`, `ggpubr`, `patchwork`, `cowplot`,
`ggrepel`, `ggseqlogo`, `ggh4x`, `ComplexUpset`, `ggVennDiagram`, `ggvenn`, `ComplexHeatmap`,
`circlize`, `gridtext`, `GetoptLong`, `digest`, `mallinfo`, `jsonlite`, and
[`piglet`](https://github.com/ayeletperes/piglet) 1.5.0.999, which provides the allele
similarity clustering, the ASC-to-IUIS vocabulary and the QTL scans. Install the exact
version used here, pinned by commit:

```r
remotes::install_github("ayeletperes/piglet", ref = "8f0b81533208fa51fad9b1cff0f0684e62d71c71")
```

The repertoire and QTL stages hold the full in-house repertoires in memory (about 8 GB).

## Data

The inputs (about 436 MB) are archived on Zenodo, DOI `TODO` (assigned at publication).
Point `fetch_data.sh` at the archive and it unpacks it into `data/`:

```bash
ZENODO_URL="<archive URL>" bash fetch_data.sh   # download and unpack into data/
bash fetch_data.sh --link-local                 # or symlink from a sibling husa_manuscript checkout
```

`data/README.md` lists every input and the manuscript run it comes from.

## Run

```bash
bash run_pipeline.sh                    # all stages, in order
bash run_pipeline.sh qtl                # one stage
bash run_pipeline.sh figures supp6      # one figure
bash run_pipeline.sh --force qtl        # redo a stage whose output exists
```

A stage is skipped when its output exists. Every path is relative to the repository root, so
run from there (the runner does).

## Stages

| stage | script | reads | writes |
|---|---|---|---|
| husa_table | `R/01_husa_table.R` | GGS (3 cohorts), genotypes, baseline ASC, IgLabel labels | `results/husa/husa.tsv`, `husa_rss_filter.tsv` |
| repertoire | `R/02_repertoire.R` | husa.tsv, genotypes, repertoires | `results/repertoire/gg_repertoire_data_<CHAIN>_genotype_corrected.csv.gz`, collapsed genotypes |
| rss_leader | `R/03_rss_leader.R` | in-house GGS, husa.tsv | `results/rss_leader/rss_leader_iuis_data.csv.gz` |
| qtl | `R/04_qtl.R` | repertoires, SNP dosage, gene and feature BEDs | `results/qtl/{source_data,reports,igkv1d13}/` |
| figures | `R/figures/<figure>.R` | the above plus the IMGT and metadata references | `results/figures/<figure>.pdf`, `results/figures/source_data/` |

Each figure script computes the tables it draws from, writes them to
`results/figures/source_data/`, and then draws. `supp4_5.R` draws both supplementary
figures 4 and 5 from one computation; `figure_cdr3_sharing.R` reads tables written by
`supp6.R`, so the runner orders them.

Shared code lives in `R/lib/`: the RSS alignment helpers, the heatmap ordering and
annotation helpers behind figure 4 and supplementary figures 4 and 5, the conditional
pairing tables and panels behind the three pairing figures, and the upset variant with a
corner panel used by figure 2.

## Layout

```
run_pipeline.sh    runner
fetch_data.sh      inputs from Zenodo (or a local husa_manuscript checkout)
R/00_setup.R       input and output paths
R/01_...04_...     one script per stage
R/figures/         one script per figure
R/lib/             helpers used by more than one script
data/              inputs (git-ignored)
results/           outputs (git-ignored)
```
