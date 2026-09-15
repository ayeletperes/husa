# HUSA: Human Unified Set of Alleles

Analysis code and figures for the HUSA manuscript: the human immunoglobulin allele
table built from the genomic germline sets, the genotype-corrected repertoires, the
RSS and leader analysis, the gene-usage and D-J pairing QTL scans, and the figures.

## Repository structure

- `scripts/` - the R scripts (`01`-`04` are the analysis stages; `figures/` are the
  figure scripts; `lib/` are shared helpers).
- `figures/` - the manuscript figures.
- `figure_data/` - the source table behind each figure.
- `data/` - inputs and derived tables (not in the repository; see **Data**).

## Setup

Install the R packages once (CRAN, Bioconductor, and piglet pinned to `v1.5.0`):

```bash
Rscript install.R
```

## Reproduce the figures

The source tables are in `figure_data/`, so the figures redraw without the input
data. Run a figure script from the repository root:

```bash
Rscript scripts/figures/figure2.R
```

Each figure script reads its table from `figure_data/` and writes the figure to
`figures/`. `supp4_5.R` draws supplementary figures 4 and 5; run `supp6.R` before
`figure_cdr3_sharing.R`, which reuses its tables.

## Run the full analysis

This needs the input data (see **Data**) unpacked into `data/`. Run the stages in
order, then the figure scripts:

```bash
Rscript scripts/01_husa_table.R
Rscript scripts/02_repertoire.R
Rscript scripts/03_rss_leader.R
Rscript scripts/04_qtl.R
```

## Data

The inputs (`husa_data.tar.gz`, about 293 MB) are archived on Zenodo, DOI
[10.5281/zenodo.22759738](https://doi.org/10.5281/zenodo.22759738) (released at
publication). Download the archive and unpack it into `data/`. The AIRR-seq and
genomic sequencing datasets are in the SRA under accession `PRJNA1274485`.

## License

CC BY 4.0 (see `LICENSE`). Please cite the manuscript.
