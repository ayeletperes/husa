# Install the R packages the HUSA scripts use.
#   Rscript install.R

cran <- c("data.table", "alakazam", "tigger", "stringdist", "stringi", "igraph",
          "ape", "dendextend", "seriation", "ggplot2", "ggpubr", "patchwork",
          "cowplot", "ggrepel", "ggseqlogo", "ggh4x", "ComplexUpset",
          "ggVennDiagram", "ggvenn", "circlize", "gridtext", "GetoptLong",
          "digest", "mallinfo", "jsonlite", "remotes", "BiocManager")
bioc <- c("Biostrings", "ComplexHeatmap")

have <- rownames(installed.packages())
if (length(setdiff(cran, have))) install.packages(setdiff(cran, have))
if (length(setdiff(bioc, have))) BiocManager::install(setdiff(bioc, have), update = FALSE, ask = FALSE)

# piglet: allele clustering, the ASC-to-IUIS vocabulary and the QTL scans (pinned).
remotes::install_github("ayeletperes/piglet", ref = "v1.5.0")
