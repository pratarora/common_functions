if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")


# List of packages to check and install
packages <- c("tidyverse","ggplot2", "dplyr", "gridGraphics", "gridExtra", "ggrepel", "DBI", "svglite", "htmltools", 
              "org.Mm.eg.db", "igraph", "psych", "org.Dr.eg.db", "clusterProfiler", "ggplot2", "gridExtra",
              "edgeR", "UpSetR", "plyr", "DESeq2", "readxl", "plotly", "data.table", "installr", "knitr",
              "gtools", "vsn", "SummarizedExperiment", "gridExtra", "sva", "BiocParallel", "biomaRt", "curl",
              "stats", "rvcheck", "tidyverse", "grid", "VennDiagram", "RColorBrewer", "ComplexHeatmap",
              "svglite", "colorspace", "circlize", "tidyr", "rlang", "sodium", "devtools", "shinymanager",
              "DT", "fontawesome", "shinydashboard", "shinycssloaders", "igraph", "rmarkdown", "markdown",
              "visNetwork", "heatmaply", "plotly", "brew", "RColorBrewer", "curl", "wesanderson", "shinyauthr",
              "Seurat",  "Seurat", "patchwork", "sva", "ggbiplot", "edgeR", "DiffBind",
              "GenomicFeatures", "ChIPseeker", "rstatix", "ggpubr", "viridis", "ReactomePA", "clusterProfiler",
              "org.Mm.eg.db", "Cardinal", "RforMassSpectrometry/RforMassSpectrometry", "MetaboAnnotation", "IRkernel", "languageserver")

# Loop through the packages and install them if they are not installed
new_pkgs <- packages[!packages %in% installed.packages()]
new_pkgs <- unique(new_pkgs)
new_pkgs
install.packages(new_pkgs)

new_pkgs <- packages[!packages %in% installed.packages()]
new_pkgs <- unique(new_pkgs)
new_pkgs
BiocManager::install(new_pkgs)

install.packages("devtools")
devtools::install_github("neurorestore/Libra")


library(tidyverse) # for data wrangling
library(Cardinal) # for peak picking and clustering of MSI data
library(MetaboAnnotation) # for mass2mz and annotation
library(MetaboCoreUtils) # for adducts
library(Spectra) # for reading mzML files
library(CompoundDb) # for reading compound databases
library(repr) 