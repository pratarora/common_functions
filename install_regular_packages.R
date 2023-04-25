
# List of packages to check and install
packages <- c("tidyverse","ggplot2", "dplyr", "gridGraphics", "gridExtra", "ggrepel", "DBI", "svglite", "htmltools", "biomaRt","BiocManager",
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

new_pkgs <- packages[!packages %in% installed.packages()]
new_pkgs <- unique(new_pkgs)
new_pkgs
install.packages(new_pkgs)

new_pkgs <- packages[!packages %in% installed.packages()]
new_pkgs <- unique(new_pkgs)
new_pkgs
BiocManager::install(new_pkgs)

ainstall.packages("devtools")
devtools::install_github("neurorestore/Libra")


