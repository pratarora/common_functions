list_of_packages <- c("ggplot2", "dplyr", "gridGraphics", "gridExtra", 
                      "ggrepel", "DBI", "svglite", "htmltools", "org.Mm.eg.db", 
                      "igraph", "psych", "org.Dr.eg.db", "clusterProfiler", 
                      "edgeR", "UpSetR", "plyr", "DESeq2", "readxl", "plotly", 
                      "data.table", "installr", "knitr", "gtools", "vsn", 
                      "RSamtools", "GenomicAlignments",
                      "SummarizedExperiment", "sva", "BiocParallel", "biomaRt", 
                      "curl", "stats", "rvcheck", "tidyverse", "grid", "brew",
                      "VennDiagram", "RColorBrewer", "ComplexHeatmap", 
                      "colorspace", "circlize", "tidyr", "matrixStats", "rlang",
                      "sodium", "devtools", "shinymanager", "DT", "fontawesome",
                      "shinydashboard", "shinycssloaders", "rmarkdown", 
                      "markdown", "visNetwork", "heatmaply", "wesanderson", 
                      "shinyauthr", "Seurat", "shinythemes","ggVennDiagram",
                      "caret", "tidymodels", "BiocManager", "languageserver", 
                      "patchwork", "ggbiplot", "DiffBind", "GenomicFeatures", 
                      "ChIPseeker", "rstatix", "ggpubr", "viridis", 
                      "ReactomePA", "Cardinal", "MetaboAnnotation", "IRkernel", 
                      "RforMassSpectrometry/RforMassSpectrometry", "openxlsx", 
                      'ggdendro', 'shinyhelper')



new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
new_packages <- unique(new_packages)
new_packages
if(length(new_packages)) install.packages(new_packages,,Ncpus=8)




new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
new_packages <- unique(new_packages)
new_packages
BiocManager::install(new_packages)

devtools::install_github("neurorestore/Libra")
devtools::install_github("SGDDNB/ShinyCell")
