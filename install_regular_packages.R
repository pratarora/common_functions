list_of_packages <- c("ggplot2", "dplyr", "gridGraphics", "gridExtra", 
                    "ggrepel", "DBI", "svglite", "htmltools", "org.Mm.eg.db", 
                    "igraph", "psych", "org.Dr.eg.db", "clusterProfiler", 
                    "edgeR", "UpSetR", "plyr", "DESeq2", "readxl", "plotly", 
                    "data.table", "installr", "knitr", "gtools", "vsn", 
                    "Rsamtools", "GenomicAlignments",
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
                    'ggdendro', 'shinyhelper', "Matrix", "glue", "hdf5r", 
                    "reticulate", "R.utils", "readr", "future", "shiny", 
                    "magrittr", "msigdb", "GSEABase", "limma", "EBImage", "paletteer", "shiny.semantic")


new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
new_packages <- unique(new_packages)
new_packages
if(length(new_packages)) install.packages(new_packages, Ncpus=8)




new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
new_packages <- unique(new_packages)
new_packages
BiocManager::install(new_packages)




devtools::install_github("neurorestore/Libra")



# for shinycell2

# run these in bash
# Install dependencies for trackplot (on command line)
# git clone https://github.com/CRG-Barcelona/libbeato.git
# cd ./libbeato
# git checkout 0c30432
# ./configure
# make
# make install
    
# cd ..
# git clone https://github.com/CRG-Barcelona/bwtool.git
# cd ./bwtool
# ./configure
# make
# make check
# make install
#  OR directly through R
# system("git clone https://github.com/CRG-Barcelona/libbeato.git")
# setwd("./libbeato")
# system("git checkout 0c30432")
# system("./configure")
# system("make")
# system("make install")
# setwd("..")
# system("git clone https://github.com/CRG-Barcelona/bwtool.git")
# setwd("./bwtool")
# system("./configure")
# system("make")
# system("make check")
# system("make install")
# setwd("..")



devtools::install_github("the-ouyang-lab/ShinyCell2")


