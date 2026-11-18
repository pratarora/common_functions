# common_functions
Commonly used R functions.

## Module files (`R/`) — `source("R/source_all.R")` loads all
| File | Description |
| --- | --- |
| viz_functions.R | `PCA()` and `draw_volcano()` (replace the old pca.R and Volcano_ggplot.R), plus shared colour helpers `lab_palette_colors()`, `lab_scale_continuous()`, `lab_scale_discrete()`. Every plotting function takes `palette =` (and `reverse_palette =`): `"red_blue"`, viridis options (`"viridis"`, `"magma"`, ...), colour-blind-safe sets (`"tol_bright"`, `"tol_muted"`, `"okabe_ito"`, `"carto_safe"`, ...), any RColorBrewer name, or a vector of colours (a named vector maps groups to colours). `PCA()` defaults to `"tol_bright"`; `draw_volcano()` has `col_up`, `col_down`, `col_ns` |
| enrichment_functions.R | `ORA_function()`, `GSEA_function()` (GO / Reactome / MSigDB), `run_ORA_loop()`, `run_GSEA_loop()` — run in parallel (BiocParallel) |
| orthology_functions.R | Alliance-based zebrafish <-> mouse/human conversion (replace the Zebrafish_to_*_gene_conversion*.R scripts) |
| meme_functions.R, run_meme_suite.sh | MEME Suite motif analysis |
| sc_functions.R | Single-cell helpers: clustering QC, annotation, pseudobulk DE, LIANA, ORA/GSEA for single-cell (replaces GO_function_sc.R) |
| utils_functions.R | `make_gct()`, `process_venn_sets()` |
| parallel_helpers.R | Worker count / BiocParallel helpers (`workers=`, `options(mercader.workers=)`, env `MERCADER_WORKERS`) |

## Standalone files kept (not covered by the modules)
GO_function.R and GSEA_function.R (custom GMT support), Gene_symbol_ensmbl.R, gene_symbol_to_entrez.R, gene_synonyms.R, granges_to_bed.R, Single_cell_bulk_sum.R, umap.R, Heatmap.R, install_regular_packages.R, Ubuntu_dependencies_for_R, r.snippets.
