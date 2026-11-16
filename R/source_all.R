# =============================================================================
# Mercader Lab — Source All Reusable Functions
# source_all.R
#
# Sources all shared function files.
# Call at the top of every analysis script / notebook cell.
#
# Usage: source("reusable_functions/source_all.R")
#
# Parallel loops (run_ORA_loop, run_GSEA_loop, run_*_sc, run_de_per_celltype,
# annotate_celltypes_enrichr, make_motif_logos, ...) pick their worker count from
# options(mercader.workers = N) / env MERCADER_WORKERS; pass workers = 1 for serial.
# =============================================================================

# Resolve the library directory from this file's own location (works from any
# working directory and inside SnowParam workers); falls back to the repo-relative path.
.LAB_SRC_DIR <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (is.character(f) && grepl("source_all", basename(f)))
      return(normalizePath(dirname(f), mustWork = FALSE))
  }
  normalizePath("reusable_functions", mustWork = FALSE)
})

suppressPackageStartupMessages({
  source(file.path(.LAB_SRC_DIR, "parallel_helpers.R"))
  source(file.path(.LAB_SRC_DIR, "viz_functions.R"))
  source(file.path(.LAB_SRC_DIR, "enrichment_functions.R"))
  source(file.path(.LAB_SRC_DIR, "utils_functions.R"))
  source(file.path(.LAB_SRC_DIR, "orthology_functions.R"))
  source(file.path(.LAB_SRC_DIR, "meme_functions.R"))
  source(file.path(.LAB_SRC_DIR, "sc_functions.R"))
})

cat("✓ All Mercader Lab reusable functions loaded.\n")
cat("  Bulk:        PCA(), draw_volcano(), make_gct(), process_venn_sets()\n")
cat("  Enrichment:  ORA_function(), GSEA_function(), run_ORA_loop(), run_GSEA_loop()\n")
cat("  Orthology:   convert_zf_to_mouse_alliance(), convert_zf_to_human_alliance(), convert_to_zebrafish_alliance()\n")
cat("  Motifs:      export_peaks_bed(), regions_to_fasta(), random_regions_fasta(), run_meme_suite(),\n")
cat("               run_meme_suite_batch(), meme_wait(), meme_run_status(), parse_xstreme(),\n")
cat("               map_tfs_to_zebrafish(), make_motif_logos()\n")
cat("  Single-cell: run_liana_zebrafish(), run_liana_loop(), plot_liana_dotplot(),\n")
cat("               run_de_per_celltype(), pseudobulk_de(), run_ORA_loop_sc(), run_GSEA_loop_sc()\n")
cat(" ", .lab_banner_parallel(), "\n")
