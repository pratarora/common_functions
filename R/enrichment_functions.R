# =============================================================================
# Mercader Lab — Pathway Enrichment Functions
# enrichment_functions.R
#
# Functions:
#   - ORA_function()     : Single ORA (GO / Reactome) enrichment call
#   - GSEA_function()    : Single GSEA call (GO, Reactome, or MSigDB)
#   - run_ORA_loop()     : Loop ORA_function() over all contrasts × thresholds
#   - run_GSEA_loop()    : Loop GSEA_function() over all contrasts × thresholds
#
# Output folder structure:
#   pathway_analysis/
#     ORA/
#       lfc2/
#         TreatmentVsControl/
#           GO/
#             BP/  CC/  MF/
#           Reactome/
#     GSEA/
#       lfc2/
#         TreatmentVsControl/
#           GO/
#             BP/  CC/  MF/
#           Reactome/
#
# Usage: source("reusable_functions/enrichment_functions.R")
# =============================================================================

# -----------------------------------------------------------------------------
# Shared helpers (check_required_packages, lab_workers, lab_bpparam, lab_run_tasks,
# .push_seed) live in parallel_helpers.R — source it if this file is sourced alone.
# -----------------------------------------------------------------------------
if (!exists("lab_run_tasks", mode = "function")) local({
  d <- "reusable_functions"
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (is.character(f) && grepl("enrichment_functions", basename(f))) { d <- dirname(f); break }
  }
  source(file.path(d, "parallel_helpers.R"))
})

# Packages needed for one ORA / GSEA call, given the ontology and organism.
.enrichment_pkgs <- function(onto, org, gsea = FALSE) {
  c("ggplot2", "dplyr", "clusterProfiler",
    if (org == "human") "org.Hs.eg.db" else "org.Mm.eg.db",
    if (onto == "reactome") c("ReactomePA", "reactome.db"),
    if (gsea) "enrichplot",
    if (gsea && onto == "msigdb") "msigdbr")
}

# MSigDB is slow to load; cache per species for the session (and for forked workers).
.msigdbr_cache <- new.env(parent = emptyenv())
.msigdbr_cached <- function(species) {
  if (is.null(.msigdbr_cache[[species]]))
    .msigdbr_cache[[species]] <- msigdbr::msigdbr(species = species)
  .msigdbr_cache[[species]]
}

# Fail early, with install instructions, if the core stack is not installed.
check_required_packages(c("ggplot2", "dplyr", "clusterProfiler", "enrichplot", "ReactomePA",
                          "reactome.db", "msigdbr", "org.Hs.eg.db", "org.Mm.eg.db",
                          "BiocParallel"), "enrichment_functions.R")

library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(msigdbr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(BiocParallel)


# -----------------------------------------------------------------------------
# .make_enrichment_dir()  [internal helper]
#
# Builds and creates the output directory path for a given analysis type.
#
# Args:
#   analysis_type : "ORA" or "GSEA"
#   lfc_label     : threshold key e.g. "lfc2", "lfc1", "lfc0.5"
#   contrast      : contrast name e.g. "TreatmentVsControl"
#   onto          : "BP", "MF", "CC", "reactome", or "msigdb"
#
# Returns: directory path string (also creates the directory)
# -----------------------------------------------------------------------------
.make_enrichment_dir <- function(analysis_type, lfc_label, contrast, onto,
                                 out_root = "pathway_analysis") {
  onto_folder <- switch(onto,
    "BP"       = "GO/BP",
    "MF"       = "GO/MF",
    "CC"       = "GO/CC",
    "reactome" = "Reactome",
    "msigdb"   = "MSigDB",
    onto
  )
  path <- file.path(out_root, analysis_type, lfc_label, contrast, onto_folder)
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  return(path)
}


# -----------------------------------------------------------------------------
# ORA_function()
#
# Runs GO (BP/MF/CC) or Reactome overrepresentation analysis on a gene list.
# Auto-detects ENSEMBL vs SYMBOL input. Saves CSV + dotplot PDF.
#
# Args:
#   gene_list    : character vector of gene identifiers
#   pval         : adjusted p-value cutoff (default 0.05)
#   onto         : "BP", "MF", "CC", or "reactome"
#   contrast     : contrast name e.g. "TreatmentVsControl" — used in folder path
#   lfc_label    : threshold key e.g. "lfc2" — used in folder path
#   gene_set     : gene set direction e.g. "up", "down", "all" — used in filename
#   org          : "human" (default) or "mouse"
#   num_pathways : max pathways to plot (default 15)
#   universe     : background gene set, same ID type as gene_list (default NULL =
#                  every annotated gene in the OrgDb). Pass the genes that could
#                  have been detected (e.g. all tested genes, converted to the
#                  same species) — required for assay-restricted data such as
#                  ChIP/CUT&RUN peaks, where the whole genome is not a valid background.
#                  gene_list is intersected with universe (clusterProfiler default).
#                  Results with a universe go to <out_root>/ORA_universe/, whole-genome
#                  results to <out_root>/ORA/, so both can sit side by side.
#   out_root     : root output folder (default "pathway_analysis", relative to the
#                  working directory) — e.g. "myproject/downstream/pathway_analysis"
#   seed         : RNG seed set at the start of the call (default 123; the caller's
#                  RNG state is restored afterwards). NULL = do not touch the RNG.
#
# Output files (in <out_root>/ORA[_universe]/<lfc_label>/<contrast>/GO/<onto>/):
#   <contrast>_<lfc_label>_<gene_set>_ORA_<onto>.csv
#   <contrast>_<lfc_label>_<gene_set>_ORA_<onto>.pdf  (+ .png, 300 dpi)
#
# Returns: enrichResult object or NULL
# -----------------------------------------------------------------------------
ORA_function <- function(gene_list, pval = 0.05, onto = "MF",
                         contrast = "contrast", lfc_label = "lfc1",
                         gene_set = "all", org = "human", num_pathways = 15,
                         universe = NULL, out_root = "pathway_analysis", seed = 123,
                         .universe_entrez = NULL) {
  # .universe_entrez: internal — Entrez ids of `universe`, precomputed once per contrast
  #                   by run_ORA_loop() so the conversion is not repeated per task.

  if (length(gene_list) == 0) {
    cat("⚠ Warning: Empty gene list provided to ORA_function\n"); return(NULL)
  }
  check_required_packages(.enrichment_pkgs(onto, org), paste0("ORA_function(onto = '", onto, "')"))
  restore_seed <- .push_seed(seed); on.exit(restore_seed(), add = TRUE)

  orgdb      <- if (org == "human") org.Hs.eg.db else org.Mm.eg.db
  is_ensembl <- grepl("^ENS[A-Z]*G\\d+", gene_list)
  key_type   <- if (sum(is_ensembl) > length(gene_list) * 0.5) "ENSEMBL" else "SYMBOL"
  compGO     <- NULL
  if (!is.null(universe)) universe <- unique(as.character(universe))

  if (onto %in% c("MF", "CC", "BP")) {
    compGO <- tryCatch(
      enrichGO(gene = gene_list, universe = universe, pvalueCutoff = pval, keyType = key_type,
               pAdjustMethod = "BH", OrgDb = orgdb, ont = onto),
      error = function(e) { cat("✗ Error in ORA enrichment:", e$message, "\n"); NULL })

  } else if (onto == "reactome") {
    gene_list_df <- tryCatch(
      bitr(gene_list, fromType = key_type, toType = "ENTREZID", OrgDb = orgdb),
      error = function(e) { cat("✗ Error in gene ID conversion:", e$message, "\n"); NULL })
    universe_entrez <- .universe_entrez
    if (!is.null(universe) && is.null(universe_entrez)) {
      universe_df <- tryCatch(
        bitr(universe, fromType = key_type, toType = "ENTREZID", OrgDb = orgdb),
        error = function(e) { cat("✗ Error in universe ID conversion:", e$message, "\n"); NULL })
      if (!is.null(universe_df)) universe_entrez <- unique(universe_df$ENTREZID)
    }

    if (!is.null(gene_list_df) && nrow(gene_list_df) > 0) {
      compGO <- tryCatch(
        enrichPathway(gene = gene_list_df$ENTREZID, universe = universe_entrez, pvalueCutoff = pval,
                      organism = if (org == "human") "human" else "mouse",
                      readable = TRUE),
        error = function(e) { cat("✗ Error in Reactome ORA:", e$message, "\n"); NULL })
    }
  }

  if (!is.null(compGO) && nrow(as.data.frame(compGO)) > 0) {
    compGO_df <- as.data.frame(compGO)
    gr <- do.call(rbind, strsplit(as.character(compGO_df$GeneRatio), "/", fixed = TRUE))
    compGO_df$GeneRatio_decimal <- as.numeric(gr[, 1]) / as.numeric(gr[, 2])

    out_dir  <- .make_enrichment_dir(if (is.null(universe)) "ORA" else "ORA_universe",
                                     lfc_label, contrast, onto, out_root = out_root)
    prefix   <- paste0(basename(contrast), "_", lfc_label, "_", gene_set)  # basename: sc wrappers pass "contrast/celltype"

    write.csv(compGO_df,
              file.path(out_dir, paste0(prefix, "_ORA_", onto, ".csv")),
              row.names = FALSE)

    pathways_selected <- unique(compGO_df$Description)[1:min(num_pathways, nrow(compGO_df))]
    compGO_df_subset  <- compGO_df[compGO_df$Description %in% pathways_selected, ]

    p <- ggplot(compGO_df_subset,
                aes(x = GeneRatio_decimal,
                    y = reorder(Description, GeneRatio_decimal),
                    color = p.adjust, size = Count)) +
      geom_point(alpha = 0.7) +
      scale_color_gradient(low = "red", high = "blue",
                           guide = guide_colourbar(reverse = TRUE)) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 9)) +
      labs(title    = paste0("ORA: ", onto),
           subtitle = paste(contrast, lfc_label, gene_set, "| background:",
                            if (is.null(universe)) "whole genome" else paste0("universe (", length(universe), " genes)")),
           y = "Pathway", x = "Gene Ratio",
           color = "Adj. P-value", size = "Count")

    ggsave(file.path(out_dir, paste0(prefix, "_ORA_", onto, ".pdf")),
           p, width = 12, height = max(6, num_pathways * 0.4))
    ggsave(file.path(out_dir, paste0(prefix, "_ORA_", onto, ".png")),
           p, width = 12, height = max(6, num_pathways * 0.4), dpi = 300)

    cat("✓ ORA completed:", onto, "-", nrow(compGO_df), "pathways,",
        nrow(compGO_df_subset), "plotted →", out_dir, "\n")
    return(compGO)
  } else {
    cat("✗ No ORA enrichment found for", onto, "\n"); return(NULL)
  }
}


# -----------------------------------------------------------------------------
# GSEA_function()
#
# Runs Gene Set Enrichment Analysis using a ranked gene list (by log2FC).
# Supports GO (BP/MF/CC), Reactome, and MSigDB collections.
# Saves CSV, dotplot PDF, and ranked enrichment PDFs for top pathways.
#
# Args:
#   df                : DESeq2 results data frame; requires columns:
#                       padj, log2FoldChange, gene
#   pval_deg          : padj cutoff for building the ranked list (default 0.05)
#   pval_enrich       : enrichment p-value cutoff (default 0.2)
#   onto              : "BP", "MF", "CC", "reactome", or "msigdb"
#   contrast          : contrast name e.g. "TreatmentVsControl" — used in folder path
#   lfc_label         : threshold key e.g. "lfc2" — used in folder path
#   organism          : "human" (default) or "mouse"
#   msigdb_collection : MSigDB category string e.g. "H" (Hallmarks);
#                       only used when onto = "msigdb"
#   num_pathways      : max pathways in dotplot (default 15)
#   top_n_ranked      : number of individual enrichment plots to save (default 3)
#   seed              : RNG seed set before the permutation test (default 123; the
#                       caller's RNG state is restored afterwards). NULL = do not
#                       touch the RNG (results then vary slightly run to run).
#
# Output files (in GSEA/<lfc_label>/<contrast>/GO/<onto>/):
#   <contrast>_<lfc_label>_GSEA_<onto>.csv
#   <contrast>_<lfc_label>_GSEA_<onto>_dotplot.pdf
#   <contrast>_<lfc_label>_GSEA_<onto>_ranked_<term>.pdf
#
# Returns: gseaResult object or NULL
# -----------------------------------------------------------------------------

# Internal helper: replace ENTREZID with SYMBOL in core_enrichment column
.entrezid_to_symbol <- function(result_df, OrgDb, id_map = NULL) {
  if (!"core_enrichment" %in% names(result_df)) return(result_df)
  ids_list <- strsplit(as.character(result_df$core_enrichment), "/", fixed = TRUE)
  if (is.null(id_map)) {
    all_ids <- unique(unlist(ids_list))
    id_map  <- tryCatch(
      bitr(all_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = OrgDb),
      error = function(e) NULL)
  }
  if (is.null(id_map) || nrow(id_map) == 0) return(result_df)
  result_df$core_enrichment <- vapply(ids_list, function(ids) {
    syms <- id_map$SYMBOL[match(ids, id_map$ENTREZID)]
    syms[is.na(syms)] <- ids[is.na(syms)]
    paste(syms, collapse = "/")
  }, character(1))
  return(result_df)
}

GSEA_function <- function(df, pval_deg = 0.05, pval_enrich = 0.2, onto = "MF",
                          contrast = "contrast", lfc_label = "lfc1",
                          lfc_threshold = 0,
                          organism = "human", msigdb_collection = NULL,
                          num_pathways = 15, top_n_ranked = 3,
                          out_root = "pathway_analysis", seed = 123,
                          .id_map = NULL) {
  # .id_map: internal — list(key_type, entrez = bitr(<key_type> -> ENTREZID),
  #          entrez2sym = bitr(ENTREZID -> SYMBOL), sym = bitr(ENSEMBL -> SYMBOL) or NULL)
  #          for ALL genes of `df`, precomputed once per contrast by run_GSEA_loop().

  check_required_packages(.enrichment_pkgs(onto, organism, gsea = TRUE),
                          paste0("GSEA_function(onto = '", onto, "')"))
  restore_seed <- .push_seed(seed); on.exit(restore_seed(), add = TRUE)

  gene_list_df <- df %>%
    filter(padj <= pval_deg, !is.na(gene)) %>%
    filter(abs(log2FoldChange) >= lfc_threshold) %>%
    arrange(desc(log2FoldChange)) %>%
    filter(!duplicated(gene))
  gene_list_df %>% summarise(total_genes = n(), upregulated = sum(log2FoldChange > 0),
                 downregulated = sum(log2FoldChange < 0)) %>%
    print()
  if (nrow(gene_list_df) == 0) {
    cat("⚠ No significant genes found for GSEA\n"); return(NULL)
  }
  gene_list <- sort(setNames(gene_list_df$log2FoldChange, gene_list_df$gene),
                   decreasing = TRUE)
  cat("  Ranked gene list:", length(gene_list), "genes\n")

  # Auto-detect input ID type (same logic as ORA_function)
  is_ensembl <- grepl("^ENS[A-Z]*G\\d+", names(gene_list))
  key_type   <- if (sum(is_ensembl) > length(gene_list) * 0.5) "ENSEMBL" else "SYMBOL"
  cat("  Detected ID type:", key_type, "\n")

  OrgDb  <- if (organism == "human") org.Hs.eg.db else org.Mm.eg.db
  compGO <- NULL
  # precomputed maps are only used when they were built for the same key type
  use_map <- !is.null(.id_map) && identical(.id_map$key_type, key_type)
  map_entrez <- if (use_map) .id_map$entrez else NULL
  map_e2s    <- if (use_map) .id_map$entrez2sym else NULL
  map_sym    <- if (use_map) .id_map$sym else NULL

  if (onto %in% c("MF", "CC", "BP")) {
    # Convert to ENTREZID — more robust than SYMBOL keyType in gseGO
    gene_entrez_go <- if (!is.null(map_entrez)) map_entrez else tryCatch(
      bitr(names(gene_list), fromType = key_type, toType = "ENTREZID", OrgDb = OrgDb),
      error = function(e) { cat("  ✗ Error in", key_type, "->ENTREZID:", e$message, "\n"); NULL })

    if (is.null(gene_entrez_go) || nrow(gene_entrez_go) == 0) {
      cat("  ✗ No genes mapped to ENTREZID for", onto, "\n"); return(NULL)
    }
    gene_list_ez <- gene_list[names(gene_list) %in% gene_entrez_go[[key_type]]]
    names(gene_list_ez) <- gene_entrez_go$ENTREZID[
      match(names(gene_list_ez), gene_entrez_go[[key_type]])]
    gene_list_ez <- sort(gene_list_ez, decreasing = TRUE)
    gene_list_ez <- gene_list_ez[!duplicated(names(gene_list_ez))]
    cat("  Mapped:", length(gene_list_ez), "genes to ENTREZID\n")

    compGO <- tryCatch(
      gseGO(geneList = gene_list_ez, pvalueCutoff = pval_enrich, keyType = "ENTREZID",
            pAdjustMethod = "BH", OrgDb = OrgDb, ont = onto),
      error = function(e) { cat("  ✗ Error in gseGO:", e$message, "\n"); NULL })

  } else if (onto == "reactome") {
    gene_entrez <- if (!is.null(map_entrez)) map_entrez else tryCatch(
      bitr(names(gene_list), fromType = key_type, toType = "ENTREZID", OrgDb = OrgDb),
      error = function(e) { cat("  ✗ Error in gene ID conversion:", e$message, "\n"); NULL })

    if (!is.null(gene_entrez) && nrow(gene_entrez) > 0) {
      gene_list_entrez <- gene_list[names(gene_list) %in% gene_entrez[[key_type]]]
      names(gene_list_entrez) <- gene_entrez$ENTREZID[
        match(names(gene_list_entrez), gene_entrez[[key_type]])]
      gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
      gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
      compGO <- tryCatch(
        gsePathway(geneList = gene_list_entrez, pvalueCutoff = pval_enrich,
                   organism = organism),
        error = function(e) { cat("  ✗ Error in gsePathway:", e$message, "\n"); NULL })
    }

  } else if (onto == "msigdb" && !is.null(msigdb_collection)) {
    # MSigDB TERM2GENE uses gene symbols — convert ENSEMBL to SYMBOL if needed
    gene_list_sym <- gene_list
    if (key_type == "ENSEMBL") {
      sym_map <- if (!is.null(map_sym)) map_sym else tryCatch(
        bitr(names(gene_list), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = OrgDb),
        error = function(e) { cat("  ✗ Error in ENSEMBL->SYMBOL:", e$message, "\n"); NULL })
      if (is.null(sym_map) || nrow(sym_map) == 0) {
        cat("  ✗ Failed to convert ENSEMBL to SYMBOL for MSigDB\n"); return(NULL)
      }
      gene_list_sym <- gene_list[names(gene_list) %in% sym_map$ENSEMBL]
      names(gene_list_sym) <- sym_map$SYMBOL[match(names(gene_list_sym), sym_map$ENSEMBL)]
      gene_list_sym <- sort(gene_list_sym, decreasing = TRUE)
      gene_list_sym <- gene_list_sym[!duplicated(names(gene_list_sym))]
    }
    msigdbr_df       <- .msigdbr_cached(if (organism == "human") "Homo sapiens" else "Mus musculus")
    cat_col          <- if ("gs_cat" %in% names(msigdbr_df)) "gs_cat" else "gs_collection"
    msig_df_filtered <- msigdbr_df %>% filter(.data[[cat_col]] == msigdb_collection)
    term2gene        <- msig_df_filtered %>% dplyr::select(gs_name, gene_symbol)
    compGO <- tryCatch(
      GSEA(geneList = gene_list_sym, TERM2GENE = term2gene,
           pvalueCutoff = pval_enrich, pAdjustMethod = "BH"),
      error = function(e) { cat("  ✗ Error in GSEA:", e$message, "\n"); NULL })
  }

  if (!is.null(compGO) && nrow(as.data.frame(compGO)) > 0) {
    compGO@result$.sign <- ifelse(compGO@result$NES > 0, "Activated", "Suppressed")

    out_dir <- .make_enrichment_dir("GSEA", lfc_label, contrast, onto, out_root = out_root)
    prefix  <- paste0(basename(contrast), "_", lfc_label)  # basename: sc wrappers pass "contrast/celltype"

    # Convert ENTREZID in core_enrichment back to gene symbols for readability
    if (onto %in% c("BP", "MF", "CC", "reactome")) {
      compGO@result <- .entrezid_to_symbol(compGO@result, OrgDb, id_map = map_e2s)
    }

    write.csv(compGO@result,
              file.path(out_dir, paste0(prefix, "_GSEA_", onto, ".csv")),
              row.names = FALSE)

    compGO_df_subset <- compGO@result %>%
      slice_max(order_by = abs(NES), n = num_pathways, with_ties = FALSE)

    if (nrow(compGO_df_subset) > 0) {
      p <- ggplot(compGO_df_subset,
                  aes(x = NES, y = reorder(Description, NES),
                      color = p.adjust, size = setSize)) +
        geom_point(alpha = 0.7) +
        scale_color_gradient(low = "red", high = "blue") +
        facet_grid(. ~ .sign, scales = "free_y", space = "free_y") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 8)) +
        labs(title    = paste0("GSEA: ", onto),
             subtitle = paste(contrast, lfc_label),
             y = "Pathway", x = "NES",
             color = "Adj. P-value", size = "Set Size")
      ggsave(file.path(out_dir, paste0(prefix, "_GSEA_", onto, "_dotplot.pdf")),
             p, width = 12, height = max(6, num_pathways * 0.35))
      ggsave(file.path(out_dir, paste0(prefix, "_GSEA_", onto, "_dotplot.png")),
             p, width = 12, height = max(6, num_pathways * 0.35), dpi = 300)
    }

    for (term in head(compGO@result$ID, top_n_ranked)) {
      term_name <- gsub("[^A-Za-z0-9_-]", "_", term)
      pdf_file  <- file.path(out_dir, paste0(prefix, "_GSEA_", onto,
                                             "_ranked_", term_name, ".pdf"))
      .plot_ranked <- function() {   # device is closed even if gseaplot2 fails
        pdf(pdf_file, width = 10, height = 6); on.exit(dev.off(), add = TRUE)
        print(enrichplot::gseaplot2(compGO, geneSetID = term, title = term))
      }
      tryCatch(.plot_ranked(),
               error = function(e) cat("  ⚠ Could not create enrichment plot for", term, "\n"))
    }

    cat("  ✓ GSEA completed:", onto, "-", nrow(compGO@result), "pathways,",
        nrow(compGO_df_subset), "plotted →", out_dir, "\n")
    return(compGO)
  } else {
    cat("  ✗ No GSEA enrichment found for", onto, "\n"); return(NULL)
  }
}


# -----------------------------------------------------------------------------
# .warm_enrichment_caches()  [internal]
#
# Populates clusterProfiler's GO table cache (ont = "ALL", so every ontology is
# served from it — get_GO_data() subsets a cached "ALL" table to BP/MF/CC) and
# ReactomePA's Reactome cache in the PARENT before the loops fork, so every
# worker inherits them copy-on-write instead of rebuilding them (tens of seconds
# per ontology switch). Guarded: if the internals change, the loops simply run
# as before.
# -----------------------------------------------------------------------------
.warm_enrichment_caches <- function(org, keytypes = character(0), go = TRUE, reactome = TRUE) {
  t0 <- Sys.time()
  orgdb <- if (org == "human") org.Hs.eg.db else org.Mm.eg.db
  if (go) for (kt in unique(keytypes))
    tryCatch(invisible(clusterProfiler:::get_GO_data(orgdb, "ALL", kt)), error = function(e) NULL)
  if (reactome)
    tryCatch(invisible(ReactomePA:::get_Reactome_DATA(if (org == "human") "human" else "mouse")),
             error = function(e) NULL)
  invisible(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
}

# Rough relative cost of one analysis, used to dispatch the longest tasks first.
.ORA_ONTO_COST  <- c(BP = 4, reactome = 2, MF = 1.5, CC = 1)
.GSEA_ONTO_COST <- c(BP = 6, reactome = 3, MF = 2, CC = 1.5, msigdb = 1)

.detect_key_type <- function(ids) {
  is_ensembl <- grepl("^ENS[A-Z]*G\\d+", ids)
  if (sum(is_ensembl) > length(ids) * 0.5) "ENSEMBL" else "SYMBOL"
}

# -----------------------------------------------------------------------------
# run_ORA_loop()
#
# Runs ORA_function() over all contrasts × thresholds × ontologies × gene sets
# (× backgrounds), in parallel with BiocParallel.
#
# Args:
#   filtered_results : named list [[contrast]][[thresh_name]], each a DEG data
#                      frame with columns gene, log2FoldChange
#   contrasts        : character vector of contrast names
#   thresholds       : named list keyed by LFC label e.g. list(lfc0.5=0.5, lfc1=1, lfc2=2)
#   org              : "human" (default) or "mouse"
#   num_pathways     : passed through to ORA_function() (default 15)
#   universe         : background gene set passed to ORA_function(). NULL (default)
#                      = whole OrgDb only; a character vector (used for every contrast);
#                      or a named list keyed by contrast name (one background each)
#   background       : "both" (default when a universe is supplied) runs every
#                      analysis twice — whole-genome background -> <out_root>/ORA/ and
#                      the universe -> <out_root>/ORA_universe/; "genome" or "universe"
#                      runs just one. Without a universe it is always "genome".
#   out_root         : root output folder (default "pathway_analysis")
#   set_labels       : names used for the all/up/down gene sets in file names and plot
#                      subtitles, e.g. c(all = "all", up = "higher_in_KO", down = "higher_in_WT")
#                      (default c(all = "all", up = "up", down = "down"))
#   workers          : number of parallel workers (default NULL = lab_workers("light"):
#                      min(32, cores - 2), or options(mercader.workers) / MERCADER_WORKERS;
#                      1 = serial, no parallelism). Capped at the number of analyses.
#                      Tasks are dispatched longest-first with one chunk per task.
#   BPPARAM          : optional BiocParallel backend (e.g. MulticoreParam(16),
#                      SnowParam(4)); overrides `workers`
#   seed             : RNG seed used inside every analysis (default 123; NULL = none).
#                      Each analysis is seeded independently, so results do not depend
#                      on the number of workers.
#
# Returns (invisibly): data frame with one row per analysis — contrast, threshold,
#   onto, gene_set, background, n_genes, status ("done" / "no_enrichment" /
#   "skipped_too_few_genes" / "error"), n_terms, message.
# -----------------------------------------------------------------------------
run_ORA_loop <- function(filtered_results, contrasts, thresholds,
                         org = "human", num_pathways = 15, universe = NULL,
                         background = c("both", "genome", "universe"),
                         out_root = "pathway_analysis",
                         set_labels = c(all = "all", up = "up", down = "down"),
                         workers = NULL, BPPARAM = NULL, seed = 123) {
  background <- match.arg(background)
  if (is.null(universe)) background <- "genome"        # nothing to compare against
  # genome  -> <out_root>/ORA/...           (whole OrgDb as background)
  # universe-> <out_root>/ORA_universe/...  (custom background)
  backgrounds <- switch(background, both = c("genome", "universe"), background)
  ontologies <- c("BP", "MF", "CC", "reactome")
  gene_sets  <- c("all", "up", "down")
  total <- length(contrasts) * length(ontologies) * length(gene_sets) * length(thresholds) *
           length(backgrounds)

  check_required_packages(unique(unlist(lapply(ontologies, .enrichment_pkgs, org = org))),
                          "run_ORA_loop()")

  cat("Starting ORA enrichment analysis...\n")
  cat("Total analyses to run:", length(contrasts), "×", length(ontologies), "×",
      length(gene_sets), "×", length(thresholds), "×", length(backgrounds),
      "background(s) [", paste(backgrounds, collapse = ", "), "] =", total, "\n\n")

  # ---- 1. build the task list (cheap, serial) --------------------------------
  tasks <- list(); skipped <- list()
  for (contrast in contrasts) {
    for (thresh_name in names(thresholds)) {
      deg_data <- filtered_results[[contrast]][[thresh_name]]

      if (is.null(deg_data) || nrow(deg_data) == 0) {
        cat("⚠ No DEGs found for", contrast, thresh_name, "- skipping\n"); next
      }

      gene_sets_list <- list(
        all  = deg_data$gene,
        up   = deg_data$gene[deg_data$log2FoldChange >  0],
        down = deg_data$gene[deg_data$log2FoldChange <  0]
      )

      for (onto in ontologies) {
        for (set_name in gene_sets) {
          genes <- gene_sets_list[[set_name]]
          for (bg in backgrounds) {
            id <- list(contrast = contrast, threshold = thresh_name, onto = onto,
                       gene_set = set_labels[[set_name]], background = bg,
                       n_genes = length(genes))
            if (length(genes) < 3) {
              skipped[[length(skipped) + 1]] <- c(id, list(status = "skipped_too_few_genes",
                                                           n_terms = NA_integer_, message = ""))
              next
            }
            tasks[[length(tasks) + 1]] <- c(id, list(
              genes = genes,
              universe = if (bg == "genome") NULL else
                         if (is.list(universe)) universe[[contrast]] else universe,
              universe_key = if (bg == "genome") NA_character_ else
                             if (is.list(universe)) contrast else ".shared"))
          }
        }
      }
    }
  }

  # ---- 2. parent-side preparation (inherited by forked workers) ----------------
  bp <- lab_bpparam(workers, length(tasks), BPPARAM, profile = "light")
  if (length(tasks) > 0) {
    orgdb     <- if (org == "human") org.Hs.eg.db else org.Mm.eg.db
    key_types <- unique(vapply(tasks, function(t) .detect_key_type(t$genes), character(1)))
    # Entrez ids of each distinct universe, converted once (Reactome tasks need them)
    uni_entrez <- list()
    for (t in tasks) if (!is.na(t$universe_key) && is.null(uni_entrez[[t$universe_key]])) {
      u  <- unique(as.character(t$universe))
      kt <- .detect_key_type(t$genes)
      df <- tryCatch(suppressMessages(bitr(u, fromType = kt, toType = "ENTREZID", OrgDb = orgdb)),
                     error = function(e) NULL)
      uni_entrez[[t$universe_key]] <- if (is.null(df)) list(NULL) else list(unique(df$ENTREZID))
    }
    warm_secs <- .warm_enrichment_caches(org, keytypes = key_types, go = TRUE,
                                         reactome = "reactome" %in% ontologies)
    cat("Caches warmed (GO/Reactome/ID maps) in", warm_secs, "s\n")
  }
  cat("Running", length(tasks), "analyses on", BiocParallel::bpnworkers(bp),
      "worker(s) [", class(bp)[1], "], seed =", if (is.null(seed)) "none" else seed, "\n")

  run_one <- function(task) {
    ue <- if (is.na(task$universe_key) || task$onto != "reactome") NULL else
          uni_entrez[[task$universe_key]][[1]]
    res <- ORA_function(gene_list = task$genes, pval = 0.05, onto = task$onto,
                        contrast = task$contrast, lfc_label = task$threshold,
                        gene_set = task$gene_set, org = org, num_pathways = num_pathways,
                        universe = task$universe, out_root = out_root, seed = seed,
                        .universe_entrez = ue)
    list(status = if (is.null(res)) "no_enrichment" else "done",
         n_terms = if (is.null(res)) 0L else nrow(as.data.frame(res)), message = "")
  }
  cost <- vapply(tasks, function(t) .ORA_ONTO_COST[[t$onto]] * log1p(t$n_genes), numeric(1))
  out <- if (length(tasks) > 0) lab_run_tasks(tasks, run_one, bp, cost = cost) else list()

  # ---- 3. report (main process, task order) ------------------------------------
  rows <- lapply(seq_along(tasks), function(i) {
    t <- tasks[[i]]; o <- out[[i]]
    r <- if (is.null(o$result)) list(status = "error", n_terms = NA_integer_,
                                     message = "worker returned nothing") else o$result
    cat("  →", t$contrast, t$threshold, t$onto, "-", t$gene_set, ":", t$n_genes,
        "genes, background =", t$background, "... ")
    if (length(o$log) > 0) cat(paste(o$log, collapse = "\n"), "\n") else cat("\n")
    if (r$status == "error") cat("  ✗ Error:", r$message, "\n")
    data.frame(t[c("contrast", "threshold", "onto", "gene_set", "background", "n_genes")],
               status = r$status, n_terms = r$n_terms, message = r$message,
               stringsAsFactors = FALSE)
  })
  for (s in skipped)
    cat("  ⊘", s$contrast, s$threshold, s$onto, "-", s$gene_set,
        ": Too few genes (", s$n_genes, ")\n")
  summary_df <- dplyr::bind_rows(c(rows, lapply(skipped, as.data.frame, stringsAsFactors = FALSE)))
  if (nrow(summary_df) > 0) {
    cat("\nSummary:", paste(names(table(summary_df$status)), table(summary_df$status),
                            sep = " = ", collapse = ", "), "\n")
  }
  cat("\n✓✓✓ All ORA enrichment analyses completed ✓✓✓\n")
  invisible(summary_df)
}


# -----------------------------------------------------------------------------
# run_GSEA_loop()
#
# Runs GSEA_function() over all contrasts × thresholds × analysis types, in
# parallel with BiocParallel.
# Analysis types: BP, MF, CC, Reactome, MSigDB Hallmarks.
#
# Args:
#   results_list  : named list [[contrast]] — each element the full (unfiltered)
#                   DESeq2 results data frame with columns padj, log2FoldChange, gene
#   contrasts     : character vector of contrast names
#   thresholds    : named list keyed by LFC label — same convention as run_ORA_loop()
#   organism      : "human" (default) or "mouse"
#   num_pathways  : passed through to GSEA_function() (default 15)
#   top_n_ranked  : number of individual enrichment plots per run (default 3)
#   workers, BPPARAM, seed : as in run_ORA_loop(). The seed matters more here —
#                   GSEA uses permutations, so without a seed results shift slightly
#                   between runs.
#
# Returns (invisibly): data frame with one row per analysis — contrast, threshold,
#   onto, status, n_terms, message.
# -----------------------------------------------------------------------------
run_GSEA_loop <- function(results_list, contrasts, thresholds,
                          organism = "human", num_pathways = 15, top_n_ranked = 3,
                          pval_deg = 1.0, out_root = "pathway_analysis",
                          workers = NULL, BPPARAM = NULL, seed = 123) {
  gsea_analyses <- list(
    list(onto = "BP",       msigdb = NULL),
    list(onto = "MF",       msigdb = NULL),
    list(onto = "CC",       msigdb = NULL),
    list(onto = "reactome", msigdb = NULL),
    list(onto = "msigdb",   msigdb = "H")
  )

  total <- length(contrasts) * length(gsea_analyses) * length(thresholds)
  check_required_packages(unique(unlist(lapply(vapply(gsea_analyses, `[[`, "", "onto"),
                                               .enrichment_pkgs, org = organism, gsea = TRUE))),
                          "run_GSEA_loop()")
  cat("Starting GSEA analysis...\n")
  cat("Total analyses:", length(contrasts), "×", length(gsea_analyses), "×",
      length(thresholds), "=", total, "\n\n")

  tasks <- list()
  for (contrast in contrasts) {
    df_full <- results_list[[contrast]]
    if (is.null(df_full)) { cat("⚠ No results for", contrast, "- skipping\n"); next }
    for (thresh_name in names(thresholds))
      for (analysis in gsea_analyses)
        tasks[[length(tasks) + 1]] <- list(contrast = contrast, threshold = thresh_name,
                                           lfc = thresholds[[thresh_name]],
                                           onto = analysis$onto, msigdb = analysis$msigdb)
        # the results table itself is looked up from `results_list` by contrast key
        # (inherited copy-on-write by forked workers; serialized once per SnowParam worker)
  }

  bp <- lab_bpparam(workers, length(tasks), BPPARAM, profile = "light")
  id_maps <- list()
  if (length(tasks) > 0) {
    # ---- parent-side preparation (inherited by forked workers) -----------------
    OrgDb <- if (organism == "human") org.Hs.eg.db else org.Mm.eg.db
    # load MSigDB once in the parent so forked workers inherit the cache
    if (any(vapply(tasks, function(t) t$onto == "msigdb", logical(1))))
      invisible(.msigdbr_cached(if (organism == "human") "Homo sapiens" else "Mus musculus"))
    # ID maps for ALL genes of each contrast (every threshold uses a subset of them)
    for (contrast in unique(vapply(tasks, `[[`, "", "contrast"))) {
      genes <- unique(as.character(results_list[[contrast]]$gene))
      genes <- genes[!is.na(genes)]
      kt    <- .detect_key_type(genes)
      ez    <- tryCatch(suppressMessages(bitr(genes, fromType = kt, toType = "ENTREZID", OrgDb = OrgDb)),
                        error = function(e) NULL)
      e2s   <- if (is.null(ez)) NULL else tryCatch(
                 suppressMessages(bitr(unique(ez$ENTREZID), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = OrgDb)),
                 error = function(e) NULL)
      sym   <- if (kt == "ENSEMBL") tryCatch(
                 suppressMessages(bitr(genes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = OrgDb)),
                 error = function(e) NULL) else NULL
      id_maps[[contrast]] <- list(key_type = kt, entrez = ez, entrez2sym = e2s, sym = sym)
    }
    warm_secs <- .warm_enrichment_caches(organism, keytypes = "ENTREZID", go = TRUE,
                                         reactome = any(vapply(tasks, function(t) t$onto == "reactome", logical(1))))
    cat("Caches warmed (GO/Reactome/ID maps) in", warm_secs, "s\n")
  }
  cat("Running", length(tasks), "analyses on", BiocParallel::bpnworkers(bp),
      "worker(s) [", class(bp)[1], "], seed =", if (is.null(seed)) "none" else seed, "\n")

  run_one <- function(task) {
    res <- GSEA_function(df = results_list[[task$contrast]], pval_deg = pval_deg, pval_enrich = 0.2,
                         onto = task$onto, contrast = task$contrast,
                         lfc_label = task$threshold, lfc_threshold = task$lfc,
                         organism = organism, msigdb_collection = task$msigdb,
                         num_pathways = num_pathways, top_n_ranked = top_n_ranked,
                         out_root = out_root, seed = seed,
                         .id_map = id_maps[[task$contrast]])
    list(status = if (is.null(res)) "no_enrichment" else "done",
         n_terms = if (is.null(res)) 0L else nrow(as.data.frame(res)), message = "")
  }
  cost <- vapply(tasks, function(t) .GSEA_ONTO_COST[[t$onto]] / (1 + t$lfc), numeric(1))
  out <- if (length(tasks) > 0) lab_run_tasks(tasks, run_one, bp, cost = cost) else list()

  rows <- lapply(seq_along(tasks), function(i) {
    t <- tasks[[i]]; o <- out[[i]]
    r <- if (is.null(o$result)) list(status = "error", n_terms = NA_integer_,
                                     message = "worker returned nothing") else o$result
    cat("\n  →", t$contrast, t$threshold, t$onto, "... ")
    if (length(o$log) > 0) cat(paste(o$log, collapse = "\n"), "\n") else cat("\n")
    if (r$status == "error") cat("  ✗ Error:", r$message, "\n")
    data.frame(contrast = t$contrast, threshold = t$threshold, onto = t$onto,
               status = r$status, n_terms = r$n_terms, message = r$message,
               stringsAsFactors = FALSE)
  })
  summary_df <- dplyr::bind_rows(rows)
  if (nrow(summary_df) > 0)
    cat("\nSummary:", paste(names(table(summary_df$status)), table(summary_df$status),
                            sep = " = ", collapse = ", "), "\n")
  cat("\n✓✓✓ All GSEA analyses completed ✓✓✓\n")
  invisible(summary_df)
}

cat("✓ enrichment_functions.R loaded: ORA_function(), GSEA_function(), run_ORA_loop(), run_GSEA_loop()\n")
cat("  Loops run in parallel (BiocParallel; workers=, BPPARAM=) and are seeded (seed = 123; seed = NULL to disable)\n")
