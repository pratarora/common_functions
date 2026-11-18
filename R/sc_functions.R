# =============================================================================
# Mercader Lab — Single-Cell Functions
# sc_functions.R
#
# Functions:
#   ── Clustering quality ──────────────────────────────────────────────────────
#   compute_clustering_qc()      : Silhouette + normalized Shannon entropy (seeded)
#   plot_clustering_qc()         : Two-panel QC plot (silhouette + entropy)
#
#   ── Cell type annotation ────────────────────────────────────────────────────
#   annotate_celltypes_enrichr() : Automated cell-type labeling via Enrichr
#                                  (zf→human, fuzzy consensus, XLSX output; clusters
#                                  are queried in parallel, a few at a time)
#
#   ── DE method selection ─────────────────────────────────────────────────────
#   compute_data_characteristics(): Sparsity & depth per cell type
#   recommend_de_method()        : Nguyen 2023 lookup → recommended DE methods
#
#   ── Nguyen 2023 DE implementations ──────────────────────────────────────────
#   run_pseudobulk_limma()       : Pseudobulk_limma  — pseudobulk + limma-trend
#   run_limmatrend_cov()         : limmatrend / limmatrend_Cov
#   run_limmavoom_cov()          : limmavoom / limmavoom_Cov
#   run_edgeR_cov()              : Pseudobulk_edgeR / ZW_edgeR_Cov
#   run_mast_cov()               : MAST / MAST_Cov
#   run_wilcox_raw()             : Raw_Wilcox
#   run_de_per_celltype()        : Dispatcher — selects & runs recommended method,
#                                  cell types in parallel (workers=)
#
#   ── LIANA ───────────────────────────────────────────────────────────────────
#   run_liana_zebrafish()        : Build zebrafish LIANA ortholog resource
#   run_liana_loop()             : LIANA per condition → CSVs + dotplots
#   plot_liana_dotplot()         : Save publication LIANA dotplot PDF + PNG
#
#   ── SC enrichment wrappers ──────────────────────────────────────────────────
#   run_ORA_loop_sc()            : ORA with <CellType> subfolder (parallel, via run_ORA_loop)
#   run_GSEA_loop_sc()           : GSEA with <CellType> subfolder (parallel, via run_GSEA_loop)
#
#   ── Legacy ──────────────────────────────────────────────────────────────────
#   pseudobulk_de()              : Alias for run_pseudobulk_limma()
#
# Performance notes
#   * Raw counts are read ONCE per top-level call (.sc_counts(), joining Seurat v5
#     split layers first) and sliced per cell type; the Seurat object is never
#     subset (subsetting copies every assay / reduction / graph).
#   * Pseudobulk is one sparse matrix product (counts %*% indicator).
#   * run_de_per_celltype() runs cell types in parallel (lab_run_tasks, "heavy"
#     profile); every DE function's numbers are unchanged.
#
# Usage: source("reusable_functions/sc_functions.R")
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tibble); library(ggplot2); library(patchwork)
})

if (!exists("lab_run_tasks", mode = "function")) local({
  d <- "reusable_functions"
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (is.character(f) && grepl("sc_functions", basename(f))) { d <- dirname(f); break }
  }
  source(file.path(d, "parallel_helpers.R"))
})
# Colour helpers (lab_scale_continuous, ...) live in viz_functions.R
if (!exists("lab_palette_colors", mode = "function")) local({
  d <- "reusable_functions"
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (is.character(f) && grepl("sc_functions", basename(f))) { d <- dirname(f); break }
  }
  source(file.path(d, "viz_functions.R"))
})


# =============================================================================
# 0. INTERNAL HELPERS (counts access, slicing, pseudobulk)
# =============================================================================

# Raw counts of `assay` as ONE dgCMatrix (genes x cells). Seurat v5 objects
# whose counts are split into layers (counts.sample1, counts.sample2, ...) are
# joined first — GetAssayData(layer = "counts") errors on them and LayerData()
# silently returns only the first layer.
.sc_counts <- function(seurat_obj, assay = "RNA") {
  a <- seurat_obj[[assay]]
  if (inherits(a, "Assay5")) {
    lyr <- SeuratObject::Layers(a)
    if (!"counts" %in% lyr || sum(grepl("^counts", lyr)) > 1) a <- SeuratObject::JoinLayers(a)
    m <- SeuratObject::LayerData(a, layer = "counts")
  } else {
    m <- SeuratObject::GetAssayData(a, layer = "counts")
  }
  if (!inherits(m, "dgCMatrix")) m <- methods::as(Matrix::Matrix(m, sparse = TRUE), "CsparseMatrix")
  m
}

# Cells of one cell type (character vector of cell names), from the metadata.
.sc_cells_of <- function(meta, cell_type, cell_type_col) {
  rownames(meta)[!is.na(meta[[cell_type_col]]) & meta[[cell_type_col]] == cell_type]
}

# Pseudobulk by sample: one sparse product instead of one rowSums() per sample.
# Samples with fewer than `min_cells` cells are dropped; sample order = first
# appearance in `sample` (as the previous per-sample loop). Returns a dense
# genes x samples matrix, or NULL if no sample passes.
.sc_pseudobulk <- function(counts, sample, min_cells) {
  f  <- factor(as.character(sample), levels = unique(as.character(sample)))
  nc <- table(f)
  keep <- names(nc)[nc >= min_cells]
  if (length(keep) == 0) return(NULL)
  agg <- counts %*% Matrix::t(Matrix::fac2sparse(f))     # genes x levels
  pb  <- as.matrix(agg[, keep, drop = FALSE])
  dimnames(pb) <- list(rownames(counts), keep)
  pb
}

# Opt-in gene filter for the single-cell (non-pseudobulk) methods: genes detected in
# >= max(3, 1 % of cells), evaluated on the sparse matrix (no dense copy).
.sc_filter_genes <- function(counts) {
  keep <- Matrix::rowSums(counts > 0) >= max(3, ceiling(0.01 * ncol(counts)))
  counts[keep, , drop = FALSE]
}

# Internal switch for MAST: use several cores in the parent (serial path),
# exactly one inside a parallel worker.
.sc_with_mc_cores <- function(expr) {
  in_worker <- isTRUE(getOption("mercader.in_worker"))
  old <- getOption("mc.cores")
  if (!in_worker) options(mc.cores = lab_workers(Inf, "heavy"))
  on.exit(options(mc.cores = old), add = TRUE)
  expr
}


# =============================================================================
# 1. CLUSTERING QUALITY METRICS
# =============================================================================

# -----------------------------------------------------------------------------
# compute_clustering_qc()
#
# Silhouette width (cluster compactness, PCA/Harmony space, subsampled) and
# normalized Shannon entropy (sample mixing per cluster = integration quality).
#
# Args:
#   seurat_obj    : Seurat object with Harmony reduction
#   resolutions   : numeric vector — must already be run via FindClusters()
#   reduction     : reduction embedding used for silhouette (default "harmony")
#   prefix        : cluster column prefix (default "SCT_snn_res.")
#   subsample_n   : max cells for silhouette distance calc (default 5000)
#   sample_col    : metadata column for sample ID used in entropy (default "sample")
#   seed          : seed for the subsample (default 123; NULL = unseeded, as before).
#                   ONE subsample / distance matrix is drawn and reused for every
#                   resolution, so silhouettes are comparable across resolutions
#                   (previously each resolution drew its own, unseeded, sample).
#   method        : "exact" (cluster::silhouette on the subsample, default) or
#                   "approx" (bluster::approxSilhouette on ALL cells, no distance matrix)
#
# Returns: data frame — resolution, n_clusters, mean_silhouette, norm_entropy
# -----------------------------------------------------------------------------
compute_clustering_qc <- function(seurat_obj,
                                   resolutions,
                                   reduction   = "harmony",
                                   prefix      = "SCT_snn_res.",
                                   subsample_n = 5000,
                                   sample_col  = "sample",
                                   seed        = 123,
                                   method      = c("exact", "approx")) {
  method <- match.arg(method)
  check_required_packages(if (method == "approx") "bluster" else "cluster", "compute_clustering_qc()")
  restore_seed <- .push_seed(seed); on.exit(restore_seed(), add = TRUE)

  emb  <- Embeddings(seurat_obj, reduction)
  meta <- seurat_obj@meta.data

  # one subsample + one distance matrix, shared by all resolutions
  idx <- if (method == "exact" && nrow(emb) > subsample_n) sample(nrow(emb), subsample_n) else seq_len(nrow(emb))
  d   <- if (method == "exact") stats::dist(emb[idx, , drop = FALSE]) else NULL

  results <- lapply(resolutions, function(res) {
    col <- paste0(prefix, res)
    if (!col %in% colnames(meta)) {
      cat("  ⚠ Column not found:", col, "— skipping\n"); return(NULL)
    }
    clust   <- as.integer(factor(meta[[col]]))
    n_clust <- length(unique(clust))

    # ── Silhouette ───────────────────────────────────────────────────────────
    mean_sil <- if (method == "exact") {
      sil <- cluster::silhouette(clust[idx], d)
      mean(sil[, 3], na.rm = TRUE)
    } else {
      mean(bluster::approxSilhouette(emb, clust)$width, na.rm = TRUE)
    }

    # ── Normalized Shannon entropy (sample mixing) ────────────────────────────
    if (!sample_col %in% colnames(meta)) {
      norm_ent <- NA_real_
    } else {
      cluster_ent <- tapply(meta[[sample_col]], clust, function(s) {
        tbl <- table(s); p <- tbl / sum(tbl); p <- p[p > 0]
        ent <- -sum(p * log2(p))
        n_s <- length(unique(s))
        if (n_s < 2) return(0)
        ent / log2(n_s)   # normalize to [0, 1]
      })
      norm_ent <- mean(unlist(cluster_ent), na.rm = TRUE)
    }

    data.frame(resolution = res, n_clusters = n_clust,
               mean_silhouette = mean_sil, norm_entropy = norm_ent)
  })
  do.call(rbind, Filter(Negate(is.null), results))
}


# -----------------------------------------------------------------------------
# plot_clustering_qc()
#
# Two-panel publication plot: silhouette width (left) + normalized entropy
# (right) across resolutions. Saves PDF + PNG.
#
# Args:
#   qc_df   : output of compute_clustering_qc()
#   out_dir : output directory (default "QC")
#   palette : bar fill colours for the two panels (silhouette, entropy). Default
#             c("steelblue", "coral3"); a palette name ("viridis", "Set1", ...) takes
#             its first two colours.
#   reverse_palette : swap the two colours (default FALSE)
# -----------------------------------------------------------------------------
plot_clustering_qc <- function(qc_df, out_dir = "QC",
                               palette = c("steelblue", "coral3"), reverse_palette = FALSE) {
  bar_cols <- unname(lab_palette_colors(palette, 2, reverse_palette))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  p_sil <- ggplot(qc_df, aes(x = factor(resolution), y = mean_silhouette)) +
    geom_col(fill = bar_cols[1], alpha = 0.8) +
    geom_point(size = 3, color = "black") +
    labs(title = "Silhouette Width", x = "Resolution", y = "Mean Silhouette") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_ent <- ggplot(qc_df, aes(x = factor(resolution), y = norm_entropy)) +
    geom_col(fill = bar_cols[2], alpha = 0.8) +
    geom_point(size = 3, color = "black") +
    labs(title = "Normalized Sample Entropy\n(Integration Quality)",
         x = "Resolution", y = "Norm. Shannon Entropy") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p <- p_sil | p_ent
  ggsave(file.path(out_dir, "clustering_qc_silhouette_entropy.pdf"),  p, width = 12, height = 5)
  ggsave(file.path(out_dir, "clustering_qc_silhouette_entropy.png"),  p, width = 12, height = 5, dpi = 300)
  cat("✓ Clustering QC plot saved →", out_dir, "\n")
  print(qc_df)
  invisible(p)
}


# =============================================================================
# 2. AUTOMATED CELL TYPE ANNOTATION VIA ENRICHR
# =============================================================================

# -----------------------------------------------------------------------------
# annotate_celltypes_enrichr()
#
# Per-cluster automated cell-type labeling:
#   1. Top markers (avg_log2FC-ranked, FDR < p_adj_cut) selected per cluster
#   2. Zebrafish genes → human orthologs (convert_zf_to_human_alliance)
#   3. 15 Enrichr databases queried
#   4. Fuzzy consensus (Jaro-Winkler >= sim_thresh, spanning >= min_libs)
#   5. Per-cluster XLSX saved to cluster_enrichr/ (color-coded by DB)
#   6. Summary saved to marker_results/enrichr_predictions.xlsx
#
# Args:
#   all_markers   : FindAllMarkers() output — needs cluster, gene, avg_log2FC, p_val_adj
#   alliance_file : path to ORTHOLOGY-ALLIANCE_COMBINED_*.tsv
#   top_n         : max markers sent per cluster (default 250)
#   n_top_terms   : top Enrichr hits collected per DB (default 10)
#   p_adj_cut     : FDR threshold for marker selection (default 0.05)
#   min_genes     : min human orthologs required (default 5)
#   sim_thresh    : Jaro-Winkler threshold for fuzzy term grouping (default 0.85)
#   min_libs      : min libraries required for confident consensus (default 2)
#   out_dir       : directory for per-cluster XLSX (default "cluster_enrichr")
#   summary_dir   : directory for summary XLSX (default "marker_results")
#   workers       : clusters queried concurrently (default 1 = one cluster at a time).
#                   Enrichr answers concurrent requests with HTTP 429 (rate limit), so
#                   values > 1 usually do NOT help; kept as an option (fresh R processes,
#                   SnowParam). Every request times out after
#                   options(mercader.enrichr_timeout = 180) s and a failed / throttled
#                   query is retried 3 times with back-off before the cluster gets NA.
#   BPPARAM       : optional BiocParallel backend (overrides workers)
#
# Returns: summary data frame (cluster, n_markers_sent, cell_type, n_libs_support)
# -----------------------------------------------------------------------------
annotate_celltypes_enrichr <- function(all_markers,
                                        alliance_file,
                                        top_n       = 250,
                                        n_top_terms = 10,
                                        p_adj_cut   = 0.05,
                                        min_genes   = 5,
                                        sim_thresh  = 0.85,
                                        min_libs    = 2,
                                        out_dir     = "cluster_enrichr",
                                        summary_dir = "marker_results",
                                        workers     = 1,
                                        BPPARAM     = NULL) {

  check_required_packages(c("enrichR", "stringdist", "openxlsx"), "annotate_celltypes_enrichr()")
  suppressPackageStartupMessages(library(enrichR))   # enrichR sets its site options in .onAttach: it must be attached

  dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)
  dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

  # ── Ortholog table ──────────────────────────────────────────────────────────
  if (!file.exists(alliance_file))
    stop("Alliance orthology file not found: ", alliance_file,
         "\nDownload from: https://www.alliancegenome.org/downloads")

  ortholog_table <- convert_zf_to_human_alliance(unique(all_markers$gene), alliance_file)

  zf_to_human <- function(genes) {
    hits <- ortholog_table[ortholog_table$zebrafish_gene_id %in% genes, "human_gene_symbol"]
    unique(hits[!is.na(hits) & hits != ""])
  }

  # ── Enrichr setup (in the parent; the package keeps global state) ───────────
  enrichR::setEnrichrSite("Enrichr")
  target_dbs <- c(
    "CellMarker_2024", "Descartes_Cell_Types_and_Tissue_2021",
    "CCLE_Proteomics_2020", "Tabula_Sapiens",
    "HuBMAP_ASCTplusB_augmented_2022", "CellMarker_Augmented_2021",
    "PanglaoDB_Augmented_2021", "Tabula_Muris", "Human_Gene_Atlas",
    "Azimuth_Cell_Types_2021", "ProteomicsDB_2020", "Mouse_Gene_Atlas",
    "ARCHS4_Tissues", "GTEx_Tissues_V8_2023", "Azimuth_2023"
  )
  available_dbs <- tryCatch(enrichR::listEnrichrDbs()$libraryName, error = function(e) character(0))
  dbs_use       <- intersect(target_dbs, available_dbs)
  cat("Querying", length(dbs_use), "Enrichr databases\n")

  # ── XLSX helper (single sheet, color-coded by DB) ───────────────────────────
  db_palette <- c("FFE4E1","E1F0FF","E1FFE8","FFF8E1","F0E1FF",
                  "E1FFFF","FFE8D6","E8FFE1","FFE1F0","E1E8FF",
                  "FFFDE1","F5E1FF","E1FFF5","FFE1E1","E1FFEB")

  .write_xlsx <- function(per_db_results, filepath) {
    combined    <- do.call(rbind, per_db_results)
    dbs_present <- unique(combined$database)
    wb          <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, "Enrichr_results")
    hdr_style   <- openxlsx::createStyle(fontColour="#000000", fgFill="#D9D9D9", halign="left",
                                         textDecoration="Bold", border="Bottom", borderColour="#888888")
    openxlsx::writeData(wb, sheet=1, combined, headerStyle=hdr_style)
    color_map <- setNames(db_palette[seq_along(dbs_present)], dbs_present)
    for (db in dbs_present) {
      rows  <- which(combined$database == db) + 1L
      style <- openxlsx::createStyle(fgFill=paste0("#", color_map[[db]]))
      openxlsx::addStyle(wb, sheet=1, style=style, rows=rows, cols=seq_len(ncol(combined)),
                         gridExpand=TRUE, stack=TRUE)
    }
    openxlsx::freezePane(wb, sheet=1, firstRow=TRUE)
    openxlsx::setColWidths(wb, sheet=1, cols=seq_len(ncol(combined)), widths="auto")
    openxlsx::saveWorkbook(wb, filepath, overwrite=TRUE)
  }

  # ── Fuzzy consensus ─────────────────────────────────────────────────────────
  .fuzzy_consensus <- function(top_terms, sim_thr=sim_thresh, min_l=min_libs) {
    all_terms <- unlist(top_terms)
    clean     <- trimws(gsub("\\s+(Human|Mouse|FACS|10x|Smart-seq2)$", "", all_terms, ignore.case=TRUE))
    if (length(clean) <= 1) return(NULL)
    sim_mat  <- as.matrix(1 - stringdist::stringdistmatrix(tolower(clean), method="jw", p=0.1))
    groups   <- integer(length(clean)); gid <- 0L
    for (i in seq_along(clean)) {
      if (groups[i] != 0L) next
      gid <- gid + 1L; groups[i] <- gid
      for (j in seq_along(clean))
        if (j != i && groups[j] == 0L && sim_mat[i,j] >= sim_thr) groups[j] <- gid
    }
    db_per_term <- rep(names(top_terms), lengths(top_terms))
    group_libs  <- tapply(db_per_term, groups, function(x) length(unique(x)))
    best_gid    <- as.integer(names(which.max(group_libs)))
    if (group_libs[as.character(best_gid)] < min_l) return(NULL)
    members <- clean[groups == best_gid]; orig <- all_terms[groups == best_gid]
    list(label=members[which.min(nchar(members))],
         n_libs=as.integer(group_libs[as.character(best_gid)]),
         supporting=unique(orig))
  }

  # ── One cluster (runs on a worker) ───────────────────────────────────────────
  clusters <- sort(unique(all_markers$cluster))
  one_cluster <- function(cl) {
    # Workers are fresh R processes (SnowParam): attach enrichR there and give every HTTP
    # request a timeout so a stalled Enrichr call fails (-> NA row) instead of hanging.
    if (is.null(getOption("enrichR.base.address"))) {
      suppressPackageStartupMessages(library(enrichR)); enrichR::setEnrichrSite("Enrichr")
    }
    if (requireNamespace("httr", quietly = TRUE)) {
      httr::set_config(httr::timeout(getOption("mercader.enrichr_timeout", 180)))
      on.exit(httr::reset_config(), add = TRUE)
    }
    cat("\n-- Cluster", cl, "--\n")
    blank_row <- data.frame(cluster=cl, n_markers_sent=0, cell_type=NA_character_,
                             n_libs_support=0L, supporting_terms=NA_character_,
                             stringsAsFactors=FALSE)
    cl_df <- all_markers[all_markers$cluster == cl & !is.na(all_markers$p_val_adj) &
                         all_markers$p_val_adj < p_adj_cut & all_markers$avg_log2FC > 0, ]

    if (nrow(cl_df) < min_genes) { cat("  Too few markers — skipping\n"); return(blank_row) }
    cl_df       <- cl_df[order(-cl_df$avg_log2FC), ]
    zf_genes    <- head(cl_df$gene, top_n)
    human_genes <- zf_to_human(zf_genes)

    if (length(human_genes) < min_genes) {
      cat("  Too few orthologs — skipping\n")
      blank_row$n_markers_sent <- length(human_genes); return(blank_row)
    }
    cat("  Sending", length(human_genes), "human genes to Enrichr\n")

    enriched <- NULL
    for (attempt in 1:3) {        # retry: Enrichr returns HTTP 429 / partial results when busy
      enriched <- tryCatch(enrichR::enrichr(human_genes, dbs_use),
                           error=function(e) { cat("  Enrichr error:", e$message, "\n"); NULL })
      ok <- is.list(enriched) && any(vapply(enriched, is.data.frame, logical(1)))
      if (ok) break
      enriched <- NULL
      if (attempt < 3) { cat("  Enrichr returned no usable result — retrying in", 15 * attempt, "s\n"); Sys.sleep(15 * attempt) }
    }
    if (is.null(enriched)) { blank_row$n_markers_sent <- length(human_genes); return(blank_row) }

    per_db <- lapply(dbs_use, function(db) {
      res <- enriched[[db]]
      if (!is.data.frame(res) || nrow(res) == 0) return(NULL)
      res$database <- db; res$zf_genes <- paste(zf_genes, collapse=";")
      res$human_genes <- paste(human_genes, collapse=";")
      head(res[order(res$Adjusted.P.value), ], n_top_terms)
    })
    names(per_db) <- dbs_use; per_db <- Filter(Negate(is.null), per_db)

    if (length(per_db) > 0)
      .write_xlsx(per_db, file.path(out_dir, paste0("cluster_", cl, "_enrichr.xlsx")))

    consensus <- tryCatch(.fuzzy_consensus(lapply(per_db, function(d) d$Term)),
                          error=function(e) NULL)

    if (!is.null(consensus)) {
      cat("  ->", consensus$label, "(", consensus$n_libs, "libraries)\n")
      data.frame(cluster=cl, n_markers_sent=length(human_genes),
                 cell_type=consensus$label, n_libs_support=consensus$n_libs,
                 supporting_terms=paste(consensus$supporting, collapse=" | "),
                 stringsAsFactors=FALSE)
    } else {
      cat("  -> No confident consensus (< 2 libraries) — NA assigned\n")
      blank_row$n_markers_sent <- length(human_genes)
      blank_row
    }
  }

  # ── Run clusters (a few at a time) ───────────────────────────────────────────
  # Fresh worker processes rather than forks (forked children share the parent's open
  # sockets and Enrichr requests were seen to hang). Default is serial: Enrichr rate-limits.
  w  <- lab_workers(length(clusters), "light", min(workers, 4L))
  bp <- if (!is.null(BPPARAM)) BPPARAM else if (w == 1L) BiocParallel::SerialParam() else
        BiocParallel::SnowParam(workers = w, tasks = length(clusters))
  out <- lab_run_tasks(as.list(clusters), one_cluster, bp)
  summary_rows <- list()
  for (i in seq_along(clusters)) {
    cl <- clusters[i]; o <- out[[i]]
    if (length(o$log) > 0) cat(paste(o$log, collapse = "\n"), "\n")
    r <- o$result
    if (is.list(r) && identical(r$status, "error")) {
      cat("  ✗ Error:", r$message, "\n")
      r <- data.frame(cluster=cl, n_markers_sent=0, cell_type=NA_character_, n_libs_support=0L,
                      supporting_terms=NA_character_, stringsAsFactors=FALSE)
    }
    summary_rows[[as.character(cl)]] <- r
  }

  # ── Summary XLSX ─────────────────────────────────────────────────────────────
  enrichr_summary <- do.call(rbind, summary_rows)
  cat("\n=== Enrichr predictions ===\n")
  print(enrichr_summary[, c("cluster","cell_type","n_libs_support")])

  wb_s <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb_s, "predictions")
  hdr  <- openxlsx::createStyle(fontColour="#000000", fgFill="#D9D9D9", halign="left",
                                textDecoration="Bold", border="Bottom", borderColour="#888888")
  openxlsx::writeData(wb_s, sheet=1, enrichr_summary, headerStyle=hdr)
  openxlsx::freezePane(wb_s, sheet=1, firstRow=TRUE)
  openxlsx::setColWidths(wb_s, sheet=1, cols=seq_len(ncol(enrichr_summary)), widths="auto")
  openxlsx::saveWorkbook(wb_s, file.path(summary_dir, "enrichr_predictions.xlsx"), overwrite=TRUE)
  cat("Summary saved:", file.path(summary_dir, "enrichr_predictions.xlsx"), "\n")
  return(invisible(enrichr_summary))
}


# =============================================================================
# 3. DATA CHARACTERISTICS & NGUYEN 2023 DE METHOD SELECTION
# =============================================================================

# -----------------------------------------------------------------------------
# compute_data_characteristics()
#
# Compute sparsity (zero rate) and median UMI depth per cell type.
# Call before run_de_per_celltype() to inform method selection.
#
# Args:
#   seurat_obj    : Seurat object (RNA assay = raw counts)
#   cell_type_col : metadata column with the cell types
#   counts        : optional raw counts matrix (from .sc_counts()); read from the
#                   object when NULL
#
# Returns: data frame — cell_type, sparsity, median_depth, n_cells
# -----------------------------------------------------------------------------
compute_data_characteristics <- function(seurat_obj, cell_type_col = "cell_type", counts = NULL) {
  if (is.null(counts)) counts <- .sc_counts(seurat_obj, "RNA")
  meta    <- seurat_obj@meta.data
  ct_list <- levels(factor(meta[[cell_type_col]]))
  results <- lapply(ct_list, function(ct) {
    cells <- .sc_cells_of(meta, ct, cell_type_col)
    if (length(cells) < 10) return(NULL)
    mat <- counts[, cells, drop = FALSE]
    data.frame(cell_type    = ct,
               sparsity     = round(1 - Matrix::nnzero(mat) / prod(dim(mat)), 3),   # no dense (mat == 0)
               median_depth = round(median(Matrix::colSums(mat)), 1),
               n_cells      = length(cells))
  })
  out <- do.call(rbind, Filter(Negate(is.null), results))
  cat("Data characteristics:\n"); print(out); invisible(out)
}


# -----------------------------------------------------------------------------
# recommend_de_method()
#
# Nguyen 2023 (doi:10.1093/bib/bbad218) lookup table.
# Categorizes sparsity (>= 0.6 = high), depth (>= 40 = high, 8-39 = medium,
# < 8 = low) and batch status to recommend DE workflows.
#
# Returns: list(primary, recommended, profile)
# -----------------------------------------------------------------------------
recommend_de_method <- function(sparsity, depth, has_batch = FALSE) {

  nguyen_table <- list(
    # ── High sparsity (≈80%) ──────────────────────────────────────────────────
    high_high_substantial = list(
      label   = "High sparsity / High depth / Substantial batch",
      methods = c("MAST_Cov","ZW_edgeR_Cov","ZW_DESeq2_Cov","scVI_limmatrend","DESeq2_FEM","limmatrend_Cov"),
      primary = "limmatrend_Cov"),
    high_high_small = list(
      label   = "High sparsity / High depth / Small batch",
      methods = c("MAST","ZW_edgeR","ZW_DESeq2","Pseudobulk_limma","DESeq2_FEM","limmatrend_Cov"),
      primary = "Pseudobulk_limma"),
    high_medium_substantial = list(
      label   = "High sparsity / Medium depth / Substantial batch",
      methods = c("DESeq2_Cov","limmatrend_Cov","DESeq2_wFisher","LogN_FEM","MAST_Cov","Raw_Wilcox"),
      primary = "limmatrend_Cov"),
    high_medium_small = list(
      label   = "High sparsity / Medium depth / Small batch",
      methods = c("DESeq2","limmatrend_Cov","LogN_FEM","Pseudobulk_edgeR","Pseudobulk_limma","Raw_Wilcox"),
      primary = "Pseudobulk_limma"),
    high_low_substantial = list(
      label   = "High sparsity / Low depth / Substantial batch",
      methods = c("LogN_FEM","limmatrend","Raw_Wilcox","RISC_QP"),
      primary = "limmatrend"),
    high_low_small = list(
      label   = "High sparsity / Low depth / Small batch",
      methods = c("LogN_FEM","limmatrend","Raw_Wilcox","RISC_QP"),
      primary = "limmatrend"),

    # ── Low sparsity (≈40%) ───────────────────────────────────────────────────
    low_high_substantial = list(
      label   = "Low sparsity / High depth / Substantial batch",
      methods = c("MAST_Cov","ZW_edgeR_Cov","ZW_DESeq2_Cov","limma_BEC_Wilcox","Scanorama_limmatrend","LogN_FEM"),
      primary = "limmatrend_Cov"),
    low_high_small = list(
      label   = "Low sparsity / High depth / Small batch",
      methods = c("Pseudobulk_limma","Raw_Wilcox","Pseudobulk_edgeR","DESeq2","LogN_FEM"),
      primary = "Pseudobulk_limma"),
    low_low_substantial = list(
      label   = "Low sparsity / Low depth / Substantial batch",
      methods = c("limmavoom_Cov","limmatrend_Cov","ZW_edgeR_Cov","ZW_DESeq2_Cov","LogN_FEM","limma_BEC_Wilcox"),
      primary = "limmatrend_Cov"),
    low_low_small = list(
      label   = "Low sparsity / Low depth / Small batch",
      methods = c("Pseudobulk_limma","limmavoom","limmatrend","Pseudobulk_edgeR","Raw_Wilcox","edgeR"),
      primary = "Pseudobulk_limma")
  )

  sp_cat <- if (sparsity >= 0.6) "high" else "low"
  dp_cat <- if (depth >= 40) "high" else if (depth >= 8) "medium" else "low"
  if (sp_cat == "low" && dp_cat == "medium") dp_cat <- "low"   # no medium row for low sparsity
  bt_cat <- if (has_batch) "substantial" else "small"
  key    <- paste(sp_cat, dp_cat, bt_cat, sep="_")

  if (!key %in% names(nguyen_table)) {
    cat("⚠ No exact match for key:", key, "— falling back to high_medium_small\n")
    key <- "high_medium_small"
  }
  row <- nguyen_table[[key]]

  cat("─────────────────────────────────────────────────────\n")
  cat("DE method recommendation (Nguyen 2023)\n")
  cat("  Profile  :", row$label, "\n")
  cat("  Sparsity :", round(sparsity*100), "% →", sp_cat, "\n")
  cat("  Depth    :", depth, "UMIs →", dp_cat, "\n")
  cat("  Batch    :", bt_cat, "\n")
  cat("  PRIMARY  :", row$primary, "\n")
  cat("  All      :", paste(row$methods, collapse=", "), "\n")
  cat("─────────────────────────────────────────────────────\n")

  list(primary=row$primary, recommended=row$methods,
       profile=list(sparsity_cat=sp_cat, depth_cat=dp_cat, batch_cat=bt_cat))
}


# =============================================================================
# 4. NGUYEN 2023 DE METHOD IMPLEMENTATIONS
# =============================================================================
# All single-cell methods expect a Seurat object with RNA assay (raw counts).
# Standard output columns: gene, log2FoldChange, pvalue, padj, method
#
# Each public run_*() slices the raw counts + metadata of one cell type
# (.sc_slice_obj) and hands them to an internal .de_*() worker that does the
# statistics on plain matrices. run_de_per_celltype() calls the .de_*() workers
# directly on slices it builds once, in parallel.

# counts + metadata of one cell type, read from the object
.sc_slice_obj <- function(seurat_obj, cell_type, cell_type_col, counts = NULL) {
  meta  <- seurat_obj@meta.data
  cells <- .sc_cells_of(meta, cell_type, cell_type_col)
  if (is.null(counts)) counts <- .sc_counts(seurat_obj, "RNA")
  list(cells = cells, counts = counts[, cells, drop = FALSE], meta = meta[cells, , drop = FALSE])
}

.de_tidy_limma <- function(fit, method, coef = NULL) {
  tt <- if (is.null(coef)) limma::topTable(fit, number = Inf, sort.by = "P") else
          limma::topTable(fit, n = Inf, adjust.method = "BH", coef = coef)
  tt %>%
    dplyr::rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val) %>%
    tibble::rownames_to_column("gene") %>% dplyr::mutate(method = method)
}

# -----------------------------------------------------------------------------
# run_pseudobulk_limma()   ← Pseudobulk_limma
# Aggregates per sample → limma-trend. Best all-around primary method.
# -----------------------------------------------------------------------------
.de_pseudobulk_limma <- function(counts, meta, cell_type, group_col, sample_col,
                                 contrast_str = NULL, min_cells = 5) {
  if (ncol(counts) < min_cells*2) { cat("  ⚠ Too few cells for", cell_type, "\n"); return(NULL) }
  pb_mat <- .sc_pseudobulk(counts, meta[[sample_col]], min_cells)
  if (is.null(pb_mat) || ncol(pb_mat) < 3) { cat("  ⚠ Too few samples for", cell_type, "\n"); return(NULL) }

  pb_meta <- meta[!duplicated(meta[[sample_col]]), ]
  rownames(pb_meta) <- pb_meta[[sample_col]]; pb_meta <- pb_meta[colnames(pb_mat),, drop=FALSE]

  y      <- edgeR::DGEList(counts=pb_mat); keep <- edgeR::filterByExpr(y, group=pb_meta[[group_col]])
  y      <- y[keep, ]
  design <- model.matrix(~ 0 + factor(pb_meta[[group_col]]))
  colnames(design) <- levels(factor(pb_meta[[group_col]]))

  if (is.null(contrast_str)) {
    gs <- levels(factor(pb_meta[[group_col]]))
    if (length(gs) != 2) stop("contrast_str required when >2 groups")
    contrast_str <- paste(rev(gs), collapse=" - ")
  }
  cont_mat <- limma::makeContrasts(contrasts=contrast_str, levels=design)
  logcpm   <- edgeR::cpm(y, log=TRUE, prior.count=3)
  fit      <- limma::lmFit(logcpm, design) %>% limma::contrasts.fit(cont_mat) %>%
              limma::eBayes(trend=TRUE, robust=TRUE)
  res <- .de_tidy_limma(fit, "Pseudobulk_limma")
  cat("  ✓ Pseudobulk_limma:", nrow(res), "genes,", sum(res$padj < 0.05, na.rm=TRUE), "sig\n")
  return(res)
}

run_pseudobulk_limma <- function(seurat_obj, cell_type,
                                  group_col="group", sample_col="sample",
                                  cell_type_col="cell_type", contrast_str=NULL,
                                  min_cells=5) {
  check_required_packages(c("limma", "edgeR"), "run_pseudobulk_limma()")
  sl <- .sc_slice_obj(seurat_obj, cell_type, cell_type_col)
  .de_pseudobulk_limma(sl$counts, sl$meta, cell_type, group_col, sample_col, contrast_str, min_cells)
}


# -----------------------------------------------------------------------------
# run_limmatrend_cov()   ← limmatrend / limmatrend_Cov
# Single-cell limma-trend + optional batch covariate (TMM-normalized logCPM).
#   filter_genes : FALSE (default) = every gene, as before. TRUE = keep only genes
#                  detected (count > 0) in >= max(3, 1 % of the cells), decided on the
#                  sparse matrix BEFORE the dense copy is made (much less memory on big
#                  cell types; changes the TMM factors / prior slightly, so it is opt-in).
# -----------------------------------------------------------------------------
.de_limmatrend_cov <- function(counts, meta, group_col, sample_col, has_batch = TRUE,
                               filter_genes = FALSE) {
  if (ncol(counts) < 10) { cat("  ⚠ Too few cells\n"); return(NULL) }
  cellinfo <- data.frame(Group=factor(meta[[group_col]]), Batch=factor(meta[[sample_col]]))
  if (filter_genes) counts <- .sc_filter_genes(counts)

  y      <- edgeR::DGEList(as.matrix(counts))
  y      <- edgeR::calcNormFactors(y)
  logCPM <- edgeR::cpm(y, log=TRUE, prior.count=3)
  design <- if (has_batch) model.matrix(~ Group + Batch, data=cellinfo) else
                           model.matrix(~ Group, data=cellinfo)
  fit    <- limma::lmFit(logCPM, design) %>% limma::eBayes(trend=TRUE, robust=TRUE)
  res    <- .de_tidy_limma(fit, ifelse(has_batch, "limmatrend_Cov", "limmatrend"), coef = 2)
  cat("  ✓", unique(res$method), ":", nrow(res), "genes,", sum(res$padj < 0.05, na.rm=TRUE), "sig\n")
  return(res)
}

run_limmatrend_cov <- function(seurat_obj, cell_type,
                                group_col="group", sample_col="sample",
                                cell_type_col="cell_type", has_batch=TRUE,
                                filter_genes=FALSE) {
  check_required_packages(c("limma", "edgeR"), "run_limmatrend_cov()")
  sl <- .sc_slice_obj(seurat_obj, cell_type, cell_type_col)
  .de_limmatrend_cov(sl$counts, sl$meta, group_col, sample_col, has_batch, filter_genes)
}


# -----------------------------------------------------------------------------
# run_limmavoom_cov()   ← limmavoom / limmavoom_Cov
# limma-voom with TMM normalization + optional batch covariate.
# -----------------------------------------------------------------------------
.de_limmavoom_cov <- function(counts, meta, group_col, sample_col, has_batch = TRUE,
                              filter_genes = FALSE) {
  if (ncol(counts) < 10) { cat("  ⚠ Too few cells\n"); return(NULL) }
  cellinfo  <- data.frame(Group=factor(meta[[group_col]]), Batch=factor(meta[[sample_col]]))
  if (filter_genes) counts <- .sc_filter_genes(counts)
  count_mat <- as.matrix(counts)
  nf        <- edgeR::calcNormFactors(count_mat)
  design    <- if (has_batch) model.matrix(~ Group + Batch, data=cellinfo) else
                              model.matrix(~ Group, data=cellinfo)
  voom_d    <- limma::voom(count_mat, design=design, lib.size=colSums(count_mat)*nf)
  fit       <- limma::lmFit(voom_d, design) %>% limma::eBayes()
  res       <- .de_tidy_limma(fit, ifelse(has_batch, "limmavoom_Cov", "limmavoom"), coef = 2)
  cat("  ✓", unique(res$method), ":", nrow(res), "genes,", sum(res$padj < 0.05, na.rm=TRUE), "sig\n")
  return(res)
}

run_limmavoom_cov <- function(seurat_obj, cell_type,
                               group_col="group", sample_col="sample",
                               cell_type_col="cell_type", has_batch=TRUE,
                               filter_genes=FALSE) {
  check_required_packages(c("limma", "edgeR"), "run_limmavoom_cov()")
  sl <- .sc_slice_obj(seurat_obj, cell_type, cell_type_col)
  .de_limmavoom_cov(sl$counts, sl$meta, group_col, sample_col, has_batch, filter_genes)
}


# -----------------------------------------------------------------------------
# run_edgeR_cov()   ← Pseudobulk_edgeR / ZW_edgeR_Cov
# glmQLFit edgeR. pseudobulk=TRUE → Pseudobulk_edgeR; FALSE → ZW_edgeR_Cov.
# -----------------------------------------------------------------------------
.de_edgeR_cov <- function(counts, meta, group_col, sample_col, has_batch = TRUE,
                          pseudobulk = FALSE, filter_genes = FALSE) {
  if (ncol(counts) < 10) { cat("  ⚠ Too few cells\n"); return(NULL) }

  if (pseudobulk) {
    count_mat <- .sc_pseudobulk(counts, meta[[sample_col]], 3)
    if (is.null(count_mat) || ncol(count_mat) < 3) { cat("  ⚠ Too few samples\n"); return(NULL) }
    pb_meta   <- meta[!duplicated(meta[[sample_col]]), ]; rownames(pb_meta) <- pb_meta[[sample_col]]
    pb_meta   <- pb_meta[colnames(count_mat),, drop=FALSE]
    cellGroup  <- factor(pb_meta[[group_col]]); cellBatch <- factor(pb_meta[[sample_col]])
  } else {
    cellGroup  <- factor(meta[[group_col]]); cellBatch <- factor(meta[[sample_col]])
    if (filter_genes) counts <- .sc_filter_genes(counts)
    count_mat <- as.matrix(counts)
  }

  y      <- edgeR::DGEList(counts=count_mat, group=cellGroup) %>% edgeR::calcNormFactors()
  design <- if (has_batch) model.matrix(~ cellGroup + cellBatch) else model.matrix(~ cellGroup)
  rownames(design) <- colnames(y)
  y      <- edgeR::estimateDisp(y, design, robust=TRUE)
  qlf    <- edgeR::glmQLFit(y, design, robust=TRUE, prior.df=0) %>% edgeR::glmQLFTest(coef=2)
  res    <- qlf$table %>%
    dplyr::mutate(padj=p.adjust(PValue, method="BH")) %>%
    dplyr::rename(log2FoldChange=logFC, pvalue=PValue) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(method=if(pseudobulk) "Pseudobulk_edgeR" else
                  ifelse(has_batch, "ZW_edgeR_Cov", "edgeR"))

  cat("  ✓", unique(res$method), ":", nrow(res), "genes,", sum(res$padj < 0.05, na.rm=TRUE), "sig\n")
  return(res)
}

run_edgeR_cov <- function(seurat_obj, cell_type,
                           group_col="group", sample_col="sample",
                           cell_type_col="cell_type", has_batch=TRUE, pseudobulk=FALSE,
                           filter_genes=FALSE) {
  check_required_packages("edgeR", "run_edgeR_cov()")
  sl <- .sc_slice_obj(seurat_obj, cell_type, cell_type_col)
  .de_edgeR_cov(sl$counts, sl$meta, group_col, sample_col, has_batch, pseudobulk, filter_genes)
}


# -----------------------------------------------------------------------------
# run_mast_cov()   ← MAST / MAST_Cov
# MAST via Seurat FindMarkers. has_batch=TRUE adds sample as latent variable.
# Runs on the full object with the two cell-name vectors as idents (no subset
# copy; identical statistics). MAST uses options(mc.cores) — set to the "heavy"
# worker count in the parent, 1 inside a parallel worker.
# -----------------------------------------------------------------------------
.de_mast_cov <- function(seurat_obj, cells, group, group_col, sample_col, has_batch = TRUE) {
  if (length(cells) < 10) { cat("  ⚠ Too few cells\n"); return(NULL) }
  groups  <- sort(levels(factor(group)))
  cells_1 <- cells[group == groups[1]]; cells_2 <- cells[group == groups[2]]
  # Seurat reads a length-1 ident as an identity-class name, not a cell: for a 1-cell
  # group fall back to subsetting + group idents (the previous approach).
  if (length(cells_1) < 2 || length(cells_2) < 2) {
    obj <- seurat_obj[, cells]; Idents(obj) <- group
    id_1 <- groups[1]; id_2 <- groups[2]
  } else { obj <- seurat_obj; id_1 <- cells_1; id_2 <- cells_2 }

  markers <- tryCatch(.sc_with_mc_cores({
    if (has_batch)
      FindMarkers(obj, ident.1=id_2, ident.2=id_1, test.use="MAST",
                  logfc.threshold=0, latent.vars=sample_col,
                  min.cells.feature=0, min.cells.group=0, min.pct=0)
    else
      FindMarkers(obj, ident.1=id_2, ident.2=id_1, test.use="MAST",
                  logfc.threshold=0, min.cells.feature=0, min.cells.group=0, min.pct=0)
  }), error=function(e) { cat("  ✗ MAST error:", e$message, "\n"); NULL })

  if (is.null(markers)) return(NULL)
  res <- markers %>%
    dplyr::rename(log2FoldChange=avg_log2FC, pvalue=p_val, padj=p_val_adj) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(method=ifelse(has_batch, "MAST_Cov", "MAST"))

  cat("  ✓", unique(res$method), ":", nrow(res), "genes,", sum(res$padj < 0.05, na.rm=TRUE), "sig\n")
  return(res)
}

run_mast_cov <- function(seurat_obj, cell_type,
                          group_col="group", sample_col="sample",
                          cell_type_col="cell_type", has_batch=TRUE) {
  check_required_packages("MAST", "run_mast_cov()")
  meta  <- seurat_obj@meta.data
  cells <- .sc_cells_of(meta, cell_type, cell_type_col)
  .de_mast_cov(seurat_obj, cells, meta[cells, group_col], group_col, sample_col, has_batch)
}


# -----------------------------------------------------------------------------
# run_wilcox_raw()   ← Raw_Wilcox
# Wilcoxon on log-normalized counts. Fallback for very low depth data only.
# WARNING: does not correct for batch effects.
# The RNA counts of the cell type are log-normalized (LogNormalize is per cell,
# so this equals normalizing the whole object) into a minimal Seurat object.
# -----------------------------------------------------------------------------
.de_wilcox_raw <- function(counts, meta, group_col) {
  if (ncol(counts) < 10) { cat("  ⚠ Too few cells\n"); return(NULL) }
  obj <- CreateSeuratObject(counts = counts, meta.data = meta)
  obj <- NormalizeData(obj, verbose=FALSE)
  Idents(obj) <- obj@meta.data[[group_col]]
  groups <- sort(levels(factor(obj@meta.data[[group_col]])))

  markers <- tryCatch(
    FindMarkers(obj, ident.1=groups[2], ident.2=groups[1], test.use="wilcox",
                logfc.threshold=0, min.cells.feature=0, min.cells.group=0, min.pct=0),
    error=function(e) { cat("  ✗ Wilcox error:", e$message, "\n"); NULL })

  if (is.null(markers)) return(NULL)
  res <- markers %>%
    dplyr::rename(log2FoldChange=avg_log2FC, pvalue=p_val, padj=p_val_adj) %>%
    tibble::rownames_to_column("gene") %>% dplyr::mutate(method="Raw_Wilcox")

  cat("  ✓ Raw_Wilcox:", nrow(res), "genes,", sum(res$padj < 0.05, na.rm=TRUE), "sig\n")
  return(res)
}

run_wilcox_raw <- function(seurat_obj, cell_type,
                            group_col="group", cell_type_col="cell_type") {
  sl <- .sc_slice_obj(seurat_obj, cell_type, cell_type_col)
  .de_wilcox_raw(sl$counts, sl$meta, group_col)
}


# -----------------------------------------------------------------------------
# run_de_per_celltype()
#
# Master dispatcher:
#   1. Runs compute_data_characteristics() per cell type (or takes `char_df`)
#   2. Calls recommend_de_method() for each
#   3. Dispatches to the recommended primary method — cell types in PARALLEL
#   4. Falls back to Pseudobulk_limma if method not implemented
#
# Args:
#   seurat_obj    : Seurat object
#   cell_types    : character vector of cell types (NULL = all in cell_type_col)
#   group_col     : metadata group column
#   sample_col    : metadata sample/batch column
#   cell_type_col : metadata cell type column
#   contrast_str  : contrast string for pseudobulk methods
#   has_batch     : logical override (NULL = auto-detect from sample/group counts)
#   char_df       : output of compute_data_characteristics() if already computed
#                   (saves recomputing it)
#   workers       : cell types run concurrently (default NULL = lab_workers("heavy"):
#                   min(12, cores/4); 1 = one cell type at a time, as before)
#   BPPARAM       : optional BiocParallel backend (overrides workers)
#
# Returns: named list [[cell_type]] of DE result data frames
# -----------------------------------------------------------------------------
run_de_per_celltype <- function(seurat_obj,
                                 cell_types    = NULL,
                                 group_col     = "group",
                                 sample_col    = "sample",
                                 cell_type_col = "cell_type",
                                 contrast_str  = NULL,
                                 has_batch     = NULL,
                                 char_df       = NULL,
                                 workers       = NULL,
                                 BPPARAM       = NULL) {
  check_required_packages(c("limma", "edgeR"), "run_de_per_celltype()")
  meta <- seurat_obj@meta.data
  if (is.null(cell_types))
    cell_types <- levels(factor(meta[[cell_type_col]]))

  if (is.null(has_batch)) {
    n_samples  <- length(unique(meta[[sample_col]]))
    n_groups   <- length(unique(meta[[group_col]]))
    has_batch  <- (n_samples > n_groups)
    cat("Auto-detected batch:", has_batch, "(", n_samples, "samples,", n_groups, "groups)\n")
  }

  counts <- .sc_counts(seurat_obj, "RNA")                       # read once, shared by every task
  if (is.null(char_df)) char_df <- compute_data_characteristics(seurat_obj, cell_type_col, counts = counts)

  # ---- method per cell type (serial, prints the recommendation as before) ----
  plan <- lapply(cell_types, function(ct) {
    cat("\n=== DE for:", ct, "===\n")
    ct_char <- char_df[char_df$cell_type == ct, ]
    rec <- if (nrow(ct_char) == 0) {
      cat("  ⚠ No characteristics — defaulting to Pseudobulk_limma\n")
      list(primary="Pseudobulk_limma")
    } else {
      recommend_de_method(ct_char$sparsity, ct_char$median_depth, has_batch)
    }
    list(cell_type = ct, method = rec$primary)
  })

  # ---- one cell type (runs on a worker; counts / seurat_obj inherited copy-on-write)
  run_one <- function(task) {
    ct <- task$cell_type; m <- task$method
    cells <- .sc_cells_of(meta, ct, cell_type_col)
    sl_counts <- counts[, cells, drop = FALSE]; sl_meta <- meta[cells, , drop = FALSE]
    switch(m,
      "Pseudobulk_limma"  = .de_pseudobulk_limma(sl_counts, sl_meta, ct, group_col, sample_col, contrast_str),
      "limmatrend_Cov"    = .de_limmatrend_cov(sl_counts, sl_meta, group_col, sample_col, has_batch=TRUE),
      "limmatrend"        = .de_limmatrend_cov(sl_counts, sl_meta, group_col, sample_col, has_batch=FALSE),
      "limmavoom_Cov"     = .de_limmavoom_cov(sl_counts, sl_meta, group_col, sample_col, has_batch=TRUE),
      "limmavoom"         = .de_limmavoom_cov(sl_counts, sl_meta, group_col, sample_col, has_batch=FALSE),
      "Pseudobulk_edgeR"  = .de_edgeR_cov(sl_counts, sl_meta, group_col, sample_col, has_batch=FALSE, pseudobulk=TRUE),
      "ZW_edgeR_Cov"      = .de_edgeR_cov(sl_counts, sl_meta, group_col, sample_col, has_batch=TRUE,  pseudobulk=FALSE),
      "MAST_Cov"          = .de_mast_cov(seurat_obj, cells, sl_meta[[group_col]], group_col, sample_col, has_batch=TRUE),
      "MAST"              = .de_mast_cov(seurat_obj, cells, sl_meta[[group_col]], group_col, sample_col, has_batch=FALSE),
      "Raw_Wilcox"        = .de_wilcox_raw(sl_counts, sl_meta, group_col),
      {
        cat("  ⚠ Method", m, "not implemented — falling back to Pseudobulk_limma\n")
        .de_pseudobulk_limma(sl_counts, sl_meta, ct, group_col, sample_col, contrast_str)
      }
    )
  }

  bp  <- lab_bpparam(workers, length(plan), BPPARAM, profile = "heavy")
  cat("\nRunning", length(plan), "cell type(s) on", BiocParallel::bpnworkers(bp), "worker(s) [", class(bp)[1], "]\n")
  cost <- vapply(plan, function(t) sum(meta[[cell_type_col]] == t$cell_type, na.rm = TRUE), numeric(1))
  out  <- lab_run_tasks(plan, run_one, bp, cost = cost)

  results <- list()
  for (i in seq_along(plan)) {
    ct <- plan[[i]]$cell_type; o <- out[[i]]
    cat("\n---", ct, "(", plan[[i]]$method, ")\n")
    if (length(o$log) > 0) cat(paste(o$log, collapse = "\n"), "\n")
    r <- o$result
    if (is.list(r) && !is.data.frame(r) && identical(r$status, "error")) {
      cat("  ✗ Error:", r$message, "\n"); r <- NULL
    }
    if (!is.null(r)) results[[ct]] <- r   # as before: cell types without a result are omitted
  }
  cat("\n✓✓✓ DE complete for all cell types ✓✓✓\n")
  return(results)
}

# Legacy alias
pseudobulk_de <- function(...) run_pseudobulk_limma(...)


# =============================================================================
# 5. LIANA — LIGAND-RECEPTOR ANALYSIS
# =============================================================================

run_liana_zebrafish <- function(target_organism=7955, resource_name="Consensus") {
  check_required_packages("liana", "run_liana_zebrafish()")
  cat("Generating zebrafish LIANA ortholog resource...\n")
  or <- liana::generate_homologs(op_resource=liana::select_resource(resource_name)[[1]],
                                  target_organism=target_organism)
  cat("✓ Resource ready:", nrow(or), "interactions\n"); return(or)
}

# run_liana_loop()
#   liana_workers : passed to liana_wrap(parallelize = TRUE, workers = ) — parallelism
#                   INSIDE one LIANA run (NULL = liana's default, serial)
#   workers       : conditions run concurrently (default 1; set e.g. 3 for 3 conditions —
#                   each worker holds one condition's subset of the object)
#   palette, reverse_palette : colour scale of the dotplots, passed to plot_liana_dotplot()
#                   (default NULL = liana's own scale; any lab palette name or colour vector)
run_liana_loop <- function(seurat_obj, conditions, ortholog_resource,
                            source_groups, target_groups,
                            condition_col="condition", cell_type_col="cell_type",
                            agg_rank_cutoff=0.05, ntop=20, out_dir="LR_analysis",
                            liana_workers=NULL, workers=1, BPPARAM=NULL,
                            palette=NULL, reverse_palette=FALSE) {
  check_required_packages("liana", "run_liana_loop()")
  if (!is.null(palette)) lab_palette_colors(palette, 2)   # fail fast on an invalid palette
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
  if (is.character(source_groups)) source_groups <- list(all=source_groups)
  if (is.character(target_groups)) target_groups <- list(all=target_groups)
  if (length(source_groups) != length(target_groups))
    stop("source_groups and target_groups must have same length")

  one_condition <- function(cond) {
    cat("\n=== LIANA: condition =", cond, "===\n")
    cells_cond <- colnames(seurat_obj)[seurat_obj@meta.data[[condition_col]] == cond]
    if (length(cells_cond) < 10) { cat("  ⚠ Too few cells\n"); return(NULL) }
    seurat_sub <- seurat_obj[, cells_cond]
    Idents(seurat_sub) <- seurat_sub@meta.data[[cell_type_col]]
    DefaultAssay(seurat_sub) <- "SCT"

    wrap_args <- list(seurat_sub, resource="custom", external_resource=ortholog_resource)
    if (!is.null(liana_workers)) wrap_args <- c(wrap_args, list(parallelize=TRUE, workers=liana_workers))
    liana_res <- tryCatch(
      do.call(liana::liana_wrap, wrap_args) %>% liana::liana_aggregate(),
      error=function(e) { cat("  ✗ LIANA error:", e$message, "\n"); NULL })
    if (is.null(liana_res) || nrow(liana_res) == 0) return(NULL)

    liana_sig <- liana_res %>% dplyr::filter(aggregate_rank < agg_rank_cutoff)
    write.csv(liana_res, file.path(out_dir, paste0("liana_", cond, ".csv")), row.names=FALSE)
    write.csv(liana_sig, file.path(out_dir, paste0("liana_", cond, "_sig_agg", agg_rank_cutoff, ".csv")), row.names=FALSE)
    cat("  ✓", nrow(liana_res), "total,", nrow(liana_sig), "significant\n")

    for (pair_name in names(source_groups))
      plot_liana_dotplot(liana_sig, source_groups[[pair_name]], target_groups[[pair_name]],
                         condition = cond, pair_name = pair_name, ntop = ntop, out_dir = out_dir,
                         palette = palette, reverse_palette = reverse_palette)
    liana_res
  }

  bp  <- lab_bpparam(workers, length(conditions), BPPARAM, profile = "heavy")
  out <- lab_run_tasks(as.list(conditions), one_condition, bp)
  results_all <- list()
  for (i in seq_along(conditions)) {
    o <- out[[i]]
    if (length(o$log) > 0) cat(paste(o$log, collapse = "\n"), "\n")
    r <- o$result
    if (is.list(r) && !is.data.frame(r) && identical(r$status, "error")) { cat("  ✗ Error:", r$message, "\n"); next }
    if (!is.null(r)) results_all[[conditions[i]]] <- r
  }
  cat("\n✓✓✓ LIANA complete ✓✓✓\n"); return(invisible(results_all))
}

# Saves the dotplot as PDF + PNG (300 dpi): liana_<condition>_<pair>_<src>__TO__<tgt>_agg05.pdf/.png
# palette = NULL keeps liana's own colour scale; otherwise any lab palette (see
# viz_functions.R) replaces the colour scale of the aggregate-rank magnitude.
# Returns the ggplot invisibly (NULL if nothing to plot / on error).
plot_liana_dotplot <- function(liana_result, source_groups, target_groups,
                                condition="cond", pair_name="pair", ntop=20,
                                out_dir="LR_analysis", width=12, height=14,
                                palette=NULL, reverse_palette=FALSE) {
  check_required_packages("liana", "plot_liana_dotplot()")
  src_str <- paste(source_groups, collapse="_"); tgt_str <- paste(target_groups, collapse="_")
  fname   <- file.path(out_dir, paste0("liana_", condition, "_", pair_name,
                                        "_", src_str, "__TO__", tgt_str, "_agg05"))
  plot_data <- liana_result %>% dplyr::filter(source %in% source_groups | target %in% target_groups)
  if (nrow(plot_data) == 0) { cat("  ⚠ No interactions to plot\n"); return(invisible(NULL)) }
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
  tryCatch({
    p <- liana::liana_dotplot(plot_data, source_groups=source_groups, target_groups=target_groups,
           specificity="natmi.edge_specificity", magnitude="aggregate_rank",
           colour.label="Aggregate Rank\n(p-val)", ntop=ntop) +
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
    if (!is.null(palette))
      p <- suppressMessages(p + lab_scale_continuous(palette, "colour", reverse_palette,
                                                     name="Aggregate Rank\n(p-val)"))
    ggsave(paste0(fname, ".pdf"), p, width=width, height=height, limitsize=FALSE)
    ggsave(paste0(fname, ".png"), p, width=width, height=height, dpi=300, limitsize=FALSE)
    cat("  ✓ Dotplot →", paste0(fname, ".pdf/.png"), "\n")
    invisible(p)
  }, error=function(e) { cat("  ✗ Dotplot error:", e$message, "\n"); invisible(NULL) })
}


# =============================================================================
# 6. SC ENRICHMENT WRAPPERS (adds <CellType> subfolder)
# =============================================================================
# Thin adapters over run_ORA_loop() / run_GSEA_loop(): every cell type x threshold
# x ontology x gene set becomes one task of ONE parallel, load-balanced pool.
# Output layout is unchanged: <out_root>/ORA/<lfc>/<contrast>/<CellType>/GO/BP/...
#
# Gene sets follow run_ORA_loop(): for threshold t, `all` = padj < 0.05 & |LFC| > t,
# `up` / `down` = its positive / negative half. The default thresholds therefore
# include lfc0, whose `all` set is every padj < 0.05 gene (the previous meaning of
# `all`, now under the lfc0 folder) and lfc0.585 (= log2(1.5), the sc default).
#
# Extra arguments (`...`) go to run_ORA_loop() / run_GSEA_loop(): out_root, universe,
# background, set_labels, workers, BPPARAM, seed, top_n_ranked, pval_deg,
# palette, reverse_palette, make_barplot (ORA + GSEA), ranked_color (GSEA).

.sc_ct_safe <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

run_ORA_loop_sc <- function(de_results_sc, contrast,
                            thresholds = list(lfc0 = 0, lfc0.585 = 0.585),
                            org = "mouse", num_pathways = 15, ...) {
  cat("Starting single-cell ORA for contrast:", contrast, "\n")
  filtered_results <- list()
  for (cell_type in names(de_results_sc)) {
    res_df <- de_results_sc[[cell_type]]
    if (is.null(res_df) || nrow(res_df) == 0) { cat("  ⚠ No results for", cell_type, "\n"); next }
    key <- file.path(contrast, .sc_ct_safe(cell_type))
    filtered_results[[key]] <- lapply(thresholds, function(thresh)
      res_df %>% dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > thresh) %>%
        dplyr::select(gene, log2FoldChange))
  }
  if (length(filtered_results) == 0) { cat("⚠ Nothing to run\n"); return(invisible(NULL)) }
  summary_df <- run_ORA_loop(filtered_results, contrasts = names(filtered_results),
                             thresholds = thresholds, org = org, num_pathways = num_pathways, ...)
  cat("\n✓✓✓ Single-cell ORA complete ✓✓✓\n")
  invisible(summary_df)
}

run_GSEA_loop_sc <- function(de_results_sc, contrast,
                             thresholds = list(lfc0 = 0, lfc0.585 = 0.585),
                             organism = "mouse", num_pathways = 15, top_n_ranked = 3, ...) {
  cat("Starting single-cell GSEA for contrast:", contrast, "\n")
  results_list <- list()
  for (cell_type in names(de_results_sc)) {
    res_df <- de_results_sc[[cell_type]]
    if (is.null(res_df) || nrow(res_df) == 0) { cat("  ⚠ No results for", cell_type, "\n"); next }
    results_list[[file.path(contrast, .sc_ct_safe(cell_type))]] <- res_df
  }
  if (length(results_list) == 0) { cat("⚠ Nothing to run\n"); return(invisible(NULL)) }
  summary_df <- run_GSEA_loop(results_list, contrasts = names(results_list), thresholds = thresholds,
                              organism = organism, num_pathways = num_pathways,
                              top_n_ranked = top_n_ranked, ...)
  cat("\n✓✓✓ Single-cell GSEA complete ✓✓✓\n")
  invisible(summary_df)
}


cat("✓ sc_functions.R loaded\n")
cat("  Clustering QC  : compute_clustering_qc(), plot_clustering_qc()\n")
cat("  Annotation     : annotate_celltypes_enrichr()\n")
cat("  DE selection   : compute_data_characteristics(), recommend_de_method()\n")
cat("  DE dispatcher  : run_de_per_celltype()  [cell types in parallel: workers=]\n")
cat("  DE methods     : run_pseudobulk_limma(), run_limmatrend_cov(),\n")
cat("                   run_limmavoom_cov(), run_edgeR_cov(),\n")
cat("                   run_mast_cov(), run_wilcox_raw()\n")
cat("  LIANA          : run_liana_zebrafish(), run_liana_loop(), plot_liana_dotplot()\n")
cat("  SC enrichment  : run_ORA_loop_sc(), run_GSEA_loop_sc()  [parallel via run_*_loop]\n")
