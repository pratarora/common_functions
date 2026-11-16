# =============================================================================
# Mercader Lab — Orthology Conversion Functions
# orthology_functions.R
#
# Functions for converting gene identifiers across species using the Alliance
# Genome orthology table (bidirectional best-score filter).
#
# Public:
#   - convert_zf_to_mouse_alliance()  : zebrafish -> mouse symbols
#   - convert_zf_to_human_alliance()  : zebrafish -> human symbols
#   - convert_to_zebrafish_alliance() : mouse / human -> zebrafish symbols
#   - alliance_cache_clear()          : drop the session / disk cache
#
# Performance: the Alliance TSV (~170 MB) is parsed ONCE per session with
# data.table::fread (only the 6 needed columns, multi-threaded), reduced to the
# high-quality Danio/Mus/Homo pairs and cached in memory. A small .rds copy of
# that reduced table is written next to the TSV (<file>.hq.rds; or under
# tools::R_user_dir("mercader", "cache") if that folder is read-only) so a new
# kernel loads it in well under a second. Both caches are invalidated when the
# TSV's modification time or size changes.
# =============================================================================

if (!exists("lab_run_tasks", mode = "function")) local({
  d <- "reusable_functions"
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (is.character(f) && grepl("orthology_functions", basename(f))) { d <- dirname(f); break }
  }
  source(file.path(d, "parallel_helpers.R"))
})

.ALLIANCE_COLS    <- c("Gene1Symbol", "Gene2Symbol", "Gene1SpeciesName", "Gene2SpeciesName",
                       "IsBestScore", "IsBestRevScore")
.ALLIANCE_SPECIES <- c("Danio rerio", "Mus musculus", "Homo sapiens")
.alliance_cache   <- new.env(parent = emptyenv())

# Cache key: path + mtime + size, so a re-downloaded file invalidates the cache.
.alliance_key <- function(path) {
    info <- file.info(path)
    paste(normalizePath(path, mustWork = FALSE), as.numeric(info$mtime), info$size, sep = "|")
}

.alliance_disk_cache_path <- function(path) {
    beside <- paste0(path, ".hq.rds")
    dir_ok <- file.access(dirname(path), 2) == 0
    if (dir_ok) return(beside)
    cache_dir <- tryCatch(tools::R_user_dir("mercader", "cache"), error = function(e) NULL)
    if (is.null(cache_dir)) return(NULL)
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    file.path(cache_dir, paste0(basename(path), ".", substr(.alliance_path_hash(path), 1, 8), ".hq.rds"))
}

.alliance_path_hash <- function(path) {
    # short stable hash of the absolute path without the digest package dependency
    x <- utf8ToInt(normalizePath(path, mustWork = FALSE))
    sprintf("%08x", sum(x * seq_along(x)) %% 4294967295)
}

# Parse the TSV (fread, 6 columns) and reduce it to the high-quality pairs of
# the three species of interest. Falls back to the original read.delim path if
# fread / the header is not found.
.alliance_parse <- function(alliance_file) {
    alliance_data <- NULL
    if (requireNamespace("data.table", quietly = TRUE)) {
        alliance_data <- tryCatch(
            data.table::setDF(data.table::fread(
                alliance_file, skip = "Gene1ID", sep = "\t", quote = "", header = TRUE,
                select = .ALLIANCE_COLS, showProgress = FALSE,
                nThread = lab_workers(Inf, "light"))),
            error = function(e) NULL)
    }
    if (is.null(alliance_data) || !all(.ALLIANCE_COLS %in% colnames(alliance_data))) {
        # original reader (slow path)
        alliance_data <- read.delim(alliance_file, header = TRUE, sep = "\t",
                                    stringsAsFactors = FALSE, quote = "", comment.char = "#")
        if (ncol(alliance_data) == 1) {
            alliance_data <- read.table(alliance_file, header = TRUE, sep = "\t",
                                        stringsAsFactors = FALSE, quote = "", fill = TRUE,
                                        comment.char = "")
        }
        if (any(!.ALLIANCE_COLS %in% colnames(alliance_data))) {
            cat("ERROR: Missing required columns in Alliance file\n")
            cat("Expected:", paste(.ALLIANCE_COLS, collapse = ", "), "\n")
            cat("Found:", paste(colnames(alliance_data), collapse = ", "), "\n")
            return(NULL)
        }
        alliance_data <- alliance_data[, .ALLIANCE_COLS]
    }
    # High-quality orthologs (bidirectional best scores) — exactly "Yes", as before
    hq <- alliance_data[alliance_data$IsBestScore == "Yes" & alliance_data$IsBestRevScore == "Yes", ]
    n_hq <- nrow(hq)
    hq <- hq[hq$Gene1SpeciesName %in% .ALLIANCE_SPECIES & hq$Gene2SpeciesName %in% .ALLIANCE_SPECIES,
             c("Gene1Symbol", "Gene2Symbol", "Gene1SpeciesName", "Gene2SpeciesName")]
    rownames(hq) <- NULL
    list(pairs = hq, n_high_quality = n_hq)
}

# Returns list(pairs = data.frame(Gene1Symbol, Gene2Symbol, Gene1SpeciesName,
# Gene2SpeciesName), n_high_quality = <count before the species filter>) or NULL.
.read_alliance <- function(alliance_file, use_cache = TRUE) {
    key <- .alliance_key(alliance_file)
    if (use_cache && !is.null(.alliance_cache[[key]])) return(.alliance_cache[[key]])

    disk <- if (use_cache) .alliance_disk_cache_path(alliance_file) else NULL
    if (!is.null(disk) && file.exists(disk)) {
        obj <- tryCatch(readRDS(disk), error = function(e) NULL)
        if (is.list(obj) && identical(obj$key, key)) {
            cat("Alliance orthology table loaded from cache:", basename(disk), "\n")
            .alliance_cache[[key]] <- obj$data
            return(obj$data)
        }
    }

    cat("Reading Alliance orthology file...\n")
    parsed <- .alliance_parse(alliance_file)
    if (is.null(parsed)) return(NULL)
    if (use_cache) {
        .alliance_cache[[key]] <- parsed
        if (!is.null(disk))
            tryCatch(saveRDS(list(key = key, data = parsed), disk),
                     error = function(e) NULL)   # read-only location: just skip the disk cache
    }
    parsed
}

# Clears the in-memory cache and (optionally) deletes the on-disk cache of one file.
alliance_cache_clear <- function(alliance_file = NULL, disk = TRUE) {
    rm(list = ls(.alliance_cache, all.names = TRUE), envir = .alliance_cache)
    if (disk && !is.null(alliance_file)) {
        for (p in unique(c(paste0(alliance_file, ".hq.rds"),
                           tryCatch(.alliance_disk_cache_path(alliance_file), error = function(e) NULL))))
            if (!is.null(p) && file.exists(p)) unlink(p)
    }
    invisible(NULL)
}

# Generic pair extractor. Builds the mapping <from species symbol> -> <to species
# symbol> from both orientations of the table, cleans/deduplicates, then filters
# for `genes` (optionally case-insensitive on the `from` side, reporting the
# symbol as spelled in the input).
.alliance_pairs <- function(pairs, from_sp, to_sp, genes, ignore_case, from_col, to_col) {
    fwd <- pairs[pairs$Gene1SpeciesName == from_sp & pairs$Gene2SpeciesName == to_sp, ]
    rev <- pairs[pairs$Gene1SpeciesName == to_sp   & pairs$Gene2SpeciesName == from_sp, ]
    mapping <- data.frame(a = c(fwd$Gene1Symbol, rev$Gene2Symbol),
                          b = c(fwd$Gene2Symbol, rev$Gene1Symbol),
                          stringsAsFactors = FALSE)
    mapping <- mapping[!is.na(mapping$a) & !is.na(mapping$b) & mapping$a != "" & mapping$b != "", ]
    mapping <- mapping[!duplicated(mapping), ]
    n_pairs <- nrow(mapping)

    if (ignore_case) {
        up_genes <- toupper(genes)
        hit <- match(toupper(mapping$a), up_genes)
        result <- mapping[!is.na(hit), ]
        result$a <- genes[hit[!is.na(hit)]]   # spelled as in the input so it joins back
        result <- result[!duplicated(result), ]
    } else {
        result <- mapping[mapping$a %in% genes, ]
    }
    rownames(result) <- NULL
    names(result) <- c(from_col, to_col)
    list(result = result, n_pairs = n_pairs)
}

.alliance_missing_msg <- function(alliance_file) {
    cat("ERROR: Alliance orthology file not found:", alliance_file, "\n")
    cat("Please download from Alliance Genome (https://www.alliancegenome.org/downloads)\n")
}

#' Convert Zebrafish Genes to Mouse Orthologs Using Alliance Database
#'
#' Converts zebrafish gene symbols to mouse orthologs using high-quality
#' Alliance Genome orthology data with bidirectional best score filtering.
#'
#' @param zf_genes Character vector of zebrafish gene symbols
#' @param alliance_file Path to Alliance orthology TSV file (.tsv or .tsv.gz)
#'        Default: "ORTHOLOGY-ALLIANCE_COMBINED_7_23052025.tsv"
#' @param ignore_case Logical (default FALSE). If TRUE, zebrafish symbols are
#'        matched case-insensitively (e.g. a GTF that writes "CDKAL1" or "Foxo1"
#'        matches Alliance "cdkal1" / "foxo1"). Only the zebrafish side is
#'        affected; mouse symbols are returned exactly as in the Alliance file
#'        so they still match org.Mm.eg.db. Use for newer Ensembl annotations
#'        (e.g. GRCz12ab) whose gene names are not all lower case.
#' @param use_cache Logical (default TRUE). Reuse the parsed table cached in the
#'        session / on disk (<alliance_file>.hq.rds). FALSE re-parses the TSV.
#'
#' @return Data frame with columns:
#'   - zebrafish_gene_id: Input zebrafish gene symbol (with ignore_case = TRUE,
#'                        spelled as in the input, so it joins back to your table)
#'   - mouse_gene_symbol: Corresponding mouse gene symbol
#'   Returns NULL if file not found or parsing fails.
#'
#' @details
#' Uses Alliance Genome orthology data filtered for high confidence:
#'   - IsBestScore = "Yes" (best zebrafish→mouse match)
#'   - IsBestRevScore = "Yes" (best mouse→zebrafish match)
#'
#' Bidirectional mapping ensures orthologs are reciprocal best hits.
#'
#' Download Alliance file from:
#' https://www.alliancegenome.org/downloads
#'
#' @examples
#' zf_genes <- c("tp53", "myca", "nkx2.5")
#' orthologs <- convert_zf_to_mouse_alliance(zf_genes)
#'
#' @export
convert_zf_to_mouse_alliance <- function(zf_genes,
                                         alliance_file = "ORTHOLOGY-ALLIANCE_COMBINED_7_23052025.tsv",
                                         ignore_case = FALSE, use_cache = TRUE) {
    if (!file.exists(alliance_file)) { .alliance_missing_msg(alliance_file); return(NULL) }
    tryCatch({
        al <- .read_alliance(alliance_file, use_cache = use_cache)
        if (is.null(al)) return(NULL)
        cat("High-quality ortholog pairs (bidirectional best scores):", al$n_high_quality, "\n")
        m <- .alliance_pairs(al$pairs, "Danio rerio", "Mus musculus", zf_genes, ignore_case,
                             "zebrafish_gene_id", "mouse_gene_symbol")
        cat("Total unique ortholog pairs:", m$n_pairs, "\n")
        cat("Mapped", length(unique(m$result$zebrafish_gene_id)), "of",
            length(unique(zf_genes)), "input genes to mouse orthologs",
            if (ignore_case) "(case-insensitive)" else "", "\n")
        m$result
    }, error = function(e) {
        cat("ERROR reading Alliance file:", e$message, "\n")
        NULL
    })
}

#' Convert Zebrafish Genes to Human Orthologs Using Alliance Database
#'
#' Converts zebrafish gene symbols to human orthologs using high-quality
#' Alliance Genome orthology data with bidirectional best score filtering.
#' Mirror of convert_zf_to_mouse_alliance() for the Danio rerio <-> Homo sapiens
#' direction; used by annotate_celltypes_enrichr() for Enrichr cell-type calls.
#'
#' @param zf_genes Character vector of zebrafish gene symbols
#' @param alliance_file Path to Alliance orthology TSV file
#'        Default: "ORTHOLOGY-ALLIANCE_COMBINED_7_23052025.tsv"
#' @param ignore_case Logical (default FALSE). If TRUE, zebrafish symbols are
#'        matched case-insensitively (see convert_zf_to_mouse_alliance()); human
#'        symbols are returned exactly as in the Alliance file.
#' @param use_cache Logical (default TRUE). See convert_zf_to_mouse_alliance().
#'
#' @return Data frame with columns:
#'   - zebrafish_gene_id: Input zebrafish gene symbol (with ignore_case = TRUE,
#'                        spelled as in the input)
#'   - human_gene_symbol: Corresponding human gene symbol
#'   Returns NULL if file not found or parsing fails.
#'
#' @export
convert_zf_to_human_alliance <- function(zf_genes,
                                         alliance_file = "ORTHOLOGY-ALLIANCE_COMBINED_7_23052025.tsv",
                                         ignore_case = FALSE, use_cache = TRUE) {
    if (!file.exists(alliance_file)) { .alliance_missing_msg(alliance_file); return(NULL) }
    tryCatch({
        al <- .read_alliance(alliance_file, use_cache = use_cache)
        if (is.null(al)) return(NULL)
        cat("High-quality ortholog pairs (bidirectional best scores):", al$n_high_quality, "\n")
        m <- .alliance_pairs(al$pairs, "Danio rerio", "Homo sapiens", zf_genes, ignore_case,
                             "zebrafish_gene_id", "human_gene_symbol")
        cat("Total unique ortholog pairs:", m$n_pairs, "\n")
        cat("Mapped", length(unique(m$result$zebrafish_gene_id)), "of",
            length(unique(zf_genes)), "input genes to human orthologs",
            if (ignore_case) "(case-insensitive)" else "", "\n")
        m$result
    }, error = function(e) {
        cat("ERROR reading Alliance file:", e$message, "\n")
        NULL
    })
}

#' Convert Mouse / Human Gene Symbols to Zebrafish Orthologs Using Alliance Database
#'
#' Reverse of convert_zf_to_mouse_alliance() / convert_zf_to_human_alliance():
#' maps mouse and/or human symbols (e.g. the TF names in MEME motif databases)
#' to zebrafish symbols with the same high-confidence Alliance filter
#' (IsBestScore = "Yes" and IsBestRevScore = "Yes", both directions combined).
#'
#' @param genes Character vector of mouse and/or human gene symbols
#' @param alliance_file Path to Alliance orthology TSV file
#'        Default: "ORTHOLOGY-ALLIANCE_COMBINED_7_23052025.tsv"
#' @param from Which source species to use: "mouse", "human" or both (default).
#'        A symbol present in both species maps through each.
#' @param ignore_case Logical (default TRUE). Match symbols case-insensitively
#'        (motif databases mix "SP1", "Zfp341", "Arnt"). Zebrafish symbols are
#'        returned exactly as in the Alliance file.
#' @param use_cache Logical (default TRUE). See convert_zf_to_mouse_alliance().
#'
#' @return Data frame with columns:
#'   - input_gene: the symbol as spelled in `genes`
#'   - zebrafish_gene_symbol: Corresponding zebrafish gene symbol
#'   - source_species: "mouse" or "human" (the species the match came through)
#'   Returns NULL if file not found or parsing fails.
#'
#' @export
convert_to_zebrafish_alliance <- function(genes,
                                          alliance_file = "ORTHOLOGY-ALLIANCE_COMBINED_7_23052025.tsv",
                                          from = c("mouse", "human"),
                                          ignore_case = TRUE, use_cache = TRUE) {
    from <- match.arg(from, several.ok = TRUE)
    if (!file.exists(alliance_file)) { .alliance_missing_msg(alliance_file); return(NULL) }
    tryCatch({
        al <- .read_alliance(alliance_file, use_cache = use_cache)
        if (is.null(al)) return(NULL)
        species_name <- c(mouse = "Mus musculus", human = "Homo sapiens")
        genes <- unique(genes[!is.na(genes) & genes != ""])
        mapping <- do.call(rbind, lapply(from, function(src) {
            sp <- species_name[[src]]
            fwd <- al$pairs[al$pairs$Gene1SpeciesName == "Danio rerio" & al$pairs$Gene2SpeciesName == sp, ]
            rev <- al$pairs[al$pairs$Gene1SpeciesName == sp & al$pairs$Gene2SpeciesName == "Danio rerio", ]
            data.frame(input_gene = c(fwd$Gene2Symbol, rev$Gene1Symbol),
                       zebrafish_gene_symbol = c(fwd$Gene1Symbol, rev$Gene2Symbol),
                       source_species = src, stringsAsFactors = FALSE)
        }))
        mapping <- mapping[!is.na(mapping$input_gene) & !is.na(mapping$zebrafish_gene_symbol) &
                           mapping$input_gene != "" & mapping$zebrafish_gene_symbol != "", ]
        mapping <- mapping[!duplicated(mapping), ]
        if (ignore_case) {
            hit <- match(toupper(mapping$input_gene), toupper(genes))
            result <- mapping[!is.na(hit), ]
            result$input_gene <- genes[hit[!is.na(hit)]]
            result <- result[!duplicated(result), ]
        } else {
            result <- mapping[mapping$input_gene %in% genes, ]
        }
        rownames(result) <- NULL
        cat("Mapped", length(unique(result$input_gene)), "of", length(genes),
            "input genes to zebrafish orthologs", if (ignore_case) "(case-insensitive)" else "", "\n")
        result
    }, error = function(e) {
        cat("ERROR reading Alliance file:", e$message, "\n")
        NULL
    })
}
