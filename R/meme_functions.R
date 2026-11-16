# =============================================================================
# Mercader Lab — MEME Suite / motif analysis functions
# meme_functions.R
#
# Functions:
#   - export_peaks_bed()      : region table (1-based) -> 0-based BED
#   - regions_to_fasta()      : region table -> FASTA from a BSgenome
#   - random_regions_fasta()  : seeded random control regions (same width) -> FASTA
#   - run_meme_suite()        : run xstreme / sea / centrimo / meme-chip on one FASTA
#   - run_meme_suite_batch()  : launch many runs in the background, N at a time
#   - meme_wait()             : block until background runs have finished
#   - meme_run_status()       : poll RUN_STATUS files of background runs
#   - parse_xstreme()         : xstreme.html -> tidy motif table (CSV [+ xlsx])
#   - map_tfs_to_zebrafish()  : motif TF symbols (mouse/human) -> zebrafish orthologs (Alliance)
#   - make_motif_logos()      : ceqlogo logos (PDF + PNG) for the motifs of a .meme file
#
# Shell side: reusable_functions/run_meme_suite.sh (MEME 5.5.9 from the conda env
# `archr_tools`; override with meme_env =). No project paths are hard-coded here.
#
# Usage: source("reusable_functions/meme_functions.R")   (loaded by source_all.R)
# =============================================================================

library(dplyr)
library(stringr)

if (!exists("lab_run_tasks", mode = "function")) local({
  d <- "reusable_functions"
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (is.character(f) && grepl("meme_functions", basename(f))) { d <- dirname(f); break }
  }
  source(file.path(d, "parallel_helpers.R"))
})


# -----------------------------------------------------------------------------
# export_peaks_bed()
#
# Writes a 0-based BED from a region table with 1-based starts (ArchR, DiffBind,
# ChIPseeker and csaw tables all use 1-based starts).
#
# Args:
#   regions : data.frame with columns seqnames, start (1-based), end
#   bed_path: output path (parent folder is created)
#   chroms  : optional vector of chromosome names to keep (e.g. as.character(1:25))
#   name    : prefix for the BED name column (default "region")
#
# Returns: the filtered region data.frame (rows written), invisibly
# -----------------------------------------------------------------------------
export_peaks_bed <- function(regions, bed_path, chroms = NULL, name = "region") {
  stopifnot(all(c("seqnames", "start", "end") %in% names(regions)))
  regions$seqnames <- as.character(regions$seqnames)
  if (!is.null(chroms)) regions <- regions[regions$seqnames %in% chroms, ]
  regions <- regions[regions$end > regions$start, ]
  dir.create(dirname(bed_path), showWarnings = FALSE, recursive = TRUE)
  bed <- data.frame(chrom = regions$seqnames, start = regions$start - 1L, end = regions$end,
                    name = paste0(name, "_", seq_len(nrow(regions))), score = 0L, strand = ".")
  write.table(bed, bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat("✓ BED:", bed_path, "-", nrow(bed), "regions\n")
  invisible(regions)
}


# -----------------------------------------------------------------------------
# regions_to_fasta()
#
# Extracts region sequences from a BSgenome (no samtools / .gzi index needed).
# FASTA headers are "<chrom>:<start>-<end>" (1-based, inclusive). Regions that
# are mostly N are dropped.
#
# Args:
#   regions    : data.frame with seqnames, start (1-based), end
#   genome     : BSgenome object, e.g. BSgenome.Drerio.Ensembl.GRCz12ab
#   fasta_path : output FASTA (parent folder is created)
#   chroms     : optional chromosomes to keep
#   max_n_frac : drop regions with more than this fraction of N (default 0.1)
#
# Returns: number of sequences written, invisibly
# -----------------------------------------------------------------------------
regions_to_fasta <- function(regions, genome, fasta_path, chroms = NULL, max_n_frac = 0.1) {
  regions$seqnames <- as.character(regions$seqnames)
  if (!is.null(chroms)) regions <- regions[regions$seqnames %in% chroms, ]
  regions <- regions[regions$end > regions$start, ]
  gr   <- GenomicRanges::GRanges(regions$seqnames, IRanges::IRanges(regions$start, regions$end))
  seqs <- Biostrings::getSeq(genome, gr)
  names(seqs) <- paste0(regions$seqnames, ":", regions$start, "-", regions$end)
  keep <- Biostrings::letterFrequency(seqs, "N", as.prob = TRUE)[, 1] <= max_n_frac
  seqs <- seqs[keep]
  dir.create(dirname(fasta_path), showWarnings = FALSE, recursive = TRUE)
  Biostrings::writeXStringSet(seqs, fasta_path)
  cat("✓ FASTA:", fasta_path, "-", length(seqs), "sequences",
      if (any(!keep)) paste0("(", sum(!keep), " dropped for N content)") else "", "\n")
  invisible(length(seqs))
}


# -----------------------------------------------------------------------------
# random_regions_fasta()
#
# Seeded random control regions of a fixed width, for xstreme --n. Chromosomes
# are sampled in proportion to their length. Regions overlapping `exclude` (the
# primary regions) or mostly N are dropped. A BED of the regions is written next
# to the FASTA so the control set is traceable.
#
# Args:
#   n          : number of control sequences wanted
#   width      : width of every control region (match the primary regions)
#   genome     : BSgenome object
#   fasta_path : output FASTA (a .bed is written beside it)
#   seed       : random seed (set here, so the same seed gives the same regions)
#   chroms     : chromosomes to sample from (e.g. as.character(1:25)); NULL = all
#   exclude    : optional GRanges of regions to avoid
#   max_n_frac : drop regions with more than this fraction of N (default 0.1)
#
# Returns: GRanges of the control regions, invisibly
# -----------------------------------------------------------------------------
random_regions_fasta <- function(n, width, genome, fasta_path, seed, chroms = NULL,
                                 exclude = NULL, max_n_frac = 0.1) {
  sl <- GenomeInfoDb::seqlengths(genome)
  if (!is.null(chroms)) sl <- sl[names(sl) %in% chroms]
  sl <- sl[sl > width]
  set.seed(seed)
  m     <- n * 4L                                       # oversample, then trim
  chr   <- sample(names(sl), m, replace = TRUE, prob = sl - width + 1)
  start <- floor(stats::runif(m) * (sl[chr] - width + 1)) + 1
  gr    <- GenomicRanges::GRanges(chr, IRanges::IRanges(start, width = width))
  if (!is.null(exclude)) gr <- gr[!IRanges::overlapsAny(gr, exclude)]
  seqs  <- Biostrings::getSeq(genome, gr)
  ok    <- Biostrings::letterFrequency(seqs, "N", as.prob = TRUE)[, 1] <= max_n_frac
  sel   <- which(ok)[seq_len(min(n, sum(ok)))]
  gr    <- gr[sel]
  if (length(gr) < n) warning("random_regions_fasta(): only ", length(gr), " of ", n, " regions kept")
  seqs  <- seqs[sel]                                   # no second getSeq()
  names(seqs) <- paste0(as.character(GenomicRanges::seqnames(gr)), ":",
                        GenomicRanges::start(gr), "-", GenomicRanges::end(gr))
  dir.create(dirname(fasta_path), showWarnings = FALSE, recursive = TRUE)
  Biostrings::writeXStringSet(seqs, fasta_path)
  bed_path <- sub("\\.(fa|fasta)$", ".bed", fasta_path)
  write.table(data.frame(as.character(GenomicRanges::seqnames(gr)), GenomicRanges::start(gr) - 1L, GenomicRanges::end(gr)),
              bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat("✓ Random control FASTA:", fasta_path, "-", length(gr), "x", width, "bp, seed", seed, "\n")
  invisible(gr)
}


# -----------------------------------------------------------------------------
# run_meme_suite()
#
# Runs MEME Suite tools on one FASTA through run_meme_suite.sh.
#
# Args:
#   fasta       : primary sequences
#   out_dir     : output folder (created); one sub-folder per tool + RUN_STATUS + log
#   motif_dbs   : character vector of .meme motif databases
#   neg_fasta   : optional control FASTA (xstreme / sea --n)
#   tools       : any of "xstreme", "sea", "centrimo", "meme-chip" (default "xstreme")
#   threads     : OpenMP threads for xstreme / meme-chip (MEME step)
#   seed        : random seed passed to the tools (default 0 = MEME's own default)
#   wait        : TRUE blocks until done; FALSE starts a detached process and returns
#                 at once — poll with meme_run_status() (preferred for long jobs)
#   script      : path to run_meme_suite.sh
#   meme_env    : conda env with MEME (default: the script's /opt/conda/envs/archr_tools)
#   xstreme_args: replace the default xstreme options (lab settings) if given
#
# Returns: for wait = TRUE the exit status (0 = OK), otherwise invisible(NULL)
# -----------------------------------------------------------------------------
run_meme_suite <- function(fasta, out_dir, motif_dbs, neg_fasta = NULL, tools = "xstreme",
                           threads = 1, seed = 0, wait = TRUE,
                           script = "reusable_functions/run_meme_suite.sh",
                           meme_env = NULL, xstreme_args = NULL) {
  stopifnot(file.exists(script), file.exists(fasta), all(file.exists(motif_dbs)))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  env <- c(sprintf("FASTA=%s", normalizePath(fasta)),
           sprintf("OUT_DIR=%s", normalizePath(out_dir)),
           sprintf("MOTIF_DBS=%s", paste(normalizePath(motif_dbs), collapse = ":")),
           sprintf("TOOLS=%s", paste(tools, collapse = ",")),
           sprintf("MEME_THREADS=%d", as.integer(threads)),
           sprintf("SEED=%d", as.integer(seed)))
  if (!is.null(neg_fasta))    env <- c(env, sprintf("NEG_FASTA=%s", normalizePath(neg_fasta)))
  if (!is.null(meme_env))     env <- c(env, sprintf("MEME_ENV=%s", meme_env))
  if (!is.null(xstreme_args)) env <- c(env, sprintf("XSTREME_ARGS=%s", xstreme_args))
  unlink(file.path(out_dir, "RUN_STATUS"))
  if (wait) {   # output goes to <out_dir>/meme_suite_run.log; the exit status is returned
    return(invisible(system2("bash", script, env = env, stdout = FALSE, stderr = FALSE)))
  }
  system2("setsid", c("nohup", "bash", script), env = env, wait = FALSE, stdout = FALSE, stderr = FALSE)
  cat("▶ started (background):", out_dir, "\n")
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# run_meme_suite_batch()
#
# Launches several run_meme_suite() runs as detached background processes,
# at most `max_concurrent` at a time (a semaphore on the RUN_STATUS files), so
# many region sets / control sets use the machine without oversubscribing it.
#
# Args:
#   runs            : named list; each element is a list of run_meme_suite()
#                     arguments (fasta, out_dir, motif_dbs, neg_fasta, tools, seed, ...)
#   max_concurrent  : runs at a time (default: light workers %/% threads_per_run)
#   threads_per_run : OpenMP threads given to every run (`threads`), default 4
#   wait            : TRUE (default) blocks until every run has finished and returns
#                     the statuses; FALSE returns as soon as the last run is launched
#   poll_sec        : seconds between checks while waiting (default 30)
#
# Returns: named character vector of statuses (as meme_run_status()), invisibly
#          when wait = FALSE
# -----------------------------------------------------------------------------
run_meme_suite_batch <- function(runs, max_concurrent = NULL, threads_per_run = 4, wait = TRUE,
                                 poll_sec = 30) {
  stopifnot(is.list(runs), length(runs) > 0)
  if (is.null(names(runs))) names(runs) <- paste0("run", seq_along(runs))
  if (is.null(max_concurrent))
    max_concurrent <- max(1L, lab_workers(Inf, "light") %/% max(1L, as.integer(threads_per_run)))
  out_dirs <- vapply(runs, function(r) r$out_dir, character(1))
  cat("▶ MEME batch:", length(runs), "runs,", max_concurrent, "at a time,", threads_per_run, "threads each\n")
  running <- character(0)
  for (nm in names(runs)) {
    while (length(running) >= max_concurrent) {
      Sys.sleep(poll_sec)
      running <- running[meme_run_status(running) == "RUNNING / NOT STARTED"]
    }
    args <- runs[[nm]]; args$wait <- FALSE; args$threads <- threads_per_run
    do.call(run_meme_suite, args)
    running <- c(running, args$out_dir)
  }
  if (!wait) return(invisible(setNames(meme_run_status(out_dirs), names(runs))))
  meme_wait(setNames(out_dirs, names(runs)), poll_sec = poll_sec)
}


# -----------------------------------------------------------------------------
# meme_wait()
#
# One blocking call that returns when every run has written its RUN_STATUS (or
# after `timeout_min`). Use it instead of polling meme_run_status() by hand.
#
# Args:  out_dirs : run output folders (named vector/list ok)
#        poll_sec : seconds between checks (default 30)
#        timeout_min : give up after this many minutes (default Inf)
# Returns: named character vector of statuses ("OK", "ERROR <tool>", or
#          "RUNNING / NOT STARTED" for runs still going at the timeout)
# -----------------------------------------------------------------------------
meme_wait <- function(out_dirs, poll_sec = 30, timeout_min = Inf) {
  t0 <- Sys.time()
  repeat {
    st <- meme_run_status(out_dirs)
    if (all(st != "RUNNING / NOT STARTED")) break
    if (as.numeric(difftime(Sys.time(), t0, units = "mins")) > timeout_min) {
      warning("meme_wait(): timeout after ", timeout_min, " min; ", sum(st == "RUNNING / NOT STARTED"), " run(s) still running")
      break
    }
    Sys.sleep(poll_sec)
  }
  cat("MEME runs:", paste(names(st) %||% seq_along(st), st, sep = " = ", collapse = ", "), "\n")
  st
}
`%||%` <- function(a, b) if (is.null(a)) b else a


# -----------------------------------------------------------------------------
# meme_run_status()
#
# Args:  out_dirs : run output folders (named vector/list ok)
# Returns: character vector, one of "OK", "ERROR <tool>", "RUNNING / NOT STARTED"
# -----------------------------------------------------------------------------
meme_run_status <- function(out_dirs) {
  vapply(out_dirs, function(d) {
    f <- file.path(d, "RUN_STATUS")
    if (file.exists(f)) trimws(readLines(f, n = 1, warn = FALSE)) else "RUNNING / NOT STARTED"
  }, character(1))
}


# --- XSTREME parsing helpers (internal) --------------------------------------
.get_all_joined <- function(match_col, field) {
  sapply(match_col, function(x) {
    if (is.null(x) || length(x) == 0L || !is.data.frame(x) || nrow(x) == 0L) return(NA_character_)
    paste(x[[field]], collapse = "; ")
  })
}

.get_formatted_matches <- function(match_col) {
  sapply(match_col, function(x) {
    if (is.null(x) || length(x) == 0L || !is.data.frame(x) || nrow(x) == 0L) return(NA_character_)
    paste(paste0(x$id, " (", x$alt, ")"), collapse = "; ")
  })
}

# Pull gene names out of motif-match strings, dropping species and DBD-score tokens:
#   "(Zfp341)_(Homo_sapiens)_(DBD_0.96)" -> "Zfp341";  "STREME-1 Patz1" -> "Patz1"
.clean_gene_names <- function(text_vector) {
  purrr::map_chr(text_vector, function(raw_str) {
    if (is.na(raw_str) || raw_str == "") return(NA_character_)
    entries   <- unlist(strsplit(raw_str, ";"))
    gene_list <- character(0)
    for (entry in entries) {
      raw_m   <- stringr::str_extract_all(entry, "\\(([^)]+)\\)")[[1]]
      cleaned <- gsub("^\\(|\\)$", "", raw_m)
      if (length(cleaned) == 0L && !grepl("^M[A-Z0-9_.]+$", trimws(entry))) cleaned <- trimws(entry)
      for (item in cleaned) {
        is_sp    <- grepl("Homo_sapiens|Mus_musculus|Danio_rerio|Drosophila|Rattus|Xenopus|Caenorhabditis|Saccharomyces|Arabidopsis|Gallus|Bos_taurus|Sus_scrofa|Oryza|Ciona|Anopheles|Aedes",
                          item, ignore.case = TRUE)
        is_score <- grepl("^DBD_[0-9.]+$", item)
        if (!is_sp && !is_score && nzchar(item)) gene_list <- c(gene_list, item)
      }
    }
    final <- unique(gene_list)
    if (length(final) > 0L) paste(final, collapse = ", ") else NA_character_
  })
}

.xstreme_motifs_from_html <- function(html_path) {
  if (!file.exists(html_path)) { message("  [SKIP] Not found: ", html_path); return(NULL) }
  lines     <- readLines(html_path, warn = FALSE)
  start_idx <- grep("var data\\s*=\\s*\\{", lines)[1]
  if (is.na(start_idx)) { message("  [SKIP] No JSON in: ", basename(html_path)); return(NULL) }
  end_offset <- grep("^\\s*\\};", lines[start_idx:length(lines)])[1]
  if (is.na(end_offset)) end_offset <- grep("</script>", lines[start_idx:length(lines)])[1]
  if (is.na(end_offset)) { message("  [SKIP] JSON end not found: ", basename(html_path)); return(NULL) }
  json_lines <- lines[start_idx:(start_idx + end_offset - 1L)]
  json_lines[1] <- sub("var data\\s*=\\s*", "", json_lines[1])
  json_lines[length(json_lines)] <- sub(";\\s*$", "", json_lines[length(json_lines)])
  data <- tryCatch(jsonlite::fromJSON(paste(json_lines, collapse = "\n")),
                   error = function(e) { message("  [ERROR] JSON parse: ", e$message); NULL })
  if (is.null(data) || is.null(data$motifs)) { message("  [SKIP] No motifs object in: ", basename(html_path)); return(NULL) }
  data$motifs
}


# -----------------------------------------------------------------------------
# parse_xstreme()
#
# Reads <xstreme_dir>/xstreme.html (embedded JSON), extracts the motif table,
# cleans TF names and writes <xstreme_dir>/xstreme_summary.csv (+ .xlsx if
# writexl is installed).
#
# Columns: Rank, analysis_algorithm (STREME / MEME / SEA = known motif),
#   Extracted_Gene_Names (TF symbols as in the motif database, human/mouse/
#   zebrafish, comma separated), Motif_ID, Alt_ID, Consensus,
#   E_value (= the "E-value" column of xstreme.html = JSON pgm_evalue: MEME / STREME
#   E-value for de novo motifs, SEA E-value for known motifs; use this for significance,
#   E_value < 0.05; MEME E-values are not comparable with STREME / SEA ones),
#   Sites, SEA_P_value / SEA_E_value / SEA_Q_value (SEA enrichment vs the control or
#   shuffled sequences; NOT accurate for de novo motifs, which were found in the same
#   sequences; XSTREME ranks by SEA_P_value), Similar_Known_Motifs, Match_IDs, Match_Evalues
#
# Returns: data.frame (NULL if nothing could be parsed)
# -----------------------------------------------------------------------------
parse_xstreme <- function(xstreme_dir, write = TRUE) {
  motifs_df <- .xstreme_motifs_from_html(file.path(xstreme_dir, "xstreme.html"))
  if (is.null(motifs_df)) return(NULL)

  if ("tomtom_matches" %in% names(motifs_df)) {
    motifs_df$Similar_Known_Motifs <- .get_formatted_matches(motifs_df$tomtom_matches)
    motifs_df$Match_Names_Raw      <- .get_all_joined(motifs_df$tomtom_matches, "alt")
    motifs_df$Match_IDs            <- .get_all_joined(motifs_df$tomtom_matches, "id")
    motifs_df$Match_Evalues        <- .get_all_joined(motifs_df$tomtom_matches, "evalue")
  } else {
    motifs_df[c("Similar_Known_Motifs", "Match_Names_Raw", "Match_IDs", "Match_Evalues")] <- NA
  }

  motifs_df <- motifs_df %>%
    dplyr::mutate(
      analysis_algorithm = dplyr::case_when(grepl("STREME", alt, ignore.case = TRUE) ~ "STREME",
                                     grepl("MEME", alt, ignore.case = TRUE) ~ "MEME",
                                     TRUE ~ "SEA"),
      Raw_Search_Text = dplyr::if_else(analysis_algorithm == "SEA",
                                dplyr::if_else(is.na(Match_Names_Raw), alt, paste(alt, Match_Names_Raw, sep = "; ")), Match_Names_Raw),
      Extracted_Gene_Names = .clean_gene_names(Raw_Search_Text),
      Extracted_Gene_Names = gsub("NA; |; NA", "", Extracted_Gene_Names),
      Extracted_Gene_Names = dplyr::if_else(Extracted_Gene_Names == "NA", NA_character_, Extracted_Gene_Names),
      dplyr::across(dplyr::any_of(c("sea_pvalue", "sea_evalue", "sea_qvalue", "pgm_evalue")), as.numeric))

  keep_cols <- c(Rank = "rank", analysis_algorithm = "analysis_algorithm",
                 Extracted_Gene_Names = "Extracted_Gene_Names", Motif_ID = "id", Alt_ID = "alt",
                 Consensus = "consensus", E_value = "pgm_evalue", Sites = "sites",
                 SEA_P_value = "sea_pvalue", SEA_E_value = "sea_evalue", SEA_Q_value = "sea_qvalue", Similar_Known_Motifs = "Similar_Known_Motifs",
                 Match_IDs = "Match_IDs", Match_Evalues = "Match_Evalues")
  present <- keep_cols[keep_cols %in% names(motifs_df)]
  out <- motifs_df[, present, drop = FALSE]
  names(out) <- names(present)
  if (nrow(out) == 0L) return(NULL)

  if (write) {
    utils::write.csv(out, file.path(xstreme_dir, "xstreme_summary.csv"), row.names = FALSE)
    if (requireNamespace("writexl", quietly = TRUE))
      writexl::write_xlsx(out, file.path(xstreme_dir, "xstreme_summary.xlsx"))
  }
  out
}


# -----------------------------------------------------------------------------
# map_tfs_to_zebrafish()
#
# Motif databases name TFs by mouse/human (sometimes zebrafish) symbols. This
# links every TF symbol of an xstreme summary to zebrafish genes in two ways:
#   - "Alliance_mouse" / "Alliance_human": bidirectional-best Alliance orthologs
#     (convert_to_zebrafish_alliance(), case-insensitive)
#   - "same_symbol": the symbol itself matches a gene in the data (case-insensitive)
# and flags whether that zebrafish gene was tested / has a peak in each region set,
# i.e. whether the TF's own gene carries H3K27ac signal.
#
# Args:
#   summary_df    : output of parse_xstreme()
#   alliance_file : Alliance TSV
#   tested_genes  : zebrafish symbols of all genes with a tested region
#   gene_sets     : named list of zebrafish symbol vectors (e.g. genes with a
#                   higher-in-KO / higher-in-WT / differential peak)
#
# Returns: long data.frame, one row per motif x TF symbol x zebrafish gene
# -----------------------------------------------------------------------------
map_tfs_to_zebrafish <- function(summary_df, alliance_file, tested_genes, gene_sets = list()) {
  tf <- summary_df %>%
    dplyr::filter(!is.na(Extracted_Gene_Names)) %>%
    tidyr::separate_rows(Extracted_Gene_Names, sep = ",\\s*") %>%
    dplyr::rename(TF_symbol = Extracted_Gene_Names) %>%
    dplyr::filter(nzchar(TF_symbol), TF_symbol != "NA") %>%
    dplyr::select(dplyr::any_of(c("Rank", "analysis_algorithm", "Motif_ID", "Alt_ID", "E_value", "TF_symbol"))) %>%
    dplyr::distinct()   # explicit dplyr:: — AnnotationDbi (org.*.eg.db) masks select()/filter() when loaded later
  if (nrow(tf) == 0L) return(tf)

  orth <- convert_to_zebrafish_alliance(unique(tf$TF_symbol), alliance_file, from = c("mouse", "human"), ignore_case = TRUE)
  orth <- if (is.null(orth) || nrow(orth) == 0L) {
    data.frame(TF_symbol = character(), zebrafish_gene = character(), link = character())
  } else {
    data.frame(TF_symbol = orth$input_gene, zebrafish_gene = orth$zebrafish_gene_symbol,
               link = paste0("Alliance_", orth$source_species), stringsAsFactors = FALSE)
  }
  same <- data.frame(TF_symbol = unique(tf$TF_symbol), stringsAsFactors = FALSE)
  same$zebrafish_gene <- tested_genes[match(toupper(same$TF_symbol), toupper(tested_genes))]
  same <- same[!is.na(same$zebrafish_gene), ]
  same$link <- rep("same_symbol", nrow(same))
  links <- rbind(orth, same)
  links <- links[!duplicated(links[, c("TF_symbol", "zebrafish_gene", "link")]), ]

  out <- dplyr::left_join(tf, links, by = "TF_symbol", relationship = "many-to-many")
  out$gene_tested <- !is.na(out$zebrafish_gene) & toupper(out$zebrafish_gene) %in% toupper(tested_genes)
  for (s in names(gene_sets))
    out[[paste0("peak_", s)]] <- !is.na(out$zebrafish_gene) & toupper(out$zebrafish_gene) %in% toupper(gene_sets[[s]])
  out
}


# -----------------------------------------------------------------------------
# make_motif_logos()
#
# ceqlogo logos for the motifs of a MEME-format file, named
# <prefix>_<gene>.{pdf,png}. ceqlogo only writes EPS/PNG, so the EPS is converted
# to PDF (ps2pdf / epstopdf, whichever exists) and removed.
#
# Args:
#   meme_file : MEME motif file (e.g. <xstreme_dir>/combined.meme)
#   out_dir   : output folder (created)
#   prefix    : file-name prefix
#   motif_ids : optional subset of motif ids (default: every MOTIF in the file)
#   reverse   : also write the reverse-complement logo (default FALSE)
#   meme_env  : conda env with MEME
#   workers   : motifs rendered concurrently (default NULL = lab_workers("light");
#               each motif is 2-3 external processes). 1 = one after the other.
#   BPPARAM   : optional BiocParallel backend (overrides workers)
#
# Returns: character vector of PNG paths written, invisibly
# -----------------------------------------------------------------------------
make_motif_logos <- function(meme_file, out_dir, prefix = "motif", motif_ids = NULL,
                             reverse = FALSE, meme_env = "/opt/conda/envs/archr_tools",
                             workers = NULL, BPPARAM = NULL) {
  ceqlogo <- Sys.glob(c(file.path(meme_env, "bin", "ceqlogo"), file.path(meme_env, "libexec", "meme-*", "ceqlogo")))[1]   # libexec in source builds
  stopifnot(!is.na(ceqlogo), file.exists(meme_file))
  # ceqlogo calls ghostscript for PNG; make the env's gs / libs visible to it
  env_vars <- c(paste0("PATH=", file.path(meme_env, "bin"), ":", Sys.getenv("PATH")),
                paste0("LD_LIBRARY_PATH=", file.path(meme_env, "lib"), ":", Sys.getenv("LD_LIBRARY_PATH")))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  motif_lines <- grep("^MOTIF ", readLines(meme_file, warn = FALSE), value = TRUE)
  rest        <- sub("^MOTIF +", "", motif_lines)
  ids         <- sub(" .*$", "", rest)
  desc        <- trimws(sub("^[^ ]+ *", "", rest))
  gene <- ifelse(desc == "", sub("[_.].*", "", ids),
          ifelse(startsWith(desc, "("), sub("^\\(([^)]*)\\).*", "\\1", desc),
                 sub(".*[ ]", "", desc)))
  if (!is.null(motif_ids)) { sel <- ids %in% motif_ids; ids <- ids[sel]; gene <- gene[sel] }
  ps2pdf <- Sys.which("ps2pdf"); if (!nzchar(ps2pdf) && file.exists(file.path(meme_env, "bin", "ps2pdf"))) ps2pdf <- file.path(meme_env, "bin", "ps2pdf")
  epstopdf <- Sys.which("epstopdf")
  gs <- Sys.which("gs"); if (!nzchar(gs) && file.exists(file.path(meme_env, "bin", "gs"))) gs <- file.path(meme_env, "bin", "gs")
  one_motif <- function(i) {
    written <- character(0)
    stem <- file.path(out_dir, paste0(prefix, "_", gsub("[^A-Za-z0-9_.-]", "_", gene[i])))
    for (rc in unique(c(FALSE, reverse))) {
      s <- if (rc) paste0(stem, "_reverse") else stem
      extra <- if (rc) c("-r") else character(0)
      system2(ceqlogo, c("-i", shQuote(meme_file), "-m", shQuote(ids[i]), extra, "-o", shQuote(paste0(s, ".eps")), "-f", "EPS"), stdout = FALSE, stderr = FALSE, env = env_vars)
      if (file.exists(paste0(s, ".eps"))) {
        if (nzchar(ps2pdf)) system2(ps2pdf, c("-dEPSCrop", shQuote(paste0(s, ".eps")), shQuote(paste0(s, ".pdf"))), env = env_vars)
        else if (nzchar(epstopdf)) system2(epstopdf, c(shQuote(paste0(s, ".eps")), paste0("--outfile=", shQuote(paste0(s, ".pdf")))))
        else warning("no ps2pdf/epstopdf found: kept the EPS, no PDF written")
        # ceqlogo's own PNG needs a ghostscript path compiled in; render the EPS ourselves (300 dpi)
        if (nzchar(gs)) system2(gs, c("-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dEPSCrop", "-sDEVICE=pngalpha", "-r300",
                                      shQuote(paste0("-sOutputFile=", s, ".png")), shQuote(paste0(s, ".eps"))), env = env_vars)
        if (file.exists(paste0(s, ".pdf"))) unlink(paste0(s, ".eps"))
      }
      if (file.exists(paste0(s, ".png"))) written <- c(written, paste0(s, ".png"))
    }
    written
  }
  bp  <- lab_bpparam(workers, length(ids), BPPARAM, profile = "light")
  out <- lab_run_tasks(as.list(seq_along(ids)), one_motif, bp)
  written <- character(0)
  for (o in out) {
    if (length(o$log) > 0) cat(paste(o$log, collapse = "\n"), "\n")
    r <- o$result
    if (is.list(r) && identical(r$status, "error")) warning("make_motif_logos(): ", r$message, call. = FALSE)
    else written <- c(written, r)
  }
  cat("✓ Logos:", length(written), "written →", out_dir, "\n")
  invisible(written)
}
