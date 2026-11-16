# =============================================================================
# Mercader Lab — Parallel / package helpers (shared by every other file)
# parallel_helpers.R
#
# Sourced FIRST by source_all.R. Every other file sources it on demand if it is
# missing (so the files can still be sourced standalone).
#
# Functions:
#   - check_required_packages() : installed-package check with install hints
#   - lab_workers()             : resolve a worker count (arg > option > env > profile default)
#   - lab_bpparam()             : build a lazily-started BiocParallel backend
#   - lab_run_tasks()           : run a task list in parallel with log capture, error
#                                 capture, longest-first ordering and worker hygiene
#   - .with_blas_threads()      : run a serial dense-linear-algebra step with N BLAS threads
#   - .push_seed()              : seed a function body, restore the caller's RNG afterwards
#
# Worker defaults (override with options(mercader.workers = N) or the env var
# MERCADER_WORKERS; per profile: options(mercader.workers.light / .heavy)):
#   light  (enrichment, Enrichr, motif logos) : min(32, cores - 2)
#   heavy  (single-cell DE, LIANA)            : min(12, cores %/% 4)
# =============================================================================

# Directory that holds the shared functions (absolute). Used to bootstrap
# SnowParam workers, which start with an empty session.
if (!exists(".LAB_SRC_DIR", inherits = FALSE)) {
  .LAB_SRC_DIR <- local({
    for (i in rev(seq_len(sys.nframe()))) {
      f <- sys.frame(i)$ofile
      if (is.character(f) && grepl("parallel_helpers|source_all", basename(f)))
        return(normalizePath(dirname(f), mustWork = FALSE))
    }
    normalizePath("reusable_functions", mustWork = FALSE)
  })
}

# Packages that live on CRAN; anything else is assumed to be Bioconductor.
.CRAN_PKGS <- c("ggplot2", "dplyr", "msigdbr", "tidyr", "tibble", "purrr", "BiocManager",
                "data.table", "RhpcBLASctl", "stringdist", "openxlsx", "enrichR", "writexl",
                "jsonlite", "stringr", "ggrepel", "patchwork", "Seurat", "cluster")

# -----------------------------------------------------------------------------
# check_required_packages()
#
# Checks that packages are installed (without loading them). If any are missing,
# prints ready-to-paste install commands.
#
# Args:
#   pkgs    : character vector of package names
#   context : short label for the message, e.g. "ORA_function(onto = 'reactome')"
#   stop_on_missing : TRUE (default) stops with an error after the message;
#                     FALSE only warns and returns FALSE
#
# Returns: TRUE (all present) or FALSE, invisibly
# -----------------------------------------------------------------------------
check_required_packages <- function(pkgs, context = "", stop_on_missing = TRUE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))

  cran <- missing[missing %in% .CRAN_PKGS]
  bioc <- setdiff(missing, cran)
  fmt  <- function(x) paste0('c("', paste(x, collapse = '", "'), '")')
  msg  <- paste0(
    "\n✗ Missing R package(s)", if (nzchar(context)) paste0(" required by ", context),
    ": ", paste(missing, collapse = ", "), "\n",
    "  Install with:\n",
    if (length(bioc) > 0 || "BiocManager" %in% missing)
      '    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")\n',
    if (length(cran) > 0) paste0("    install.packages(", fmt(setdiff(cran, "BiocManager")), ")\n"),
    if (length(bioc) > 0) paste0("    BiocManager::install(", fmt(bioc), ")\n"),
    "  (inside the Docker image: re-run with the lab image or install into the mounted R library)\n")
  if (stop_on_missing) stop(msg, call. = FALSE)
  cat(msg); warning("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(FALSE)
}

# Sets the RNG seed for the duration of the calling function and restores the
# caller's RNG state afterwards (the user's session is not reseeded).
# Usage inside a function:  restore <- .push_seed(seed); on.exit(restore(), add = TRUE)
# The RNG KIND is pinned too (Mersenne-Twister / Inversion / Rejection): BiocParallel
# switches forked workers to L'Ecuyer-CMRG, so the same set.seed() gave different
# permutation results (GSEA) inside a parallel loop than in a direct call.
.push_seed <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  old_kind <- RNGkind()
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  function() {
    if (had) assign(".Random.seed", old, envir = globalenv())     # restores state AND kind
    else {
      suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    }
    invisible(NULL)
  }
}

.lab_cores <- function() {
  n <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_)
  if (is.na(n) || n < 1L) 2L else as.integer(n)
}

# -----------------------------------------------------------------------------
# lab_workers()
#
# Resolves the number of workers for a loop.
#   workers > options(mercader.workers.<profile>) > options(mercader.workers)
#           > env MERCADER_WORKERS > profile default; always capped at n_tasks.
#
# Args:
#   n_tasks : number of independent tasks (cap); Inf = no cap
#   profile : "light" (cheap, memory-light tasks) or "heavy" (sc DE, big matrices)
#   workers : explicit request (NULL = resolve)
#
# Returns: integer >= 1
# -----------------------------------------------------------------------------
lab_workers <- function(n_tasks = Inf, profile = c("light", "heavy"), workers = NULL) {
  profile <- match.arg(profile)
  cores <- .lab_cores()
  w <- workers
  if (is.null(w)) w <- getOption(paste0("mercader.workers.", profile))
  if (is.null(w)) w <- getOption("mercader.workers")
  if (is.null(w)) { e <- Sys.getenv("MERCADER_WORKERS", ""); if (nzchar(e)) w <- suppressWarnings(as.integer(e)) }
  if (is.null(w) || is.na(w))
    w <- if (profile == "light") min(32L, cores - 2L) else min(12L, cores %/% 4L)
  w <- as.integer(w)
  cap <- if (is.finite(n_tasks)) as.integer(max(1, n_tasks)) else .Machine$integer.max
  max(1L, min(w, cap))
}

# -----------------------------------------------------------------------------
# lab_bpparam()
#
# Builds the BiocParallel backend for a loop. The backend is returned UNSTARTED
# on purpose: with MulticoreParam the workers fork inside bplapply(), after the
# parent has warmed its caches, so the children inherit them copy-on-write.
#   tasks = n_tasks gives one chunk per task (dynamic load balancing) instead of
#   BiocParallel's default of one contiguous chunk per worker.
#
# Args:
#   workers : number of workers (NULL = lab_workers(); 1 = serial)
#   n_tasks : number of tasks (cap and chunking)
#   BPPARAM : a user-supplied backend; wins over `workers`
#   profile : "light" / "heavy" (see lab_workers)
#
# Returns: a BiocParallelParam
# -----------------------------------------------------------------------------
lab_bpparam <- function(workers = NULL, n_tasks = Inf, BPPARAM = NULL, profile = "light") {
  check_required_packages("BiocParallel", "parallel loops")
  if (!is.null(BPPARAM)) {
    if (isTRUE(tryCatch(BiocParallel::bpisup(BPPARAM), error = function(e) FALSE)))
      warning("BPPARAM is already started: parent-side caches will be serialized to the workers ",
              "instead of being inherited. Pass an unstarted param (e.g. MulticoreParam(16)) for speed.",
              call. = FALSE)
    return(BPPARAM)
  }
  w <- lab_workers(n_tasks, profile, workers)
  if (w == 1L) return(BiocParallel::SerialParam())
  n_chunks <- if (is.finite(n_tasks)) as.integer(n_tasks) else 0L
  pb <- isTRUE(getOption("mercader.progressbar", FALSE))
  if (.Platform$OS.type == "windows")
    BiocParallel::SnowParam(workers = w, tasks = n_chunks, progressbar = pb)
  else
    BiocParallel::MulticoreParam(workers = w, tasks = n_chunks, progressbar = pb)
}

# Worker hygiene: one thread per worker, and never let a library fan out again
# inside a worker (fgsea / clusterProfiler via BiocParallel, MAST via mc.cores,
# OpenBLAS via its own pool).
.lab_worker_prologue <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(1L), silent = TRUE)
  }
  options(mc.cores = 1L, mercader.in_worker = TRUE)
  if (requireNamespace("BiocParallel", quietly = TRUE))
    BiocParallel::register(BiocParallel::SerialParam())
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# .with_blas_threads()
#
# Evaluates `expr` with the BLAS/OpenMP thread count set to `n` (restored
# afterwards). For SERIAL dense linear algebra in the parent (lmFit, voom,
# prcomp, dist). Inside parallel workers the prologue already pins to 1.
#
# Args:
#   n    : threads (NULL = min(16, cores))
#   expr : expression to evaluate
# -----------------------------------------------------------------------------
.with_blas_threads <- function(n = NULL, expr) {
  if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) return(expr)
  if (is.null(n)) n <- min(16L, .lab_cores())
  old_blas <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NULL)
  old_omp  <- tryCatch(RhpcBLASctl::omp_get_max_threads(), error = function(e) NULL)
  if (!is.null(old_blas)) {
    try(RhpcBLASctl::blas_set_num_threads(as.integer(n)), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(as.integer(n)), silent = TRUE)
    on.exit({
      try(RhpcBLASctl::blas_set_num_threads(old_blas), silent = TRUE)
      if (!is.null(old_omp)) try(RhpcBLASctl::omp_set_num_threads(old_omp), silent = TRUE)
    }, add = TRUE)
  }
  expr
}

# -----------------------------------------------------------------------------
# lab_run_tasks()
#
# Runs `tasks` (a list) through FUN on a BiocParallel backend.
#   - console output of every task is captured on the worker and returned, so
#     the log reads in task order instead of interleaving;
#   - errors are caught per task and returned as list(status = "error", message = ...);
#   - with `cost`, tasks are dispatched longest-first (results come back in the
#     original order);
#   - the parent pins BLAS to 1 thread around the fork (the parent's thread pool
#     is what fork duplicates) and restores it afterwards;
#   - every worker runs .lab_worker_prologue(); the serial path registers
#     SerialParam() so nothing (fgsea, ...) forks behind our back;
#   - SnowParam workers (empty sessions) bootstrap by sourcing source_all.R.
#
# Args:
#   tasks   : list of task objects
#   FUN     : function(task) -> result
#   BPPARAM : backend from lab_bpparam()
#   cost    : optional numeric vector (length(tasks)), higher = run earlier
#
# Returns: list (same length/order as tasks) of list(result = <FUN output or
#   error list>, log = <character vector of captured output>)
# -----------------------------------------------------------------------------
lab_run_tasks <- function(tasks, FUN, BPPARAM, cost = NULL) {
  n <- length(tasks)
  if (n == 0L) return(list())
  check_required_packages("BiocParallel", "lab_run_tasks()")
  serial <- inherits(BPPARAM, "SerialParam")

  worker_fun <- function(task, task_fun, src_dir, serial) {
    log <- character(); out <- NULL
    log <- utils::capture.output({
      out <- tryCatch({
        if (!exists("lab_run_tasks", mode = "function")) {       # SnowParam worker: empty session
          f <- file.path(src_dir, "source_all.R")
          if (!file.exists(f))
            stop("worker cannot find '", f, "' - run from the repo root or use the default MulticoreParam")
          suppressMessages(source(f))
        }
        if (!serial) .lab_worker_prologue()
        task_fun(task)
      }, error = function(e) list(status = "error", n_terms = NA_integer_,
                                  message = conditionMessage(e)))
    })
    list(result = out, log = log)
  }

  # parent-side guards (restored on exit)
  old_bp <- tryCatch(BiocParallel::bpparam(), error = function(e) NULL)
  BiocParallel::register(BiocParallel::SerialParam())
  on.exit(if (!is.null(old_bp)) BiocParallel::register(old_bp), add = TRUE)
  if (!serial && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_blas <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NULL)
    old_omp  <- tryCatch(RhpcBLASctl::omp_get_max_threads(), error = function(e) NULL)
    try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(1L), silent = TRUE)
    on.exit({
      if (!is.null(old_blas)) try(RhpcBLASctl::blas_set_num_threads(old_blas), silent = TRUE)
      if (!is.null(old_omp))  try(RhpcBLASctl::omp_set_num_threads(old_omp), silent = TRUE)
    }, add = TRUE)
  }

  # Force every promise in FUN's enclosing frame (lazily evaluated loop arguments such as
  # out_root = file.path(OUT, "x")): a forked worker could resolve them from inherited memory,
  # a SnowParam worker cannot ("object 'OUT' not found"). Forcing is harmless otherwise.
  fenv <- environment(FUN)
  if (!is.null(fenv) && !identical(fenv, globalenv()) && !identical(fenv, baseenv()))
    for (nm in ls(fenv, all.names = TRUE))
      tryCatch(get(nm, envir = fenv, inherits = FALSE), error = function(e) NULL)

  ord <- if (!is.null(cost) && !serial && length(cost) == n) order(cost, decreasing = TRUE) else seq_len(n)
  out <- BiocParallel::bplapply(tasks[ord], worker_fun, task_fun = FUN, src_dir = .LAB_SRC_DIR,
                                serial = serial, BPPARAM = BPPARAM)
  res <- vector("list", n); res[ord] <- out
  res
}

# Backward-compatible aliases (old internal names used by enrichment_functions.R)
.enrichment_bpparam   <- function(BPPARAM, workers, n_tasks) lab_bpparam(workers, n_tasks, BPPARAM, "light")
.run_enrichment_tasks <- function(tasks, FUN, BPPARAM) lab_run_tasks(tasks, FUN, BPPARAM)

.lab_banner_parallel <- function() {
  cores <- .lab_cores()
  sprintf("Parallel: %d cores -> light = %d, heavy = %d workers (options(mercader.workers=) / MERCADER_WORKERS)",
          cores, lab_workers(Inf, "light"), lab_workers(Inf, "heavy"))
}
