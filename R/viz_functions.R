# =============================================================================
# Mercader Lab — Visualization Functions
# viz_functions.R
#
# Functions:
#   - PCA()          : PCA plot from a count/expression matrix
#   - draw_volcano() : Volcano plot from DESeq2-style results
#   - lab_palette_colors(), lab_scale_continuous(), lab_scale_discrete() : shared colour helpers
#
# Usage: source("R/viz_functions.R")
# =============================================================================

library(ggplot2)
library(ggrepel)

if (!exists("lab_run_tasks", mode = "function")) local({
  d <- "reusable_functions"
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (is.character(f) && grepl("viz_functions", basename(f))) { d <- dirname(f); break }
  }
  source(file.path(d, "parallel_helpers.R"))
})

options(ggrepel.max.overlaps = 20)

# =============================================================================
# Shared color helpers — used by every plotting function in the library
# =============================================================================
# `palette` can be given in any of these forms, everywhere in the library:
#   "red_blue"                       : red (low) -> blue (high); lab default for p.adjust
#   "viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo"
#                                    : viridisLite colour maps
#   any RColorBrewer name            : e.g. "RdBu", "RdYlBu", "Spectral", "Set1", "Dark2"
#   colour-blind-safe categorical sets (see .LAB_PALETTES below):
#     "tol_bright" (7, PCA default), "tol_muted" (9), "tol_vibrant" (7),
#     "tol_high_contrast" (3), "carto_safe" (10), "tableau_cb" (10), "ibm" (5),
#     "okabe_ito" (8); "ggplot" = ggplot2's default hues
#     (more groups than colours -> colours are interpolated, with a warning)
#   a character vector of colours    : e.g. c("#FF0000", "white", "#0000FF") or
#                                      c("darkred", "gold"); a NAMED vector maps colours
#                                      to group levels in discrete plots (e.g. PCA)
# reverse_palette = TRUE flips the colour order (e.g. so that significant = bright
# yellow instead of dark purple with viridis).
# Functions: lab_palette_colors(), lab_scale_continuous(), lab_scale_discrete()
# -----------------------------------------------------------------------------
.LAB_VIRIDIS <- c("viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo")
# Colour-blind-safe categorical palettes (Paul Tol, CARTO, Tableau, IBM, Okabe & Ito)
.LAB_PALETTES <- list(
  tol_bright        = c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"),
  tol_muted         = c("#332288", "#CC6677", "#117733", "#DDCC77", "#88CCEE", "#882255", "#44AA99",
                        "#999933", "#AA4499"),
  tol_vibrant       = c("#0077BB", "#EE7733", "#009988", "#EE3377", "#33BBEE", "#CC3311", "#BBBBBB"),
  tol_high_contrast = c("#004488", "#BB5566", "#DDAA33"),
  carto_safe        = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99",
                        "#999933", "#882255", "#661100"),
  tableau_cb        = c("#1170AA", "#FC7D0B", "#57606C", "#5FA2CE", "#C85200", "#A3ACB9", "#7B848F",
                        "#A3CCE9", "#FFBC79", "#C8D0D9"),
  ibm               = c("#648FFF", "#DC267F", "#FFB000", "#785EF0", "#FE6100"),
  okabe_ito         = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                        "#999999"))

# Returns n colours for `palette` (see above). For named custom vectors returns them as given.
lab_palette_colors <- function(palette, n, reverse_palette = FALSE) {
  if (is.null(palette) || length(palette) == 0) stop("palette is empty")
  cols <- if (length(palette) == 1 && identical(palette, "red_blue")) {
    grDevices::colorRampPalette(c("red", "blue"))(n)
  } else if (length(palette) == 1 && palette %in% .LAB_VIRIDIS) {
    viridisLite::viridis(n, option = palette)
  } else if (length(palette) == 1 && palette == "ggplot") {
    scales::hue_pal()(n)
  } else if (length(palette) == 1 && palette %in% names(.LAB_PALETTES)) {
    base <- .LAB_PALETTES[[palette]]
    if (n <= length(base)) base[seq_len(n)] else {
      warning("palette \"", palette, "\" has ", length(base), " colours but ", n,
              " are needed; interpolating (colours get harder to tell apart). ",
              "Consider \"tol_muted\" (9) / \"carto_safe\" (10) or add shape_pca.")
      grDevices::colorRampPalette(base)(n)
    }
  } else if (length(palette) == 1 && palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    maxc <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    base <- RColorBrewer::brewer.pal(max(3, min(n, maxc)), palette)
    if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
  } else {
    bad <- palette[!vapply(palette, function(x)
      tryCatch({ grDevices::col2rgb(x); TRUE }, error = function(e) FALSE), logical(1))]
    if (length(bad) > 0)
      stop("Unknown palette / colour(s): ", paste(bad, collapse = ", "),
           ". Use \"red_blue\", a viridis option (", paste(.LAB_VIRIDIS, collapse = ", "),
           "), a lab palette (", paste(c(names(.LAB_PALETTES), "ggplot"), collapse = ", "),
           "), an RColorBrewer name, or a vector of valid colours.")
    if (!is.null(names(palette))) return(if (reverse_palette) rev(palette) else palette)
    if (n <= length(palette)) palette[seq_len(n)] else {
      if (length(palette) > 1)
        warning("palette has ", length(palette), " colours but ", n, " are needed; interpolating.")
      grDevices::colorRampPalette(palette)(n)
    }
  }
  if (reverse_palette) rev(cols) else cols
}

# Continuous colour/fill scale. "red_blue" reproduces scale_*_gradient(low = "red", high = "blue").
lab_scale_continuous <- function(palette = "red_blue", aesthetic = "colour",
                                 reverse_palette = FALSE, ...) {
  aesthetic <- if (aesthetic %in% c("color", "colour")) "colour" else aesthetic
  if (length(palette) == 1 && palette %in% .LAB_VIRIDIS) {
    f <- if (aesthetic == "fill") ggplot2::scale_fill_viridis_c else ggplot2::scale_colour_viridis_c
    return(f(option = palette, direction = if (reverse_palette) -1 else 1, ...))
  }
  cols <- if (length(palette) == 1 && palette %in% names(.LAB_PALETTES)) .LAB_PALETTES[[palette]]
          else if (length(palette) == 1) lab_palette_colors(palette, if (palette == "red_blue") 2 else 9)
          else unname(lab_palette_colors(palette, length(palette)))
  if (reverse_palette) cols <- rev(cols)
  if (length(cols) == 2) {
    f <- if (aesthetic == "fill") ggplot2::scale_fill_gradient else ggplot2::scale_colour_gradient
    f(low = cols[1], high = cols[2], ...)
  } else {
    f <- if (aesthetic == "fill") ggplot2::scale_fill_gradientn else ggplot2::scale_colour_gradientn
    f(colours = cols, ...)
  }
}

# Discrete colour/fill scale. palette = NULL -> no scale (ggplot2 default hues).
lab_scale_discrete <- function(palette = NULL, aesthetic = "colour", reverse_palette = FALSE, ...) {
  if (is.null(palette)) return(NULL)
  aesthetic <- if (aesthetic %in% c("color", "colour")) "colour" else aesthetic
  if (!is.null(names(palette))) {
    lab_palette_colors(palette, length(palette))                          # validates
    f <- if (aesthetic == "fill") ggplot2::scale_fill_manual else ggplot2::scale_colour_manual
    return(f(values = if (reverse_palette) rev(palette) else palette, ...))
  }
  lab_palette_colors(palette, 2)                                          # validate early
  ggplot2::discrete_scale(aesthetic,
    palette = function(n) unname(lab_palette_colors(palette, n, reverse_palette)), ...)
}

# -----------------------------------------------------------------------------
# PCA()
#
# Args:
#   mat          : numeric matrix, samples as columns
#   color_pca    : vector mapping samples to color groups
#   shape_pca    : vector mapping samples to shape groups
#   label_pca    : vector of sample labels for repel text
#   save_plot    : "yes" saves a PDF and a 300-dpi PNG; "no" (default) prints and returns the plot
#   name_of_plot : base filename (no extension) when save_plot = "yes"
#   comp1, comp2 : which PCs to plot (default 1, 2)
#   pdf_height, pdf_width : dimensions in inches when saving
#   palette      : colours for the color_pca groups. Default "tol_bright" (Paul Tol,
#                  colour-blind safe, 7 colours). Any palette name ("tol_muted",
#                  "okabe_ito", "viridis", "Set1", "ggplot" = old ggplot2 hues, ...),
#                  any vector of custom colours, e.g. c("#1B9E77", "orange", "grey30"),
#                  or a NAMED vector mapping groups to colours,
#                  e.g. c(WT = "grey40", KO = "firebrick"). NULL = ggplot2 default hues.
#   reverse_palette : reverse the palette order (default FALSE)
#
# Returns: ggplot object
# -----------------------------------------------------------------------------
PCA <- function(mat, color_pca = "", shape_pca = "", label_pca = "",
                save_plot = "no", name_of_plot = "PCA",
                comp1 = 1, comp2 = 2, pdf_height = 12, pdf_width = 12,
                palette = "tol_bright", reverse_palette = FALSE) {
  dt <- mat
  pca_dt <- .with_blas_threads(NULL, prcomp(t(dt)))   # multi-threaded LAPACK for big matrices
  cat("PCA running...\n")
  percentVar_dt <- pca_dt$sdev^2 / sum(pca_dt$sdev^2)
  cat("Percents calculated...\n")
  dt_f <- data.frame(
    PC1       = pca_dt$x[, comp1],
    PC2       = pca_dt$x[, comp2],
    color_pca = color_pca,
    shape_pca = shape_pca,
    label_pca = label_pca
  )
  cat("Data frame built...\n")
  cat("Plotting...\n")

  pca_p <- ggplot(data = dt_f,
                  aes(x = .data[["PC1"]], y = .data[["PC2"]],
                      color = .data[["color_pca"]], shape = .data[["shape_pca"]],
                      label = .data[["label_pca"]])) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = 50,
                    box.padding = 1.5, point.padding = 0.5, force = 50) +
    xlab(paste0("PC", comp1, ": ", round(percentVar_dt[comp1] * 100), "% variance")) +
    ylab(paste0("PC", comp2, ": ", round(percentVar_dt[comp2] * 100), "% variance")) +
    theme_bw()
  if (!is.null(palette))
    pca_p <- pca_p + lab_scale_discrete(palette, "colour", reverse_palette)

  if (save_plot == "yes") {
    pdf(paste0(name_of_plot, ".pdf"), width = pdf_width, height = pdf_height)
    cat("Saving plot as:", paste0(name_of_plot, ".pdf\n"))
    print(pca_p)
    dev.off()
    ggsave(paste0(name_of_plot, ".png"), pca_p, width = pdf_width, height = pdf_height, dpi = 300)
  }

  cat("Done\n")
  print(pca_p)
  invisible(pca_p)
}


# -----------------------------------------------------------------------------
# draw_volcano()
#
# Args:
#   fileinput : data frame with columns log2FoldChange and padj
#   title     : string in the form "GroupA vs GroupB" — group names are parsed
#               automatically to annotate axes (red = higher in GroupA,
#               blue = higher in GroupB)
#   lfcthres  : absolute log2FC threshold for coloring points (default 1)
#   col_up, col_down, col_ns : point colours for significant up / significant down /
#               other genes (defaults "#CD0000" red, "#3A5FCD" blue, "grey51")
#   label_col_up, label_col_down : colours of the "Higher in ..." axis annotations.
#               Default: "red3" / "royalblue3" while col_up / col_down are left at their
#               defaults, otherwise the same colour as col_up / col_down
#
# Returns: ggplot object (does NOT save — wrap in ggsave() or pdf() at call site)
# -----------------------------------------------------------------------------
draw_volcano <- function(fileinput, title, lfcthres = 1,
                         col_up = "#CD0000", col_down = "#3A5FCD", col_ns = "grey51",
                         label_col_up   = if (identical(col_up,   "#CD0000")) "red3"       else col_up,
                         label_col_down = if (identical(col_down, "#3A5FCD")) "royalblue3" else col_down) {
  # Clean data
  fileinput <- fileinput[is.finite(fileinput$log2FoldChange), ]
  fileinput$padj[is.na(fileinput$padj) | !is.finite(fileinput$padj) | fileinput$padj < 0] <- 1
  fileinput$padj[fileinput$padj < 1e-100] <- 1e-110
  fileinput$log2FoldChange[fileinput$log2FoldChange >  7] <-  7
  fileinput$log2FoldChange[fileinput$log2FoldChange < -7] <- -7

  xlim <- c(round(min(fileinput$log2FoldChange, na.rm = TRUE) - 1),
            round(max(fileinput$log2FoldChange, na.rm = TRUE) + 1))
  ylim <- c(-1, round(-log10(min(fileinput$padj, na.rm = TRUE)) + 1))

  sig_down <- fileinput[fileinput$padj < 0.05 & fileinput$log2FoldChange < -lfcthres, ]
  sig_up   <- fileinput[fileinput$padj < 0.05 & fileinput$log2FoldChange >  lfcthres, ]

  p <- ggplot(data = fileinput, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", col = "darkgoldenrod") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.5, size = 2, color = col_ns) +
    geom_point(data = sig_down,
               aes(x = log2FoldChange, y = -log10(padj)),
               shape = 21, color = col_down, fill = col_down,
               alpha = 0.5, size = 2, inherit.aes = FALSE) +
    geom_point(data = sig_up,
               aes(x = log2FoldChange, y = -log10(padj)),
               shape = 21, color = col_up, fill = col_up,
               alpha = 0.5, size = 2, inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq(xlim[1], xlim[2], by = 1), limits = xlim) +
    scale_y_continuous(breaks = seq(0, ylim[2], by = 4), limits = ylim) +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title   = element_text(hjust = 0.5),
          axis.text    = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))

  # Auto-annotate axes from "GroupA vs GroupB" title
  match <- regexec("^(.*)\\s+[Vv][Ss]\\s+(.*)$", title)
  parts <- regmatches(title, match)[[1]]
  if (length(parts) == 3) {
    group1  <- trimws(parts[2])
    group2  <- trimws(parts[3])
    y_annot <- -0.75
    p <- p +
      annotate("text", x = xlim[2] - 0.5, y = y_annot,
               label = paste("Higher in", group1),
               hjust = 1, size = 3.5, fontface = "italic", color = label_col_up) +
      annotate("text", x = xlim[1] + 0.5, y = y_annot,
               label = paste("Higher in", group2),
               hjust = 0, size = 3.5, fontface = "italic", color = label_col_down)
  }
  return(p)
}

cat("✓ viz_functions.R loaded: PCA(), draw_volcano(), lab_palette_colors(), lab_scale_continuous(), lab_scale_discrete()\n")
