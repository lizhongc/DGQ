# ============================================================================
#  DGQ S3 methods.
# ============================================================================

# print(DGQ): Print a concise summary of a fitted DGQ.
#
# Reports the analyzed array dimensions, selected computational options, dynamic
# medoid, constrained empirical selections, and number of unconstrained
# time-wise trajectories. The fitted object is returned invisibly and its
# trajectory arrays are not printed.
#
# Arguments:
#   x   : An object of class "DGQ".
#   ... : Further arguments (currently unused).
print.DGQ <- function(x, ...) {
  ## Summarize dimensions, computation choices, and selected constrained paths
  ## without printing the potentially large N-by-T-by-d analyzed data array.
  cat("Dynamic Geometric Quantile (DGQ)\n")
  cat(sprintf("  collection : N = %d series, T = %d times, d = %d variables\n",
              x$N, x$T, x$d))
  cat(sprintf("  metric     : %s    depth: %s%s    cumulative: %s\n",
              x$metric, x$method,
              if (!is.null(x$sample.size)) sprintf(" (s = %d)", x$sample.size) else "",
              if (isTRUE(x$cumulative)) "yes" else "no"))
  if (!is.null(x$medoid))
    cat(sprintf("  medoid     : series #%d (dynamic median, u = 0)\n", x$medoid))
  if (!is.null(x$empirical)) {
    cat("  empirical DGQ (constrained, observed series) by direction:\n")
    idx <- x$empirical$index
    for (nm in names(idx))
      cat(sprintf("    %-14s -> series #%d\n", nm, idx[[nm]]))
  }
  if (!is.null(x$timewise))
    cat(sprintf("  time-wise geometric quantile (unconstrained): %d direction(s)\n",
                dim(x$timewise)[1L]))
  invisible(x)
}

# plot(DGQ): Visualize a fitted DGQ, one panel per requested variable.
#
# type selects the view:
#   "trajectory" (default) draws all N analyzed series in grey with the
#     empirical (solid) and time-wise (dashed) DGQ trajectories overlaid per
#     direction in a shared direction color.
#   "rainbow" colors every series by its trajectory depth C(j) (deep = central),
#     shades the central-50% band, and overlays the dynamic medoid (and the
#     empirical median when a tau = 0.5 direction is present).
#   "fan" draws the time-wise quantile fan (nested bands between symmetric tau
#     pairs, plus the median) with the empirical per-tau selections overlaid as
#     lines; requires a fit with at least two tau levels.
# The fitted object is returned invisibly.
#
# Arguments:
#   components : Integer vector of variable indices to plot (default all).
#   type       : One of "trajectory", "rainbow", "fan".
plot.DGQ <- function(x, components = NULL,
                     type = c("trajectory", "rainbow", "fan"), ...) {
  type <- match.arg(type)
  comps <- if (is.null(components)) seq_len(x$d) else components
  ## Restore the caller's graphics parameters even if plotting exits early.
  op <- graphics::par(mfrow = c(1, length(comps)), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(op))
  switch(type,
         trajectory = .dgq_plot_trajectory(x, comps),
         rainbow    = .dgq_plot_rainbow(x, comps),
         fan        = .dgq_plot_fan(x, comps))
  invisible(x)
}

## Data with empirical (solid) and time-wise (dashed) DGQ trajectories overlaid.
.dgq_plot_trajectory <- function(x, comps) {
  X <- x$X
  m <- length(x$directions)
  cols <- grDevices::hcl.colors(max(m, 1L), "Dark3")
  for (cc in comps) {
    graphics::matplot(t(X[, , cc]), type = "l", lty = 1, col = "grey85",
                      xlab = "time", ylab = paste("variable", cc),
                      main = paste("DGQ - variable", cc))
    for (g in seq_len(m)) {
      if (!is.null(x$empirical))
        graphics::lines(x$empirical$trajectory[g, , cc], col = cols[g], lwd = 2.2)
      if (!is.null(x$timewise))
        graphics::lines(x$timewise[g, , cc], col = cols[g], lwd = 1.2, lty = 3)
    }
  }
}

## Depth-ranked "rainbow": series colored by depth C(j) with the central-50%
## band, the dynamic medoid, and (when available) the empirical median.
.dgq_plot_rainbow <- function(x, comps) {
  N <- x$N
  tt <- seq_len(x$T)
  ord <- order(x$depth)                       # deepest (smallest C) first
  rank_of <- integer(N); rank_of[ord] <- seq_len(N)
  cols <- grDevices::hcl.colors(N, "viridis") # deep = dark, shallow = bright
  draw_order <- order(rank_of, decreasing = TRUE)  # shallow first, deep on top
  central <- ord[seq_len(ceiling(N / 2))]
  g_med <- match("tau=0.500", names(x$directions))

  leg_txt <- c("central 50% band", "medoid (deepest)")
  leg_col <- c(grDevices::adjustcolor("grey60", 0.5), "black")
  leg_lty <- c(1, 1); leg_lwd <- c(8, 3)
  if (!is.na(g_med)) {
    leg_txt <- c(leg_txt, "empirical median")
    leg_col <- c(leg_col, "#D55E00"); leg_lty <- c(leg_lty, 2); leg_lwd <- c(leg_lwd, 2)
  }

  for (cc in comps) {
    graphics::matplot(tt, t(x$X[, , cc]), type = "n",
                      xlab = "time", ylab = paste("variable", cc),
                      main = paste("DGQ rainbow - variable", cc))
    band <- apply(matrix(x$X[central, , cc], nrow = length(central)), 2, range)
    graphics::polygon(c(tt, rev(tt)), c(band[1, ], rev(band[2, ])),
                      col = grDevices::adjustcolor("grey60", 0.25), border = NA)
    for (i in draw_order)
      graphics::lines(tt, x$X[i, , cc], col = cols[rank_of[i]], lwd = 0.8)
    graphics::lines(tt, x$X[x$medoid, , cc], col = "black", lwd = 3)
    if (!is.na(g_med))
      graphics::lines(tt, x$empirical$trajectory[g_med, , cc],
                      col = "#D55E00", lwd = 2, lty = 2)
    if (cc == comps[1])
      graphics::legend("topleft", legend = leg_txt, col = leg_col,
                       lty = leg_lty, lwd = leg_lwd, bty = "n", cex = 0.8)
  }
}

## Quantile "fan": nested time-wise bands between symmetric tau pairs plus the
## median, with the empirical per-tau selections overlaid as lines.
.dgq_plot_fan <- function(x, comps) {
  m <- length(x$directions)
  if (m < 2L)
    stop("type = \"fan\" needs at least 2 tau levels; refit DGQ() with a vector tau")
  tt <- seq_len(x$T)
  taus <- as.numeric(sub("tau=", "", names(x$directions)))
  ord <- order(taus)                          # ascending tau
  line_cols <- grDevices::hcl.colors(m, "Dark3")[ord]
  npair <- floor(m / 2)
  band_cols <- grDevices::colorRampPalette(c("#DEEBF7", "#3182BD"))(npair)
  for (cc in comps) {
    graphics::matplot(tt, t(x$X[, , cc]), type = "n",
                      xlab = "time", ylab = paste("variable", cc),
                      main = paste("DGQ fan - variable", cc))
    for (i in seq_len(x$N))
      graphics::lines(tt, x$X[i, , cc], col = "grey90")
    for (p in seq_len(npair))
      graphics::polygon(c(tt, rev(tt)),
                        c(x$timewise[ord[p], , cc],
                          rev(x$timewise[ord[m - p + 1], , cc])),
                        col = band_cols[p], border = NA)
    if (m %% 2L == 1L)
      graphics::lines(tt, x$timewise[ord[(m + 1) / 2], , cc], col = "black", lwd = 2)
    for (j in seq_along(ord))
      graphics::lines(tt, x$empirical$trajectory[ord[j], , cc],
                      col = line_cols[j], lwd = 1.6)
    if (cc == comps[1])
      graphics::legend("topleft", title = "empirical lines / time-wise bands",
                       legend = names(x$directions)[ord], col = line_cols,
                       lty = 1, lwd = 1.6, bty = "n", cex = 0.7)
  }
}
