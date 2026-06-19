# ============================================================================
#  DGQ S3 methods.
# ============================================================================

#' @describeIn DGQ Print a concise summary of a fitted DGQ.
#'
#' Reports the analyzed array dimensions, selected computational options,
#' dynamic medoid, constrained empirical selections, and number of
#' unconstrained time-wise trajectories. The fitted object is returned
#' invisibly and its trajectory arrays are not printed.
#' @param x An object of class `"DGQ"`.
#' @param ... Further arguments (currently unused).
#' @export
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

#' @describeIn DGQ Plot the data with the empirical (solid) and time-wise
#'   (dashed) DGQ trajectories, one panel per variable.
#'
#' For each requested variable, all `N` analyzed trajectories are drawn in
#' grey. Direction-specific empirical trajectories are overlaid as solid lines
#' and time-wise trajectories as dashed lines using a shared direction color.
#' The fitted object is returned invisibly.
#' @param components Integer vector of variable indices to plot (default all).
#' @export
plot.DGQ <- function(x, components = NULL, ...) {
  ## Draw one panel per requested variable. Grey lines are all analyzed paths;
  ## colors identify directions, with solid constrained empirical trajectories
  ## and dashed unconstrained time-wise trajectories.
  X <- x$X; d <- x$d
  comps <- if (is.null(components)) seq_len(d) else components
  m <- length(x$directions)
  cols <- grDevices::hcl.colors(max(m, 1L), "Dark3")
  ## Restore the caller's graphics parameters even if plotting exits early.
  op <- graphics::par(mfrow = c(1, length(comps)), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(op))
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
  invisible(x)
}
