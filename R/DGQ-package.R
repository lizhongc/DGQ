#' DGQ: Dynamic Geometric Quantiles for Multivariate Time Series
#'
#' The package computes the Dynamic Geometric Quantile (DGQ) for a collection of
#' multivariate time series. The single entry point [DGQ()] returns both the
#' unconstrained *time-wise geometric quantile* (a free trajectory, the
#' population target) and the constrained *empirical DGQ* (a real observed
#' series) for one or more tilt directions.
#'
#' The method extends the empirical dynamic quantile of Pena, Tsay and Zamar
#' (2019) to the multivariate setting via Chaudhuri's (1996) geometric quantile.
#' The numerically heavy core (the pairwise depth and the Weiszfeld iteration) is
#' implemented in C and reached through the R C interface.
#'
#' @references
#' Chaudhuri, P. (1996). On a geometric notion of quantiles for multivariate
#' data. *J. Amer. Statist. Assoc.* 91(434), 862--872.
#'
#' Pena, D., Tsay, R. S., & Zamar, R. (2019). Empirical dynamic quantiles for
#' visualization of high-dimensional time series. *Technometrics* 61(4),
#' 429--444.
#'
#' Kaufman, L. & Rousseeuw, P. J. (1990). *Finding Groups in Data*. Wiley
#' (program CLARA).
#'
#' @useDynLib DGQ, .registration = TRUE
#' @importFrom stats cov
#' @keywords internal
"_PACKAGE"
