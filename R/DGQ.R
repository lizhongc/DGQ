# ============================================================================
#  DGQ: main user-facing function.
#  Heavy numerics (depth C(j)/V(j) and the time-wise Weiszfeld iteration) are
#  done in C; see src/empirical.c, src/timewise.c, and the R/SEXP boundary in
#  src/R_export.c.
# ============================================================================

# Dynamic Geometric Quantile of a collection of multivariate time series.
#
# Computes the Dynamic Geometric Quantile (DGQ) for a collection of N
# multivariate trajectories. For each quantile level in `tau`, the function
# builds a tilt trajectory in the open unit ball and returns two objects:
#
#   * the empirical DGQ (the constrained estimator): the real observed series
#     minimising C(j) + N <u, V(j)>; and
#   * the time-wise geometric quantile (the unconstrained population target): a
#     free trajectory whose value at each time t is Chaudhuri's geometric
#     u-quantile of the time-t marginal.
#
# Both are returned.
#
# Let Z_i(t) denote trajectory i at time t after applying the selected metric
# transformation, and let Zbar(t) be its cross-sectional mean. The empirical
# branch computes
#
#   C(j) = sum_i sum_t ||Z_i(t) - Z_j(t)||
#
# and
#
#   V(j) = sum_t {Zbar(t) - Z_j(t)}.
#
# For each direction u, it returns the observed trajectory whose index minimizes
# C(j) + N <u, V(j)>. At u = 0, this is the dynamic medoid that minimizes total
# trajectory distance. With method = "clara", the sum over i in C(j) is estimated
# from reference subsamples and scaled by N/s; V(j) remains exact. Empirical
# candidate calculations can be distributed over OpenMP threads with n.cores > 1;
# builds without OpenMP support ignore the request and execute serially.
#
# Independently at each time t, the time-wise branch minimizes
#
#   sum_i ||Z_i(t) - q|| - N <u, q>
#
# over unrestricted q. Its first-order equation is
# sum_i (q - Z_i(t)) / ||q - Z_i(t)|| = N u, which is solved by a modified
# Weiszfeld iteration. The resulting standardized trajectory is transformed back
# to the analyzed data units before it is returned.
#
# Geometry is applied before either branch. "pooled" uses one mean and covariance
# over all series-time observations; "timewise" uses a separate cross-sectional
# mean and covariance at each time; "none" leaves the data unchanged. If
# cumulative = TRUE, cumulative paths are formed before this metric transformation.
#
# Arguments:
#   X           : A numeric array of dimension c(N, T, d) (series, time,
#                 variables), or a list of N matrices each of dimension c(T, d).
#   u           : A length-d numeric vector giving one constant direction over
#                 time, or a T x d numeric matrix giving a time-varying
#                 direction. Required; no time point may have zero norm.
#                 Directions are normalized time-wise before computing the tilt
#                 trajectory (2 * tau - 1) * u_norm.
#   tau         : A scalar or numeric vector of quantile levels strictly inside
#                 (0, 1) (default 0.5, the median). Multiple values produce
#                 multiple returned directions, each using the same normalized u.
#   metric      : Underlying geometry: "pooled" (default, a single
#                 pooled-whitening Mahalanobis metric giving exact affine
#                 invariance), "timewise" (Mahalanobis using the per-time
#                 cross-sectional covariance), or "none" (raw Euclidean).
#   method      : Depth computation for the empirical branch: "exact" (default,
#                 exact C(j)) or "clara" (a CLARA subsample estimate, cheaper for
#                 large N).
#   sample.size : CLARA subsample size s. Defaults to
#                 min(N, ceiling(40 + 4 * sqrt(N))).
#   clara.draws : Number of CLARA subsamples to draw; the draw whose medoid has
#                 the smallest depth is kept (default 5).
#   iter, eps   : Maximum iterations and tolerance for the time-wise Weiszfeld
#                 iteration.
#   tol         : Numerical tolerance for tie-breaking.
#   cumulative  : Logical; if TRUE, replace every series-variable path with its
#                 cumulative sum over time before computing the metric and DGQs.
#                 Returned data and trajectories are also cumulative. Default FALSE.
#   n.cores     : Positive integer giving the requested number of OpenMP
#                 cores/threads for empirical candidate evaluation in both exact
#                 and CLARA modes. The default 1L preserves serial execution.
#                 Requests are limited to N cores, and builds without OpenMP
#                 support run serially.
#
# Value:
#   An object of class "DGQ": a list with components empirical (indices and an
#   m x T x d array of observed trajectories), timewise (an m x T x d array of
#   trajectories), directions, medoid, depth (the C(j) vector), metric, method,
#   sample.size and the analyzed data X (cumulative when cumulative = TRUE). See
#   print.DGQ() and plot.DGQ().
#
# Examples:
#   set.seed(1)
#   N <- 60; T <- 20; d <- 2
#   X <- array(rnorm(N * T * d), dim = c(N, T, d))
#   fit <- DGQ(X, u = c(1, 0))                        # dynamic median (tau = 0.5)
#   fit
#   fit2 <- DGQ(X, u = c(1, 0), tau = c(0.25, 0.75))  # tilted DGQs at multiple tau
#   fit2$empirical$index
#
# References:
#   Chaudhuri, P. (1996). On a geometric notion of quantiles for multivariate
#   data. J. Amer. Statist. Assoc. 91(434), 862--872.
#
#   Pena, D., Tsay, R. S., & Zamar, R. (2019). Empirical dynamic quantiles for
#   visualization of high-dimensional time series. Technometrics 61(4), 429--444.
DGQ <- function(X,
                u,
                tau         = 0.5,
                metric      = c("pooled", "timewise", "none"),
                method      = c("exact", "clara"),
                sample.size = NULL,
                clara.draws = 5L,
                iter        = 200L,
                eps         = 1e-8,
                tol         = 1e-9,
                cumulative  = FALSE,
                n.cores     = 1L) {
  ## Resolve the requested metric and depth implementations before
  ## constructing any derived arrays.
  metric <- match.arg(metric)
  method <- match.arg(method)

  if (!is.logical(cumulative) || length(cumulative) != 1L || is.na(cumulative))
    stop("'cumulative' must be a single non-missing logical value")
  if (!(is.integer(n.cores) || is.double(n.cores)) ||
      length(n.cores) != 1L || is.na(n.cores) ||
      !is.finite(n.cores) || n.cores < 1 || n.cores != floor(n.cores))
    stop("'n.cores' must be a single positive integer")

  ## Normalize all accepted inputs to X[i,t,k], then optionally replace each
  ## path X_i(1:T)[k] by its cumulative levels before computing geometry.
  X <- .as_dgq_array(X)
  NN <- dim(X)[1L]
  TT <- dim(X)[2L]
  dd <- dim(X)[3L]
  if (NN < 2L) stop("need at least N = 2 series")
  n.cores <- as.integer(min(n.cores, NN))
  if (cumulative) X <- .cumulative_series(X)

  ## U is a named list of T x d matrices of tilt trajectories (one per tau level).
  U <- .build_directions(u, tau, dd, TT, tol)
  m <- length(U)

  ## Standardization makes the C core's Euclidean norm represent the selected
  ## raw, pooled-Mahalanobis, or time-wise-Mahalanobis geometry.
  std <- .apply_metric(X, metric)
  Z   <- std$Z

  ## ---- empirical depth and tilt terms ------------------------------------
  ## Compute C(j) and V(j) once because only the linear tilt N<v,V(j)> changes
  ## across directions. Exact mode references every series. CLARA estimates C
  ## from several reference subsamples and retains the draw with the smallest
  ## candidate-medoid depth.
  C       <- NULL
  V       <- NULL
  medoid  <- NULL
  smp     <- NULL
  if (method == "exact") {
    ref <- seq_len(NN) - 1L
    pr  <- .Call("r_dgq_empirical", Z, as.integer(ref), n.cores,
                 PACKAGE = "DGQ")
    C <- pr$C
    V <- pr$V
  } else {
    s <- if (is.null(sample.size))
      min(NN, ceiling(40 + 4 * sqrt(NN))) else min(NN, as.integer(sample.size))
    smp      <- s
    best_obj <- Inf
    for (draw in seq_len(clara.draws)) {
      ref <- sort(sample.int(NN, s)) - 1L
      pr  <- .Call("r_dgq_empirical", Z, as.integer(ref), n.cores,
                   PACKAGE = "DGQ")
      obj <- min(pr$C)            # depth of the candidate medoid
      if (obj < best_obj) {
        best_obj <- obj
        C        <- pr$C
        V        <- pr$V
      }
    }
  }
  ## u = 0 removes the tilt, so this index minimizes C(j).
  medoid <- .emp_index(C, V, matrix(0, TT, dd), NN, tol)

  ## ---- direction-specific constrained and unconstrained results -----------
  ## The empirical branch selects an observed path in original analyzed units.
  ## The time-wise branch solves in standardized space and is re-transformed.
  emp_idx  <- integer(m)
  emp_traj <- array(NA_real_, dim = c(m, TT, dd))
  tw_traj  <- array(NA_real_, dim = c(m, TT, dd))

  for (g in seq_len(m)) {
    ug <- U[[g]]
    idx <- .emp_index(C, V, ug, NN, tol)
    emp_idx[g]      <- idx
    emp_traj[g, , ] <- X[idx, , ]
    out <- .Call("r_dgq_timewise", Z, as.double(ug),
                 as.integer(iter), as.double(eps), PACKAGE = "DGQ")
    tw_traj[g, , ] <- .retransform(out, metric, std$mu, std$Winv)
  }

  ## Package direction-indexed outputs into consistently named m-by-T-by-d
  ## arrays, then retain the data and metadata needed by S3 methods.
  dirnames <- names(U)
  dimnames(emp_traj) <- list(dirnames, NULL, NULL)
  empirical <- list(
    index      = stats::setNames(emp_idx, dirnames), 
    trajectory = emp_traj
  )
  dimnames(tw_traj) <- list(dirnames, NULL, NULL)

  structure(
    list(
      call        = match.call(),
      metric      = metric,
      method      = method,
      sample.size = smp,
      cumulative  = cumulative,
      directions  = U,
      N           = NN,
      T           = TT,
      d           = dd,
      empirical   = empirical, 
      timewise    = tw_traj,
      medoid      = medoid, 
      depth       = C, 
      X           = X),
    class = "DGQ")
}
