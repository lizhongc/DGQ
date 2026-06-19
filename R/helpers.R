# ============================================================================
#  DGQ internal helpers.
# ============================================================================

## Return the symmetric inverse square root of a covariance matrix.
##
## For symmetric S = Q diag(lambda) Q', this computes
##   S^(-1/2) = Q diag(max(lambda, 1e-12)^(-1/2)) Q'.
## Multiplying a centered row vector by this matrix whitens it. Flooring the
## eigenvalues regularizes singular or nearly singular sample covariances and
## prevents division by zero.
.inv_sqrt <- function(S) {
  e <- eigen(S, symmetric = TRUE)
  V <- e$vectors
  dd <- pmax(e$values, 1e-12)
  V %*% diag(1 / sqrt(dd), length(dd)) %*% t(V)
}

## Return the regularized symmetric square root of a covariance matrix.
##
## Using the same eigensystem and eigenvalue floor as .inv_sqrt(), this computes
##   S^(1/2) = Q diag(max(lambda, 1e-12)^(1/2)) Q'.
## It is the inverse-whitening map used to return time-wise quantiles from
## standardized coordinates to analyzed data units.
.sqrt_mat <- function(S) {
  e <- eigen(S, symmetric = TRUE)
  V <- e$vectors
  dd <- pmax(e$values, 1e-12)
  V %*% diag(sqrt(dd), length(dd)) %*% t(V)
}

## Normalize user input to the package's canonical numeric array representation.
##
## The output always has dimensions (N,T,d), where X[i,t,k] is variable k of
## series i at time t, and double storage compatible with the native routines.
## A list input is interpreted as N equally sized T-by-d trajectory matrices;
## assigning each matrix to A[i,,] preserves the same mathematical indexing.
.as_dgq_array <- function(X) {
  if (is.array(X) && length(dim(X)) == 3L) {
    storage.mode(X) <- "double"
    return(X)
  }
  if (is.list(X)) {
    if (length(X) < 1L) stop("'X' is an empty list")
    mats <- lapply(X, as.matrix)
    Td <- dim(mats[[1L]])
    ## The anonymous predicate verifies the shared T-by-d trajectory shape.
    ok <- vapply(mats, function(m) identical(dim(m), Td), logical(1))
    if (!all(ok))
      stop("all matrices in 'X' must share the same dimensions (T x d)")
    NN <- length(mats); TT <- Td[1L]; dd <- Td[2L]
    A <- array(0, dim = c(NN, TT, dd))
    for (i in seq_len(NN)) A[i, , ] <- mats[[i]]
    storage.mode(A) <- "double"
    return(A)
  }
  stop("'X' must be a 3-d array [N, T, d] or a list of T x d matrices")
}

## Convert increment/flow paths into cumulative-level paths.
##
## Independently for every series i and variable k, the returned array satisfies
##   out[i,t,k] = sum_{s=1}^t X[i,s,k].
## Dimensions and storage layout are unchanged, so all later metric and DGQ
## calculations operate on cumulative levels without special cases.
.cumulative_series <- function(X) {
  NN <- dim(X)[1L]; dd <- dim(X)[3L]
  out <- X
  for (i in seq_len(NN))
    for (k in seq_len(dd))
      out[i, , k] <- cumsum(X[i, , k])
  out
}

## Map analyzed data into the Euclidean coordinates used by the C core.
##
## For metric = "none", Z_i(t) = X_i(t). For metric = "pooled", one mean mu and
## covariance Sigma are estimated from all N*T observations and
##   Z_i(t) = {X_i(t) - mu} Sigma^(-1/2).
## For metric = "timewise", separate mu_t and Sigma_t are estimated from the N
## observations at each time and
##   Z_i(t) = {X_i(t) - mu_t} Sigma_t^(-1/2).
##
## The return value contains Z plus mu and Winv = Sigma^(1/2), in pooled or
## time-indexed form, so .retransform() can invert the map for unconstrained
## time-wise trajectories. Empirical trajectories are selected directly from X.
.apply_metric <- function(X, metric) {
  NN <- dim(X)[1L]; TT <- dim(X)[2L]; dd <- dim(X)[3L]
  if (metric == "none") {
    return(list(Z = X, mu = NULL, Winv = NULL))
  }
  if (metric == "pooled") {
    ## Flatten all series-time pairs into rows without changing variable order.
    flat <- matrix(aperm(X, c(2, 1, 3)), NN * TT, dd)
    mu  <- colMeans(flat)
    Sig <- stats::cov(flat)
    W    <- .inv_sqrt(Sig)
    Winv <- .sqrt_mat(Sig)
    Z <- array(0, dim = c(NN, TT, dd))
    for (t in seq_len(TT))
      Z[, t, ] <- sweep(matrix(X[, t, ], NN, dd), 2, mu, "-") %*% W
    storage.mode(Z) <- "double"
    return(list(Z = Z, mu = mu, Winv = Winv))
  }
  if (metric == "timewise") {
    ## Each time point receives its own cross-sectional affine transformation.
    mu   <- matrix(0, TT, dd)
    Winv <- array(0, dim = c(dd, dd, TT))
    Z    <- array(0, dim = c(NN, TT, dd))
    for (t in seq_len(TT)) {
      Pt   <- matrix(X[, t, ], NN, dd)
      mut  <- colMeans(Pt)
      Sigt <- stats::cov(Pt)
      Z[, t, ]   <- sweep(Pt, 2, mut, "-") %*% .inv_sqrt(Sigt)
      mu[t, ]    <- mut
      Winv[, , t] <- .sqrt_mat(Sigt)
    }
    storage.mode(Z) <- "double"
    return(list(Z = Z, mu = mu, Winv = Winv))
  }
  stop("unknown metric")
}

## Invert the metric transformation for a standardized T-by-d trajectory.
##
## The inverse maps are q(t) Sigma^(1/2) + mu for a pooled metric and
## q(t) Sigma_t^(1/2) + mu_t for a time-wise metric. With metric = "none", the
## C result already uses analyzed data units and is returned unchanged.
.retransform <- function(out, metric, mu, Winv) {
  if (metric == "none") return(out)
  TT <- nrow(out); dd <- ncol(out)
  if (metric == "pooled")
    return(sweep(out %*% Winv, 2, mu, "+"))
  ## timewise
  res <- matrix(0, TT, dd)
  for (t in seq_len(TT)) res[t, ] <- as.vector(out[t, ] %*% Winv[, , t]) + mu[t, ]
  res
}

## Assemble and validate geometric-quantile tilt trajectories.
##
## A vector u supplies one fixed direction that is replicated over all TT time
## points. A matrix u supplies a time-varying direction and must already have
## shape TT-by-dd. Each time row is normalized separately, then each tau value
## scales the normalized trajectory by (2*tau - 1). List names label outputs.
.build_directions <- function(u, tau, dd, TT, tol) {
  if (is.null(tau) || !is.numeric(tau) || any(is.na(tau)) || any(tau <= 0) || any(tau >= 1)) {
    stop("tau must be numeric with values strictly in the open interval (0, 1)")
  }

  if (is.null(u)) {
    stop("u cannot be NULL")
  }

  if (is.null(dim(u))) {
    if (length(u) != dd) {
      stop(sprintf("tilt direction u vector must have length %d, got %d", dd, length(u)))
    }
    U_raw <- matrix(u, TT, dd, byrow = TRUE)
  } else {
    U_raw <- as.matrix(u)
    if (nrow(U_raw) != TT || ncol(U_raw) != dd) {
      stop(sprintf("tilt direction u matrix must have size %d x %d, got %d x %d",
                   TT, dd, nrow(U_raw), ncol(U_raw)))
    }
  }
  nrms <- sqrt(rowSums(U_raw^2))
  if (any(nrms < tol)) {
    stop("tilt direction u cannot be a zero vector at any time point (norm must be > 0)")
  }
  ## Always normalize to unit length time-wise.
  U_norm <- U_raw / nrms

  ## Construct one tilt trajectory per quantile level.
  U <- list()
  for (ti in tau) {
    U_val <- (2 * ti - 1) * U_norm
    rownames(U_val) <- rownames(u)
    label <- sprintf("tau=%.3f", ti)
    U[[label]] <- U_val
  }
  U
}

## Select the constrained empirical DGQ for one tilt direction.
##
## Candidate observed trajectory j has criterion
##   L_u(j) = C(j) + N <u, V(j)>.
## V is stored as an N-by-T-by-d array. We flatten V to N-by-(T*d) and
## u to a vector of length T*d, so V %*% u evaluates every inner product.
.emp_index <- function(C, V, u, N, tol) {
  dim(V) <- c(N, length(V) / N)
  crit <- C + N * as.vector(V %*% as.vector(u))
  mn   <- min(crit)
  best <- which(crit <= mn + tol * (abs(mn) + 1))
  if (length(best) > 1L) best[which.min(C[best])] else best[1L]
}
