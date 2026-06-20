# Pure-R reference implementations used to cross-check the C core ------------

ref_depth <- function(Z, ref0) {            # ref0: 0-based reference indices
  N <- dim(Z)[1]; T <- dim(Z)[2]; d <- dim(Z)[3]
  refi <- ref0 + 1L
  C <- numeric(N)
  for (j in seq_len(N)) {
    cj <- 0
    for (i in refi)
      for (t in seq_len(T)) cj <- cj + sqrt(sum((Z[i, t, ] - Z[j, t, ])^2))
    C[j] <- (N / length(refi)) * cj
  }
  Zbar <- apply(Z, c(2, 3), mean)
  V <- array(0, dim = c(N, T, d))
  for (j in seq_len(N)) {
    for (t in seq_len(T)) {
      V[j, t, ] <- Zbar[t, ] - Z[j, t, ]
    }
  }
  list(C = C, V = V)
}

ref_tw <- function(Z, u, iter = 200, eps = 1e-8) {
  N <- dim(Z)[1]; T <- dim(Z)[2]; d <- dim(Z)[3]
  out <- matrix(0, T, d)
  for (t in seq_len(T)) {
    ut <- u[t, ]
    P <- matrix(Z[, t, ], N, d); Q <- colMeans(P)
    for (k in seq_len(iter)) {
      dv <- sqrt(rowSums(sweep(P, 2, Q, "-")^2)); dv[dv < eps] <- eps
      Qn <- (colSums(P * (1 / dv)) + N * ut) / sum(1 / dv)
      if (sqrt(sum((Qn - Q)^2)) < eps) { Q <- Qn; break }
      Q <- Qn
    }
    out[t, ] <- Q
  }
  out
}

# 1. d = 1 reduction to the EDQ ---------------------------------------------

test_that("d=1: empirical DGQ at u = 2p-1 is monotone in p (EDQ reduction)", {
  set.seed(99)
  levels_i <- sort(rnorm(120))
  X1 <- array(0, dim = c(120, 40, 1))
  for (i in seq_len(120)) X1[i, , 1] <- levels_i[i] + rnorm(40, sd = 0.3)
  ps  <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  sel <- vapply(ps, function(p)
    DGQ(X1, u = 1, tau = p, metric = "none")$empirical$index,
    integer(1))
  meanlev <- vapply(sel, function(j) mean(X1[j, , 1]), numeric(1))
  expect_true(all(diff(meanlev) > 0))
})

# 2. medoid equals the depth minimiser --------------------------------------

test_that("u = 0 selects the depth-minimising medoid", {
  set.seed(1)
  X <- array(rnorm(50 * 10 * 2), dim = c(50, 10, 2))
  f <- DGQ(X, u = c(1, 0), metric = "none")
  expect_equal(f$medoid, which.min(f$depth))
  expect_equal(unname(DGQ(X, u = c(1, 0), tau = 0.5,
                           metric = "none")$empirical$index), f$medoid)
})

# 3. C core matches the pure-R reference ------------------------------------

test_that("r_dgq_empirical and r_dgq_timewise match the R reference", {
  set.seed(2)
  N <- 30; T <- 8; d <- 2
  Z <- array(rnorm(N * T * d), dim = c(N, T, d))
  ref0 <- 0:(N - 1)
  pr <- .Call("r_dgq_empirical", Z, as.integer(ref0), 1L, PACKAGE = "DGQ")
  rr <- ref_depth(Z, ref0)
  expect_equal(pr$C, rr$C, tolerance = 1e-10)
  expect_equal(pr$V, rr$V, tolerance = 1e-10)

  u  <- matrix(rnorm(T * d), T, d)
  tw <- .Call("r_dgq_timewise", Z, as.double(u), 200L, 1e-8, PACKAGE = "DGQ")
  expect_equal(tw, ref_tw(Z, u), tolerance = 1e-6)

  ## CLARA scaling: a subsample reference reproduces ref_depth with N/s factor
  set.seed(10); ref0s <- sort(sample.int(N, 12)) - 1L
  prs <- .Call("r_dgq_empirical", Z, as.integer(ref0s), 1L, PACKAGE = "DGQ")
  expect_equal(prs$C, ref_depth(Z, ref0s)$C, tolerance = 1e-10)
})

# 4. invariance / equivariance ----------------------------------------------

test_that("medoid is affine-invariant under pooled and timewise whitening", {
  set.seed(3)
  N <- 60; T <- 12; d <- 2
  X <- array(rnorm(N * T * d), dim = c(N, T, d))

  A <- matrix(c(2, 0.5, -0.3, 1.5), 2, 2); b <- c(1, -2)
  Xa <- X
  for (t in seq_len(T)) Xa[, t, ] <- X[, t, ] %*% t(A) + matrix(b, N, d, byrow = TRUE)
  expect_equal(DGQ(Xa, u = c(1, 0), metric = "pooled")$medoid, DGQ(X, u = c(1, 0), metric = "pooled")$medoid)

  ## per-time affine maps -> timewise Mahalanobis must be invariant
  Xt <- X
  for (t in seq_len(T)) {
    At <- matrix(rnorm(4), 2, 2); At <- At + diag(2) * 2   # well-conditioned
    Xt[, t, ] <- X[, t, ] %*% t(At) + matrix(rnorm(d), N, d, byrow = TRUE)
  }
  expect_equal(DGQ(Xt, u = c(1, 0), metric = "timewise")$medoid,
               DGQ(X, u = c(1, 0), metric = "timewise")$medoid)
})

test_that("orthogonal equivariance: direction Au gives the same index", {
  set.seed(4)
  N <- 60; T <- 10; d <- 2
  X <- array(rnorm(N * T * d), dim = c(N, T, d))
  th <- 0.7
  A <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2, 2)  # rotation
  Xr <- X
  for (t in seq_len(T)) Xr[, t, ] <- X[, t, ] %*% t(A)
  u <- c(0.3, -0.4)
  i1 <- DGQ(X,  u = u,                metric = "none")$empirical$index
  i2 <- DGQ(Xr, u = as.vector(A %*% u), metric = "none")$empirical$index
  expect_equal(unname(i1), unname(i2))
})

# 5. array and list inputs agree --------------------------------------------

test_that("array and list-of-matrices inputs give identical results", {
  set.seed(5)
  N <- 40; T <- 9; d <- 2
  X <- array(rnorm(N * T * d), dim = c(N, T, d))
  Xl <- lapply(seq_len(N), function(i) X[i, , ])
  dirs <- c(0.4, 0.2)
  fa <- DGQ(X,  u = dirs, tau = 0.7)
  fl <- DGQ(Xl, u = dirs, tau = 0.7)
  expect_equal(fa$empirical$index, fl$empirical$index)
  expect_equal(fa$timewise, fl$timewise, tolerance = 1e-10)
})

# 6. cumulative series -------------------------------------------------------

test_that("cumulative mode matches manually cumulative input", {
  set.seed(6)
  N <- 35; T <- 7; d <- 2
  X <- array(rnorm(N * T * d), dim = c(N, T, d))
  Xcum <- X
  for (i in seq_len(N))
    for (k in seq_len(d))
      Xcum[i, , k] <- cumsum(X[i, , k])

  dirs <- c(0.3, 0.2)
  fc <- DGQ(X, u = dirs, tau = 0.7, metric = "none", cumulative = TRUE)
  fm <- DGQ(Xcum, u = dirs, tau = 0.7, metric = "none", cumulative = FALSE)

  expect_true(fc$cumulative)
  expect_equal(fc$X, Xcum)
  expect_equal(fc$empirical$index, fm$empirical$index)
  expect_equal(fc$depth, fm$depth, tolerance = 1e-10)
  expect_equal(fc$empirical$trajectory, fm$empirical$trajectory)
  expect_equal(fc$timewise, fm$timewise, tolerance = 1e-10)
  expect_output(print(fc), "cumulative: yes")

  f_default <- DGQ(X, u = dirs, tau = 0.7, metric = "none")
  f_false <- DGQ(X, u = dirs, tau = 0.7, metric = "none", cumulative = FALSE)
  expect_false(f_default$cumulative)
  expect_equal(f_default$empirical$index, f_false$empirical$index)
  expect_equal(f_default$depth, f_false$depth)
  expect_equal(f_default$timewise, f_false$timewise)
})

# 7. constrained -> unconstrained as N grows (targeting, section 2) ----------

test_that("empirical trajectory approaches the time-wise target as N grows", {
  gap <- function(N, seed) {
    set.seed(seed)
    X <- array(rnorm(N * 6 * 2), dim = c(N, 6, 2))
    f <- DGQ(X, u = c(0.3, 0.1))
    emp <- f$empirical$trajectory[1, , ]; tw <- f$timewise[1, , ]
    mean(sqrt(rowSums((emp - tw)^2)))
  }
  gsmall <- mean(vapply(1:6, function(s) gap(25,  s), numeric(1)))
  glarge <- mean(vapply(1:6, function(s) gap(800, s), numeric(1)))
  expect_lt(glarge, gsmall)
})

# 8. CLARA approximation -----------------------------------------------------

test_that("CLARA recovers exact depth/medoid for large subsamples", {
  set.seed(7)
  N <- 150; T <- 10; d <- 2
  X <- array(rnorm(N * T * d), dim = c(N, T, d))
  fe <- DGQ(X, u = c(1, 0), metric = "none")

  set.seed(123)
  fc <- DGQ(X, u = c(1, 0), method = "clara", sample.size = N, metric = "none")  # s = N => exact
  expect_equal(fc$medoid, fe$medoid)
  expect_equal(fc$depth, fe$depth, tolerance = 1e-8)

  set.seed(5)
  fc2 <- DGQ(X, u = c(1, 0), method = "clara", sample.size = 120, clara.draws = 8, metric = "none")
  expect_equal(fc2$medoid, fe$medoid)
})

# 9. empirical OpenMP cores preserve exact and CLARA results --------------

test_that("empirical core counts preserve exact and CLARA results", {
  set.seed(81)
  N <- 70; T <- 12; d <- 3
  X <- array(rnorm(N * T * d), dim = c(N, T, d))
  dirs <- c(0.3, 0.1, 0.2)

  f_default <- DGQ(X, u = dirs, tau = 0.7, metric = "none")
  f1 <- DGQ(X, u = dirs, tau = 0.7, metric = "none",
             n.cores = 1L)
  for (nth in c(2L, 4L, N + 10L)) {
    fp <- DGQ(X, u = dirs, tau = 0.7, metric = "none",
               n.cores = nth)
    expect_identical(fp$depth, f1$depth)
    expect_identical(fp$empirical$index, f1$empirical$index)
    expect_identical(fp$empirical$trajectory, f1$empirical$trajectory)
  }
  expect_identical(f_default$depth, f1$depth)
  expect_identical(f_default$empirical$index, f1$empirical$index)

  for (nth in c(1L, 2L, 4L, N + 10L)) {
    set.seed(812)
    fc <- DGQ(X, u = dirs, tau = 0.7, metric = "none",
               method = "clara", sample.size = 35L, clara.draws = 4L,
               n.cores = nth)
    if (nth == 1L) {
      fc1 <- fc
    } else {
      expect_identical(fc$depth, fc1$depth)
      expect_identical(fc$empirical$index, fc1$empirical$index)
      expect_identical(fc$empirical$trajectory, fc1$empirical$trajectory)
    }
  }
})

test_that("native empirical wrapper validates and uses core count", {
  set.seed(82)
  N <- 25; T <- 6; d <- 2
  Z <- array(rnorm(N * T * d), dim = c(N, T, d))
  ref0 <- 0:(N - 1)
  p1 <- .Call("r_dgq_empirical", Z, as.integer(ref0), 1L, PACKAGE = "DGQ")
  p4 <- .Call("r_dgq_empirical", Z, as.integer(ref0), 4L, PACKAGE = "DGQ")
  expect_identical(p4, p1)
  expect_error(.Call("r_dgq_empirical", Z, as.integer(ref0), 0L,
                     PACKAGE = "DGQ"), "positive integer")
  expect_error(.Call("r_dgq_empirical", Z, as.integer(ref0), 2,
                     PACKAGE = "DGQ"), "positive integer")
})

# 10. argument validation and S3 methods ------------------------------------

test_that("invalid directions and basic S3 methods behave", {
  set.seed(8)
  X <- array(rnorm(40 * 8 * 2), dim = c(40, 8, 2))
  expect_error(DGQ(X, u = c(0, 0)), "zero vector")
  expect_error(DGQ(X, u = NULL, tau = 0.7), "cannot be NULL")
  expect_error(DGQ(X),                      "argument \"u\" is missing")
  expect_error(DGQ(X, tau = 0), "open interval")
  expect_error(DGQ(X, tau = 1.2), "open interval")
  expect_error(DGQ(X, u = c(0.5, 0.5, 0.5)), "length")
  expect_error(DGQ(X, cumulative = NA), "single non-missing logical")
  expect_error(DGQ(X, cumulative = c(TRUE, FALSE)), "single non-missing logical")
  expect_error(DGQ(X, cumulative = 1), "single non-missing logical")
  expect_error(DGQ(X, n.cores = 0), "single positive integer")
  expect_error(DGQ(X, n.cores = -1), "single positive integer")
  expect_error(DGQ(X, n.cores = 1.5), "single positive integer")
  expect_error(DGQ(X, n.cores = NA), "single positive integer")
  expect_error(DGQ(X, n.cores = "2"), "single positive integer")
  expect_error(DGQ(X, n.cores = 1 + 0i), "single positive integer")

  tw_default <- DGQ(X, u = c(1, 0), metric = "none")
  tw_cores <- DGQ(X, u = c(1, 0), metric = "none", n.cores = 4L)
  expect_identical(tw_cores$timewise, tw_default$timewise)

  f <- DGQ(X, u = c(1, 0), tau = 0.7)
  expect_output(print(f), "Dynamic Geometric Quantile")
  pf <- tempfile(fileext = ".pdf"); pdf(pf)
  expect_invisible(plot(f)); dev.off(); unlink(pf)
})

test_that("plot() rainbow and fan types render", {
  set.seed(11)
  X <- array(rnorm(40 * 8 * 2), dim = c(40, 8, 2))
  f1 <- DGQ(X, u = c(1, 0), tau = 0.5)
  f5 <- DGQ(X, u = c(1, 0), tau = c(0.1, 0.25, 0.5, 0.75, 0.9))
  pf <- tempfile(fileext = ".pdf"); pdf(pf)
  expect_invisible(plot(f1, type = "trajectory"))
  expect_invisible(plot(f1, type = "rainbow"))
  expect_invisible(plot(f5, type = "fan"))
  expect_invisible(plot(f5, type = "fan", components = 1))
  dev.off(); unlink(pf)
  expect_error(plot(f1, type = "fan"), "at least 2 tau")
  expect_error(plot(f1, type = "nope"))   # match.arg rejects unknown type
})

# 11. time-varying directions verification -----------------------------------

test_that("T x d matrix directions work and equal replicated vector u", {
  set.seed(99)
  N <- 50; T <- 8; d <- 2
  X <- array(rnorm(N * T * d), dim = c(N, T, d))
  u_vec <- c(0.6, -0.8)
  u_mat <- matrix(u_vec, T, d, byrow = TRUE)

  f_vec <- DGQ(X, u = u_vec, tau = 0.75, metric = "pooled")
  f_mat <- DGQ(X, u = u_mat, tau = 0.75, metric = "pooled")

  expect_equal(f_vec$empirical$index, f_mat$empirical$index)
  expect_equal(f_vec$timewise, f_mat$timewise, tolerance = 1e-10)

  # Test time-varying directions that change over time
  u_var <- matrix(rnorm(T * d), T, d)
  f_var <- DGQ(X, u = u_var, tau = 0.8)
  expect_equal(length(f_var$directions), 1L)
  expect_equal(dim(f_var$directions[[1]]), c(T, d))

  # Test error when matrix has wrong dimensions
  expect_error(DGQ(X, u = matrix(1, T - 1, d)), "must have size")
  expect_error(DGQ(X, u = matrix(1, T, d + 1)), "must have size")
})
