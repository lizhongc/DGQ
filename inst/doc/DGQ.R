## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(comment = "#>", fig.width = 7, fig.height = 4)

## ----library------------------------------------------------------------------
library(DGQ)

## ----data---------------------------------------------------------------------
set.seed(1)
N <- 200; T <- 60; d <- 2
tt <- seq_len(T) / T
offsets <- cbind(rnorm(N, sd = 2), rnorm(N, sd = 1.4))
X <- array(0, dim = c(N, T, d))
for (i in seq_len(N)) {
  f1 <- runif(1, 1, 3); f2 <- runif(1, 1, 3)
  a  <- runif(1, 0.8, 1.6)
  osc <- cbind(a * sin(2 * pi * f1 * tt), a * cos(2 * pi * f2 * tt))
  X[i, , ] <- sweep(osc + matrix(rnorm(T * d, sd = 0.15), T, d),
                     2, offsets[i, ], "+")
}

## ----median-------------------------------------------------------------------
fit0 <- DGQ(X, u = c(1, 0))  # u = c(1, 0), pooled-whitening metric (the default)
fit0

## ----directions---------------------------------------------------------------
fit <- DGQ(X, u = c(1, 0), tau = c(0.15, 0.5, 0.85))
fit$empirical$index          # observed series representing each direction

## ----plot, fig.width=7, fig.height=4------------------------------------------
plot(fit)

## ----clara--------------------------------------------------------------------
fit_clara <- DGQ(X, u = c(1, 0), tau = c(0.15, 0.5, 0.85), method = "clara", sample.size = 80)
fit_clara$empirical$index

## ----input-formats------------------------------------------------------------
Xl <- lapply(seq_len(N), function(i) X[i, , ])
identical(DGQ(Xl, u = c(1, 0))$medoid, DGQ(X, u = c(1, 0))$medoid)

## ----cumulative---------------------------------------------------------------
fit_cumulative <- DGQ(X, u = c(1, 0), cumulative = TRUE)
fit_cumulative$cumulative

