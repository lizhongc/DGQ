# DGQ

Dynamic Geometric Quantiles for multivariate time series.

`DGQ` computes a quantile-like trajectory for a collection of multivariate time
series. The main function, `DGQ()`, returns both a constrained empirical DGQ
that is an observed series and an unconstrained time-wise geometric quantile.

## Installation

Install from a local source checkout with:

```r
install.packages("DGQ", repos = NULL, type = "source")
```

After CRAN acceptance, install with:

```r
install.packages("DGQ")
```

For Git hosting, replace `lizhongc/DGQ` with the repository location:

```r
# install.packages("remotes")
remotes::install_github("lizhongc/DGQ")
```

## Example

```r
library(DGQ)

set.seed(1)
N  <- 60
TT <- 20
d  <- 2
X  <- array(rnorm(N * TT * d), dim = c(N, TT, d))

fit <- DGQ(X, u = c(1, 0), tau = c(0.25, 0.5, 0.75))
fit
fit$empirical$index

plot(fit)
```

## Main Options

- `metric = "pooled"` uses one pooled-whitening Mahalanobis metric. This is the default.
- `metric = "timewise"` uses a separate cross-sectional Mahalanobis metric at each time point.
- `metric = "none"` uses raw Euclidean geometry.
- `method = "exact"` computes exact empirical depth.
- `method = "clara"` uses CLARA-style reference subsamples for cheaper approximate empirical depth when `N` is large.
- `cumulative = TRUE` analyzes cumulative paths, useful for increment or flow data.
- `n.cores` requests OpenMP threads for empirical candidate evaluation. Builds without OpenMP support run serially.

## References

Chaudhuri, P. (1996). On a geometric notion of quantiles for multivariate data. *Journal of the American Statistical Association*, 91(434), 862-872.
https://doi.org/10.2307/2291681

Pena, D., Tsay, R. S., & Zamar, R. (2019). Empirical dynamic quantiles for visualization of high-dimensional time series. *Technometrics*, 61(4), 429-444. 
https://doi.org/10.1080/00401706.2019.1575285
