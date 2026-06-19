# DGQ 1.0.0

## First public release

- Added `DGQ()`, the main entry point for Dynamic Geometric Quantiles of
  multivariate time series.
- Added the constrained empirical DGQ, which selects an observed trajectory for
  each requested quantile level.
- Added the unconstrained time-wise geometric quantile, computed by a modified
  Weiszfeld iteration.
- Added raw Euclidean, pooled-whitening Mahalanobis, and time-wise Mahalanobis
  metrics.
- Added exact empirical depth and CLARA-subsampled empirical depth modes.
- Added cumulative-series analysis for increment or flow data.
- Added optional OpenMP candidate-level parallelism for the empirical branch.
- Added S3 `print()` and `plot()` methods.
- Added tests comparing the C core with pure-R reference implementations and
  checking invariance, input formats, cumulative mode, CLARA mode, and OpenMP
  consistency.
