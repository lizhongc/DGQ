\name{DGQ}
\alias{DGQ}
\alias{print.DGQ}
\alias{plot.DGQ}
\title{Dynamic Geometric Quantile of a collection of multivariate time series}
\usage{
DGQ(
  X,
  u,
  tau = 0.5,
  metric = c("pooled", "timewise", "none"),
  method = c("exact", "clara"),
  sample.size = NULL,
  clara.draws = 5L,
  iter = 200L,
  eps = 1e-08,
  tol = 1e-09,
  cumulative = FALSE,
  n.cores = 1L
)

\method{print}{DGQ}(x, ...)

\method{plot}{DGQ}(x, components = NULL, type = c("trajectory", "rainbow", "fan"), ...)
}
\arguments{
\item{X}{A numeric array of dimension \code{c(N, T, d)} (series, time, variables),
or a list of \code{N} matrices each of dimension \code{c(T, d)}.}

\item{u}{A length-\code{d} numeric vector giving one constant direction over time,
or a \verb{T x d} numeric matrix giving a time-varying direction. Required; no
time point may have zero norm. Directions are normalized time-wise before
computing the tilt trajectory \code{(2 * tau - 1) * u_norm}.}

\item{tau}{A scalar or numeric vector of quantile levels strictly inside
\verb{(0, 1)} (default 0.5, representing the median). Multiple values produce
multiple returned directions, each using the same normalized \code{u}.}

\item{metric}{Metric used for the underlying geometry: \code{"pooled"} (default,
a single pooled-whitening Mahalanobis metric giving exact affine
invariance), \code{"timewise"} (Mahalanobis using the per-time cross-sectional
covariance), or \code{"none"} (raw Euclidean).}

\item{method}{Depth computation for the empirical branch: \code{"exact"}
(default, exact \code{C(j)}) or \code{"clara"} (a CLARA subsample estimate, cheaper
for large \code{N}).}

\item{sample.size}{CLARA subsample size \code{s}. Defaults to
\code{min(N, ceiling(40 + 4 * sqrt(N)))}.}

\item{clara.draws}{Number of CLARA subsamples to draw; the draw whose medoid
has the smallest depth is kept (default 5).}

\item{iter, eps}{Maximum iterations and tolerance for the time-wise Weiszfeld
iteration.}

\item{tol}{Numerical tolerance for tie-breaking.}

\item{cumulative}{Logical; if \code{TRUE}, replace every series-variable path with
its cumulative sum over time before computing the metric and DGQs. Returned
data and trajectories are also cumulative. Defaults to \code{FALSE}.}

\item{n.cores}{Positive integer giving the requested number of OpenMP cores/threads
for empirical candidate evaluation in both exact and CLARA modes. The
default \code{1L} preserves serial execution. Requests are limited to \code{N}
cores, and builds without OpenMP support run serially.}

\item{x}{An object of class \code{"DGQ"}.}

\item{...}{Further arguments (currently unused).}

\item{components}{Integer vector of variable indices to plot (default all).}

\item{type}{Character; the plot view, one of \code{"trajectory"} (default),
\code{"rainbow"}, or \code{"fan"}. See \emph{Methods (by generic)} below.}
}
\value{
An object of class \code{"DGQ"}: a list with components \code{empirical}
(indices and an \verb{m x T x d} array of observed trajectories),
\code{timewise} (an \verb{m x T x d} array of trajectories), \code{directions},
\code{medoid}, \code{depth} (the \code{C(j)} vector), \code{metric}, \code{method}, \code{sample.size}
and the analyzed data \code{X} (cumulative when \code{cumulative = TRUE}). See
\code{\link[=print.DGQ]{print.DGQ()}} and \code{\link[=plot.DGQ]{plot.DGQ()}}.
}
\description{
Computes the Dynamic Geometric Quantile (DGQ) for a collection of \code{N}
multivariate trajectories. For each quantile level in \code{tau}, the function
builds a tilt trajectory in the open unit ball and returns two objects:
}
\details{
\itemize{
\item the \strong{empirical DGQ} (the \emph{constrained} estimator):
the real observed series minimising \verb{C(j) + N <u, V(j)>}; and
\item the \strong{time-wise geometric quantile} (the
\emph{unconstrained} population target): a free trajectory whose value at each
time \code{t} is Chaudhuri's geometric \code{u}-quantile of the time-\code{t} marginal.
}

Both are returned.

Let \code{Z_i(t)} denote trajectory \code{i} at time \code{t} after applying the selected
metric transformation, and let \code{Zbar(t)} be its cross-sectional mean. The
empirical branch computes

\verb{C(j) = sum_i sum_t ||Z_i(t) - Z_j(t)||}

and

\verb{V(j) = sum_t \{Zbar(t) - Z_j(t)\}}.

For each direction \code{u}, it returns the observed trajectory whose index
minimizes \verb{C(j) + N <u, V(j)>}. At \code{u = 0}, this is the dynamic medoid that
minimizes total trajectory distance. With \code{method = "clara"}, the sum over
\code{i} in \code{C(j)} is estimated from reference subsamples and scaled by \code{N/s};
\code{V(j)} remains exact. Empirical candidate calculations can be distributed
over OpenMP threads with \code{n.cores > 1}; builds without OpenMP support ignore
the request and execute serially.

Independently at each time \code{t}, the time-wise branch minimizes

\verb{sum_i ||Z_i(t) - q|| - N <u, q>}

over unrestricted \code{q}. Its first-order equation is
\verb{sum_i (q - Z_i(t)) / ||q - Z_i(t)|| = N u}, which is solved by a modified
Weiszfeld iteration. The resulting standardized trajectory is transformed
back to the analyzed data units before it is returned.

Geometry is applied before either branch. \code{"pooled"} uses one mean and
covariance over all series-time observations; \code{"timewise"} uses a separate
cross-sectional mean and covariance at each time; \code{"none"} leaves the data
unchanged. If \code{cumulative = TRUE}, cumulative paths are formed before this
metric transformation.
}
\section{Methods (by generic)}{
\itemize{
\item \code{print(DGQ)}: Print a concise summary of a fitted DGQ.

Reports the analyzed array dimensions, selected computational options,
dynamic medoid, constrained empirical selections, and number of
unconstrained time-wise trajectories. The fitted object is returned
invisibly and its trajectory arrays are not printed.

\item \code{plot(DGQ)}: Visualize a fitted DGQ, one panel per requested
variable, with the view chosen by \code{type}:
\itemize{
\item \code{"trajectory"} (default): all \code{N} analyzed trajectories in grey,
with direction-specific empirical trajectories overlaid as solid lines and
time-wise trajectories as dashed lines in a shared direction color.
\item \code{"rainbow"}: every series colored by its trajectory depth \code{C(j)}
(deep = central), the central-50\% band shaded, and the dynamic medoid
overlaid (plus the empirical median when a \code{tau = 0.5} direction is present).
\item \code{"fan"}: the time-wise quantile fan (nested bands between symmetric
\code{tau} pairs plus the median) with the empirical per-\code{tau} selections
overlaid as lines; requires a fit with at least two \code{tau} levels.
}
The fitted object is returned invisibly.

}}
\examples{
set.seed(1)
N <- 60; T <- 20; d <- 2
X <- array(rnorm(N * T * d), dim = c(N, T, d))
fit <- DGQ(X, u = c(1, 0))                   # dynamic median (tau = 0.5)
fit
fit2 <- DGQ(X, u = c(1, 0), tau = c(0.25, 0.75)) # tilted DGQs at multiple tau levels
fit2$empirical$index

}
\references{
Chaudhuri, P. (1996). On a geometric notion of quantiles for multivariate
data. \emph{J. Amer. Statist. Assoc.} 91(434), 862--872.

Pena, D., Tsay, R. S., & Zamar, R. (2019). Empirical dynamic quantiles for
visualization of high-dimensional time series. \emph{Technometrics} 61(4),
429--444.
}
