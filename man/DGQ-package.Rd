\docType{package}
\name{DGQ-package}
\alias{DGQ-package}
\title{DGQ: Dynamic Geometric Quantiles for Multivariate Time Series}
\description{
The package computes the Dynamic Geometric Quantile (DGQ) for a collection of
multivariate time series. The single entry point \code{\link[=DGQ]{DGQ()}} returns both the
unconstrained \emph{time-wise geometric quantile} (a free trajectory, the
population target) and the constrained \emph{empirical DGQ} (a real observed
series) for one or more tilt directions.
}
\details{
The method extends the empirical dynamic quantile of Pena, Tsay and Zamar
(2019) to the multivariate setting via Chaudhuri's (1996) geometric quantile.
The numerically heavy core (the pairwise depth and the Weiszfeld iteration) is
implemented in C and reached through the R C interface.
}
\references{
Chaudhuri, P. (1996). On a geometric notion of quantiles for multivariate
data. \emph{J. Amer. Statist. Assoc.} 91(434), 862--872.

Pena, D., Tsay, R. S., & Zamar, R. (2019). Empirical dynamic quantiles for
visualization of high-dimensional time series. \emph{Technometrics} 61(4),
429--444.

Kaufman, L. & Rousseeuw, P. J. (1990). \emph{Finding Groups in Data}. Wiley
(program CLARA).
}
\author{
\strong{Maintainer}: Lizhong Chen \email{nealforeverfriend@gmail.com}

Authors:
\itemize{
  \item Lizhong Chen \email{nealforeverfriend@gmail.com}
  \item Xueming Li
}

}
\keyword{internal}
