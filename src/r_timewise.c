/* r_timewise.c -- R/SEXP boundary for the time-wise DGQ computation.
 *
 * This wrapper validates R-level dimensions, coerces inputs to the primitive
 * storage types expected by dgq.h, allocates the R-owned output, and protects
 * every allocated/coerced SEXP across subsequent allocations.
 */
#include <R.h>
#include <Rinternals.h>

#include "dgq.h"

/* R wrapper for compute_timewise().
 *
 * .Call contract:
 *   Z      numeric-coercible standardized array with dim c(N,T,d);
 *   u      numeric-coercible length-d tilt direction;
 *   iter   integer-coercible maximum Weiszfeld iteration count;
 *   eps    numeric-coercible convergence tolerance and distance guard.
 * Returns the unconstrained time-wise trajectory as a numeric T-by-d matrix.
 *
 * R-side DGQ() validates the direction norm and supplies normal positive
 * iteration controls. This boundary additionally checks the shape needed to
 * prevent an out-of-bounds read of u.
 */
SEXP r_dgq_timewise(SEXP Z, SEXP u, SEXP iter, SEXP eps)
{
    SEXP dim = getAttrib(Z, R_DimSymbol);
    if (dim == R_NilValue || LENGTH(dim) != 3)
        error("'Z' must be a 3-dimensional array with dim (N, T, d)");
    int N = INTEGER(dim)[0];
    int T = INTEGER(dim)[1];
    int d = INTEGER(dim)[2];

    SEXP Zr = PROTECT(coerceVector(Z, REALSXP));
    SEXP ur = PROTECT(coerceVector(u, REALSXP));
    if (LENGTH(ur) != T * d)
        error("length of 'u' (%d) must equal T * d (%d)", LENGTH(ur), T * d);
    int it = asInteger(iter);
    double ep = asReal(eps);

    /* The pure C output layout is already R's column-major T-by-d layout. */
    SEXP out = PROTECT(allocMatrix(REALSXP, T, d));
    compute_timewise(REAL(Zr), N, T, d, REAL(ur), it, ep, REAL(out));

    UNPROTECT(3);
    return out;
}
