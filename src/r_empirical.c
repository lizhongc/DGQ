/* r_empirical.c -- R/SEXP boundary for the empirical DGQ computation.
 *
 * This wrapper validates R-level dimensions, coerces inputs to the primitive
 * storage types expected by dgq.h, allocates R-owned outputs, and protects
 * every allocated/coerced SEXP across subsequent allocations.
 */
#include <R.h>
#include <Rinternals.h>

#include "dgq.h"

/* R wrapper for compute_empirical().
 *
 * .Call contract:
 *   Z     numeric-coercible standardized array with dim c(N,T,d);
 *   ref      integer-coercible vector of zero-based reference-series indices;
 *   n_cores  positive integer requested for candidate-level OpenMP execution.
 * Returns list(C = numeric(N), V = array(N,T,d)).
 *
 * The wrapper reads dimensions before coercion because coerceVector preserves
 * data but need not be relied upon to preserve array attributes. C and V are
 * allocated in the same layouts consumed by compute_empirical(), so no output
 * reshaping or copy is required. The T-by-d mean scratch buffer is allocated
 * on R's transient main-thread heap before entering any OpenMP region.
 */
SEXP r_dgq_empirical(SEXP Z, SEXP ref, SEXP n_cores)
{
    SEXP dim = getAttrib(Z, R_DimSymbol);
    if (dim == R_NilValue || LENGTH(dim) != 3)
        error("'Z' must be a 3-dimensional array with dim (N, T, d)");
    int N = INTEGER(dim)[0];
    int T = INTEGER(dim)[1];
    int d = INTEGER(dim)[2];

    /* Coercions may allocate; protect both converted inputs until return. */
    SEXP Zr = PROTECT(coerceVector(Z, REALSXP));
    SEXP refi = PROTECT(coerceVector(ref, INTSXP));
    int nref = LENGTH(refi);
    if (TYPEOF(n_cores) != INTSXP || LENGTH(n_cores) != 1 ||
        INTEGER(n_cores)[0] == NA_INTEGER || INTEGER(n_cores)[0] < 1)
        error("'n.cores' must be a single positive integer");
    int nthreads = INTEGER(n_cores)[0] < N ? INTEGER(n_cores)[0] : N;

    /* R owns these output buffers while the pure C routine fills them. */
    SEXP C = PROTECT(allocVector(REALSXP, N));
    SEXP V = PROTECT(alloc3DArray(REALSXP, N, T, d));
    double *mean_tk = R_Calloc((size_t) T * (size_t) d, double);

    compute_empirical(REAL(Zr), N, T, d, INTEGER(refi), nref,
                      mean_tk, nthreads, REAL(C), REAL(V));

    R_Free(mean_tk);

    SEXP out = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, C);
    SET_VECTOR_ELT(out, 1, V);
    SEXP nm = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(nm, 0, mkChar("C"));
    SET_STRING_ELT(nm, 1, mkChar("V"));
    setAttrib(out, R_NamesSymbol, nm);

    UNPROTECT(6);
    return out;
}
