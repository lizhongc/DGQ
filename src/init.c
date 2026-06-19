/* init.c -- registration-only entry point for the DGQ shared library.
 *
 * Keeping registration separate from r_empirical.c and r_timewise.c makes the
 * SEXP wrappers easy to identify while exposing only the listed interfaces.
 */
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP r_dgq_empirical(SEXP Z, SEXP ref, SEXP threads);
extern SEXP r_dgq_timewise(SEXP Z, SEXP u, SEXP iter, SEXP eps);

/* Map each R-visible native name to its wrapper address and fixed arity. The
 * NULL sentinel terminates the table as required by R's registration API.
 */
static const R_CallMethodDef CallEntries[] = {
    {"r_dgq_empirical", (DL_FUNC) &r_dgq_empirical, 3},
    {"r_dgq_timewise", (DL_FUNC) &r_dgq_timewise, 4},
    {NULL, NULL, 0}
};

/* Register the package's fixed-arity .Call interface when R loads DGQ. */
void R_init_DGQ(DllInfo *dll)
{
    /* R calls R_init_DGQ when loading the DLL. Disabling dynamic lookup means
     * unregistered symbols cannot be called accidentally by name.
     */
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
