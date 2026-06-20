/* empirical.c -- empirical DGQ depth and tilt implementation.
 *
 * See dgq.h for the complete mathematical contract. This file maps the depth
 * C(j) and tilt V(j) definitions directly to loops over the R-compatible
 * column-major buffers.
 */
#include "dgq.h"
#include <math.h>

/* Convert the mathematical subscript Z_i(t)[k] to its linear buffer offset. */
#define ZIDX(i, t, k, N, T) ((size_t)(i) + (size_t)(N) * (t) + (size_t)(N) * (T) * (k))

/* Evaluate C(j) and V(j) for every observed candidate trajectory j.
 *
 * C uses either all series (exact method) or the supplied CLARA references.
 * V always uses the complete cross-section. Output buffers are allocated and
 * validated by the R boundary wrapper; this pure numeric routine does not
 * allocate memory or report R errors.
 */
void compute_empirical(const double *Z, int N, int T, int d,
                       const int *ref, int nref,
                       double *mean_tk, int nthreads,
                       double *C, double *V)
{
    int i, j, k, t;
    double scale = (nref > 0) ? ((double) N / (double) nref) : 0.0;
    #ifndef _OPENMP
    (void) nthreads;
    #endif

    /* Precompute Zbar(t)[k] once. The i-loop preserves the original
     * cross-sectional summation order used by the serial implementation.
     */
    for (k = 0; k < d; k++) {
        for (t = 0; t < T; t++) {
            double mean = 0.0;
            for (i = 0; i < N; i++)
                mean += Z[ZIDX(i, t, k, N, T)];
            mean_tk[(size_t) t + (size_t) T * k] = mean / (double) N;
        }
    }

    /* Estimate the full depth by scaling the reference-sample distance sum:
     * C(j) = (N/nref) sum_r sum_t ||Z_ref[r](t) - Z_j(t)||.
     * The innermost coordinate loop forms one squared Euclidean norm. Each
     * candidate j is independent and writes only C[j] and row j of V, making
     * this outer loop safe for static OpenMP parallelism. Builds without
     * OpenMP ignore the pragma and execute the same loop serially.
     */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(nthreads) \
        default(none) shared(Z, N, T, d, ref, nref, scale, mean_tk, C, V)
    #endif
    for (j = 0; j < N; j++) {
        double cj = 0.0;
        for (int rr = 0; rr < nref; rr++) {
            int ri = ref[rr];
            for (int tt = 0; tt < T; tt++) {
                double ss = 0.0;
                for (int kk = 0; kk < d; kk++) {
                    double diff = Z[ZIDX(ri, tt, kk, N, T)] - Z[ZIDX(j, tt, kk, N, T)];
                    ss += diff * diff;
                }
                cj += sqrt(ss);
            }
        }
        C[j] = scale * cj;

        /* Compute the exact tilt coordinates in the original time-variable order:
         * V[j, t, k] = Zbar(t)[k] - Z_j(t)[k].
         */
        for (int kk = 0; kk < d; kk++) {
            for (int tt = 0; tt < T; tt++) {
                V[ZIDX(j, tt, kk, N, T)] = mean_tk[(size_t) tt + (size_t) T * kk] - Z[ZIDX(j, tt, kk, N, T)];
            }
        }
    }
}
