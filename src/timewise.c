/* timewise.c -- time-wise geometric quantile implementation.
 *
 * See dgq.h for the objective, estimating equation, modified Weiszfeld
 * update, convergence rule, and array-layout contract.
 */
#include "dgq.h"
#include <math.h>

/* Convert the mathematical subscript Z_i(t)[k] to its linear buffer offset. */
#define ZIDX(i, t, k, N, T) ((size_t)(i) + (size_t)(N) * (t) + (size_t)(N) * (T) * (k))

/* Solve the geometric-u-quantile problem independently for every time point.
 *
 * q and qn are caller-owned length-d scratch buffers holding the current and
 * next d-dimensional iterates; they are reused across time points.
 */
void compute_timewise(const double *Z, int N, int T, int d,
                      const double *u, int iter, double eps,
                      double *q, double *qn, double *out)
{
    int i, k, t, it;

    for (t = 0; t < T; t++) {
        /* Initialize q^(0) at the time-t cross-sectional mean, a stable
         * central starting point for the iterative solve.
         */
        for (k = 0; k < d; k++) {
            double m = 0.0;
            for (i = 0; i < N; i++) m += Z[ZIDX(i, t, k, N, T)];
            q[k] = m / (double) N;
        }

        for (it = 0; it < iter; it++) {
            double wsum = 0.0;
            double num_extra;            /* directional term N*u[k] */
            for (k = 0; k < d; k++) qn[k] = 0.0;

            /* Form inverse-distance weights w_i and the weighted numerator
             * sum_i w_i Z_i(t). max(distance, eps) prevents division by zero.
             */
            for (i = 0; i < N; i++) {
                double ss = 0.0;
                for (k = 0; k < d; k++) {
                    double diff = Z[ZIDX(i, t, k, N, T)] - q[k];
                    ss += diff * diff;
                }
                double dist = sqrt(ss);
                if (dist < eps) dist = eps;
                double w = 1.0 / dist;
                wsum += w;
                for (k = 0; k < d; k++)
                    qn[k] += w * Z[ZIDX(i, t, k, N, T)];
            }

            /* Complete q^(a+1) = (sum_i w_i Z_i + N*u) / sum_i w_i
             * and accumulate the squared Euclidean update length.
             */
            double step = 0.0;
            for (k = 0; k < d; k++) {
                num_extra = (double) N * u[(size_t) t + (size_t) T * k];
                double newk = (qn[k] + num_extra) / wsum;
                double dq = newk - q[k];
                step += dq * dq;
                qn[k] = newk;
            }
            for (k = 0; k < d; k++) q[k] = qn[k];
            if (sqrt(step) < eps) break;
        }

        /* R expects the result as a column-major T-by-d matrix. */
        for (k = 0; k < d; k++)
            out[(size_t) t + (size_t) T * k] = q[k];
    }
}
