/* dgq.h -- mathematical and memory-layout contract for the pure C core.
 *
 * The routines declared here contain no R API calls. They receive plain C
 * buffers from R_export.c and compute in the standardized coordinate system.
 *
 * DATA AND ARRAY LAYOUT
 * ---------------------
 * There are N observed trajectories, T time points, and d variables. Write
 * Z_i(t) in R^d for standardized series i at time t and Z_i(t)[k] for its
 * k-th coordinate. R stores Z as a column-major array with dim = c(N, T, d):
 *
 *   Z_i(t)[k] = Z[i + N*t + N*T*k],     i,t,k are zero based.
 *
 * A matrix with N rows and d columns similarly stores entry (j,k) at j + N*k.
 * The C routines preserve these layouts so their output can be returned to R
 * without copying or transposition.
 *
 * MATHEMATICAL CONVENTIONS
 * ------------------------
 * All norms and inner products below are Euclidean in standardized space.
 * Whitening in R turns the requested Mahalanobis geometry into this Euclidean
 * geometry before the C routines run. A tilt direction u must lie in the open
 * unit ball, ||u|| < 1, for the geometric-quantile objective to be coercive.
 */
#ifndef DGQ_H
#define DGQ_H

#include <stddef.h>

/* Compute the depth and tilt terms used by the empirical DGQ.
 *
 * For candidate observed trajectory Z_j, define
 *
 *   C(j) = (N / nref) sum_r sum_t ||Z_ref[r](t) - Z_j(t)||,
 *   V(j, t) = Zbar(t) - Z_j(t),
 *   Zbar(t) = (1/N) sum_i Z_i(t).
 *
 * The empirical DGQ with tilt trajectory v(t) is the observed-series index minimizing
 *
 *   L_v(j) = C(j) + N sum_t <v(t), V(j, t)>.
 *
 * If ref contains all N indices, C(j) is the exact sum of pairwise trajectory
 * distances. For a size-nref reference sample, N/nref scales the sampled sum
 * to the CLARA estimate of the full sum. V does not depend on ref and is exact.
 *
 * Inputs:
 *   Z       standardized (N,T,d) array;
 *   ref     nref zero-based series indices;
 *   N,T,d   dimensions, assumed positive;
 *   nref    number of reference indices.
 *   mean_tk caller-owned length-T*d scratch buffer;
 *   nthreads requested candidate-level thread count, assumed positive.
 * Outputs:
 *   C       length-N vector, C[j] = C(j);
 *   V       column-major (N,T,d) array, V[j + N*t + N*T*k] = Zbar(t)[k] - Z_j(t)[k].
 *
 * mean_tk is caller-provided T-by-d scratch space. Cross-sectional means are
 * computed once in O(N*T*d), after which independent candidate evaluations
 * may run in parallel. The total work is O(N*nref*T*d + N*T*d), and no N-by-N
 * matrix is formed. With OpenMP support, nthreads controls the candidate-level
 * parallel region; otherwise execution is serial.
 */
void compute_empirical(const double *, int, int, int,
                       const int *, int,
                       double *, int,
                       double *, double *);

/* Compute the time-wise geometric u-quantile, the unconstrained DGQ target.
 *
 * Independently at each time t, q_u(t) minimizes
 *
 *   Q_t(q) = sum_i ||Z_i(t) - q|| - N <u(t), q>,
 *
 * up to the additive constant in the equivalent tilted-distance formulation.
 * Away from observations, its first-order equation is
 *
 *   sum_i {q - Z_i(t)} / ||q - Z_i(t)|| = N u(t).
 *
 * Given iterate q^(a), the implemented modified Weiszfeld update is
 *
 *   w_i^(a) = 1 / max(||Z_i(t) - q^(a)||, eps),
 *   q^(a+1) = {sum_i w_i^(a) Z_i(t) + N u(t)} / sum_i w_i^(a).
 *
 * Iteration starts at the cross-sectional mean and stops when
 * ||q^(a+1)-q^(a)|| < eps or after iter updates. The max-with-eps guard avoids
 * division by zero when an iterate coincides with an observation.
 *
 * Inputs:
 *   Z       standardized (N,T,d) array;
 *   u       column-major (T,d) matrix of tilt directions;
 *   iter    maximum updates per time point;
 *   eps     positive distance guard and convergence tolerance.
 *   q, qn   caller-owned length-d scratch buffers for the current and next
 *           iterates; reused across time points.
 * Output:
 *   out     column-major (T,d) matrix, out[t + T*k] = q_u(t)[k].
 *
 * Time complexity is O(T*iter*N*d); scratch (q, qn) is caller-provided and O(d).
 */
void compute_timewise(const double *, int, int, int,
                      const double *, int, double,
                      double *, double *, double *);

#endif /* DGQ_H */
