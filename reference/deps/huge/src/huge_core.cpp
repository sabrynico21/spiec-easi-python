// huge_core.cpp — Standalone core implementations.
// No Rcpp, no pybind11, no Eigen. Uses BLAS for hot-path linear algebra.
#include "huge/huge_core.h"
#include "huge/blas_config.h"
#include <cstring>
#include <stdexcept>

// BLAS constants reused throughout
static const char   BLAS_N   = 'N';
static const char   BLAS_T   = 'T';
static const int    BLAS_1   = 1;
static const double BLAS_ONE = 1.0, BLAS_ZERO = 0.0;

namespace huge {

static inline double threshold_l1(double x, double thr) {
    if (x > thr) return x - thr;
    if (x < -thr) return x + thr;
    return 0.0;
}

// Column-major accessor for raw const pointer
static inline double cm(const double* data, int nrows, int r, int c) {
    return data[static_cast<size_t>(c) * nrows + r];
}

// Emit one lambda's active coefficients in ascending coordinate order.
// Both language wrappers build compressed sparse columns, which require
// ascending row indices; sorting here removes the reorder pass the R layer
// used to do. Sorts a scratch copy — never the solver's own active list,
// whose order determines subsequent sweep order.
static void collect_sorted(ColResult& col, int lambda_idx, int d,
                           const double* w, const int* active, int size_active,
                           std::vector<int>& scratch) {
    scratch.assign(active, active + size_active);
    std::sort(scratch.begin(), scratch.end());
    for (int j = 0; j < size_active; j++) {
        col.vals.push_back(w[scratch[j]]);
        col.indices.push_back(lambda_idx * d + scratch[j]);
    }
}

// =========================================================================
// Glasso: inner coordinate-descent solver
// =========================================================================

static constexpr double GLASSO_DEFAULT_CONVERGENCE_TOL = 1e-4;
static constexpr double GLASSO_REFINEMENT_CONVERGENCE_TOL = 1e-8;

static void glasso_sub(Matrix& S, Matrix& W, Matrix& T, int d,
                       double ilambda, int& df, bool scr, bool& hit_max_iter,
                       double convergence_tolerance)
{
    const double thol_act = convergence_tolerance;
    const double thol_ext = convergence_tolerance;
    const int MAX_ITER_EXT = 100;
    const int MAX_ITER_INT = 10000;
    const int MAX_ITER_ACT = 10000;
    const double neg_ilambda = -ilambda;

    // idx_a, idx_i stored as d columns of max-d rows each
    std::vector<int> idx_a(static_cast<size_t>(d) * d, 0);
    std::vector<int> idx_i(static_cast<size_t>(d) * d, 0);
    std::vector<int> size_a(static_cast<size_t>(d), 0);
    std::vector<double> w1(static_cast<size_t>(d), 0.0);
    std::vector<double> ww(static_cast<size_t>(d), 0.0);
    std::vector<double> temp(static_cast<size_t>(d), 0.0);  // reused in W-update

    auto idx_a_col = [&](int col) -> int* { return idx_a.data() + static_cast<size_t>(col) * d; };
    auto idx_i_col = [&](int col) -> int* { return idx_i.data() + static_cast<size_t>(col) * d; };

    // Recover initial solution for each individual lasso
    for (int i = 0; i < d; i++) {
        int* ia = idx_a_col(i);
        int* ii = idx_i_col(i);
        double* T_col = T.col_ptr(i);
        const double* S_col = S.col_ptr(i);

        W(i, i) = S_col[i] + ilambda;
        size_a[i] = 0;
        double diag_T = T_col[i];
        T_col[i] = 0;

        for (int j = 0; j < d; j++) {
            if (scr && std::fabs(S_col[j]) <= ilambda) {
                ii[j] = -1;
                T_col[j] = 0;
                continue;
            }
            if (T_col[j] != 0) {
                ia[size_a[i]++] = j;
                ii[j] = -1;
                T_col[j] = -T_col[j] / diag_T;
            } else {
                ii[j] = 1;
            }
        }
        ii[i] = -1;
    }

    double gap_ext = 1;
    int iter_ext = 0;
    // Convergence bookkeeping: coef_abs_sum accumulates |w| over the last
    // active-set pass of a column; the outer gap is (total change in T
    // columns) / (total coefficient mass).
    double coef_abs_sum = 0, ext_change_sum = 0, ext_coef_sum = 0;
    while (gap_ext > thol_ext && iter_ext < MAX_ITER_EXT) {
        ext_coef_sum = 0;
        ext_change_sum = 0;
        for (int i = 0; i < d; i++) {
            int active_size = size_a[i];
            int* ia = idx_a_col(i);
            int* ii = idx_i_col(i);
            double* T_col = T.col_ptr(i);
            const double* S_col = S.col_ptr(i);

            int gap_int = 1;
            int iter_int = 0;
            for (int j = 0; j < d; j++) ww[j] = T_col[j];

            while (gap_int != 0 && iter_int < MAX_ITER_INT) {
                int size_a_prev = active_size;
                for (int j = 0; j < d; j++) {
                    if (ii[j] != -1) {
                        double r = S_col[j];
                        for (int k = 0; k < active_size; k++)
                            r -= W(ia[k], j) * T_col[ia[k]];

                        double w_new = 0.0;
                        if (r > ilambda) {
                            w_new = (r - ilambda) / W(j, j);
                            ia[active_size++] = j;
                            ii[j] = -1;
                        } else if (r < neg_ilambda) {
                            w_new = (r + ilambda) / W(j, j);
                            ia[active_size++] = j;
                            ii[j] = -1;
                        }
                        w1[j] = w_new;
                        T_col[j] = w_new;
                    }
                }
                gap_int = active_size - size_a_prev;

                double gap_act = 1;
                int iter_act = 0;
                while (gap_act > thol_act && iter_act < MAX_ITER_ACT) {
                    double act_change_sum = 0;
                    coef_abs_sum = 0;
                    for (int j = 0; j < active_size; j++) {
                        int w_idx = ia[j];
                        double r = S_col[w_idx] + T_col[w_idx] * W(w_idx, w_idx);
                        for (int k = 0; k < active_size; k++)
                            r -= W(ia[k], w_idx) * T_col[ia[k]];

                        double w_new = 0.0;
                        if (r > ilambda)
                            w_new = (r - ilambda) / W(w_idx, w_idx);
                        else if (r < neg_ilambda)
                            w_new = (r + ilambda) / W(w_idx, w_idx);

                        coef_abs_sum += std::fabs(w_new);
                        act_change_sum += std::fabs(w_new - T_col[w_idx]);
                        w1[w_idx] = w_new;
                        T_col[w_idx] = w_new;
                    }
                    gap_act = coef_abs_sum > 0 ? act_change_sum / coef_abs_sum : 0;
                    iter_act++;
                }

                // Move false active variables back
                int junk_a = 0;
                for (int j = 0; j < active_size; j++) {
                    int w_idx = ia[j];
                    if (w1[w_idx] == 0) {
                        junk_a++;
                        ii[w_idx] = 1;
                    } else {
                        ia[j - junk_a] = w_idx;
                    }
                }
                active_size -= junk_a;
                iter_int++;
            }
            size_a[i] = active_size;

            // Update W from current T column: W[:,i] = W * T[:,i].
            // T_col is nonzero only on the active set, so accumulating those
            // columns with daxpy does exactly the useful work; a fused dgemv
            // processes all d columns including zeros. Calibrated 2026-07-20
            // at d=2000 (10-lambda path, interleaved, 2 rounds): always-daxpy
            // 132.0s/132.0s vs always-dgemv 184.8s/189.0s, with threshold
            // variants in between — so no dense fallback.
            std::fill(temp.begin(), temp.end(), 0.0);
            for (int k = 0; k < active_size; k++) {
                int col = ia[k];
                double coef = T_col[col];
                if (coef != 0.0)
                    daxpy_(&d, &coef, W.col_ptr(col), &BLAS_1,
                           temp.data(), &BLAS_1);
            }
            for (int j = 0; j < d; j++) {
                if (j != i) { W(j, i) = temp[j]; W(i, j) = temp[j]; }
            }

            for (int j = 0; j < d; j++)
                ext_change_sum += std::fabs(ww[j] - T_col[j]);
            ext_coef_sum += coef_abs_sum;
        }
        gap_ext = ext_coef_sum > 0 ? ext_change_sum / ext_coef_sum : 0;
        iter_ext++;
    }
    if (gap_ext > thol_ext && iter_ext >= MAX_ITER_EXT) hit_max_iter = true;

    // Compute final T (precision matrix)
    for (int i = 0; i < d; i++) {
        double* T_col = T.col_ptr(i);
        double quad_form = ddot_(&d, W.col_ptr(i), &BLAS_1, T_col, &BLAS_1);
        quad_form -= W(i, i) * T_col[i];
        double inv_diag = 1.0 / (W(i, i) - quad_form);
        double neg_inv_diag = -inv_diag;
        dscal_(&d, &neg_inv_diag, T_col, &BLAS_1);
        T_col[i] = inv_diag;
    }
    for (int i = 0; i < d; i++) df += size_a[i];
}

// Scale-equilibrated Cholesky for the precision log-determinant. Scaling each
// diagonal to one avoids an absolute pivot cutoff across tiny/huge matrices,
// while a non-positive pivot correctly rejects a non-SPD precision estimate.
static double spd_log_det_colmajor(const Matrix& m) {
    const int n = m.rows;
    if (n <= 0 || n != m.cols) return -std::numeric_limits<double>::infinity();
    Matrix L(n, n);
    std::vector<double> inv_scale(n);
    double ldet = 0.0;
    for (int i = 0; i < n; i++) {
        double diagonal = m(i, i);
        if (!(diagonal > 0.0) || !std::isfinite(diagonal))
            return -std::numeric_limits<double>::infinity();
        inv_scale[i] = 1.0 / std::sqrt(diagonal);
        if (!std::isfinite(inv_scale[i]))
            return -std::numeric_limits<double>::infinity();
        ldet += std::log(diagonal);
    }
    for (int j = 0; j < n; j++) {
        double pivot = 1.0;
        for (int k = 0; k < j; k++) pivot -= L(j, k) * L(j, k);
        if (!(pivot > 0.0) || !std::isfinite(pivot))
            return -std::numeric_limits<double>::infinity();
        double diagonal = std::sqrt(pivot);
        L(j, j) = diagonal;
        ldet += 2.0 * std::log(diagonal);
        for (int i = j + 1; i < n; i++) {
            double value = m(i, j) * inv_scale[i] * inv_scale[j];
            for (int k = 0; k < j; k++) value -= L(i, k) * L(j, k);
            value /= diagonal;
            if (!std::isfinite(value))
                return -std::numeric_limits<double>::infinity();
            L(i, j) = value;
        }
    }
    return ldet;
}

// Infinity norm of covariance * precision - I. Both matrices are already
// available even when covariance output is not requested by the caller.
static double inverse_residual_inf(const Matrix& covariance,
                                   const Matrix& precision) {
    const int n = covariance.rows;
    if (n <= 0 || covariance.cols != n ||
            precision.rows != n || precision.cols != n)
        return std::numeric_limits<double>::infinity();
    Matrix product(n, n);
    dgemm_(&BLAS_N, &BLAS_N, &n, &n, &n, &BLAS_ONE,
           covariance.v.data(), &n, precision.v.data(), &n,
           &BLAS_ZERO, product.v.data(), &n);
    double residual = 0.0;
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) {
            double value = product(i, j) - (i == j ? 1.0 : 0.0);
            if (!std::isfinite(value))
                return std::numeric_limits<double>::infinity();
            row_sum += std::fabs(value);
        }
        if (!std::isfinite(row_sum))
            return std::numeric_limits<double>::infinity();
        residual = std::max(residual, row_sum);
    }
    return residual;
}

static constexpr double GLASSO_INVERSE_RESIDUAL_TOL = 1e-2;
static constexpr double GLASSO_REFINEMENT_RESIDUAL_TRIGGER =
    0.5 * GLASSO_INVERSE_RESIDUAL_TOL;

// trace(A*B) = sum_k dot(A[:,k], B[k,:]).  B is symmetric here (sub_S), so
// B[k,:] == B[:,k] and both dot operands are contiguous columns — bitwise
// the same result as the strided form, without the stride-d walk.
static double trace_matmul(const Matrix& a, const Matrix& b) {
    const int d = a.rows;
    double tr = 0;
    for (int k = 0; k < d; k++)
        tr += ddot_(&d, a.col_ptr(k), &BLAS_1, b.col_ptr(k), &BLAS_1);
    return tr;
}

static bool matrix_is_finite(const Matrix& matrix) {
    return std::all_of(
        matrix.v.begin(), matrix.v.end(),
        [](double value) { return std::isfinite(value); }
    );
}

static void symmetrize_square_in_place(Matrix& matrix) {
    for (int i = 0; i < matrix.rows; i++) {
        for (int j = i + 1; j < matrix.cols; j++) {
            double average = 0.5 * matrix(i, j) + 0.5 * matrix(j, i);
            matrix(i, j) = average;
            matrix(j, i) = average;
        }
    }
}

static void validate_regularization_inputs(const double* matrix, int d,
                                           const double* lambda,
                                           int nlambda) {
    if (d <= 0)
        throw std::invalid_argument(
            "regularization matrix dimension must be positive.");
    if (matrix == nullptr)
        throw std::invalid_argument(
            "regularization matrix must not be null.");
    if (nlambda <= 0)
        throw std::invalid_argument(
            "regularization nlambda must be positive.");
    if (lambda == nullptr)
        throw std::invalid_argument(
            "regularization lambda must not be null.");
    for (int i = 0; i < nlambda; i++) {
        if (!std::isfinite(lambda[i]) || lambda[i] <= 0.0)
            throw std::invalid_argument(
                "regularization lambda values must be positive and finite.");
    }
    for (int i = 1; i < nlambda; i++) {
        if (lambda[i] > lambda[i - 1])
            throw std::invalid_argument(
                "regularization lambda values must be non-increasing "
                "(ties are allowed).");
    }
}

// =========================================================================
// Glasso: outer driver with pre-screening
// =========================================================================

static GlassoResult glasso_impl(const double* S_data, int d,
                               const double* lambda, int nlambda,
                               bool scr, bool cov_output,
                               bool include_path)
{
    validate_regularization_inputs(S_data, d, lambda, nlambda);
    GlassoResult res;
    res.loglik.assign(nlambda, -static_cast<double>(d));
    res.sparsity.assign(nlambda, 0.0);
    res.df.assign(nlambda, 0);
    res.icov.resize(nlambda);
    if (cov_output) res.cov.resize(nlambda);

    // Build S matrix
    Matrix S(d, d);
    std::memcpy(S.v.data(), S_data, static_cast<size_t>(d) * d * sizeof(double));

    std::vector<double> s_diag(d);
    for (int i = 0; i < d; i++) s_diag[i] = S(i, i);

    const double sparsity_denom = d > 1 ? static_cast<double>(d) * (d - 1) : 1.0;

    std::vector<Matrix> tmp_icov(nlambda), tmp_cov;
    if (cov_output) tmp_cov.resize(nlambda);
    const Matrix* prev_icov_ptr = nullptr;
    const Matrix* prev_cov_ptr = nullptr;
    // Two W buffers, alternated per lambda: cur is written while prev is
    // still being read for warm starts (a single reused buffer would zero
    // the previous solution before it is consumed).
    Matrix cov_pingpong[2];
    int cov_flip = 0;

    // Union-find scratch for the exact connected-component decomposition
    // (Witten/Friedman & Mazumder/Hastie): the glasso solution is block
    // diagonal over the connected components of {(j,k): |S_jk| > lambda},
    // so each component can be solved independently — and cross-component
    // entries of icov AND cov are exactly zero at the optimum.
    std::vector<int> uf_parent(d);
    std::vector<int> comp_of(d);
    auto uf_find = [&](int a) {
        while (uf_parent[a] != a) {
            uf_parent[a] = uf_parent[uf_parent[a]];
            a = uf_parent[a];
        }
        return a;
    };

    for (int i = nlambda - 1; i >= 0; i--) {
        double lambda_i = lambda[i];

        // Components of the thresholded adjacency (O(d^2) scan, the same
        // cost the previous single-block screening paid).
        for (int j = 0; j < d; j++) uf_parent[j] = j;
        for (int k = 1; k < d; k++) {
            const double* S_col = S.col_ptr(k);
            for (int j = 0; j < k; j++) {
                if (std::fabs(S_col[j]) > lambda_i) {
                    int ra = uf_find(j), rb = uf_find(k);
                    if (ra != rb) uf_parent[rb] = ra;
                }
            }
        }
        // Order components by root; singletons are handled by the diagonal
        // initialization below and never enter a solver.
        std::vector<std::vector<int>> comps;
        std::fill(comp_of.begin(), comp_of.end(), -1);
        for (int j = 0; j < d; j++) {
            int r = uf_find(j);
            if (r == j) continue;              // roots resolved when a member arrives
            if (comp_of[r] == -1) {
                comp_of[r] = static_cast<int>(comps.size());
                comps.emplace_back();
                comps.back().push_back(r);
            }
            comps[comp_of[r]].push_back(j);
        }
        for (auto& c : comps) {
            std::sort(c.begin(), c.end());
            // Reuse comp_of as a membership marker below.  Nodes left at -1
            // are singleton components and contribute their exact diagonal
            // log-likelihood without entering glasso_sub().
            for (int node : c) comp_of[node] = 0;
        }

        // Build full-size output matrices
        Matrix& cur_icov = tmp_icov[i];
        cur_icov.resize(d, d);

        Matrix* cur_cov_ptr;
        if (cov_output) {
            tmp_cov[i].resize(d, d);
            cur_cov_ptr = &tmp_cov[i];
        } else {
            cov_flip ^= 1;
            cov_pingpong[cov_flip].resize(d, d);   // resize zero-fills
            cur_cov_ptr = &cov_pingpong[cov_flip];
        }

        int total_edges = 0;
        double total_ldet = 0.0;
        double total_tr = 0.0;
        bool ldet_finite = true;

        // Diagonal initialization.  For singleton components this is the
        // complete solution, so include its exact contribution to
        // log(det(Theta)) - trace(S * Theta).  The previous fixed -1 term was
        // only correct when lambda == 0 and S(j,j) == 1.
        for (int j = 0; j < d; j++) {
            double diagonal_cov = s_diag[j] + lambda_i;
            cur_icov(j, j) = 1.0 / diagonal_cov;
            (*cur_cov_ptr)(j, j) = diagonal_cov;
            if (comp_of[j] == -1) {
                if (diagonal_cov > 0.0 && std::isfinite(diagonal_cov)) {
                    total_ldet -= std::log(diagonal_cov);
                    total_tr += s_diag[j] / diagonal_cov;
                } else {
                    ldet_finite = false;
                }
            }
        }

        // Components are solved serially. Parallelizing this loop was
        // measured 2026-07-20 and REMOVED: typical correlation structures
        // form one giant component at every useful lambda (zero gain), while
        // OpenMP-nesting Accelerate BLAS calls adds platform-dependent risk.
        Matrix sub_S, sub_W, sub_T;
        for (const auto& z : comps) {
            int q = static_cast<int>(z.size());
            sub_S.resize(q, q); sub_W.resize(q, q); sub_T.resize(q, q);
            for (int ii = 0; ii < q; ii++) {
                for (int jj = 0; jj < q; jj++) {
                    sub_S(ii, jj) = S(z[ii], z[jj]);
                    if (prev_cov_ptr == nullptr || prev_icov_ptr == nullptr) {
                        sub_W(ii, jj) = S(z[ii], z[jj]);
                        sub_T(ii, jj) = (ii == jj) ? 1.0 : 0.0;
                    } else {
                        sub_W(ii, jj) = (*prev_cov_ptr)(z[ii], z[jj]);
                        sub_T(ii, jj) = (*prev_icov_ptr)(z[ii], z[jj]);
                    }
                }
            }
            int solver_directed_df = 0;
            bool component_hit_max_iter = false;
            glasso_sub(sub_S, sub_W, sub_T, q, lambda_i,
                       solver_directed_df, scr,
                       component_hit_max_iter,
                       GLASSO_DEFAULT_CONVERGENCE_TOL);
            double raw_inverse_residual =
                inverse_residual_inf(sub_W, sub_T);

            // Column-wise coordinate solves are only approximately symmetric
            // at finite tolerance. Project onto the symmetric matrices before
            // exposing the precision estimate, deriving its graph, computing
            // likelihood metadata, or using it as the next warm start.
            symmetrize_square_in_place(sub_T);
            double projected_ldet = spd_log_det_colmajor(sub_T);
            double projected_residual =
                inverse_residual_inf(sub_W, sub_T);

            // Symmetrization can amplify the inverse residual of a column-wise
            // solution. Refine a finite SPD candidate once at a tighter
            // tolerance, including a useful iterate that reached the first
            // solve's limit; the independent checks below still reject an
            // inconsistent refined result.
            if (matrix_is_finite(sub_T) && matrix_is_finite(sub_W)) {
                if (std::isfinite(projected_ldet) &&
                        std::isfinite(projected_residual) &&
                        std::isfinite(raw_inverse_residual) &&
                        (projected_residual >
                             GLASSO_REFINEMENT_RESIDUAL_TRIGGER ||
                         raw_inverse_residual >
                             GLASSO_INVERSE_RESIDUAL_TOL)) {
                    bool refinement_hit_max_iter = false;
                    solver_directed_df = 0;
                    glasso_sub(
                        sub_S, sub_W, sub_T, q, lambda_i,
                        solver_directed_df, scr, refinement_hit_max_iter,
                        GLASSO_REFINEMENT_CONVERGENCE_TOL
                    );
                    component_hit_max_iter =
                        component_hit_max_iter || refinement_hit_max_iter;
                    raw_inverse_residual =
                        inverse_residual_inf(sub_W, sub_T);
                    symmetrize_square_in_place(sub_T);
                    projected_ldet = spd_log_det_colmajor(sub_T);
                    projected_residual =
                        inverse_residual_inf(sub_W, sub_T);
                }
            }
            res.hit_max_iter =
                res.hit_max_iter || component_hit_max_iter;

            int component_edges = 0;
            for (int ii = 0; ii < q; ii++) {
                for (int jj = ii + 1; jj < q; jj++) {
                    if (sub_T(ii, jj) != 0.0) component_edges++;
                }
            }
            total_edges += component_edges;

            if (!matrix_is_finite(sub_T) || !matrix_is_finite(sub_W))
                throw std::runtime_error("glasso produced non-finite estimates.");
            if (!std::isfinite(projected_ldet))
                throw std::runtime_error(
                    "glasso produced a precision estimate that is not positive definite.");
            if (!std::isfinite(raw_inverse_residual) ||
                    raw_inverse_residual > GLASSO_INVERSE_RESIDUAL_TOL)
                throw std::runtime_error(
                    "glasso produced inconsistent precision and covariance estimates.");
            if (!std::isfinite(projected_residual) ||
                    projected_residual > GLASSO_INVERSE_RESIDUAL_TOL)
                throw std::runtime_error(
                    "glasso produced inconsistent precision and covariance estimates.");

            for (int ii = 0; ii < q; ii++) {
                for (int jj = 0; jj < q; jj++) {
                    cur_icov(z[ii], z[jj]) = sub_T(ii, jj);
                    (*cur_cov_ptr)(z[ii], z[jj]) = sub_W(ii, jj);
                }
            }
            total_ldet += projected_ldet;
            total_tr += trace_matmul(sub_T, sub_S);
        }

        res.sparsity[i] = 2.0 * total_edges / sparsity_denom;
        res.df[i] = total_edges;
        res.loglik[i] = ldet_finite
            ? (total_ldet - total_tr)
            : -std::numeric_limits<double>::infinity();
        if (!matrix_is_finite(cur_icov) ||
            !matrix_is_finite(*cur_cov_ptr) ||
            !std::isfinite(res.loglik[i])) {
            throw std::runtime_error("glasso produced non-finite estimates.");
        }
        prev_icov_ptr = &cur_icov;
        prev_cov_ptr = cur_cov_ptr;
    }

    res.icov = std::move(tmp_icov);
    if (include_path) {
        res.path.resize(nlambda);
        for (int k = 0; k < nlambda; k++) {
            Matrix& path = res.path[k];
            path.resize(d, d);
            for (int j = 0; j < d; j++) {
                for (int i = 0; i < d; i++) {
                    path(i, j) =
                        (i != j && res.icov[k](i, j) != 0.0) ? 1.0 : 0.0;
                }
            }
        }
    }
    if (cov_output) res.cov = std::move(tmp_cov);
    return res;
}

GlassoResult glasso(const double* S_data, int d,
                    const double* lambda, int nlambda,
                    bool scr, bool cov_output)
{
    return glasso_impl(
        S_data, d, lambda, nlambda, scr, cov_output, true);
}

GlassoResult glasso_compact(const double* S_data, int d,
                            const double* lambda, int nlambda,
                            bool scr, bool cov_output)
{
    return glasso_impl(
        S_data, d, lambda, nlambda, scr, cov_output, false);
}

// =========================================================================
// MB graph estimation — shared column-solver pieces
// =========================================================================

// Residual gradient of coordinate j for the m-th nodewise lasso:
//   r_j = S[m,j] - sum_k S[j, a_k] * w[a_k]  over the active set.
// NOTE: S(j, idx_a[k]) reads look stride-d, but with j advancing in the
// caller's outer loop each fixed k walks its column sequentially — size_a
// parallel prefetch streams. Rewriting to stay within column j (using
// symmetry) measured 2x SLOWER at d=2000; do not "fix" this access pattern.
static inline double mb_partial_residual(const double* S_data, int d, int m,
                                         int j, const int* idx_a, int size_a,
                                         const double* w) {
    double r = cm(S_data, d, m, j);
    for (int k = 0; k < size_a; k++)
        r -= cm(S_data, d, j, idx_a[k]) * w[idx_a[k]];
    return r;
}

// Cycle coordinate descent over the current active set until the relative
// coefficient change drops below thol. Shared verbatim by mb() and mb_scr().
static void mb_refine_active(const double* S_data, int d, int m,
                             double ilambda, double thol, int max_iter,
                             const int* idx_a, int size_a,
                             double* w0, double* w1) {
    double gap_int = 1;
    int iter_int = 0;
    while (gap_int > thol && iter_int < max_iter) {
        double change_sum = 0, coef_sum = 0;
        for (int j = 0; j < size_a; j++) {
            int w_idx = idx_a[j];
            double r = cm(S_data, d, m, w_idx) + w0[w_idx];
            for (int k = 0; k < size_a; k++)
                r -= cm(S_data, d, w_idx, idx_a[k]) * w0[idx_a[k]];

            w1[w_idx] = threshold_l1(r, ilambda);
            coef_sum += std::fabs(w1[w_idx]);
            change_sum += std::fabs(w1[w_idx] - w0[w_idx]);
            w0[w_idx] = w1[w_idx];
        }
        gap_int = (coef_sum > 0) ? change_sum / coef_sum : 0;
        iter_int++;
    }
}

// =========================================================================
// MB graph estimation (lossless screening)
// =========================================================================

MBResult mb(const double* S_data, int d,
            const double* lambda, int nlambda)
{
    validate_regularization_inputs(S_data, d, lambda, nlambda);
    MBResult res;
    res.columns.resize(d);

    const double thol = 1e-4;
    const int MAX_ITER = 10000;
    // Written under "#pragma omp atomic write" from worker threads; only
    // ever transitions false -> true.
    bool shared_hit_max_iter = false;

    #ifdef _OPENMP
    #pragma omp parallel
    {
    #endif
    std::vector<double> w0(d, 0.0), w1(d, 0.0);
    std::vector<int> idx_a(d), idx_i(d);
    // Sequential strong rule state (picasso actgd pattern): grad_abs[j] holds
    // |gradient_j| refreshed by the previous lambda's certification pass; per
    // lambda only coordinates with grad_abs > 2*lambda[i] - lambda[i-1] are
    // swept, then one full KKT pass certifies the screen and repairs rare
    // violations (re-running the strong loop when any are found).
    std::vector<unsigned char> strong(d, 0);
    std::vector<double> grad_abs(d, 0.0);
    std::vector<int> sort_scratch;

    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for (int m = 0; m < d; m++) {
        bool hit_any_max_iter = false;
        ColResult& col = res.columns[m];
        idx_i[m] = 0;
        for (int j = 0; j < m; j++) idx_i[j] = 1;
        for (int j = m + 1; j < d; j++) idx_i[j] = 1;

        int size_a = 0;
        std::fill(w0.begin(), w0.end(), 0.0);
        std::fill(w1.begin(), w1.end(), 0.0);

        // Gradient at w = 0 is the m-th row of S.
        for (int j = 0; j < d; j++) grad_abs[j] = std::fabs(cm(S_data, d, m, j));

        for (int i = 0; i < nlambda; i++) {
            double ilambda = lambda[i];
            double strong_thr = (i > 0) ? 2.0 * lambda[i] - lambda[i - 1]
                                        : 2.0 * lambda[i];
            for (int j = 0; j < d; j++)
                strong[j] = (idx_i[j] == 1 && grad_abs[j] > strong_thr) ? 1 : 0;

            int iter_ext = 0;
            bool certified = false;
            while (!certified && iter_ext < MAX_ITER) {
              int gap_ext = 1;
              while (gap_ext != 0 && iter_ext < MAX_ITER) {
                int size_a_prev = size_a;
                for (int j = 0; j < d; j++) {
                    if (idx_i[j] == 1 && strong[j]) {
                        double r = mb_partial_residual(S_data, d, m, j,
                                                       idx_a.data(), size_a,
                                                       w0.data());
                        grad_abs[j] = std::fabs(r);
                        w1[j] = threshold_l1(r, ilambda);
                        if (w1[j] != 0) {
                            idx_a[size_a++] = j;
                            idx_i[j] = 0;
                        }
                        w0[j] = w1[j];
                    }
                }
                gap_ext = size_a - size_a_prev;

                mb_refine_active(S_data, d, m, ilambda, thol, MAX_ITER,
                                 idx_a.data(), size_a, w0.data(), w1.data());
                int junk_a = 0;
                for (int j = 0; j < size_a; j++) {
                    int w_idx = idx_a[j];
                    if (w1[w_idx] == 0) {
                        junk_a++;
                        idx_i[w_idx] = 1;
                        strong[w_idx] = 1;  // recently active: keep sweeping
                    }
                    else idx_a[j - junk_a] = w_idx;
                }
                size_a -= junk_a;
                iter_ext++;
            }

            // Certification: one full pass over every inactive coordinate.
            // Refreshes grad_abs for the next lambda's screen and catches
            // strong-rule violations; violators activate here and the strong
            // loop re-runs to converge them.
            int violations = 0;
            for (int j = 0; j < d; j++) {
                if (idx_i[j] == 1) {
                    double r = mb_partial_residual(S_data, d, m, j,
                                                   idx_a.data(), size_a,
                                                   w0.data());
                    grad_abs[j] = std::fabs(r);
                    w1[j] = threshold_l1(r, ilambda);
                    if (w1[j] != 0) {
                        idx_a[size_a++] = j;
                        idx_i[j] = 0;
                        strong[j] = 1;
                        violations++;
                    }
                    w0[j] = w1[j];
                }
            }
            certified = (violations == 0);
            }
            if (!certified) hit_any_max_iter = true;

            collect_sorted(col, i, d, w1.data(), idx_a.data(), size_a,
                           sort_scratch);
        }
        if (hit_any_max_iter) {
            #ifdef _OPENMP
            #pragma omp atomic write
            #endif
            shared_hit_max_iter = true;
        }
    }
    #ifdef _OPENMP
    }
    #endif
    res.hit_max_iter = shared_hit_max_iter;
    return res;
}

// =========================================================================
// MB graph estimation (lossy screening)
// =========================================================================

MBResult mb_scr(const double* S_data, int d,
                const double* lambda, int nlambda,
                const int* idx_scr, int nscr)
{
    validate_regularization_inputs(S_data, d, lambda, nlambda);
    MBResult res;
    if (nscr <= 0 || nscr >= d)
        throw std::invalid_argument(
            "idx_scr must have between 1 and d - 1 rows.");
    if (idx_scr == nullptr)
        throw std::invalid_argument("idx_scr must not be null.");

    // Validate before entering the OpenMP region. The public matrix contains
    // only real zero-based predictors; -1 is reserved for the private
    // mutable copy below after a predictor enters the active set.
    std::vector<int> seen(d, -1);
    for (int m = 0; m < d; m++) {
        for (int j = 0; j < nscr; j++) {
            int idx = idx_scr[static_cast<size_t>(m) * nscr + j];
            if (idx < 0 || idx >= d)
                throw std::invalid_argument(
                    "idx_scr entries must be zero-based indices in [0, d).");
            if (idx == m)
                throw std::invalid_argument(
                    "idx_scr columns must exclude their response index.");
            if (seen[idx] == m)
                throw std::invalid_argument(
                    "idx_scr columns must contain distinct indices.");
            seen[idx] = m;
        }
    }
    res.columns.resize(d);

    const double thol = 1e-4;
    const int MAX_ITER = 10000;
    bool shared_hit_max_iter = false;  // omp atomic write; false -> true only

    #ifdef _OPENMP
    #pragma omp parallel
    {
    #endif
    std::vector<double> w0(d, 0.0), w1(d, 0.0);
    std::vector<int> idx_a(nscr), idx_i_local(nscr);
    std::vector<int> sort_scratch;

    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for (int m = 0; m < d; m++) {
        bool hit_any_max_iter = false;
        ColResult& col = res.columns[m];
        int size_a = 0;

        for (int j = 0; j < nscr; j++)
            idx_i_local[j] = idx_scr[static_cast<size_t>(m) * nscr + j];

        std::fill(w0.begin(), w0.end(), 0.0);
        std::fill(w1.begin(), w1.end(), 0.0);

        for (int i = 0; i < nlambda; i++) {
            double ilambda = lambda[i];
            int gap_ext = 1, iter_ext = 0;
            while (iter_ext < MAX_ITER && gap_ext > 0) {
                int size_a_prev = size_a;
                for (int j = 0; j < nscr; j++) {
                    int w_idx = idx_i_local[j];
                    if (w_idx != -1) {
                        double r = mb_partial_residual(S_data, d, m, w_idx,
                                                       idx_a.data(), size_a,
                                                       w0.data());
                        w1[w_idx] = threshold_l1(r, ilambda);
                        if (w1[w_idx] != 0) {
                            idx_a[size_a++] = w_idx;
                            idx_i_local[j] = -1;
                        }
                        w0[w_idx] = w1[w_idx];
                    }
                }
                gap_ext = size_a - size_a_prev;

                mb_refine_active(S_data, d, m, ilambda, thol, MAX_ITER,
                                 idx_a.data(), size_a, w0.data(), w1.data());
                iter_ext++;
            }
            if (gap_ext > 0 && iter_ext >= MAX_ITER) hit_any_max_iter = true;
            collect_sorted(col, i, d, w1.data(), idx_a.data(), size_a,
                           sort_scratch);
        }
        if (hit_any_max_iter) {
            #ifdef _OPENMP
            #pragma omp atomic write
            #endif
            shared_hit_max_iter = true;
        }
    }
    #ifdef _OPENMP
    }
    #endif
    res.hit_max_iter = shared_hit_max_iter;
    return res;
}

// =========================================================================
// TIGER (sqrt-lasso graph estimation)
// =========================================================================

static void validate_tiger_lambda_path(const double* lambda, int nlambda)
{
    if (lambda == nullptr)
        throw std::invalid_argument("tiger lambda path must not be null.");
    for (int i = 0; i < nlambda; i++) {
        if (!std::isfinite(lambda[i]) || lambda[i] <= 0.0)
            throw std::invalid_argument(
                "tiger lambda values must be positive and finite.");
    }
    for (int i = 1; i < nlambda; i++) {
        if (lambda[i] > lambda[i - 1])
            throw std::invalid_argument(
                "tiger lambda values must be non-increasing (ties are allowed).");
    }
}

TigerResult tiger(const double* data_colmajor, int n, int d,
                  const double* lambda, int nlambda)
{
    TigerResult res;
    if (n <= 0 || d <= 0 || nlambda <= 0) {
        res.columns.resize(d > 0 ? d : 0);
        return res;
    }
    validate_tiger_lambda_path(lambda, nlambda);
    res.columns.resize(d);

    // Pre-initialize icov matrices
    res.icov.resize(nlambda);
    for (int i = 0; i < nlambda; i++) res.icov[i].resize(d, d);

    const double prec = 1e-4;
    const int max_iter = 1000;
    const int num_relaxation_round = 3;
    const double eps = 1e-12;
    bool shared_hit_max_iter = false;  // omp atomic write; false -> true only

    // Precompute squared column norms: xx_dot[j] = sum_t X[t,j]^2
    // Immutable; shared safely across OMP threads.
    std::vector<double> xx_dot(d);
    for (int j = 0; j < d; j++)
        xx_dot[j] = ddot_(&n, data_colmajor + static_cast<size_t>(j) * n, &BLAS_1,
                               data_colmajor + static_cast<size_t>(j) * n, &BLAS_1);

    #ifdef _OPENMP
    #pragma omp parallel
    {
    #endif
    // Thread-local workspaces
    std::vector<double> Xb(n, 0.0), r_vec(n, 0.0), grad(d, 0.0), w1(d, 0.0);
    std::vector<double> Y(n, 0.0), gr(d, 0.0), rx(n, 0.0);
    std::vector<double> Xb_master(n, 0.0), w1_master(d, 0.0);
    std::vector<int> actset_indcat(d, 0), actset_indcat_master(d, 0);
    std::vector<int> actset_idx;
    std::vector<double> old_coef(d, 0.0), grad_master(d, 0.0);
    std::vector<int> sort_scratch;

    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for (int m = 0; m < d; m++) {
        bool hit_any_max_iter = false;
        ColResult& col = res.columns[m];

        std::fill(Xb.begin(), Xb.end(), 0.0);
        std::fill(w1.begin(), w1.end(), 0.0);
        std::fill(grad.begin(), grad.end(), 0.0);
        std::fill(gr.begin(), gr.end(), 0.0);
        std::fill(w1_master.begin(), w1_master.end(), 0.0);
        std::fill(Xb_master.begin(), Xb_master.end(), 0.0);
        std::fill(actset_indcat.begin(), actset_indcat.end(), 0);
        std::fill(actset_indcat_master.begin(), actset_indcat_master.end(), 0);
        std::fill(old_coef.begin(), old_coef.end(), 0.0);
        std::fill(grad_master.begin(), grad_master.end(), 0.0);

        const double* Y_col = data_colmajor + static_cast<size_t>(m) * n;
        std::memcpy(Y.data(), Y_col, n * sizeof(double));

        double L = 0, sum_r2 = 0;
        auto refresh_residual = [&]() {
            sum_r2 = ddot_(&n, r_vec.data(), &BLAS_1, r_vec.data(), &BLAS_1);
            if (sum_r2 < eps) sum_r2 = eps;
            L = std::sqrt(sum_r2 / n);
            if (L < eps) L = eps;
        };

        for (int t = 0; t < n; t++) r_vec[t] = Y[t] - Xb[t];
        refresh_residual();
        double dev_thr = std::fabs(L) * prec;

        for (int j = 0; j < d; j++) {
            const double* x_col = data_colmajor + static_cast<size_t>(j) * n;
            double s = ddot_(&n, r_vec.data(), &BLAS_1, x_col, &BLAS_1);
            grad[j] = s / (n * L);
            gr[j] = std::fabs(grad[j]);
            grad_master[j] = gr[j];
            w1_master[j] = w1[j];
        }
        std::memcpy(Xb_master.data(), Xb.data(), n * sizeof(double));

        for (int i = 0; i < nlambda; i++) {
            double stage_lambda = lambda[i];
            w1 = w1_master;
            Xb = Xb_master;
            for (int j = 0; j < d; j++) {
                gr[j] = grad_master[j];
                actset_indcat[j] = actset_indcat_master[j];
            }

            double threshold = (i > 0) ? 2 * lambda[i] - lambda[i - 1] : 2 * lambda[i];
            for (int j = 0; j < d; j++)
                if (j != m && gr[j] > threshold) actset_indcat[j] = 1;

            for (int t = 0; t < n; t++) r_vec[t] = Y[t] - Xb[t];
            refresh_residual();

            double tmp_change = 0, local_change = 0;

            // update_coordinate: uses precomputed xx_dot[ci] and BLAS ddot for
            // speed; rx is a thread-local temp buffer (size n). NOTE: a fused
            // single-pass scalar loop for both reductions measured 1.5x SLOWER
            // than buffer + two Accelerate ddots at d=1000 — keep the BLAS form.
            auto update_coordinate = [&](int ci) {
                const double* x_col = data_colmajor + static_cast<size_t>(ci) * n;
                for (int t = 0; t < n; t++) rx[t] = r_vec[t] * x_col[t];
                double dot_rxrx = ddot_(&n, rx.data(),     &BLAS_1, rx.data(),     &BLAS_1);
                double dot_rx   = ddot_(&n, r_vec.data(),  &BLAS_1, x_col,         &BLAS_1);
                double sum_wxx  = xx_dot[ci] - dot_rxrx / sum_r2;
                double a = sum_wxx / (n * L);
                double g = (sum_wxx * w1[ci] + dot_rx) / (n * L);
                double oldv = w1[ci];
                w1[ci] = (std::fabs(a) > eps) ? threshold_l1(g, stage_lambda) / a : 0.0;
                double delta = w1[ci] - oldv;
                if (delta != 0) {
                    daxpy_(&n, &delta, x_col, &BLAS_1, Xb.data(), &BLAS_1);
                    double neg_delta = -delta;
                    daxpy_(&n, &neg_delta, x_col, &BLAS_1, r_vec.data(), &BLAS_1);
                    // O(1) residual-norm update in place of an O(n) ddot:
                    // ||r - d*x||^2 = ||r||^2 - 2*d*(r.x) + d^2*||x||^2, with
                    // dot_rx taken against the pre-update residual. The full
                    // recomputes at sweep/stage boundaries reset accumulated
                    // rounding; clamps mirror refresh_residual().
                    sum_r2 += -2.0 * delta * dot_rx + delta * delta * xx_dot[ci];
                    if (sum_r2 < eps) sum_r2 = eps;
                    L = std::sqrt(sum_r2 / n);
                    if (L < eps) L = eps;
                }
            };

            int loopcnt_level_0 = 0;
            while (loopcnt_level_0 < num_relaxation_round) {
                loopcnt_level_0++;

                int loopcnt_level_1 = 0;
                bool level1_converged = false;
                while (loopcnt_level_1 < max_iter) {
                    loopcnt_level_1++;
                    for (int j = 0; j < d; j++) old_coef[j] = w1[j];
                    refresh_residual();

                    actset_idx.clear();
                    for (int j = 0; j < d; j++) {
                        if (j == m || !actset_indcat[j]) continue;
                        update_coordinate(j);
                        if (std::fabs(w1[j]) > 0) actset_idx.push_back(j);
                    }

                    // Level 2: proximal newton on active set
                    int loopcnt_level_2 = 0;
                    while (loopcnt_level_2 < max_iter) {
                        loopcnt_level_2++;
                        bool term2 = true;
                        for (int k = 0; k < static_cast<int>(actset_idx.size()); k++) {
                            int idx = actset_idx[k];
                            double old_w1 = w1[idx];
                            update_coordinate(idx);
                            tmp_change = old_w1 - w1[idx];
                            // local_change = h * tmp_change^2 / (2Ln) is 0 when the
                            // coordinate did not move; skip the O(n) h computation
                            // (dev_thr > 0, so 0 can never exceed it).
                            if (tmp_change == 0) continue;
                            // O(1) conservative bound: drxrx >= 0 gives
                            // |hsum| <= xx_dot, so local_change <=
                            // xx_dot * tmp_change^2 / (2*L^2*n^2). When even the
                            // bound is below dev_thr the exact value cannot flip
                            // term2; skip the O(n) reduction (picasso sqrtmse).
                            if (xx_dot[idx] * tmp_change * tmp_change
                                    <= dev_thr * (2.0 * L * L * n * n))
                                continue;
                            const double* xc = data_colmajor + static_cast<size_t>(idx) * n;
                            for (int t = 0; t < n; t++) rx[t] = r_vec[t] * xc[t];
                            double drxrx = ddot_(&n, rx.data(), &BLAS_1, rx.data(), &BLAS_1);
                            double hsum = xx_dot[idx] - drxrx / sum_r2;
                            double h = std::fabs(hsum / (n * L));
                            local_change = h * tmp_change * tmp_change / (2.0 * L * n);
                            if (local_change > dev_thr) term2 = false;
                        }
                        if (term2) break;
                    }

                    // Check stopping criterion 1
                    bool term1 = true;
                    for (size_t k = 0; k < actset_idx.size(); k++) {
                        int idx = actset_idx[k];
                        tmp_change = old_coef[idx] - w1[idx];
                        if (tmp_change == 0) continue;  // same skips as level 2
                        if (xx_dot[idx] * tmp_change * tmp_change
                                <= dev_thr * (2.0 * L * L * n * n))
                            continue;
                        const double* xc = data_colmajor + static_cast<size_t>(idx) * n;
                        for (int t = 0; t < n; t++) rx[t] = r_vec[t] * xc[t];
                        double drxrx = ddot_(&n, rx.data(), &BLAS_1, rx.data(), &BLAS_1);
                        double hsum = xx_dot[idx] - drxrx / sum_r2;
                        double h = std::fabs(hsum / (n * L));
                        local_change = h * tmp_change * tmp_change / (2.0 * L * n);
                        if (local_change > dev_thr) term1 = false;
                    }
                    for (int t = 0; t < n; t++) r_vec[t] = Y[t] - Xb[t];
                    refresh_residual();
                    if (term1) { level1_converged = true; break; }

                    // Check stopping criterion 2: active set change
                    bool new_active = false;
                    for (int k = 0; k < d; k++) {
                        if (k == m || actset_indcat[k] != 0) continue;
                        const double* x_col_k = data_colmajor + static_cast<size_t>(k) * n;
                        double s = ddot_(&n, r_vec.data(), &BLAS_1, x_col_k, &BLAS_1);
                        grad[k] = s / (n * L);
                        gr[k] = std::fabs(grad[k]);
                        if (gr[k] > stage_lambda) {
                            actset_indcat[k] = 1;
                            new_active = true;
                        }
                    }
                    if (!new_active) { level1_converged = true; break; }
                }
                if (!level1_converged) hit_any_max_iter = true;

                if (loopcnt_level_0 == 1) {
                    w1_master = w1;
                    Xb_master = Xb;
                    for (int j = 0; j < d; j++) {
                        grad_master[j] = gr[j];
                        actset_indcat_master[j] = actset_indcat[j];
                    }
                }
            }

            // Collect results
            collect_sorted(col, i, d, w1.data(), actset_idx.data(),
                           static_cast<int>(actset_idx.size()), sort_scratch);

            for (int t = 0; t < n; t++) r_vec[t] = Y[t] - Xb[t];
            refresh_residual();
            double tal = L;

            // Write icov column m — each thread writes different column, no race
            Matrix& icov_ref = res.icov[i];
            icov_ref(m, m) = (tal > 0) ? 1.0 / (tal * tal) : 0.0;
            for (int j = 0; j < d; j++)
                if (j != m) icov_ref(j, m) = -icov_ref(m, m) * w1[j];
        }
        if (hit_any_max_iter) {
            #ifdef _OPENMP
            #pragma omp atomic write
            #endif
            shared_hit_max_iter = true;
        }
    }
    #ifdef _OPENMP
    }
    #endif
    res.hit_max_iter = shared_hit_max_iter;

    // Symmetrize icov (upper triangle only: writing both sides in a full
    // sweep would average against already-averaged values)
    for (int i = 0; i < nlambda; i++) {
        Matrix& ic = res.icov[i];
        for (int c0 = 0; c0 < d; c0++)
            for (int r0 = 0; r0 < c0; r0++) {
                double avg = 0.5 * (ic(r0, c0) + ic(c0, r0));
                ic(r0, c0) = avg;
                ic(c0, r0) = avg;
            }
    }
    return res;
}

// Build the sample correlation matrix in native code.  Raw observations are
// centered and normalized column-by-column before one BLAS cross-product;
// covariance input is converted with D^{-1/2} S D^{-1/2}.  In both cases the
// returned matrix has an exact unit diagonal and is explicitly symmetric.
static bool correlation_is_positive_semidefinite(const Matrix& corr)
{
    const int d = corr.rows;
    if (d <= 0 || corr.cols != d) return false;

    double spectral_bound = 1.0;
    for (int r = 0; r < d; r++) {
        double row_sum = 0.0;
        for (int c = 0; c < d; c++) row_sum += std::fabs(corr(r, c));
        spectral_bound = std::max(spectral_bound, row_sum);
    }
    const double tolerance = 100.0 * std::numeric_limits<double>::epsilon()
        * static_cast<double>(std::max(1, d)) * spectral_bound;

    // Cholesky of corr + tolerance * I accepts singular PSD input while
    // rejecting negative eigenvalues beyond floating-point roundoff.
    Matrix lower(d, d);
    for (int j = 0; j < d; j++) {
        double pivot = corr(j, j) + tolerance;
        for (int k = 0; k < j; k++) pivot -= lower(j, k) * lower(j, k);
        if (!(pivot > 0.0) || !std::isfinite(pivot)) return false;
        lower(j, j) = std::sqrt(pivot);

        for (int i = j + 1; i < d; i++) {
            double value = corr(i, j);
            for (int k = 0; k < j; k++)
                value -= lower(i, k) * lower(j, k);
            value /= lower(j, j);
            if (!std::isfinite(value)) return false;
            lower(i, j) = value;
        }
    }
    return true;
}

static Matrix tiger_build_correlation(const double* input, int n, int d,
                                      bool covariance_input)
{
    if (input == nullptr || d <= 0)
        throw std::invalid_argument("tiger input must be a non-empty numeric matrix.");

    Matrix corr(d, d);
    if (covariance_input) {
        if (n != d)
            throw std::invalid_argument("tiger covariance input must be square.");

        std::vector<double> variance(d);
        std::vector<double> sd(d);
        std::vector<double> inv_sd(d);
        for (int j = 0; j < d; j++) {
            variance[j] = cm(input, d, j, j);
            if (!std::isfinite(variance[j]) || variance[j] <= 0.0)
                throw std::invalid_argument(
                    "tiger covariance input must have a positive finite diagonal.");
            sd[j] = std::sqrt(variance[j]);
            inv_sd[j] = 1.0 / sd[j];
        }

        const double symmetry_tolerance =
            100.0 * std::numeric_limits<double>::epsilon();
        for (int c = 0; c < d; c++) {
            corr(c, c) = 1.0;
            for (int r = 0; r < c; r++) {
                double src_rc = cm(input, d, r, c);
                double src_cr = cm(input, d, c, r);
                if (!std::isfinite(src_rc) || !std::isfinite(src_cr))
                    throw std::invalid_argument(
                        "tiger covariance input must contain only finite values.");
                bool equal = src_rc == src_cr;
                if (!equal) {
                    double covariance_scale = variance[r] == variance[c]
                        ? variance[c] : sd[r] * sd[c];
                    double reference = std::max(
                        std::max(std::fabs(src_rc), std::fabs(src_cr)),
                        covariance_scale);
                    double normalized_rc = src_rc / reference;
                    double normalized_cr = src_cr / reference;
                    double scale = std::max(
                        std::fabs(normalized_rc), std::fabs(normalized_cr));
                    double threshold = scale <= symmetry_tolerance
                        ? symmetry_tolerance : symmetry_tolerance * scale;
                    double difference = std::fabs(
                        normalized_rc - normalized_cr);
                    if (difference > threshold)
                        throw std::invalid_argument(
                            "tiger covariance input must be symmetric within numeric tolerance.");
                }
                double symmetric = equal
                    ? src_rc : 0.5 * src_rc + 0.5 * src_cr;
                double inv_large = std::max(inv_sd[r], inv_sd[c]);
                double inv_small = std::min(inv_sd[r], inv_sd[c]);
                double value = (symmetric * inv_large) * inv_small;
                if (!std::isfinite(value) ||
                        std::fabs(value) > 1.0 + 1e-8)
                    throw std::invalid_argument(
                        "tiger covariance input is not a valid covariance matrix.");
                value = std::max(-1.0, std::min(1.0, value));
                corr(r, c) = value;
                corr(c, r) = value;
            }
        }
        if (!correlation_is_positive_semidefinite(corr))
            throw std::invalid_argument(
                "tiger covariance input must be positive semidefinite.");
        return corr;
    }

    if (n < 2)
        throw std::invalid_argument("tiger raw input requires at least two observations.");

    Matrix standardized(n, d);
    auto add_compensated = [](double value, double& sum,
                              double& correction) {
        double adjusted = value - correction;
        double updated = sum + adjusted;
        correction = (updated - sum) - adjusted;
        sum = updated;
    };
    for (int j = 0; j < d; j++) {
        const double* source = input + static_cast<size_t>(j) * n;
        double column_scale = 0.0;
        for (int i = 0; i < n; i++) {
            if (!std::isfinite(source[i]))
                throw std::invalid_argument(
                    "tiger raw input must contain only finite values.");
            column_scale = std::max(column_scale, std::fabs(source[i]));
        }
        if (column_scale == 0.0)
            throw std::invalid_argument("tiger raw input contains a constant column.");

        double* target = standardized.col_ptr(j);
        // Exact power-of-two scaling bounds every value without erasing
        // representable ULP differences near the floating-point maximum.
        int column_exponent = 0;
        std::frexp(column_scale, &column_exponent);
        const int scale_exponent = -column_exponent;
        const double scaled_reference = std::scalbn(source[0], scale_exponent);
        double delta_sum = 0.0;
        double delta_correction = 0.0;
        for (int i = 0; i < n; i++) {
            double scaled = std::scalbn(source[i], scale_exponent);
            double delta = scaled - scaled_reference;
            target[i] = delta;
            add_compensated(delta, delta_sum, delta_correction);
        }

        const double mean_delta = delta_sum / static_cast<double>(n);
        double sum_squares = 0.0;
        double square_correction = 0.0;
        for (int i = 0; i < n; i++) {
            target[i] -= mean_delta;
            add_compensated(target[i] * target[i], sum_squares,
                            square_correction);
        }
        if (!std::isfinite(sum_squares) || sum_squares <= 0.0)
            throw std::invalid_argument("tiger raw input contains a constant column.");

        double inv_norm = 1.0 / std::sqrt(sum_squares);
        dscal_(&n, &inv_norm, target, &BLAS_1);
    }

    dgemm_(&BLAS_T, &BLAS_N, &d, &d, &n, &BLAS_ONE,
           standardized.v.data(), &n, standardized.v.data(), &n,
           &BLAS_ZERO, corr.v.data(), &d);
    for (int c = 0; c < d; c++) {
        corr(c, c) = 1.0;
        for (int r = 0; r < c; r++) {
            double value = 0.5 * (corr(r, c) + corr(c, r));
            value = std::max(-1.0, std::min(1.0, value));
            corr(r, c) = value;
            corr(c, r) = value;
        }
    }
    return corr;
}

static std::vector<double> tiger_prepare_lambda(const Matrix& corr,
                                                const double* lambda,
                                                int nlambda,
                                                double lambda_min_ratio)
{
    if (nlambda <= 0)
        throw std::invalid_argument("tiger nlambda must be positive.");

    std::vector<double> path(nlambda);
    if (lambda != nullptr) {
        validate_tiger_lambda_path(lambda, nlambda);
        for (int i = 0; i < nlambda; i++) path[i] = lambda[i];
        return path;
    }

    if (!std::isfinite(lambda_min_ratio) || lambda_min_ratio <= 0.0 ||
            lambda_min_ratio > 1.0)
        throw std::invalid_argument("tiger lambda_min_ratio must lie in (0, 1].");

    double lambda_max = 0.0;
    for (int c = 1; c < corr.cols; c++)
        for (int r = 0; r < c; r++)
            lambda_max = std::max(lambda_max, std::fabs(corr(r, c)));
    // Keep the path finite and strictly positive for identity/one-variable
    // inputs, matching the native Python front end's established floor.
    if (lambda_max == 0.0) lambda_max = 1e-3;

    if (nlambda == 1) {
        path[0] = lambda_max;
        return path;
    }
    double log_max = std::log(lambda_max);
    double lambda_min = lambda_max * lambda_min_ratio;
    if (lambda_min > 0.0) {
        // Preserve the established path bit-for-bit when its endpoint is
        // representable.
        double log_min = std::log(lambda_min);
        for (int i = 0; i < nlambda; i++) {
            double fraction = static_cast<double>(i) / (nlambda - 1);
            path[i] = std::exp(
                log_max + fraction * (log_min - log_max));
        }
    } else {
        // Retain representable interior points in the requested geometric
        // grid, then saturate only its unrepresentable tail at the smallest
        // positive double.  This avoids log(0), zero lambda values, and NaNs.
        double smallest_positive =
            std::numeric_limits<double>::denorm_min();
        if (!(smallest_positive > 0.0))
            smallest_positive = std::numeric_limits<double>::min();
        double log_ratio = std::log(lambda_min_ratio);
        path[0] = lambda_max;
        for (int i = 1; i < nlambda; i++) {
            double fraction = static_cast<double>(i) / (nlambda - 1);
            double candidate = std::exp(log_max + fraction * log_ratio);
            if (!(candidate > 0.0)) candidate = smallest_positive;
            path[i] = std::min(path[i - 1], candidate);
        }
    }
    return path;
}

// Solve the fixed-tau Lasso subproblem from the variational SQRT-Lasso
// formulation.  residual_gradient is s - R*beta and is updated in O(d) per
// coordinate, so inactive KKT certification is exact and inexpensive.
static bool tiger_refine_lasso(const Matrix& corr, int response,
                               double penalty,
                               std::vector<double>& beta,
                               std::vector<double>& residual_gradient,
                               std::vector<int>& active,
                               std::vector<unsigned char>& is_active)
{
    const int d = corr.rows;
    // Correlation matrices are singular whenever d > n, so cap each fixed-tau
    // solve.  A separate normalized KKT check below decides whether an iterate
    // is safe to expose; exhausting this budget never masquerades as success.
    const int max_sweeps = 100;
    const double solver_tol = 1e-8;
    const double kkt_tol = 1e-10 * std::max(1.0, penalty);

    for (int sweep = 0; sweep < max_sweeps; sweep++) {
        for (int position = 0; position < static_cast<int>(active.size()); position++) {
            int j = active[position];
            double diagonal = corr(j, j);
            double partial = residual_gradient[j] + diagonal * beta[j];
            double updated = threshold_l1(partial, penalty) / diagonal;
            double delta = updated - beta[j];
            if (delta != 0.0) {
                beta[j] = updated;
                const double* corr_col = corr.col_ptr(j);
                for (int k = 0; k < d; k++)
                    residual_gradient[k] -= corr_col[k] * delta;
            }
        }

        int violations = 0;
        for (int j = 0; j < d; j++) {
            if (j == response || is_active[j]) continue;
            if (std::fabs(residual_gradient[j]) > penalty + kkt_tol) {
                is_active[j] = 1;
                active.push_back(j);
                violations++;
            }
        }

        if (violations == 0) {
            double max_kkt_error = 0.0;
            for (int position = 0; position < static_cast<int>(active.size()); position++) {
                int j = active[position];
                double error = beta[j] == 0.0
                    ? std::max(std::fabs(residual_gradient[j]) - penalty, 0.0)
                    : std::fabs(residual_gradient[j] -
                                penalty * (beta[j] > 0.0 ? 1.0 : -1.0));
                max_kkt_error = std::max(max_kkt_error, error);
            }
            if (max_kkt_error <= solver_tol * std::max(1.0, penalty))
                return true;
        }
    }
    return false;
}

static bool tiger_kkt_certified(const Matrix& corr, int response,
                                double lambda, double tau,
                                const std::vector<double>& beta,
                                const std::vector<double>& residual_gradient)
{
    const double public_kkt_tol = 8e-7;
    if (!std::isfinite(tau) || tau <= 0.0) return false;

    for (int j = 0; j < corr.rows; j++) {
        if (j == response) continue;
        double score = residual_gradient[j] / tau;
        double error = beta[j] == 0.0
            ? std::max(std::fabs(score) - lambda, 0.0)
            : std::fabs(score - lambda * (beta[j] > 0.0 ? 1.0 : -1.0));
        if (!std::isfinite(error) || error > public_kkt_tol) return false;
    }
    return true;
}

static TigerResult tiger_from_correlation(const Matrix& corr,
                                           std::vector<double> lambda,
                                           bool generated_path)
{
    const int d = corr.rows;
    const int nlambda = static_cast<int>(lambda.size());
    const double tau_floor = 1e-12;
    const double tau_degenerate = 1e-8;
    const double outer_tol = 2.5e-7;
    const int max_outer = 50;

    TigerResult res;
    res.lambda = std::move(lambda);
    res.columns.resize(d);
    res.icov.resize(nlambda);
    for (int i = 0; i < nlambda; i++) res.icov[i].resize(d, d);
    std::vector<int> valid_prefix(d, nlambda);

    #ifdef _OPENMP
    #pragma omp parallel
    {
    #endif
    std::vector<double> beta(d, 0.0), residual_gradient(d, 0.0);
    std::vector<unsigned char> is_active(d, 0);
    std::vector<int> active, nonzero, sort_scratch;
    active.reserve(d); nonzero.reserve(d);

    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for (int m = 0; m < d; m++) {
        int certified_count = 0;
        std::fill(beta.begin(), beta.end(), 0.0);
        std::fill(is_active.begin(), is_active.end(), 0);
        active.clear();
        for (int j = 0; j < d; j++) residual_gradient[j] = corr(j, m);
        residual_gradient[m] = 0.0;
        is_active[m] = 1;
        double tau = std::sqrt(std::max(corr(m, m), tau_floor * tau_floor));

        for (int i = 0; i < nlambda; i++) {
            bool outer_converged = false;
            bool numerically_degenerate = false;
            for (int outer = 0; outer < max_outer; outer++) {
                double penalty = res.lambda[i] * tau;
                bool lasso_converged = tiger_refine_lasso(
                    corr, m, penalty, beta, residual_gradient, active,
                    is_active);

                // q = R_mm - 2*s'*beta + beta'*R*beta.  Because
                // residual_gradient = s - R*beta, this is the equivalent
                // O(d) expression R_mm - beta'*(s + residual_gradient).
                double q = corr(m, m);
                for (int j = 0; j < d; j++) {
                    if (j == m || beta[j] == 0.0) continue;
                    q -= beta[j] * (corr(j, m) + residual_gradient[j]);
                }
                if (!std::isfinite(q) || q < -1e-10 ||
                        q <= tau_degenerate * tau_degenerate) {
                    numerically_degenerate = true;
                    break;
                }
                double tau_new = std::sqrt(std::max(q, tau_floor * tau_floor));
                double tau_change = std::fabs(tau_new - tau);
                tau = tau_new;
                if (lasso_converged &&
                        tau_change <= outer_tol * tau &&
                        tiger_kkt_certified(corr, m, res.lambda[i], tau,
                                            beta, residual_gradient)) {
                    outer_converged = true;
                    break;
                }

                // Repeating a fully exhausted coordinate-descent solve at a
                // new tau multiplies two iteration limits and can explode
                // runtime for singular d > n correlations.  Return the best
                // iterate and report non-convergence instead.
                if (!lasso_converged)
                    break;
            }

            bool certified = outer_converged && !numerically_degenerate &&
                tiger_kkt_certified(corr, m, res.lambda[i], tau, beta,
                                    residual_gradient);
            if (!certified) break;

            nonzero.clear();
            for (int j = 0; j < d; j++)
                if (j != m && beta[j] != 0.0) nonzero.push_back(j);
            collect_sorted(res.columns[m], i, d, beta.data(), nonzero.data(),
                           static_cast<int>(nonzero.size()), sort_scratch);

            Matrix& icov = res.icov[i];
            double inverse_variance = 1.0 / (tau * tau);
            icov(m, m) = inverse_variance;
            for (int j = 0; j < d; j++)
                if (j != m) icov(j, m) = -inverse_variance * beta[j];
            certified_count = i + 1;
        }
        valid_prefix[m] = certified_count;
    }
    #ifdef _OPENMP
    }
    #endif
    int certified_nlambda = nlambda;
    for (int m = 0; m < d; m++)
        certified_nlambda = std::min(certified_nlambda, valid_prefix[m]);

    if (certified_nlambda < nlambda) {
        res.hit_max_iter = true;
        if (!generated_path)
            throw std::runtime_error(
                "tiger solver could not certify a supplied lambda; use a larger lambda.");
        if (certified_nlambda <= 0)
            throw std::runtime_error(
                "tiger solver could not certify the generated lambda path.");

        res.path_truncated = true;
        res.lambda.resize(certified_nlambda);
        res.icov.resize(certified_nlambda);
        int encoded_limit = certified_nlambda * d;
        for (int m = 0; m < d; m++) {
            ColResult& col = res.columns[m];
            auto first_suffix = std::lower_bound(
                col.indices.begin(), col.indices.end(), encoded_limit);
            size_t keep = static_cast<size_t>(first_suffix - col.indices.begin());
            col.indices.resize(keep);
            col.vals.resize(keep);
        }
    }

    for (int i = 0; i < static_cast<int>(res.lambda.size()); i++) {
        Matrix& icov = res.icov[i];
        for (int c = 0; c < d; c++)
            for (int r = 0; r < c; r++) {
                double average = 0.5 * (icov(r, c) + icov(c, r));
                icov(r, c) = average;
                icov(c, r) = average;
            }
    }
    return res;
}

TigerResult tiger_fit(const double* input_colmajor, int n, int d,
                      bool covariance_input,
                      const double* lambda, int nlambda,
                      double lambda_min_ratio)
{
    Matrix corr = tiger_build_correlation(input_colmajor, n, d,
                                          covariance_input);
    std::vector<double> path = tiger_prepare_lambda(
        corr, lambda, nlambda, lambda_min_ratio);
    return tiger_from_correlation(corr, std::move(path), lambda == nullptr);
}

// =========================================================================
// RIC (Rotation Information Criterion)
// =========================================================================

double ric(const double* X_data, int n, int d, const int* r, int t)
{
    if (d <= 1 || n <= 0 || t <= 0) return 0.0;

    // Column 2-norms, used to certify numerical zeros below.  A rotation
    // permutes rows, so the rotated column keeps its norm.
    std::vector<double> col_norm(static_cast<size_t>(d));
    for (int j = 0; j < d; j++)
        col_norm[j] = std::sqrt(ddot_(&n, X_data + static_cast<size_t>(j) * n,
                                     &BLAS_1,
                                     X_data + static_cast<size_t>(j) * n,
                                     &BLAS_1));

    // Standard forward-error bound for a length-n dot product:
    //   |computed - exact| <= gamma_n * sum_i |u_i v_i| <= gamma_n * ||u|| ||v||
    // with gamma_n = n*eps / (1 - n*eps).  A BLAS is free to reassociate,
    // vectorize, block, and use fused multiply-add, so a mathematically zero
    // inner product may come back as a tiny nonzero value whose magnitude
    // differs between implementations (this is exactly what ATLAS does on
    // rank-deficient input).  Any |C[j,k]| within the pair's bound is
    // therefore indistinguishable from zero and is certified to zero, which
    // makes the selected lambda reproducible across BLAS implementations.
    // The bound is pair-specific and scale-aware: it never erases a
    // correlation that the working precision can actually represent.
    const double eps = std::numeric_limits<double>::epsilon();
    const double scaled_eps = static_cast<double>(n) * eps;
    const double dot_gamma = scaled_eps < 1.0
        ? scaled_eps / (1.0 - scaled_eps)
        : std::numeric_limits<double>::infinity();

    double lambda_min = std::numeric_limits<double>::infinity();

    #ifdef _OPENMP
    // Each worker owns a d x d scratch matrix.  Never start more workers than
    // rotations, or idle workers can dominate RIC's peak memory.
    const int worker_count = std::min(t, omp_get_max_threads());
    #pragma omp parallel num_threads(worker_count)
    {
    #endif
    // Per-thread d x d buffer for the rotated cross-product C = Xrot^T * X
    std::vector<double> C(static_cast<size_t>(d) * d);

    #ifdef _OPENMP
    #pragma omp for schedule(dynamic) reduction(min:lambda_min)
    #endif
    for (int i = 0; i < t; i++) {
        int tmp_r = r[i];
        if (tmp_r < 0) tmp_r = 0;
        if (tmp_r > n) tmp_r = n;
        int split = n - tmp_r;

        // Row-rotating X by tmp_r makes C[j,k] = dot(X[(.+tmp_r) mod n, j], X[., k]),
        // which splits into two contiguous row-block products:
        //   C = X[tmp_r:n, :]^T * X[0:split, :]  +  X[0:tmp_r, :]^T * X[split:n, :]
        // When split == 0 the first GEMM has k = 0 and beta = 0, which zeroes C.
        dgemm_(&BLAS_T, &BLAS_N, &d, &d, &split, &BLAS_ONE,
               X_data + tmp_r, &n, X_data, &n, &BLAS_ZERO, C.data(), &d);
        dgemm_(&BLAS_T, &BLAS_N, &d, &d, &tmp_r, &BLAS_ONE,
               X_data, &n, X_data + split, &n, &BLAS_ONE, C.data(), &d);

        // Max |C[j,k]| over strictly upper-triangular pairs (j < k), matching
        // the pair set of the original scalar loops.  Splitting the rotation
        // into two GEMMs adds one rounding of the two partial sums, so allow
        // one extra eps on top of the dot-product bound.
        double lambda_max = 0;
        for (int k = 1; k < d; k++) {
            const double* col = C.data() + static_cast<size_t>(k) * d;
            for (int j = 0; j < k; j++) {
                double tmp = std::fabs(col[j]);
                double bound = dot_gamma * col_norm[j] * col_norm[k];
                bound += eps * bound + eps * tmp;
                if (tmp <= bound) continue;
                if (tmp > lambda_max) lambda_max = tmp;
            }
        }
        if (lambda_max < lambda_min) lambda_min = lambda_max;
    }
    #ifdef _OPENMP
    }
    #endif

    if (!std::isfinite(lambda_min)) return 0.0;
    return lambda_min;
}

// =========================================================================
// Scale-free graph generator
// =========================================================================

void sfgen(int d0, int d, int* G_out, const double* rands)
{
    // G_out: d*d column-major, pre-zeroed by caller
    std::memset(G_out, 0, static_cast<size_t>(d) * d * sizeof(int));
    std::vector<int> degree(d, 0);

    // Initial cycle of d0 nodes
    for (int i = 0; i < d0 - 1; i++) {
        G_out[static_cast<size_t>(i) * d + (i + 1)] = 1;
        G_out[static_cast<size_t>(i + 1) * d + i] = 1;
    }
    G_out[static_cast<size_t>(0) * d + (d0 - 1)] = 1;
    G_out[static_cast<size_t>(d0 - 1) * d + 0] = 1;

    for (int i = 0; i < d0; i++) degree[i] = 2;
    int total = 2 * d0;

    for (int i = d0; i < d; i++) {
        double x = static_cast<double>(total) * rands[i - d0];
        int tmp = 0, j = 0;
        while (tmp < x && j < i) { tmp += degree[j]; j++; }
        if (j > 0) j--;
        G_out[static_cast<size_t>(i) * d + j] = 1;
        G_out[static_cast<size_t>(j) * d + i] = 1;
        total += 2;
        degree[j]++;
        degree[i]++;
    }
}

} // namespace huge
