// huge_core.h — Standalone C++ core for huge graph estimation algorithms.
// No dependency on Rcpp, pybind11, or Eigen. Pure standard C++.
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace huge {

// ---- Column-major matrix ------------------------------------------------

struct Matrix {
    std::vector<double> v;
    int rows = 0, cols = 0;

    Matrix() = default;
    Matrix(int r, int c) : v(static_cast<size_t>(r) * c, 0.0), rows(r), cols(c) {}

    void resize(int r, int c) { rows = r; cols = c; v.assign(static_cast<size_t>(r) * c, 0.0); }
    void set_zero() { std::fill(v.begin(), v.end(), 0.0); }

    double& operator()(int r, int c)             { return v[static_cast<size_t>(c) * rows + r]; }
    const double& operator()(int r, int c) const { return v[static_cast<size_t>(c) * rows + r]; }
    double* col_ptr(int c)             { return v.data() + static_cast<size_t>(c) * rows; }
    const double* col_ptr(int c) const { return v.data() + static_cast<size_t>(c) * rows; }
};

// ---- MB / TIGER per-column results --------------------------------------

struct ColResult {
    std::vector<double> vals;
    std::vector<int> indices; // encoded as lambda_idx * d + var_idx
};

// ---- Glasso -------------------------------------------------------------

struct GlassoResult {
    std::vector<double> loglik;   // nlambda
    std::vector<double> sparsity; // nlambda
    std::vector<int>    df;       // nlambda
    std::vector<Matrix> path;     // nlambda x (d,d)
    std::vector<Matrix> icov;     // nlambda x (d,d)
    std::vector<Matrix> cov;      // nlambda x (d,d); empty if !cov_output
    // True when any solver loop exhausted its iteration budget anywhere on
    // the path: results are then best-effort, not converged to tolerance.
    bool hit_max_iter = false;
};

GlassoResult glasso(const double* S_colmajor, int d,
                    const double* lambda, int nlambda,
                    bool scr, bool cov_output);

// Adapter-oriented variant that leaves GlassoResult::path empty.  The graph
// support is exactly the off-diagonal nonzero pattern of icov, so language
// bindings can derive their compact path output without retaining a second
// dense double matrix path inside the core result.
GlassoResult glasso_compact(const double* S_colmajor, int d,
                            const double* lambda, int nlambda,
                            bool scr, bool cov_output);

// ---- MB graph -----------------------------------------------------------

struct MBResult {
    std::vector<ColResult> columns; // size d
    bool hit_max_iter = false;      // any requested point was not certified
};

MBResult mb(const double* S_colmajor, int d,
            const double* lambda, int nlambda);

MBResult mb_scr(const double* S_colmajor, int d,
                const double* lambda, int nlambda,
                // Column-major nscr-by-d matrix of distinct zero-based
                // predictors; each column excludes its response index.
                const int* idx_scr, int nscr);

// ---- TIGER (sqrt-lasso) -------------------------------------------------

struct TigerResult {
    std::vector<ColResult> columns; // size d
    std::vector<Matrix>    icov;    // nlambda x (d,d)
    std::vector<double>    lambda;  // actual path used by the native solver
    bool hit_max_iter = false;      // any column exhausted an iteration budget
    bool path_truncated = false;    // generated path kept only its certified prefix
};

// Legacy raw-sample entry point retained while language adapters migrate to
// the correlation-domain implementation below.
TigerResult tiger(const double* data_colmajor, int n, int d,
                  const double* lambda, int nlambda);

// Build a correlation matrix from either raw n-by-d observations or a d-by-d
// covariance/correlation matrix, generate the default lambda path when
// lambda == nullptr, and solve TIGER entirely in the correlation domain.
TigerResult tiger_fit(const double* input_colmajor, int n, int d,
                      bool covariance_input,
                      const double* lambda, int nlambda,
                      double lambda_min_ratio);

// ---- RIC ----------------------------------------------------------------

double ric(const double* X_colmajor, int n, int d,
           const int* r, int t);

// ---- Scale-free graph generator -----------------------------------------
// G_out: pre-allocated d*d array (column-major), written as adjacency matrix.
// rands: array of (d - d0) uniform(0,1) random values, supplied by caller.

void sfgen(int d0, int d, int* G_out, const double* rands);

} // namespace huge
