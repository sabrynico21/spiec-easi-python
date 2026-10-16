// Rcpp thin wrapper for RIC — delegates to huge::ric()
#include <Rcpp.h>
#include <cmath>
#include <limits>
#include "huge/huge_core.h"
#include "huge_r_threads.h"
using namespace Rcpp;

//[[Rcpp::export]]
double RIC(NumericMatrix &X, int d, int n, NumericVector &r, int t)
{
    if (n <= 0 || d <= 0)
        stop("RIC n and d must be positive.");
    if (X.nrow() != n || X.ncol() != d)
        stop("RIC X dimensions must exactly match n by d.");
    for (R_xlen_t i = 0; i < X.size(); i++) {
        if (!std::isfinite(X[i]))
            stop("RIC X must contain only finite values.");
    }
    if (t <= 0)
        stop("RIC t must be positive.");
    if (static_cast<R_xlen_t>(t) > r.size())
        stop("RIC t must not exceed the length of r.");

    // Convert only the selected rotations after validating that conversion
    // cannot truncate, overflow, or trigger the core's legacy clamping.
    std::vector<int> r_int(t);
    for (int i = 0; i < t; i++) {
        double value = r[i];
        if (!std::isfinite(value) || value != std::floor(value))
            stop("RIC r values must be finite integers.");
        if (value < static_cast<double>(std::numeric_limits<int>::min()) ||
                value > static_cast<double>(std::numeric_limits<int>::max()))
            stop("RIC r values must lie within the C++ integer range.");
        int rotation = static_cast<int>(value);
        if (rotation < 0 || rotation > n)
            stop("RIC r values must be rotation indices in [0, n].");
        r_int[i] = rotation;
    }

    if (d <= 1)
        return 0.0;
    // X is column-major in R (NumericMatrix)
    huge::RThreadLimitGuard thread_limit;
    return huge::ric(X.begin(), n, d, r_int.data(), t);
}
