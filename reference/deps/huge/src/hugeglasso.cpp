// Rcpp thin wrapper for graphical lasso — delegates to huge::glasso()
//
// Memory strategy: NumericMatrix is an SEXP wrapper whose .begin() points
// directly into R's heap.  By memcpy-ing core results straight into
// NumericMatrix we avoid the extra Eigen::MatrixXd intermediate that Rcpp
// would then have to serialize again — cutting one full d×d copy per matrix.
#include <Rcpp.h>
#include "huge/huge_core.h"
using namespace Rcpp;

static void validate_lambda_dimensions(const NumericVector& lambda) {
    if (lambda.hasAttribute("dim")) {
        IntegerVector dimensions = lambda.attr("dim");
        if (dimensions.size() > 1)
            stop("lambda must be a one-dimensional numeric vector.");
    }
}

//[[Rcpp::export]]
List hugeglasso(NumericMatrix S, NumericVector lambda, bool scr, bool verbose, bool cov_output)
{
    int d = S.nrow();
    int nlambda = lambda.size();
    if (d <= 0 || S.ncol() != d)
        stop("S must be a non-empty square matrix.");
    validate_lambda_dimensions(lambda);
    if (nlambda <= 0)
        stop("lambda must contain at least one value.");

    if (verbose) {
        if (scr) Rcout << "Conducting the graphical lasso (glasso) with lossy screening...";
        else     Rcout << "Conducting the graphical lasso (glasso) with lossless screening...";
    }

    // S is column-major in R; .begin() gives a direct const-free pointer.
    huge::GlassoResult res = huge::glasso_compact(
        S.begin(), d, lambda.begin(), nlambda, scr, cov_output);

    if (verbose) Rcout << "done.\n";

    // Scalar vectors — single allocation each, no copies
    NumericVector loglik(nlambda), sparsity(nlambda);
    IntegerVector df(nlambda);
    for (int i = 0; i < nlambda; i++) {
        loglik[i]   = res.loglik[i];
        sparsity[i] = res.sparsity[i];
        df[i]       = res.df[i];
    }

    // Write directly into R SEXP memory. Precision/covariance use one memcpy
    // each; path is derived from precision without a duplicate native cube.
    List path(nlambda), icov(nlambda), cov;
    if (cov_output) cov = List(nlambda);
    const size_t mat_elements = static_cast<size_t>(d) * d;
    const size_t mat_bytes = sizeof(double) * mat_elements;
    for (int i = 0; i < nlambda; i++) {
        NumericMatrix p(d, d);
        double* path_data = p.begin();
        const double* precision = res.icov[i].v.data();
        for (size_t entry = 0; entry < mat_elements; entry++)
            path_data[entry] = precision[entry] != 0.0 ? 1.0 : 0.0;
        for (int j = 0; j < d; j++)
            p(j, j) = 0.0;
        path[i] = p;

        NumericMatrix ic(d, d);
        std::memcpy(ic.begin(), res.icov[i].v.data(), mat_bytes);
        icov[i] = ic;

        if (cov_output) {
            NumericMatrix cv(d, d);
            std::memcpy(cv.begin(), res.cov[i].v.data(), mat_bytes);
            cov[i] = cv;
        }
    }

    List result;
    result["loglik"]   = loglik;
    result["sparsity"] = sparsity;
    result["df"]       = df;
    result["path"]     = path;
    result["icov"]     = icov;
    if (cov_output) result["cov"] = cov;
    result["hit_max_iter"] = res.hit_max_iter;
    return result;
}
