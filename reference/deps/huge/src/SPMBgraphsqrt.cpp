// Rcpp thin wrapper for TIGER estimation — delegates to huge::tiger()
#include <Rcpp.h>
#include "huge/huge_core.h"
#include "huge_r_threads.h"
using namespace Rcpp;

static void validate_lambda_dimensions(const NumericVector& lambda) {
    if (lambda.hasAttribute("dim")) {
        IntegerVector dimensions = lambda.attr("dim");
        if (dimensions.size() > 1)
            stop("lambda must be a one-dimensional numeric vector.");
    }
}

static List tiger_result_to_list(const huge::TigerResult& res, int d)
{
    int nlambda = static_cast<int>(res.icov.size());
    int total_nnz = 0;
    for (int m = 0; m < d; m++) total_nnz += res.columns[m].vals.size();
    NumericVector x(total_nnz);
    IntegerVector col_cnz(d + 1), row_idx(total_nnz);
    col_cnz[0] = 0;
    int cnz = 0;
    for (int m = 0; m < d; m++) {
        const huge::ColResult& col = res.columns[m];
        for (size_t j = 0; j < col.vals.size(); j++) {
            x[cnz] = col.vals[j];
            row_idx[cnz] = col.indices[j];
            cnz++;
        }
        col_cnz[m + 1] = cnz;
    }

    List icov_list(nlambda);
    const size_t mat_bytes = sizeof(double) * static_cast<size_t>(d) * d;
    for (int i = 0; i < nlambda; i++) {
        NumericMatrix ic(d, d);
        std::memcpy(ic.begin(), res.icov[i].v.data(), mat_bytes);
        icov_list[i] = ic;
    }

    return List::create(
        _["col_cnz"] = col_cnz,
        _["row_idx"] = row_idx,
        _["x"] = x,
        _["icov"] = icov_list,
        _["lambda"] = wrap(res.lambda),
        _["hit_max_iter"] = res.hit_max_iter,
        _["path_truncated"] = res.path_truncated
    );
}

//[[Rcpp::export]]
List SPMBgraphsqrt(NumericMatrix data, NumericVector lambda, int nlambda, int d)
{
    validate_lambda_dimensions(lambda);
    int n = data.nrow();
    if (n <= 0 || d <= 0 || nlambda <= 0) {
        int d_safe = d > 0 ? d : 0;
        return List::create(
            _["col_cnz"] = IntegerVector(d_safe + 1),
            _["row_idx"] = IntegerVector(0),
            _["x"] = NumericVector(0),
            _["icov"] = List(0)
        );
    }
    if (data.ncol() != d)
        stop("tiger input dimensions are invalid.");
    if (lambda.size() != nlambda)
        stop("tiger lambda length must match nlambda.");

    huge::RThreadLimitGuard thread_limit;
    huge::TigerResult res = huge::tiger(data.begin(), n, d, lambda.begin(), nlambda);
    return tiger_result_to_list(res, d);
}

// Native correlation/covariance-first TIGER entry point.  A NULL lambda asks
// the shared core to generate the default path from its correlation matrix.
//[[Rcpp::export]]
List SPMBgraphsqrtFit(NumericMatrix input, Nullable<NumericVector> lambda,
                      int nlambda, int d, bool covariance_input,
                      double lambda_min_ratio)
{
    int n = input.nrow();
    if (n <= 0 || d <= 0 || input.ncol() != d)
        stop("tiger input dimensions are invalid.");

    std::vector<double> supplied_lambda;
    const double* lambda_ptr = nullptr;
    if (lambda.isNotNull()) {
        NumericVector values(lambda);
        validate_lambda_dimensions(values);
        supplied_lambda.assign(values.begin(), values.end());
        nlambda = static_cast<int>(supplied_lambda.size());
        lambda_ptr = supplied_lambda.data();
    }

    huge::RThreadLimitGuard thread_limit;
    huge::TigerResult res = huge::tiger_fit(
        input.begin(), n, d, covariance_input, lambda_ptr, nlambda,
        lambda_min_ratio);
    return tiger_result_to_list(res, d);
}
