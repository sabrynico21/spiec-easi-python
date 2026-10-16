// Rcpp thin wrapper for MB graph estimation — delegates to huge::mb() / huge::mb_scr()
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

// Helper: convert core ColResult vector to R sparse-like output
static List mb_result_to_list(const huge::MBResult& res, int d) {
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
    return List::create(
        _["col_cnz"] = col_cnz,
        _["row_idx"] = row_idx,
        _["x"] = x,
        _["hit_max_iter"] = res.hit_max_iter
    );
}

//[[Rcpp::export]]
List SPMBscr(NumericMatrix S, NumericVector lambda, int nlambda, int d, IntegerMatrix idx_scr, int nscr)
{
    if (d <= 0)
        stop("d must be positive.");
    if (nlambda <= 0)
        stop("nlambda must be positive.");
    if (S.nrow() != d || S.ncol() != d)
        stop("S must be a d by d matrix.");
    validate_lambda_dimensions(lambda);
    if (lambda.size() != nlambda)
        stop("lambda length must equal nlambda.");
    if (nscr <= 0 || nscr >= d)
        stop("idx_scr must have between 1 and d - 1 rows.");
    if (idx_scr.nrow() != nscr || idx_scr.ncol() != d)
        stop("idx_scr must be an nscr by d matrix.");

    huge::RThreadLimitGuard thread_limit;
    huge::MBResult res = huge::mb_scr(S.begin(), d, lambda.begin(), nlambda,
                                      idx_scr.begin(), nscr);
    return mb_result_to_list(res, d);
}

//[[Rcpp::export]]
List SPMBgraph(NumericMatrix S, NumericVector lambda, int nlambda, int d)
{
    if (d <= 0)
        stop("d must be positive.");
    if (nlambda <= 0)
        stop("nlambda must be positive.");
    if (S.nrow() != d || S.ncol() != d)
        stop("S must be a d by d matrix.");
    validate_lambda_dimensions(lambda);
    if (lambda.size() != nlambda)
        stop("lambda length must equal nlambda.");

    huge::RThreadLimitGuard thread_limit;
    huge::MBResult res = huge::mb(S.begin(), d, lambda.begin(), nlambda);
    return mb_result_to_list(res, d);
}
