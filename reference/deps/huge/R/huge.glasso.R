#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# glasso(): The graphical lasso (glasso) using sparse matrix output     #
#-----------------------------------------------------------------------#

#' The graphical lasso (glasso) using sparse matrix output
#'
#' See more details in \code{\link{huge}}
#' @param x There are 2 options: (1) \code{x} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension).
#' @param lambda A numeric scalar or non-empty one-dimensional numeric input
#'   defining a finite, strictly positive, non-increasing regularization path;
#'   tied values are allowed. Leave \code{lambda = NULL}
#'   to generate a path from \code{nlambda} and \code{lambda.min.ratio}.
#' @param nlambda The number of regularization/thresholding parameters. The default value is \code{20} for \code{method = "ct"} and \code{10} for \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}.
#' @param lambda.min.ratio If \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, it is the smallest value for \code{lambda}, as a fraction of the upperbound (\code{MAX}) of the regularization/thresholding parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. If \code{method = "ct"}, it is the largest sparsity level for estimated graphs. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda}, which makes the sparsity level of the graph path increases from \code{0} to \code{lambda.min.ratio} evenly.The default value is \code{0.1} when \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, and 0.05 when \code{method = "ct"}.
#' @param scr If \code{scr = TRUE}, the lossy screening rule is applied to preselect the neighborhood before the graph estimation. The default value is \code{FALSE}.
#' @param cov.output If \code{cov.output = TRUE}, the output will include a path of estimated covariance matrices. ONLY applicable when \code{method = "glasso"}. Since the estimated covariance matrices are generally not sparse, please use it with care, or it may take much memory under high-dimensional setting. The default value is \code{FALSE}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @param input.type How to interpret \code{x}: \code{"auto"} preserves
#'   symmetry-based detection, \code{"data"} forces an observation matrix,
#'   and \code{"covariance"} requires a square covariance or correlation
#'   matrix. For covariance input with \code{lambda = NULL}, \code{"auto"}
#'   retains the historical diagonal-sensitive lambda scale; use
#'   \code{"covariance"} for a scale based only on off-diagonal entries.
#'   Indefinite pairwise covariance estimates are accepted only when
#'   regularization produces a certified positive-definite result.
#' @seealso \code{\link{huge}}, and \code{\link{huge-package}}.
#' @export
huge.glasso = function(x, lambda = NULL, lambda.min.ratio = NULL, nlambda = NULL, scr = NULL, cov.output = FALSE, verbose = TRUE, input.type = "auto"){

  inp = .huge_validate_estimation_input(
    x, input.type = input.type, require.psd = FALSE
  )
  x = inp$input
  n = inp$n
  d = inp$d
  cov.input = inp$cov.input
  if(cov.input)
  {
    if(verbose) cat("The input is identified as the covariance matrix.\n")
    S = inp$covariance
  }
  else
  {
    S = .huge_fast_cor(x)
    if(any(!is.finite(S)))
      stop("Raw data x cannot produce a finite correlation matrix.")
  }
  if(is.null(scr)) scr = FALSE
  lam = .huge_default_lambda(
    S, d, nlambda, lambda.min.ratio, lambda,
    legacy.glasso.covariance =
      cov.input && input.type == "auto"
  )
  lambda = lam$lambda; nlambda = lam$nlambda

  fit = .Call("_huge_hugeglasso",S,lambda,scr,verbose,cov.output,PACKAGE="huge")

  if (isTRUE(fit$hit_max_iter))
    warning("glasso solver reached its iteration limit; estimates may not be fully converged.")
  fit$hit_max_iter = NULL

  fit$scr = scr
  fit$lambda = lambda
  fit$cov.input = cov.input
  fit$cov.output = cov.output

  if(verbose){
       cat("\nConducting the graphical lasso (glasso)....done.                                          \r")
       cat("\n")
      flush.console()
  }
  return(fit)
}
