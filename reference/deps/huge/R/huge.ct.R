#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.gect(): graph estimation via correlation thresholding (ct)       #
#-----------------------------------------------------------------------#

#' Graph estimation via correlation thresholding (ct) 
#' 
#' See more details in \code{\link{huge}}
#' @details The default path targets increasing numbers of undirected edges
#'   and defines each graph by the strict rule
#'   \code{abs(correlation) > lambda}. Equal-weight edges are never split, so
#'   ties can make the realized sparsity smaller than the nominal target.
#'   Reusing the returned \code{lambda} values therefore reconstructs the
#'   same path. When \code{lambda = NULL}, \code{nlambda} must be a positive
#'   integer and \code{lambda.min.ratio} must lie in \code{(0, 1]}. Supplying
#'   \code{lambda} overrides those two arguments.
#' @param x There are 2 options: (1) \code{x} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension).
#' @param lambda A numeric scalar or non-empty one-dimensional numeric input
#'   of finite, non-negative thresholds. Values are applied in the supplied
#'   order, and zero is allowed. Leave
#'   \code{lambda = NULL} to generate a path from \code{nlambda} and
#'   \code{lambda.min.ratio}.
#' @param nlambda The number of regularization/thresholding parameters. The default value is \code{20} for \code{method = "ct"} and \code{10} for \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}.
#' @param lambda.min.ratio If \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, it is the smallest value for \code{lambda}, as a fraction of the upperbound (\code{MAX}) of the regularization/thresholding parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. If \code{method = "ct"}, it is the largest sparsity level for estimated graphs. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda}, which makes the sparsity level of the graph path increases from \code{0} to \code{lambda.min.ratio} evenly.The default value is \code{0.1} when \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, and 0.05 when \code{method = "ct"}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @param input.type How to interpret \code{x}: \code{"auto"} preserves
#'   symmetry-based detection, \code{"data"} forces an observation matrix,
#'   and \code{"covariance"} requires a square covariance or correlation
#'   matrix.
#' @seealso \code{\link{huge}}, and \code{\link{huge-package}}.
#' @export
huge.ct = function(x, nlambda = NULL, lambda.min.ratio = NULL, lambda = NULL, verbose = TRUE, input.type = "auto")
{
  inp = .huge_preprocess(x, verbose, input.type = input.type)
  S = inp$S; d = inp$d
  fit = list()
  fit$cov.input = inp$cov.input

  diag(S) = 0
  S = abs(S)

  ct.progress = function(i)
  {
    cat(paste(c("Conducting the graph estimation via correlation thresholding (ct)....in progress:", floor(100*i/nlambda), "%"), collapse=""), "\r")
    flush.console()
  }

  if(is.null(lambda))
  {
    if(is.null(nlambda))
      nlambda = 20
    else
      nlambda = .huge_validate_positive_integer(nlambda, "nlambda")
    if(is.null(lambda.min.ratio))
      lambda.min.ratio = 0.05
    else
      lambda.min.ratio = .huge_validate_ratio(lambda.min.ratio)

    edge.weights = sort(S[upper.tri(S)], decreasing = TRUE)
    edge.total = length(edge.weights)

    if(edge.total == 0)
      lambda = rep(0, nlambda)
    else
    {
      target.edges = ceiling(seq(
        1, lambda.min.ratio * edge.total, length.out = nlambda
      ))
      target.edges = pmax(0L, pmin(edge.total, as.integer(target.edges)))
      lambda = vapply(target.edges, function(target) {
        if(target < edge.total)
        {
          next.edge = max(target + 1L, 1L)
          return(edge.weights[next.edge])
        }
        0
      }, numeric(1))
    }
  }
  else
    lambda = .huge_validate_lambda(lambda, allow.zero = TRUE)

  nlambda = length(lambda)
  fit$path = list()
  fit$sparsity = rep(0,nlambda)
  for(i in seq_len(nlambda))
  {
    fit$path[[i]] = Matrix(0, d, d)
    fit$path[[i]][S > lambda[i]] = 1
    if(d <= 1)
      fit$sparsity[i] = 0
    else
      fit$sparsity[i] = sum(fit$path[[i]])/d/(d-1)
    if(verbose) ct.progress(i)
  }
  fit$lambda = lambda

  if(verbose)
  {
        cat("Conducting the graph estimation via correlation thresholding (ct)....done.             \r\n")
        flush.console()
    }
  return(fit)
}
