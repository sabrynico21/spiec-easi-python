#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.mb(): Meinshausen & Buhlmann graph estimation (mb)               #
#-----------------------------------------------------------------------#

#' Meinshausen & Buhlmann graph estimation
#'
#' See more details in \code{\link{huge}}
#' @param x There are 2 options: (1) \code{x} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension).
#' @param lambda A numeric scalar or non-empty one-dimensional numeric input
#'   defining a finite, strictly positive, non-increasing regularization path;
#'   tied values are allowed. Leave \code{lambda = NULL}
#'   to generate a path from \code{nlambda} and \code{lambda.min.ratio}.
#' @param nlambda The number of regularization/thresholding parameters. The default value is \code{20} for \code{method = "ct"} and \code{10} for \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}.
#' @param lambda.min.ratio If \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, it is the smallest value for \code{lambda}, as a fraction of the upperbound (\code{MAX}) of the regularization/thresholding parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. If \code{method = "ct"}, it is the largest sparsity level for estimated graphs. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda}, which makes the sparsity level of the graph path increases from \code{0} to \code{lambda.min.ratio} evenly.The default value is \code{0.1} when \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, and 0.05 when \code{method = "ct"}.
#' @param scr If \code{scr = TRUE}, the lossy screening rule is applied to preselect the neighborhood before the graph estimation. The default value is  \code{FALSE}.
#' @param scr.num The neighborhood size after the lossy screening rule (the number of remaining neighbors per node). It must be an integer between \code{1} and \code{d - 1}. ONLY applicable when \code{scr = TRUE}. The default value is \code{n-1}. An alternative value is \code{n/log(n)}. ONLY applicable when \code{scr = TRUE} and \code{method = "mb"}.
#' @param idx.mat A \code{scr.num} by \code{d} index matrix for screening. Column \code{m} must contain distinct, zero-based indices in \code{0:(d - 1)} and must exclude the response index \code{m - 1}.
#' @param sym Symmetrize the output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as the neighbor for each other. The default value is \code{"or"}. ONLY applicable when \code{method = "mb"} or \code{"tiger"}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @param input.type How to interpret \code{x}: \code{"auto"} preserves
#'   symmetry-based detection, \code{"data"} forces an observation matrix,
#'   and \code{"covariance"} requires a square covariance or correlation
#'   matrix.
#' @seealso \code{\link{huge}}, and \code{\link{huge-package}}.
#' @export
huge.mb = function(x, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, scr = NULL, scr.num = NULL, idx.mat = NULL, sym = "or", verbose = TRUE, input.type = "auto")
{
  sym = .huge_validate_sym(sym)
  inp = .huge_preprocess(x, verbose, input.type = input.type)
  S = inp$S; n = inp$n; d = inp$d
  fit = list()
  fit$cov.input = inp$cov.input

  if(is.null(idx.mat))
  {
    if(is.null(scr))
      scr = FALSE

    if(scr)
    {
      if(is.null(scr.num))
      {
        if(n<d)
          scr.num = n-1
        if(n>=d)
        {
          if(verbose) cat("lossy screening is skipped without specifying scr.num.\n")
          scr = FALSE
        }
      }
    }
    fit$scr = scr
  }

  if(!is.null(idx.mat))
  {
    if(!is.matrix(idx.mat))
      stop("idx.mat must be a matrix.")
    if(ncol(idx.mat) != d)
      stop("idx.mat must have exactly d columns.")
    if(nrow(idx.mat) < 1 || nrow(idx.mat) >= d)
      stop("idx.mat must have between 1 and d - 1 rows.")
    if(!is.numeric(idx.mat) || any(!is.finite(idx.mat)) ||
       any(idx.mat != floor(idx.mat)))
      stop("idx.mat must contain finite integer-valued indices.")
    if(any(idx.mat < 0 | idx.mat >= d))
      stop("idx.mat must contain zero-based indices in 0:(d - 1).")

    idx.mat = matrix(as.integer(idx.mat), nrow = nrow(idx.mat), ncol = d)
    for(m in seq_len(d))
    {
      if(anyDuplicated(idx.mat[,m]))
        stop("Each column of idx.mat must contain distinct indices.")
      if((m - 1L) %in% idx.mat[,m])
        stop("Each column of idx.mat must exclude its response index.")
    }

    scr = TRUE
    fit$scr = scr
    scr.num = nrow(idx.mat)
  }

  lam = .huge_default_lambda(S, d, nlambda, lambda.min.ratio, lambda)
  lambda = lam$lambda; nlambda = lam$nlambda

  if(scr)
  {
    if(length(scr.num) != 1 || !is.numeric(scr.num) ||
       !is.finite(scr.num) || scr.num != floor(scr.num) ||
       scr.num < 1 || scr.num >= d)
      stop("scr.num must be an integer between 1 and d - 1.")
    scr.num = as.integer(scr.num)

    if(verbose)
    {
      cat("Conducting Meinshausen & Buhlmann graph estimation (mb) with lossy screening....")
      flush.console()
    }

    if(is.null(idx.mat))
    {
      idx.mat = vapply(seq_len(d), function(m) {
        candidates = seq_len(d)[-m]
        ord = order(-abs(S[candidates,m]), candidates)
        candidates[ord[seq_len(scr.num)]] - 1L
      }, integer(scr.num))
      idx.mat = matrix(idx.mat, nrow = scr.num, ncol = d)
    }

    fit$idx_mat = idx.mat
    out = .Call("_huge_SPMBscr", S, lambda, nlambda, d, idx.mat, scr.num, PACKAGE="huge")
  }
  if(!scr)
  {
    if(verbose)
    {
      cat("Conducting Meinshausen & Buhlmann graph estimation (mb)....")
      flush.console()
    }
    fit$idx_mat = NULL
    out = .Call("_huge_SPMBgraph", S, lambda, nlambda, d, PACKAGE= "huge")
  }
  if (isTRUE(out$hit_max_iter))
    warning("mb solver reached its iteration limit; estimates may not be fully converged.")

  # The core emits each column's row indices in ascending order (see
  # collect_sorted in huge_core.cpp), as dgCMatrix requires.
  nnz = out$col_cnz[d + 1]
  x_values = if(nnz > 0) out$x[seq_len(nnz)] else numeric(0)
  i_values = if(nnz > 0) out$row_idx[seq_len(nnz)] else integer(0)
  G = new("dgCMatrix", Dim = as.integer(c(d*nlambda,d)),
          x = as.vector(x_values), p = as.integer(out$col_cnz),
          i = as.integer(i_values))

  fit$beta = list()
  fit$path = list()
  fit$df = matrix(0,d,nlambda)
  fit$sparsity = rep(0,nlambda)
  for(i in 1:nlambda)
  {
    fit$beta[[i]] = G[((i-1)*d+1):(i*d),,drop = FALSE]
    fit$path[[i]] = abs(fit$beta[[i]])
    fit$df[,i] = Matrix::colSums(fit$path[[i]] != 0)

    # Matrix::t on the sparse matrix (not as.matrix) keeps the update sparse
    if(sym == "or")
      fit$path[[i]] = sign(fit$path[[i]] + Matrix::t(fit$path[[i]]))
    if(sym == "and")
      fit$path[[i]] = sign(fit$path[[i]] * Matrix::t(fit$path[[i]]))
    fit$sparsity[i] = if(d <= 1) 0 else sum(fit$path[[i]])/d/(d-1)
  }
  fit$lambda = lambda

  if(verbose)
  {
     cat("done\n")
      flush.console()
  }

  return(fit)
}
