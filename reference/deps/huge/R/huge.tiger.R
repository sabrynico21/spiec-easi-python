#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.tiger(): Tuning-insensitive graph estimation                     #
#-----------------------------------------------------------------------#

#' Tuning-insensitive graph estimation
#'
#' See more details in \code{\link{huge}}
#' @details Raw observations are centered and normalized, while covariance
#'   input is converted to a correlation matrix, inside the shared C++ core.
#'   When \code{lambda = NULL}, the same native correlation matrix determines
#'   the returned default lambda path. A user-supplied path must be finite,
#'   strictly positive, and non-increasing; tied values are allowed and are
#'   used unchanged. If smaller generated values cannot be certified to the solver's
#'   KKT tolerance, the function warns and returns the longest certified path
#'   prefix. A user-supplied value that cannot be certified raises an error.
#' @param x There are 2 options: (1) \code{x} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension).
#' @param lambda A numeric scalar or non-empty one-dimensional numeric input
#'   defining a finite, strictly positive, non-increasing regularization path.
#'   Tied values are allowed. Leave \code{lambda = NULL} to generate
#'   the path from \code{nlambda} and \code{lambda.min.ratio} in C++.
#' @param nlambda The number of regularization/thresholding parameters. The default value is \code{20} for \code{method = "ct"} and \code{10} for \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}.
#' @param lambda.min.ratio If \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, it is the smallest value for \code{lambda}, as a fraction of the upperbound (\code{MAX}) of the regularization/thresholding parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. If \code{method = "ct"}, it is the largest sparsity level for estimated graphs. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda}, which makes the sparsity level of the graph path increases from \code{0} to \code{lambda.min.ratio} evenly.The default value is \code{0.1} when \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, and 0.05 when \code{method = "ct"}.
#' @param sym Symmetrize the output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as the neighbor for each other. The default value is \code{"or"}. ONLY applicable when \code{method = "mb"} or \code{"tiger"}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @param input.type How to interpret \code{x}: \code{"auto"} preserves
#'   symmetry-based detection, \code{"data"} forces an observation matrix,
#'   and \code{"covariance"} requires a square covariance or correlation
#'   matrix. Correlation construction, covariance validation, and automatic
#'   lambda selection then occur together in C++.
#' @seealso \code{\link{huge}}, and \code{\link{huge-package}}.
#' @export
huge.tiger = function(x, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, sym = "or", verbose = TRUE, input.type = "auto")
{
		sym = .huge_validate_sym(sym)
		inp = .huge_validate_estimation_input(
			x, input.type = input.type, prepare.covariance = FALSE
		)
		d = inp$d
		cov.input = inp$cov.input
		if(cov.input && verbose)
			cat("The input is identified as the covariance matrix.\n")
		fit = list()
		fit$cov.input = cov.input

		if(!is.null(lambda))
		{
			lambda = .huge_validate_lambda(lambda)
			if(length(lambda) > 1L &&
			   any(lambda[-1L] > lambda[-length(lambda)]))
				stop("lambda must be non-increasing for method = \"tiger\"; tied values are allowed.")
			nlambda = length(lambda)
			# The native fixed-lambda branch does not use this argument.
			lambda.min.ratio = 0.1
		}
		else
		{
			if(is.null(nlambda))
				nlambda = 10
			else
				nlambda = .huge_validate_positive_integer(nlambda, "nlambda")
			if(is.null(lambda.min.ratio))
				lambda.min.ratio = 0.1
			else
				lambda.min.ratio = .huge_validate_ratio(lambda.min.ratio)
		}

	if(verbose)
	{
	  cat("Conducting graph estimation through a tuning-insensitive approach (tiger)....")
	  flush.console()
		}
		fit$idx_mat = NULL
		fit$scr = FALSE
		requested.nlambda = nlambda
		out = .Call("_huge_SPMBgraphsqrtFit", inp$input, lambda, as.integer(nlambda), d,
		            cov.input, lambda.min.ratio, PACKAGE = "huge")
		lambda = out$lambda
		nlambda = length(lambda)

	if (isTRUE(out$path_truncated))
	  warning(sprintf(
	    "tiger returned the %d-value certified prefix of the %d-value native lambda path; smaller lambda values did not converge or were numerically degenerate.",
	    nlambda, requested.nlambda))
	else if (isTRUE(out$hit_max_iter))
	  warning("tiger solver reached its iteration limit; estimates may not be fully converged.")

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
	fit$icov = out$icov
 	fit$lambda = lambda

	if(verbose)
 	{
 		cat("done\n")
      flush.console()
  }

 	return(fit)
}
