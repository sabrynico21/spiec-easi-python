#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.npn(): nonparanormal(npn) transformation                         #
#-----------------------------------------------------------------------#

#' Nonparanormal(npn) transformation
#'
#' Implements the Gaussianization to help relax the assumption of normality.
#'
#' The nonparanormal extends Gaussian graphical models to semiparametric Gaussian copula models.Motivated by sparse additive models, the nonparanormal method estimates the Gaussian copula by marginally transforming the variables using smooth functions.Computationally, the estimation of a nonparanormal transformation is very efficient and only requires one pass of the data matrix.
#'
#' @param x The \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} dimensions
#' @param npn.func The transformation function used in the npn transformation. If \code{npn.func = "truncation"}, the truncated ECDF is applied. If \code{npn.func = "shrinkage"}, the shrunken ECDF is applied. The default is \code{"shrinkage"}. If \code{npn.func = "skeptic"}, the nonparanormal skeptic is applied.
#' @param npn.thresh The truncation threshold used in nonparanormal transformation, ONLY applicable when \code{npn.func = "truncation"}. The default value is \code{1/(4*(n^0.25)*} \code{sqrt(pi*log(n)))}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @param na.last for controlling the treatment of NAs. If TRUE, missing values in the data are put last; if FALSE, they are put first; if NA, they are removed; if "keep" they are kept with rank NA. See also \code{\link{rank}}.
#' @return
#' \item{data}{
#' A \code{d} by \code{d} nonparanormal correlation matrix if \code{npn.func = "skeptic"}, and A \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} transformed dimensions otherwise.
#' }
#' @seealso \code{\link{huge}} and \code{\link{huge-package}}.
#' @examples
#' # generate nonparanormal data
#' L = huge.generator(graph = "cluster", g = 5)
#' L$data = L$data^5
#'
#' # transform the data using the shrunken ECDF
#' Q = huge.npn(L$data)
#'
#' # transform the non-Gaussian data using the truncated ECDF
#' Q = huge.npn(L$data, npn.func = "truncation")
#'
#' # estimate the nonparanormal correlation matrix using the skeptic estimator
#' Q = huge.npn(L$data, npn.func = "skeptic")
#' @export
huge.npn = function(x, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE, na.last = "keep"){
  n = nrow(x)
  x.col = colnames(x)
  x.row = rownames(x)

  if(!npn.func %in% c("shrinkage", "truncation", "skeptic"))
    stop("npn.func must be one of \"shrinkage\", \"truncation\", \"skeptic\".")

  normalize.columns = function(z){
    z.sd = apply(z, 2, sd, na.rm = TRUE)
    z.sd[!is.finite(z.sd) | z.sd == 0] = 1
    sweep(z, 2, z.sd, "/")
  }

  if(verbose) cat(sprintf("Conducting the nonparanormal (npn) transformation via %s....",
    switch(npn.func, shrinkage = "shrunken ECDF", truncation = "truncated ECDF", skeptic = "skeptic")))

  if(npn.func == "shrinkage"){
    # Shrunken ECDF: rank/(n+1) keeps the transform inside qnorm's open support
    x = normalize.columns(qnorm(apply(x,2,rank,na.last=na.last)/(n+1)))
    dimnames(x) = list(x.row, x.col)
  }

  if(npn.func == "truncation"){
    if(is.null(npn.thresh)) npn.thresh = 1/(4*(n^0.25)*sqrt(pi*log(n)))
    x = normalize.columns(qnorm(pmin(pmax(apply(x,2,rank,na.last=na.last)/n, npn.thresh), 1-npn.thresh)))
    dimnames(x) = list(x.row, x.col)
  }

  if(npn.func == "skeptic"){
    # Spearman correlation = Pearson correlation of column ranks; computing it
    # as a BLAS crossprod of standardized ranks is ~5x faster than
    # cor(method = "spearman") at d = 2000. Fall back for NA input.
    if(anyNA(x)){
      rho = cor(x, method = "spearman", use = "pairwise.complete.obs")
    } else {
      rk = scale(apply(x, 2, rank))
      rho = crossprod(rk) / (n - 1)
      diag(rho) = 1
    }
    x = 2*sin(pi/6*rho)
    # d x d correlation matrix: variable names on both dimensions
    dimnames(x) = list(x.col, x.col)
  }

  if(verbose) cat("done.\n")
  return(x)
}
