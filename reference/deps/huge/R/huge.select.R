#-------------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                   #
# huge.select(): Model selection using:                                   #
#                1.rotation information criterion (ric)                   #
#                2.stability approach to regularization selection (stars) #
#                3.extended Bayesian information criterion (ebic)        #
#-------------------------------------------------------------------------#

#' Model selection for high-dimensional undirected graph estimation
#'
#' Implements the regularization parameter selection for high dimensional undirected graph estimation. The optional approaches are rotation information criterion (ric), stability approach to regularization selection (stars) and extended Bayesian information criterion (ebic).
#'
#' Stability approach to regularization selection (stars) selects the optimal graph by variability of subsamplings and tends to overselect edges in Gaussian graphical models. It is available for \code{"mb"}, \code{"ct"}, and \code{"glasso"}. TIGER can certify different lambda-path prefixes on different subsamples, so TIGER with stars is rejected until a common certified-prefix protocol is available; use \code{"ric"} for TIGER. Besides selecting the regularization parameters, stars can also provide an additional estimated graph by merging the corresponding subsampled graphs using the frequency counts. The subsampling procedure in stars may NOT be very efficient, we also provide the recent developed highly efficient, rotation information criterion approach (ric). Instead of tuning over a grid by cross-validation or subsampling, we directly estimate the optimal regularization parameter based on random Rotations. However, ric usually has very good empirical performances but suffers from underselections sometimes. Therefore, we suggest if user are sensitive of false negative rates, they should either consider increasing \code{rep.num} or applying the stars to model selection where supported. Extended Bayesian information criterion (ebic) is another competitive approach, but the \code{ebic.gamma} can only be tuned by experience.
#'
#' For \code{criterion = "stars"}, \code{est$lambda} must be non-increasing; tied values are allowed.
#'
#' @param est An object with S3 class \code{"huge"}.
#' @param criterion Model selection criterion. \code{"ric"} is available for all four estimation methods, \code{"stars"} for \code{"mb"}, \code{"ct"}, and \code{"glasso"}, and \code{"ebic"} only for \code{"glasso"}. Defaults are \code{"ric"} for \code{"mb"} and \code{"tiger"}, \code{"stars"} for \code{"ct"}, and \code{"ebic"} for \code{"glasso"}.
#' @param ebic.gamma The tuning parameter for ebic. The default value is 0.5. Only applicable when \code{est$method = "glasso"} and \code{criterion = "ebic"}.
#' @param stars.thresh The variability threshold in stars. The default value is \code{0.1}. An alternative value is \code{0.05}. Only applicable when \code{criterion = "stars"}.
#' @param stars.subsample.ratio The subsampling ratio. The default value is \code{10*sqrt(n)/n} when \code{n>144} and \code{0.8} when \code{n<=144}, where \code{n} is the sample size. Only applicable when \code{criterion = "stars"}.
#' @param rep.num The number of subsamplings when \code{criterion = "stars"} or rotations when \code{criterion = "ric"}. The default value is \code{20}. NOT applicable when \code{criterion = "ebic"}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @param num.cores The number of CPU cores used to fit the stars subsamplings in parallel via \code{parallel::mclapply}. The default value is \code{1} (serial). At most two forked workers are used, and this argument is ignored on Windows. Each forked worker limits huge's native OpenMP code to one thread, but an external threaded BLAS may still create additional threads; \code{num.cores = 1} is the portable choice when a bounded thread budget or fork safety matters. Results are identical to the serial path for the same random seed. Only applicable when \code{criterion = "stars"}.
#' @return
#' An object with S3 class "select" is returned:
#'   \item{refit}{
#'     The optimal graph selected from the graph path
#'   }
#' \item{opt.icov}{
#'   The optimal precision matrix from the path only applicable when \code{method = "glasso"}
#' }
#' \item{opt.cov}{
#'   The optimal covariance matrix from the path only applicable when \code{method = "glasso"} and \code{est$cov} is available.
#' }
#' \item{merge}{
#'   The graph path estimated by merging the subsampling paths. Only applicable when the input \code{criterion = "stars"}.
#' }
#' \item{variability}{
#'   The variability along the subsampling paths. Only applicable when the input \code{criterion = "stars"}.
#' }
#' \item{ebic.score}{
#'   Extended BIC scores for regularization parameter selection. Only applicable when \code{criterion = "ebic"}.
#' }
#' \item{opt.index}{
#'   The index of the selected regularization parameter. NOT applicable when the input \code{criterion = "ric"}
#' }
#' \item{opt.lambda}{
#'   The selected regularization/thresholding parameter.
#' }
#' \item{opt.sparsity}{
#'   The sparsity level of \code{"refit"}.
#' }
#'
#' and anything else included in the input \code{est}
#'
#' @note The model selection is NOT available when the data input is the sample covariance matrix.
#' @seealso \code{\link{huge}} and \code{\link{huge-package}}.
#' @examples
#' #generate data
#' L = huge.generator(d = 20, graph="hub")
#' out.mb = huge(L$data)
#' out.ct = huge(L$data, method = "ct")
#' out.glasso = huge(L$data, method = "glasso")
#'
#' #model selection using ric
#' out.select = huge.select(out.mb)
#' plot(out.select)
#'
#' #model selection using stars
#' #out.select = huge.select(out.ct, criterion = "stars", stars.thresh = 0.05,rep.num=10)
#' #plot(out.select)
#'
#' #model selection using ebic
#' out.select = huge.select(out.glasso,criterion = "ebic")
#' plot(out.select)
#' @export
huge.select = function(est, criterion = NULL, ebic.gamma = 0.5, stars.thresh = 0.1, stars.subsample.ratio = NULL, rep.num = 20, verbose = TRUE, num.cores = 1){

  if(!is.list(est) || !inherits(est, "huge"))
    stop("est must be an object returned by huge().")

  required = c("method", "cov.input", "data", "lambda", "path", "sparsity")
  missing = required[vapply(required, function(name)
    is.null(est[[name]]), logical(1))]
  if(length(missing) > 0)
    stop(sprintf(
      "est is incomplete; missing field(s): %s.",
      paste(missing, collapse = ", ")
    ))

  if(length(est$method) != 1 || !is.character(est$method) ||
     !(est$method %in% c("mb", "ct", "glasso", "tiger")))
    stop("est has an invalid method field.")
  if(length(est$cov.input) != 1 || !is.logical(est$cov.input) ||
     is.na(est$cov.input))
    stop("est has an invalid cov.input field.")
  if(!is.matrix(est$data) || !is.numeric(est$data) ||
     nrow(est$data) < 1 || ncol(est$data) < 1)
    stop("est has an invalid data matrix.")
  if(!is.numeric(est$lambda) || length(est$lambda) < 1 ||
     any(!is.finite(est$lambda)))
    stop("est has an invalid lambda path.")

  nlambda = length(est$lambda)
  d = ncol(est$data)
  if(!is.list(est$path) || length(est$path) != nlambda ||
     any(!vapply(est$path, function(path)
       length(dim(path)) == 2 && all(dim(path) == c(d, d)), logical(1))))
    stop("est has an invalid graph path.")
  if(!is.numeric(est$sparsity) || length(est$sparsity) != nlambda ||
     any(!is.finite(est$sparsity)))
    stop("est has an invalid sparsity path.")
  if(est$method %in% c("mb", "tiger"))
  {
    if(length(est$sym) != 1 || !is.character(est$sym) || is.na(est$sym) ||
       !(est$sym %in% c("or", "and")))
      stop("est is missing valid symmetrization metadata.")
  }
  if(est$method %in% c("mb", "glasso"))
  {
    if(length(est$scr) != 1 || !is.logical(est$scr) || is.na(est$scr))
      stop("est is missing screening metadata.")
  }

  if(est$cov.input)
    stop("Model selection is not available when using the covariance matrix as input.")

  if(!est$cov.input)
  {
    if(is.null(criterion))
      criterion = switch(est$method,
        mb = "ric", ct = "stars", glasso = "ebic", tiger = "ric")
    else if(length(criterion) != 1 || !is.character(criterion) ||
            is.na(criterion) ||
            !(criterion %in% c("ric", "stars", "ebic")))
      stop("criterion must be exactly one of \"ric\", \"stars\", or \"ebic\".")

    if(criterion == "ebic" && est$method != "glasso")
      stop("criterion = \"ebic\" is available only for method = \"glasso\".")
    if(criterion == "stars" && est$method == "tiger")
      stop(paste(
        "TIGER with StARS is unavailable until subsample fits can share a",
        "common certified prefix; use criterion = \"ric\" for TIGER."
      ))
    if(criterion == "stars" && nlambda > 1 &&
       any(diff(est$lambda) > 0))
      stop(paste(
        "StARS requires est$lambda to be non-increasing;",
        "refit with lambda in decreasing order."
      ))

    if(criterion %in% c("ric", "stars"))
      rep.num = .huge_validate_positive_integer(rep.num, "rep.num")
    if(criterion == "stars")
    {
      num.cores = .huge_validate_positive_integer(num.cores, "num.cores")
      stars.thresh = .huge_validate_ratio(stars.thresh, "stars.thresh")
      if(!is.null(stars.subsample.ratio))
        stars.subsample.ratio = .huge_validate_ratio(
          stars.subsample.ratio, "stars.subsample.ratio"
        )
    }
    if(criterion == "ebic")
    {
      if(length(ebic.gamma) != 1 || !is.numeric(ebic.gamma) ||
         !is.finite(ebic.gamma))
        stop("ebic.gamma must be a finite numeric value.")
      ebic.gamma = as.numeric(ebic.gamma)
      if(!is.numeric(est$loglik) || length(est$loglik) != nlambda ||
         any(!is.finite(est$loglik)) ||
         !is.numeric(est$df) || length(est$df) != nlambda ||
         any(!is.finite(est$df)) ||
         !is.list(est$icov) || length(est$icov) != nlambda)
        stop("est is incomplete or invalid for criterion = \"ebic\".")
    }
    if(criterion == "stars" && d < 2)
      stop("StARS requires at least two variables.")
    if(criterion == "stars" && .Platform$OS.type != "windows" &&
       min(num.cores, rep.num, 2L) > 1)
      warning(paste(
        "num.cores > 1 uses at most two forked StARS workers. Each worker",
        "limits huge's native OpenMP code to one thread, but an external",
        "threaded BLAS may still create additional threads. Use",
        "num.cores = 1 for the portable single-worker path."
      ), call. = FALSE)

    n = nrow(est$data)

    # One fitting closure shared by the RIC refit and the stars subsampling:
    # refits est$method on `data` over `lambda` with the original fit's
    # screening/symmetrization settings. (This also fixes a historical
    # inconsistency where the RIC glasso refit dropped `scr` unless
    # cov.output was requested.) Selection is unavailable for covariance
    # input, so force raw-data routing when a square symmetric sample matrix
    # is refitted.
    refit.fn = function(data, lambda, cov.output = FALSE) {
      switch(est$method,
        mb     = huge.mb(data, lambda = lambda, scr = est$scr, idx.mat = est$idx_mat, sym = est$sym, verbose = FALSE, input.type = "data"),
        ct     = huge.ct(data, lambda = lambda, verbose = FALSE, input.type = "data"),
        glasso = huge.glasso(data, lambda = lambda, scr = est$scr, cov.output = cov.output, verbose = FALSE, input.type = "data"),
        tiger  = huge.tiger(data, lambda = lambda, sym = est$sym, verbose = FALSE, input.type = "data"))
    }

    if(criterion == "ric")
    {
      if(verbose)
      {
        cat("Conducting rotation information criterion (ric) selection....")
        flush.console()
      }

      if(n>rep.num){
        nr = rep.num
        r = sample(n,rep.num)
      }
      if(n<=rep.num){
        nr = n
        r = 1:n
      }

      # RIC must see standardized data: the lambda path is defined on the
      # correlation scale, and the rotated inner products are otherwise
      # scale-dependent (multiplying x by c scales opt.lambda by c^2).
      standardized.data = .huge_standardize(est$data)
      est$opt.lambda = .Call(
        "_huge_RIC", standardized.data, d, n, r, nr
      ) * 1.0 / n
      if(verbose){
        cat("done\n")
        flush.console()
      }

      if(verbose)
      {
        cat("Computing the optimal graph....")
        flush.console()
      }

      # A BLAS dot product can leave a tiny residual for mathematically
      # orthogonal columns.  Suppress only values within a conservative,
      # pair-specific dot-product roundoff bound; a fixed tolerance would
      # erase genuine weak correlations.
      denominator = max(n - 1, 1)
      correlation = crossprod(standardized.data) / denominator
      absolute.cross = crossprod(abs(standardized.data)) / denominator
      scaled.eps = n * .Machine$double.eps
      dot.gamma = if(scaled.eps < 1) scaled.eps / (1 - scaled.eps) else Inf
      cor.offdiag = abs(correlation)
      roundoff.bound = dot.gamma * absolute.cross +
        .Machine$double.eps * cor.offdiag
      cor.offdiag[cor.offdiag <= roundoff.bound] = 0
      diag(cor.offdiag) = 0
      if(est$opt.lambda>=max(cor.offdiag))
        est$refit = Matrix(0,d,d)
      else{
        # [[ ]] for exact-name access: est$cov would partial-match cov.input
        # on fits that lack a cov field (all methods except glasso+cov.output).
        refit.lambda = est$opt.lambda
        # Route on "too small for the solver to certify", not on a bitwise
        # zero.  A conforming BLAS may return any value inside the rotated
        # inner product's roundoff interval, so RIC's optimum can be a tiny
        # positive residual rather than exactly 0 (observed with ATLAS).
        # Comparing against the same roundoff scale keeps the fallback, and
        # therefore the selected graph, identical across BLAS implementations.
        zero.proxy = est$opt.lambda <= max(cor.offdiag) * dot.gamma &&
          est$method != "ct"
        if(zero.proxy)
          refit.lambda = max(est$opt.lambda, .Machine$double.xmin)

        tmp = if(zero.proxy) {
          tryCatch(
            refit.fn(
              est$data, refit.lambda,
              cov.output = !is.null(est[["cov"]])
            ),
            error = function(error) NULL
          )
        } else {
          refit.fn(
            est$data, refit.lambda,
            cov.output = !is.null(est[["cov"]])
          )
        }

        if(is.null(tmp)) {
          nearest = which.min(abs(est$lambda - est$opt.lambda))
          est$refit = est$path[[nearest]]
          if(!is.null(est[["icov"]]))
            est$opt.icov = est[["icov"]][[nearest]]
          if(!is.null(est[["cov"]]))
            est$opt.cov = est[["cov"]][[nearest]]
          warning(sprintf(
            paste(
              "RIC selected lambda = 0 (within dot-product roundoff), but the",
              "method could not certify the smallest positive proxy; the",
              "original fitted path at lambda %.17g was used."
            ),
            est$lambda[[nearest]]
          ), call. = FALSE)
        } else {
          est$refit = tmp$path[[1]]
          if(!is.null(tmp[["icov"]]))
            est$opt.icov = tmp[["icov"]][[1]]
          if(!is.null(est[["cov"]]) && !is.null(tmp[["cov"]]))
            est$opt.cov = tmp[["cov"]][[1]]
        }
      }
      est$opt.sparsity=if(d <= 1) 0 else sum(est$refit)/d/(d-1)

      if(verbose){
        cat("done\n")
        flush.console()
      }
    }

    if(criterion == "ebic"&&est$method == "glasso")
    {
      if(verbose)
      {
        cat("Conducting extended Bayesian information criterion (ebic) selection....")
        flush.console()
      }
      est$ebic.score = -n*est$loglik + log(n)*est$df + 4*ebic.gamma*log(d)*est$df
      est$opt.index = which.min(est$ebic.score)
      est$refit = est$path[[est$opt.index]]
      est$opt.icov = est$icov[[est$opt.index]]
      if(est$cov.output)
        est$opt.cov = est$cov[[est$opt.index]]
      est$opt.lambda = est$lambda[est$opt.index]
      est$opt.sparsity = est$sparsity[est$opt.index]
      if(verbose)
      {
        cat("done\n")
        flush.console()
      }
    }

    if(criterion == "stars"){
      if(is.null(stars.subsample.ratio))
      {
        if(n>144) stars.subsample.ratio = 10*sqrt(n)/n
        if(n<=144) stars.subsample.ratio = 0.8
      }
      subsample.size = floor(n * stars.subsample.ratio)
      if(subsample.size < 2)
        stop("stars.subsample.ratio must select at least two observations.")

      # Draw all subsample index sets up front: the fitting functions consume
      # no RNG, so this yields the same subsamples as drawing inside the loop
      # (bitwise-identical results), and makes the fits order-independent so
      # they can run in parallel.
      ind.list = lapply(1:rep.num, function(i)
        sample(c(1:n), subsample.size, replace=FALSE))

      subsample.fit = function(ind.sample)
        refit.fn(est$data[ind.sample,], est$lambda)$path

      # This wrapper runs only in fork children. The native R adapter reads
      # the private marker and temporarily limits huge's OpenMP regions to
      # one thread, avoiding nested package-owned parallelism.
      subsample.fit.fork = function(ind.sample) {
        previous = Sys.getenv("HUGE_R_FORK_WORKER", unset = NA_character_)
        on.exit({
          if(is.na(previous))
            Sys.unsetenv("HUGE_R_FORK_WORKER")
          else
            Sys.setenv(HUGE_R_FORK_WORKER = previous)
        }, add = TRUE)
        Sys.setenv(HUGE_R_FORK_WORKER = "1")
        subsample.fit(ind.sample)
      }

      # Fork-based parallelism (not on Windows); never start more workers
      # than there are independent subsample fits or CRAN's two-core limit.
      # num.cores = 1 keeps the serial path with per-replication progress.
      use.cores = min(num.cores, rep.num, 2L)
      if(.Platform$OS.type == "windows") use.cores = 1

      est$merge = vector("list", nlambda)
      for(lambda_idx in 1:nlambda)
        est$merge[[lambda_idx]] = Matrix(0,d,d)

      if(use.cores > 1) {
        if(verbose)
        {
          cat(sprintf("Conducting Subsampling....on %d cores", use.cores), "\r")
          flush.console()
        }
        paths = parallel::mclapply(
          ind.list, subsample.fit.fork, mc.cores = use.cores,
          mc.cleanup = TRUE, mc.allow.recursive = FALSE
        )
        failed = vapply(paths, function(p) inherits(p, "try-error") || is.null(p), logical(1))
        if(any(failed))
          stop("stars subsampling failed in a parallel worker; retry with num.cores = 1")
        for(lambda_idx in 1:nlambda)
          for(rep_idx in 1:rep.num)
            est$merge[[lambda_idx]] = est$merge[[lambda_idx]] + paths[[rep_idx]][[lambda_idx]]
      } else {
        for(rep_idx in 1:rep.num)
        {
          if(verbose)
          {
            mes <- paste(c("Conducting Subsampling....in progress:", floor(100*rep_idx/rep.num), "%"), collapse="")
            cat(mes, "\r")
            flush.console()
          }
          current.path = subsample.fit(ind.list[[rep_idx]])
          for(lambda_idx in 1:nlambda)
            est$merge[[lambda_idx]] = est$merge[[lambda_idx]] + current.path[[lambda_idx]]
          rm(current.path)
        }
      }

      if(verbose)
      {
        mes = "Conducting Subsampling....done.                 "
        cat(mes, "\r")
        cat("\n")
        flush.console()
      }

      est$variability = rep(0,nlambda)
      for(lambda_idx in 1:nlambda){
        est$merge[[lambda_idx]] = est$merge[[lambda_idx]]/rep.num
        # sum(m) - sum(m^2) == sum(m*(1-m)) but avoids densifying (1-m)
        m.tmp = est$merge[[lambda_idx]]
        est$variability[lambda_idx] = 4*(sum(m.tmp) - sum(m.tmp^2))/(d*(d-1))
      }

      stars.cross = which(est$variability >= stars.thresh)
      if(length(stars.cross) == 0)
        est$opt.index = nlambda
      else
        est$opt.index = max(stars.cross[1]-1,1)
      est$refit = est$path[[est$opt.index]]
      est$opt.lambda = est$lambda[est$opt.index]
      est$opt.sparsity = est$sparsity[est$opt.index]
      if(est$method == "glasso")
      {
        est$opt.icov = est$icov[[est$opt.index]]
        # [[ ]] avoids est$cov partial-matching cov.input/cov.output
        if(!is.null(est[["cov"]]))
          est$opt.cov = est[["cov"]][[est$opt.index]]
      }
      if(est$method == "tiger")
        est$opt.icov = est$icov[[est$opt.index]]
    }
    est$criterion = criterion
    class(est) = "select"
    return(est)
  }
}

#-----------------------------------------------------------------------#
# default printing function for class "select"                          #
#-----------------------------------------------------------------------#

#' Print function for S3 class "select"
#'
#' Print the information about the model usage, graph dimension, model selection criterion, sparsity level of the optimal graph.
#'
#' @param x An object with S3 class \code{"select"}.
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{huge.select}}
#' @export
print.select = function(x, ...)
{
  if(x$cov.input){
    cat("Model selection is not available when using the covariance matrix as input.")
  }
  if(!x$cov.input)
  {
    if(x$method == "ct")
      cat("Model: graph estimation via correlation thresholding (ct)\n")
    if(x$method == "glasso")
      cat("Model: graphical lasso (glasso)\n")
    if(x$method == "mb")
      cat("Model: Meinshausen & Buhlmann Graph Estimation (mb)\n")
    if(x$method == "tiger")
      cat("Model: tuning-insensitive approach (tiger)\n")

    cat("selection criterion:",x$criterion,"\n")
    if((x$method != "ct")&&x$scr)
      cat("lossy screening: on\n")
    cat("Graph dimension:",ncol(x$data),"\n")
    cat("sparsity level", x$opt.sparsity,"\n")
  }
}

#' Plot function for S3 class "select"
#'
#' Plot the optimal graph by model selection.
#'
#' @param x An object with S3 class \code{"select"}
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{huge.select}}
#' @export
plot.select = function(x, ...){
  if(x$cov.input){
    cat("Model selection is not available when using the covariance matrix as input.")
  }
  if(!x$cov.input)
  {
    old.par = .huge_graphics_state()
    on.exit(.huge_restore_graphics_state(old.par), add = TRUE)
    par(mfrow=c(1,2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))

      g = graph_from_adjacency_matrix(as.matrix(x$refit), mode="undirected", diag=FALSE)
    layout.grid = layout_with_fr(g)

    plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA)
      plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Solution path sparsity levels")
      lines(x$opt.lambda,x$opt.sparsity,type = "p")
    }
}
