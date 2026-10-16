# Internal helper functions (not exported)

# Base graphics parameters have dependencies: layouts reset text scaling,
# outer regions reset figures, and margins reset plot regions. Reapply linked
# groups in dependency order so callers receive their writable state back.
.huge_graphics_state = function() {
  par(no.readonly = TRUE)
}

.huge_restore_graphics_state = function(state) {
  # `new = TRUE` warns when an error occurred before the first valid plot.
  # Restore it separately at the end so that warning cannot mask the error.
  bulk.state = state
  bulk.state$new = NULL
  par(bulk.state)

  restore.group = function(fields, setters) {
    matches = function() {
      current = par(no.readonly = TRUE)
      all(vapply(fields, function(name)
        identical(current[[name]], state[[name]]), logical(1)))
    }
    if(matches()) return(invisible(NULL))
    for(name in setters) {
      do.call(par, stats::setNames(list(state[[name]]), name))
      if(matches()) break
    }
    invisible(NULL)
  }

  multi.panel = prod(state$mfrow) > 1L
  if(multi.panel) par(mfrow = state$mfrow)
  par(pty = state$pty, cex = state$cex, mex = state$mex)
  restore.group(c("oma", "omi", "omd"), c("omi", "omd", "oma"))
  if(multi.panel) par(mfg = state$mfg)
  if(!multi.panel) restore.group(c("fig", "fin"), c("fig", "fin"))
  restore.group(c("mar", "mai"), c("mai", "mar"))
  restore.group(c("pin", "plt"), c("pin", "plt"))
  par(usr = state$usr, xaxp = state$xaxp, yaxp = state$yaxp)
  suppressWarnings(par(new = state$new))
  invisible(NULL)
}

# Unit-variance centered columns with bounded intermediate arithmetic.
# Reference shifting before max-absolute scaling preserves representable
# differences at extreme finite magnitudes.
.huge_standardize = function(x) {
  n = nrow(x)
  d = ncol(x)
  xs = matrix(0, n, d, dimnames = dimnames(x))

  for(j in seq_len(d)) {
    source = x[,j] * 1.0
    column.scale = max(abs(source))
    if(!is.finite(column.scale) || column.scale == 0)
      stop("Data must have finite, positive sample standard deviations.")

    reference = source[[1]]
    direct = source == 0 | reference == 0 |
      ((source < 0) == (reference < 0))
    delta = source / column.scale - reference / column.scale
    delta[direct] = (source[direct] - reference) / column.scale
    delta = delta - mean(delta)

    centered.scale = max(abs(delta))
    if(!is.finite(centered.scale) || centered.scale == 0)
      stop("Data must have finite, positive sample standard deviations.")
    delta = delta / centered.scale
    column.sd = sd(delta)
    if(!is.finite(column.sd) || column.sd <= 0)
      stop("Data must have finite, positive sample standard deviations.")
    xs[,j] = delta / column.sd
  }
  xs
}

# Correlation matrix via BLAS crossprod on stable standardized columns.
.huge_fast_cor = function(x) {
  xs = tryCatch(.huge_standardize(x), error = function(e) NULL)
  if(is.null(xs))
    return(cor(x))
  n = nrow(xs)
  S = crossprod(xs) / (n - 1)
  S[S > 1] = 1
  S[S < -1] = -1
  diag(S) = 1
  S
}

.huge_is_covariance_input = function(x) {
  if(nrow(x) != ncol(x)) return(FALSE)

  values = x * 1.0
  if(any(!is.finite(values))) return(FALSE)

  tolerance = 100 * .Machine$double.eps
  dimension = nrow(values)
  if(dimension <= 1) return(TRUE)

  # Compare off-diagonal entries in implied-correlation units.  This keeps
  # covariance routing invariant to a finite uniform rescaling without
  # subtracting large values before they have been normalized.
  diagonal = abs(diag(values))
  diagonal.root = sqrt(diagonal)
  for(column in 2:dimension) {
    rows = seq_len(column - 1L)
    left = values[rows, column]
    right = values[column, rows]
    different = left != right
    if(!any(different)) next

    rows = rows[different]
    left = left[different]
    right = right[different]
    covariance.scale = diagonal.root[rows] * diagonal.root[column]
    equal.diagonal = diagonal[rows] == diagonal[column]
    covariance.scale[equal.diagonal] = diagonal[column]
    reference = pmax(abs(left), abs(right), covariance.scale)
    left = left / reference
    right = right / reference
    scale = pmax(abs(left), abs(right))
    threshold = ifelse(
      scale <= tolerance, tolerance, tolerance * scale
    )
    if(any(abs(left - right) > threshold)) return(FALSE)
  }
  TRUE
}

.huge_validate_correlation_psd = function(correlation) {
  d = nrow(correlation)
  spectral.bound = max(1, max(rowSums(abs(correlation))))
  tolerance = 100 * .Machine$double.eps * max(1, d) * spectral.bound
  shifted = correlation
  diag(shifted) = diag(shifted) + tolerance
  factor = tryCatch(chol(shifted), error = function(e) NULL)
  if(is.null(factor))
    stop("Covariance input must be positive semidefinite.")
  invisible(NULL)
}

.huge_validate_input_type = function(input.type) {
  allowed = c("auto", "data", "covariance")
  if(length(input.type) != 1L || !is.character(input.type) ||
     is.na(input.type) || !(input.type %in% allowed))
    stop('input.type must be exactly one of "auto", "data", or "covariance".')
  input.type
}

.huge_validate_estimation_input = function(
  x, input.type = "auto", prepare.covariance = TRUE, require.psd = TRUE
) {
  if(!is.matrix(x) || !is.numeric(x))
    stop("x must be a numeric matrix.")

  n = nrow(x)
  d = ncol(x)
  if(n < 1 || d < 1)
    stop("x must be a non-empty numeric matrix.")
  if(any(!is.finite(x)))
    stop("x must contain only finite values.")

  input.type = .huge_validate_input_type(input.type)
  if(input.type == "covariance") {
    if(n != d)
      stop('input.type = "covariance" requires a square matrix.')
    if(!.huge_is_covariance_input(x))
      stop("Covariance input must be symmetric within numeric tolerance.")
  }
  cov.input = switch(
    input.type,
    auto = .huge_is_covariance_input(x),
    data = FALSE,
    covariance = TRUE
  )

  covariance = NULL
  correlation = NULL
  if(cov.input && prepare.covariance)
  {
    transpose = t(x)
    covariance = x * 1.0
    different = x != transpose
    covariance[different] = x[different] / 2 + transpose[different] / 2
    diagonal = diag(covariance)
    if(any(!is.finite(diagonal)) || any(diagonal <= 0))
      stop("Covariance input must have a positive finite diagonal.")

    inv.sd = 1 / sqrt(diagonal)
    correlation = diag(1, d)
    dimnames(correlation) = dimnames(covariance)
    if(d > 1) {
      for(column in 2:d) {
        for(row in seq_len(column - 1L)) {
          inv.large = max(inv.sd[[row]], inv.sd[[column]])
          inv.small = min(inv.sd[[row]], inv.sd[[column]])
          value = (covariance[[row, column]] * inv.large) * inv.small
          if(!is.finite(value))
            stop("Covariance input cannot produce a finite correlation matrix.")
          if(abs(value) > 1 + 1e-8)
            stop("Covariance input is not a valid covariance matrix.")
          value = max(-1, min(1, value))
          correlation[[row, column]] = value
          correlation[[column, row]] = value
        }
      }
    }
    if(require.psd)
      .huge_validate_correlation_psd(correlation)
  }
  else if(!cov.input)
  {
    if(n < 2)
      stop("Raw data x must contain at least two observations.")
    constant = vapply(seq_len(d), function(j) {
      all(x[,j] == x[1,j])
    }, logical(1))
    if(any(constant))
      stop("Raw data x contains a constant column.")
  }

  list(
    input = x, covariance = covariance, correlation = correlation,
    n = n, d = d, cov.input = cov.input
  )
}

.huge_preprocess = function(x, verbose = TRUE, input.type = "auto") {
  inp = .huge_validate_estimation_input(x, input.type = input.type)
  x = inp$input
  n = inp$n
  d = inp$d
  cov.input = inp$cov.input
  if(cov.input) {
    if(verbose) cat("The input is identified as the covariance matrix.\n")
    S = inp$correlation
  } else {
    S = .huge_fast_cor(x)
    if(any(!is.finite(S)))
      stop("Raw data x cannot produce a finite correlation matrix.")
  }
  list(x = x, S = S, n = n, d = d, cov.input = cov.input)
}

.huge_validate_positive_integer = function(value, name) {
  if(length(value) != 1 || !is.numeric(value) || !is.finite(value) ||
     value != floor(value) || value < 1)
    stop(sprintf("%s must be a finite positive integer.", name))
  as.integer(value)
}

.huge_validate_ratio = function(value, name = "lambda.min.ratio") {
  if(length(value) != 1 || !is.numeric(value) || !is.finite(value) ||
     value <= 0 || value > 1)
    stop(sprintf("%s must be a finite number in (0, 1].", name))
  as.numeric(value)
}

.huge_validate_lambda = function(lambda, allow.zero = FALSE) {
  if(!is.numeric(lambda))
    stop("lambda must be a numeric scalar or one-dimensional sequence.")
  lambda.dim = dim(lambda)
  if(!is.null(lambda.dim) && length(lambda.dim) > 1L)
    stop("lambda must be a numeric scalar or one-dimensional sequence.")
  if(length(lambda) < 1 || any(!is.finite(lambda)))
    stop("lambda must be a non-empty finite numeric sequence.")
  if(allow.zero)
  {
    if(any(lambda < 0))
      stop("lambda must contain non-negative values for method = \"ct\".")
  }
  else if(any(lambda <= 0))
    stop("lambda must contain strictly positive values.")
  as.numeric(lambda)
}

.huge_validate_sym = function(sym) {
  if(length(sym) != 1 || !is.character(sym) || is.na(sym) ||
     !(sym %in% c("or", "and")))
    stop("sym must be exactly one of \"or\" or \"and\".")
  sym
}

.huge_geometric_lambda = function(lambda.max, lambda.min.ratio, nlambda) {
  if(nlambda == 1L)
    return(lambda.max)

  lambda.min = lambda.max * lambda.min.ratio
  if(lambda.min > 0) {
    lambda = exp(seq(
      log(lambda.max), log(lambda.min), length.out = nlambda
    ))
    lambda[[1L]] = lambda.max
    return(cummin(lambda))
  }

  # Preserve representable interior points even when the requested endpoint
  # underflows, then saturate only the unrepresentable tail.
  fraction = (seq_len(nlambda) - 1) / (nlambda - 1)
  lambda = exp(log(lambda.max) + fraction * log(lambda.min.ratio))
  smallest.positive = .Machine$double.xmin * .Machine$double.eps
  lambda[!is.finite(lambda) | lambda <= 0] = smallest.positive
  lambda[[1L]] = lambda.max
  cummin(lambda)
}

.huge_default_lambda = function(
  S, d, nlambda = NULL, lambda.min.ratio = NULL, lambda = NULL,
  legacy.glasso.covariance = FALSE
) {
  if(!is.null(lambda)) {
    lambda = .huge_validate_lambda(lambda)
    if(length(lambda) > 1L &&
       any(lambda[-1L] > lambda[-length(lambda)]))
      stop(paste(
        "lambda must be non-increasing for method = \"mb\" or \"glasso\";",
        "tied values are allowed."
      ))
    return(list(lambda = lambda, nlambda = length(lambda)))
  }
  if(is.null(nlambda)) nlambda = 10
  else nlambda = .huge_validate_positive_integer(nlambda, "nlambda")
  if(is.null(lambda.min.ratio)) lambda.min.ratio = 0.1
  else lambda.min.ratio = .huge_validate_ratio(lambda.min.ratio)
  if(legacy.glasso.covariance)
    lambda.max = max(abs(S - diag(d)))
  else {
    offdiag = S
    diag(offdiag) = 0
    lambda.max = max(abs(offdiag))
  }
  if(lambda.max == 0) lambda.max = 1e-3
  lambda = .huge_geometric_lambda(
    lambda.max, lambda.min.ratio, nlambda
  )
  list(lambda = lambda, nlambda = nlambda)
}
