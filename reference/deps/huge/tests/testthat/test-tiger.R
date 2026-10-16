test_that("tiger returns valid structure", {
  set.seed(30)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "tiger", verbose = FALSE)

  expect_s3_class(fit, "huge")
  expect_equal(fit$method, "tiger")
  expect_equal(length(fit$path), length(fit$lambda))
  expect_equal(length(fit$icov), length(fit$lambda))
  expect_equal(length(fit$sparsity), length(fit$lambda))
})

test_that("tiger handles single-variable raw and covariance inputs", {
  raw <- matrix(c(-1, 0, 1), ncol = 1)
  covariance <- matrix(1, nrow = 1, ncol = 1)
  cases <- list(
    raw_auto = list(x = raw, nlambda = 3, lambda.min.ratio = 0.2),
    covariance_auto = list(
      x = covariance, nlambda = 3, lambda.min.ratio = 0.2
    ),
    raw_explicit = list(x = raw, lambda = c(0.5, 0.2)),
    covariance_explicit = list(x = covariance, lambda = c(0.5, 0.2))
  )

  fits <- lapply(cases, function(args) {
    do.call(huge.tiger, c(args, list(verbose = FALSE)))
  })
  for (case in names(fits)) {
    fit <- fits[[case]]
    path_length <- if (grepl("auto$", case)) 3L else 2L
    expect_identical(
      fit$cov.input, grepl("^covariance", case), info = case
    )
    expect_length(fit$lambda, path_length)
    expect_true(all(is.finite(fit$lambda) & fit$lambda > 0), info = case)
    expect_identical(fit$sparsity, rep(0, path_length), info = case)
    expect_equal(dim(fit$df), c(1L, path_length), info = case)
    expect_identical(as.numeric(fit$df), rep(0, path_length), info = case)
    expect_length(fit$beta, path_length)
    expect_length(fit$path, path_length)
    expect_length(fit$icov, path_length)
    for (index in seq_len(path_length)) {
      expect_true(inherits(fit$beta[[index]], "sparseMatrix"), info = case)
      expect_true(inherits(fit$path[[index]], "sparseMatrix"), info = case)
      expect_silent(validObject(fit$beta[[index]]))
      expect_silent(validObject(fit$path[[index]]))
      expect_equal(dim(fit$beta[[index]]), c(1L, 1L), info = case)
      expect_equal(dim(fit$path[[index]]), c(1L, 1L), info = case)
      expect_equal(sum(fit$beta[[index]]), 0, info = case)
      expect_equal(sum(fit$path[[index]]), 0, info = case)
      expect_equal(dim(fit$icov[[index]]), c(1L, 1L), info = case)
      expect_true(all(is.finite(fit$icov[[index]])), info = case)
      expect_equal(
        fit$icov[[index]], matrix(1, 1, 1), tolerance = 1e-15,
        info = case
      )
    }
  }

  expected_auto <- exp(seq(log(1e-3), log(2e-4), length.out = 3))
  expect_identical(fits$raw_auto$lambda, expected_auto)
  expect_identical(fits$raw_auto$lambda, fits$covariance_auto$lambda)
  expect_identical(fits$raw_explicit$lambda, c(0.5, 0.2))
  expect_identical(fits$raw_explicit$lambda, fits$covariance_explicit$lambda)
})

test_that("tiger sparsity is non-decreasing", {
  set.seed(31)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "tiger", verbose = FALSE)

  expect_true(all(diff(fit$lambda) < 0))
  expect_true(all(diff(fit$sparsity) >= -1e-10))
})

test_that("tiger path matrices are symmetric and binary", {
  set.seed(32)
  L <- huge.generator(n = 80, d = 30, graph = "band", verbose = FALSE)
  fit <- huge(L$data, method = "tiger", verbose = FALSE)

  for (k in seq_along(fit$path)) {
    p <- as.matrix(fit$path[[k]])
    expect_equal(p, t(p), info = paste("path asymmetric at k =", k))
    expect_true(all(p %in% c(0, 1)))
    expect_true(all(diag(p) == 0))
  }
})

test_that("tiger works across graph types", {
  set.seed(33)
  for (g in c("hub", "band", "cluster")) {
    L <- huge.generator(n = 60, d = 20, graph = g, verbose = FALSE)
    fit <- huge(L$data, method = "tiger", verbose = FALSE)
    expect_true(all(diff(fit$sparsity) >= -1e-10),
                info = paste("non-monotone sparsity for graph =", g))
  }
})

test_that("tiger raw and covariance inputs use the same correlation problem", {
  set.seed(34)
  x = matrix(rnorm(120 * 12), 120, 12)
  x = sweep(x, 2, seq_len(ncol(x)), "*")
  lambda = c(0.35, 0.2, 0.1)

  raw_fit = huge.tiger(x, lambda = lambda, verbose = FALSE)
  cov_fit = huge.tiger(cov(x), lambda = lambda, verbose = FALSE)

  expect_identical(cov_fit$lambda, lambda)
  for (i in seq_along(lambda)) {
    expect_equal(as.matrix(raw_fit$beta[[i]]),
                 as.matrix(cov_fit$beta[[i]]), tolerance = 1e-8)
    expect_identical(as.matrix(raw_fit$path[[i]]),
                     as.matrix(cov_fit$path[[i]]))
    expect_equal(raw_fit$icov[[i]], cov_fit$icov[[i]], tolerance = 1e-8)
  }
})

test_that("tiger returns the native default lambda path", {
  set.seed(35)
  x = matrix(rnorm(100 * 10), 100, 10)
  ratio = 0.2
  nlambda = 6

  raw_fit = huge.tiger(x, nlambda = nlambda,
                       lambda.min.ratio = ratio, verbose = FALSE)
  cov_fit = huge.tiger(cov(x), nlambda = nlambda,
                       lambda.min.ratio = ratio, verbose = FALSE)
  lambda_max = max(abs(cor(x)[upper.tri(cor(x))]))

  expect_length(raw_fit$lambda, nlambda)
  expect_equal(raw_fit$lambda[[1]], lambda_max, tolerance = 1e-12)
  expect_equal(raw_fit$lambda[[nlambda]], ratio * lambda_max,
               tolerance = 1e-12)
  expect_equal(raw_fit$lambda, cov_fit$lambda, tolerance = 1e-12)
})

test_that("tiger native lambda preserves weak nonzero correlations", {
  rho = 5e-4
  a = c(-1, -1, 1, 1)
  b = c(-1, 1, -1, 1)
  x = cbind(a, rho * a + sqrt(1 - rho^2) * b)
  ratio = 0.2
  expected = exp(seq(log(rho), log(rho * ratio), length.out = 3))

  raw_fit = huge.tiger(
    x, nlambda = 3, lambda.min.ratio = ratio, verbose = FALSE
  )
  covariance_fit = huge.tiger(
    matrix(c(1, rho, rho, 1), 2, 2),
    nlambda = 3, lambda.min.ratio = ratio, verbose = FALSE
  )

  expect_equal(cor(x)[1, 2], rho, tolerance = 1e-15)
  expect_equal(raw_fit$lambda, expected, tolerance = 1e-12)
  expect_equal(covariance_fit$lambda, expected, tolerance = 1e-12)
  expect_equal(raw_fit$lambda, covariance_fit$lambda, tolerance = 1e-12)
})

test_that("tiger native lambda paths stay positive at the subnormal limit", {
  minimum = .Machine$double.xmin * .Machine$double.eps
  covariance = matrix(c(1, minimum, minimum, 1), 2, 2)

  fit = huge.tiger(
    covariance, nlambda = 3, lambda.min.ratio = .1, verbose = FALSE
  )
  native = .Call(
    "_huge_SPMBgraphsqrtFit", covariance, NULL, 3L, 2L, TRUE, .1,
    PACKAGE = "huge"
  )
  expect_identical(fit$lambda, rep(minimum, 3))
  expect_identical(native$lambda, fit$lambda)
  expect_false(native$path_truncated)
  expect_true(all(vapply(fit$icov, function(x) all(is.finite(x)), logical(1))))

  identity = huge.tiger(
    diag(2), nlambda = 3, lambda.min.ratio = minimum,
    verbose = FALSE
  )
  native.identity = .Call(
    "_huge_SPMBgraphsqrtFit", diag(2), NULL, 3L, 2L, TRUE, minimum,
    PACKAGE = "huge"
  )
  expect_equal(identity$lambda[[1]], 1e-3, tolerance = 1e-15)
  expect_identical(identity$lambda[[3]], minimum)
  expect_identical(native.identity$lambda, identity$lambda)
  expect_false(native.identity$path_truncated)
  expect_true(all(is.finite(identity$lambda) & identity$lambda > 0))
  expect_true(all(diff(identity$lambda) <= 0))
  expect_equal(
    log(identity$lambda[[2]]),
    log(1e-3) + log(minimum) / 2,
    tolerance = 1e-12
  )

  singleton = huge.tiger(
    matrix(1, 1, 1), nlambda = 3, lambda.min.ratio = minimum,
    verbose = FALSE
  )
  native.singleton = .Call(
    "_huge_SPMBgraphsqrtFit", matrix(1, 1, 1), NULL, 3L, 1L, TRUE,
    minimum, PACKAGE = "huge"
  )
  expect_identical(singleton$lambda, identity$lambda)
  expect_identical(native.singleton$lambda, singleton$lambda)
  expect_false(native.singleton$path_truncated)
})

test_that("tiger explicit lambda paths allow ties and reject increases", {
  covariance = matrix(c(1, .3, .3, 1), 2, 2)
  automatic = huge.tiger(
    covariance, nlambda = 3, lambda.min.ratio = 1, verbose = FALSE
  )
  replay = huge.tiger(
    covariance, lambda = automatic$lambda, verbose = FALSE
  )
  generic.replay = huge(
    covariance, method = "tiger", lambda = automatic$lambda,
    verbose = FALSE
  )
  native.replay = .Call(
    "_huge_SPMBgraphsqrtFit", covariance, automatic$lambda, 3L, 2L,
    TRUE, .1, PACKAGE = "huge"
  )

  expect_true(all(diff(automatic$lambda) == 0))
  expect_identical(replay$lambda, automatic$lambda)
  expect_identical(generic.replay$lambda, automatic$lambda)
  expect_identical(native.replay$lambda, automatic$lambda)
  for(index in seq_along(automatic$lambda)) {
    expect_identical(
      as.matrix(replay$beta[[index]]),
      as.matrix(automatic$beta[[index]])
    )
    expect_identical(replay$icov[[index]], automatic$icov[[index]])
  }

  expect_error(
    huge.tiger(covariance, lambda = c(.1, .2), verbose = FALSE),
    "non-increasing"
  )
  expect_error(
    huge(
      covariance, method = "tiger", lambda = c(.1, .2),
      verbose = FALSE
    ),
    "non-increasing"
  )
  expect_error(
    .Call(
      "_huge_SPMBgraphsqrtFit", covariance, c(.1, .2), 2L, 2L,
      TRUE, .1, PACKAGE = "huge"
    ),
    "non-increasing"
  )
  expect_error(
    .Call(
      "_huge_SPMBgraphsqrtFit", covariance, c(.1, .2, NaN), 3L, 2L,
      TRUE, .1, PACKAGE = "huge"
    ),
    "positive and finite"
  )

  raw = cbind(c(-1, -1, 1, 1), c(-1, 1, -1, 1))
  expect_no_error(
    .Call(
      "_huge_SPMBgraphsqrt", raw, c(.5, .5), 2L, 2L,
      PACKAGE = "huge"
    )
  )
  expect_error(
    .Call(
      "_huge_SPMBgraphsqrt", raw, c(.1, .2), 2L, 2L,
      PACKAGE = "huge"
    ),
    "non-increasing"
  )
  expect_error(
    .Call(
      "_huge_SPMBgraphsqrt", raw, c(.1, .2, NaN), 3L, 2L,
      PACKAGE = "huge"
    ),
    "positive and finite"
  )
  expect_error(
    .Call(
      "_huge_SPMBgraphsqrt", raw, .5, 2L, 2L,
      PACKAGE = "huge"
    ),
    "length must match"
  )
  expect_error(
    .Call(
      "_huge_SPMBgraphsqrt", raw, c(.5, .5), 2L, 3L,
      PACKAGE = "huge"
    ),
    "dimensions"
  )
})

test_that("tiger raw correlation is stable across finite column scales", {
  x = cbind(c(-1, 0, 1), c(-1, 1, 0))
  ratio = 0.5
  reference = huge.tiger(
    x, nlambda = 3, lambda.min.ratio = ratio, verbose = FALSE
  )
  native.reference = .Call(
    "_huge_SPMBgraphsqrtFit", x, NULL, 3L, 2L, FALSE, ratio,
    PACKAGE = "huge"
  )
  explicit.reference = .Call(
    "_huge_SPMBgraphsqrtFit", x, .25, 1L, 2L, FALSE, ratio,
    PACKAGE = "huge"
  )
  expect_equal(reference$lambda[[1]], 0.5, tolerance = 1e-15)

  minimum = .Machine$double.xmin * .Machine$double.eps
  inputs = list(
    x * 1e308,
    x * 1e-200,
    x * minimum,
    sweep(x, 2, c(1e308, 1e-200), "*")
  )
  for (current in inputs) {
    expect_no_warning({
      fit = huge.tiger(
        current, nlambda = 3, lambda.min.ratio = ratio, verbose = FALSE
      )
      native = .Call(
        "_huge_SPMBgraphsqrtFit", current, NULL, 3L, 2L, FALSE, ratio,
        PACKAGE = "huge"
      )
      explicit = .Call(
        "_huge_SPMBgraphsqrtFit", current, .25, 1L, 2L, FALSE, ratio,
        PACKAGE = "huge"
      )
    })
    expect_equal(fit$lambda, reference$lambda, tolerance = 1e-14)
    for (field in c("col_cnz", "row_idx")) {
      expect_identical(native[[field]], native.reference[[field]])
      expect_identical(explicit[[field]], explicit.reference[[field]])
    }
    for (field in c("x", "icov", "lambda")) {
      expect_equal(native[[field]], native.reference[[field]], tolerance = 1e-12)
      expect_equal(
        explicit[[field]], explicit.reference[[field]], tolerance = 1e-12
      )
    }
    for (field in c("hit_max_iter", "path_truncated")) {
      expect_identical(native[[field]], native.reference[[field]])
      expect_identical(explicit[[field]], explicit.reference[[field]])
    }
    for (i in seq_along(reference$lambda)) {
      expect_equal(
        as.matrix(fit$beta[[i]]), as.matrix(reference$beta[[i]]),
        tolerance = 1e-12
      )
      expect_identical(
        as.matrix(fit$path[[i]]), as.matrix(reference$path[[i]])
      )
      expect_equal(fit$icov[[i]], reference$icov[[i]], tolerance = 1e-12)
    }
  }
})

test_that("tiger correlation solver satisfies square-root lasso KKT conditions", {
  set.seed(9)
  n = 80
  d = 15
  lambda = 0.2
  x = matrix(rnorm(n * d), n, d)
  fit = huge.tiger(x, lambda = lambda, verbose = FALSE)
  beta = as.matrix(fit$beta[[1]])
  correlation = cor(x)

  kkt_error = vapply(seq_len(d), function(response) {
    coefficient = beta[, response]
    residual_variance = correlation[response, response] -
      2 * sum(correlation[, response] * coefficient) +
      sum(coefficient * (correlation %*% coefficient))
    tau = sqrt(max(residual_variance, .Machine$double.eps))
    score = as.numeric((correlation[, response] -
                          correlation %*% coefficient) / tau)
    score[response] = 0

    active = abs(coefficient) > 1e-8
    active[response] = FALSE
    active_error = if(any(active))
      max(abs(score[active] - lambda * sign(coefficient[active]))) else 0
    inactive = !active
    inactive[response] = FALSE
    inactive_error = if(any(inactive))
      max(pmax(abs(score[inactive]) - lambda, 0)) else 0
    max(active_error, inactive_error)
  }, numeric(1))

  expect_gt(sum(abs(beta) > 1e-8), 0)
  expect_lte(max(kkt_error), 1e-6)
})

test_that("tiger matches the two-variable analytic solution", {
  rho = 0.3
  lambda = 0.2
  correlation = matrix(c(1, rho, rho, 1), 2, 2)
  expected = rho - lambda * sqrt((1 - rho^2) / (1 - lambda^2))

  fit = huge.tiger(correlation, lambda = lambda, verbose = FALSE)
  coefficient = as.matrix(fit$beta[[1]])

  expect_equal(coefficient[1, 2], expected, tolerance = 1e-7)
  expect_equal(coefficient[2, 1], expected, tolerance = 1e-7)
})

test_that("tiger handles identity, singular, and rescaled covariance inputs", {
  identity_fit = huge.tiger(diag(4), nlambda = 3, verbose = FALSE)
  expect_equal(identity_fit$lambda, c(1e-3, sqrt(1e-7), 1e-4),
               tolerance = 1e-15)
  expect_identical(identity_fit$sparsity, rep(0, 3))
  for (icov in identity_fit$icov)
    expect_true(all(is.finite(icov)))

  set.seed(36)
  x = matrix(rnorm(100 * 6), 100, 6)
  covariance = cov(x)
  scale = diag(seq_len(ncol(x)))
  fit = huge.tiger(covariance, lambda = c(0.4, 0.2), verbose = FALSE)
  rescaled = huge.tiger(scale %*% covariance %*% scale,
                        lambda = c(0.4, 0.2), verbose = FALSE)

  for (i in 1:2) {
    expect_equal(as.matrix(fit$beta[[i]]), as.matrix(rescaled$beta[[i]]),
                 tolerance = 1e-8)
    expect_equal(fit$icov[[i]], rescaled$icov[[i]], tolerance = 1e-8)
  }

  correlation = matrix(c(1, .6, .6, 1), 2, 2)
  large = .Machine$double.xmax * correlation
  large_fit = huge.tiger(large, lambda = .25, verbose = FALSE)
  correlation_fit = huge.tiger(
    correlation, lambda = .25, verbose = FALSE
  )
  expect_equal(
    as.matrix(large_fit$beta[[1]]),
    as.matrix(correlation_fit$beta[[1]]), tolerance = 1e-12
  )
  expect_equal(
    large_fit$icov[[1]], correlation_fit$icov[[1]], tolerance = 1e-12
  )

  large_auto_fit = huge.tiger(
    large, nlambda = 3, lambda.min.ratio = .5, verbose = FALSE
  )
  correlation_auto_fit = huge.tiger(
    correlation, nlambda = 3, lambda.min.ratio = .5, verbose = FALSE
  )
  expect_equal(
    large_auto_fit$lambda, correlation_auto_fit$lambda, tolerance = 1e-14
  )
  expect_equal(large_auto_fit$lambda[[1]], .6, tolerance = 1e-14)

  nearly.symmetric = large
  nearly.symmetric[2, 1] = nearly.symmetric[1, 2] *
    (1 - .Machine$double.eps)
  expect_false(identical(nearly.symmetric[1, 2], nearly.symmetric[2, 1]))
  projected = nearly.symmetric / 2 + t(nearly.symmetric) / 2
  near_fit = huge.tiger(
    nearly.symmetric, nlambda = 3, lambda.min.ratio = .5, verbose = FALSE
  )
  projected_fit = huge.tiger(
    projected, nlambda = 3, lambda.min.ratio = .5, verbose = FALSE
  )
  expect_equal(near_fit$lambda, projected_fit$lambda, tolerance = 1e-14)
  for (i in seq_along(near_fit$lambda)) {
    expect_equal(
      as.matrix(near_fit$beta[[i]]),
      as.matrix(projected_fit$beta[[i]]), tolerance = 1e-12
    )
    expect_equal(
      near_fit$icov[[i]], projected_fit$icov[[i]], tolerance = 1e-12
    )
  }

  x[, 6] = x[, 1]
  expect_warning(
    singular_fit <- huge.tiger(cov(x), nlambda = 2,
                               lambda.min.ratio = 0.9, verbose = FALSE),
    "certified prefix"
  )
  expect_length(singular_fit$lambda, 1)
  expect_true(all(is.finite(as.matrix(singular_fit$beta[[1]]))))
  expect_true(all(is.finite(singular_fit$icov[[1]])))
})

test_that("tiger preserves weak correlations across extreme variances", {
  rho = 1e-200
  high = .Machine$double.xmax
  low = 1e-300
  off.diagonal = rho * sqrt(high) * sqrt(low)
  covariance = matrix(
    c(high, off.diagonal, off.diagonal, low), nrow = 2
  )
  correlation = matrix(c(1, rho, rho, 1), nrow = 2)

  reference = huge.tiger(
    correlation, nlambda = 3, lambda.min.ratio = .5, verbose = FALSE
  )

  for (current in list(covariance, covariance[2:1, 2:1])) {
    fit = huge.tiger(
      current, nlambda = 3, lambda.min.ratio = .5, verbose = FALSE
    )
    native = .Call(
      "_huge_SPMBgraphsqrtFit", current, NULL, 3L, 2L, TRUE, .5,
      PACKAGE = "huge"
    )

    expect_equal(
      fit$lambda / rho, reference$lambda / rho, tolerance = 1e-12
    )
    expect_equal(fit$lambda[[1]] / rho, 1, tolerance = 1e-12)
    expect_equal(
      native$lambda / rho, reference$lambda / rho, tolerance = 1e-12
    )
  }
})

test_that("tiger stops generated paths before an interpolation degeneracy", {
  signal = seq_len(30)
  duplicated_data = cbind(signal, signal)

  expect_warning(
    fit <- huge.tiger(duplicated_data, nlambda = 2,
                      lambda.min.ratio = 0.99, verbose = FALSE),
    "certified prefix"
  )
  expect_equal(fit$lambda, 1, tolerance = 1e-12)
  returned_nlambda = length(fit$lambda)
  expect_length(fit$beta, returned_nlambda)
  expect_length(fit$path, returned_nlambda)
  expect_length(fit$icov, returned_nlambda)
  expect_length(fit$sparsity, returned_nlambda)
  expect_equal(ncol(fit$df), returned_nlambda)
  expect_true(all(is.finite(fit$icov[[1]])))
  expect_equal(fit$icov[[1]], diag(2), tolerance = 1e-12)

  expect_warning(
    covariance_fit <- huge.tiger(cov(duplicated_data), nlambda = 2,
                                 lambda.min.ratio = 0.99, verbose = FALSE),
    "certified prefix"
  )
  expect_identical(covariance_fit$lambda, fit$lambda)
  expect_equal(as.matrix(covariance_fit$beta[[1]]),
               as.matrix(fit$beta[[1]]), tolerance = 1e-12)
  expect_equal(covariance_fit$icov[[1]], fit$icov[[1]], tolerance = 1e-12)

  expect_error(
    huge.tiger(duplicated_data, lambda = 0.99, verbose = FALSE),
    "could not certify a supplied lambda"
  )
})

test_that("tiger rejects invalid native preprocessing inputs", {
  constant = cbind(matrix(rnorm(60 * 3), 60, 3), 1)
  expect_error(huge.tiger(constant, verbose = FALSE), "constant column")

  invalid_covariance = diag(c(1, 0, 1))
  expect_error(huge.tiger(invalid_covariance, verbose = FALSE),
               "positive finite diagonal")

  nonfinite = matrix(rnorm(60 * 4), 60, 4)
  nonfinite[1, 1] = NA_real_
  expect_error(huge.tiger(nonfinite, verbose = FALSE), "finite values")
})

test_that("tiger native covariance input validates symmetry", {
  asymmetric = matrix(c(1, 1, -1, 1), 2, 2, byrow = TRUE)
  expect_error(
    .Call(
      "_huge_SPMBgraphsqrtFit", asymmetric, .25, 1L, 2L, TRUE, .1,
      PACKAGE = "huge"
    ),
    "symmetric"
  )

  near = matrix(c(1, .5, .5 + 1e-15, 1), 2, 2, byrow = TRUE)
  expect_no_error(.Call(
    "_huge_SPMBgraphsqrtFit", near, .25, 1L, 2L, TRUE, .1,
    PACKAGE = "huge"
  ))
})
