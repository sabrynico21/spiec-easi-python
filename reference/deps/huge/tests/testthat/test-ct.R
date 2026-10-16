test_that("ct returns valid structure", {
  set.seed(40)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "ct", verbose = FALSE)

  expect_s3_class(fit, "huge")
  expect_equal(fit$method, "ct")
  expect_equal(length(fit$path), length(fit$lambda))
  expect_equal(length(fit$sparsity), length(fit$lambda))
})

test_that("ct sparsity is non-decreasing", {
  set.seed(41)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "ct", verbose = FALSE)
  expect_true(all(diff(fit$sparsity) >= -1e-10))
})

test_that("ct path matrices are symmetric and binary", {
  set.seed(42)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "ct", verbose = FALSE)

  for (k in seq_along(fit$path)) {
    p <- as.matrix(fit$path[[k]])
    expect_equal(p, t(p))
    expect_true(all(p %in% c(0, 1)))
    expect_true(all(diag(p) == 0))
  }
})

.ct_covariance <- function(edge.weights) {
  S <- diag(4)
  S[upper.tri(S)] <- edge.weights
  S[lower.tri(S)] <- t(S)[lower.tri(S)]
  S
}

.expect_ct_refittable <- function(S, fit) {
  refit <- huge.ct(S, lambda = fit$lambda, verbose = FALSE)
  d <- ncol(S)

  expect_equal(length(fit$path), length(fit$lambda))
  expect_equal(length(fit$sparsity), length(fit$lambda))
  for (k in seq_along(fit$path)) {
    path <- as.matrix(fit$path[[k]])
    expected.sparsity <- if (d <= 1) 0 else sum(path) / d / (d - 1)

    expect_equal(path, t(path))
    expect_true(all(diag(path) == 0))
    expect_true(all(path %in% c(0, 1)))
    expect_equal(fit$sparsity[k], expected.sparsity)
    expect_equal(path, as.matrix(refit$path[[k]]))
    expect_equal(fit$sparsity[k], refit$sparsity[k])
    if (k > 1)
      expect_true(all(as.matrix(fit$path[[k - 1]]) <= path))
  }
}

test_that("ct default path handles undirected edges and ties", {
  cases <- list(
    unique = list(
      S = .ct_covariance(c(.30, .25, .20, .15, .10, .05)),
      nlambda = 3, ratio = .5,
      lambda = c(.25, .20, .15), edges = c(1, 2, 3)
    ),
    boundary_ties = list(
      S = .ct_covariance(c(.30, .20, .20, .10, .08, .05)),
      nlambda = 3, ratio = .5,
      lambda = c(.20, .20, .10), edges = c(1, 1, 3)
    ),
    all_tied = list(
      S = .ct_covariance(rep(.10, 6)),
      nlambda = 2, ratio = 1,
      lambda = c(.10, 0), edges = c(0, 6)
    ),
    identity = list(
      S = diag(4),
      nlambda = 3, ratio = .5,
      lambda = rep(0, 3), edges = c(0, 0, 0)
    )
  )

  for (case in cases) {
    fit <- huge.ct(
      case$S, nlambda = case$nlambda,
      lambda.min.ratio = case$ratio, verbose = FALSE
    )
    edge.count <- vapply(
      fit$path,
      function(path) sum(as.matrix(path)[upper.tri(path)]),
      numeric(1)
    )

    expect_equal(fit$lambda, case$lambda)
    expect_true(all(fit$lambda >= 0))
    expect_equal(edge.count, case$edges)
    .expect_ct_refittable(case$S, fit)
  }
})

test_that("ct handles one- and two-variable default paths", {
  singleton <- huge.ct(
    matrix(1, nrow = 1, ncol = 1),
    nlambda = 3, lambda.min.ratio = .5, verbose = FALSE
  )
  expect_equal(singleton$lambda, rep(0, 3))
  expect_equal(singleton$sparsity, rep(0, 3))
  .expect_ct_refittable(matrix(1, nrow = 1, ncol = 1), singleton)

  pair.S <- matrix(c(1, .3, .3, 1), nrow = 2)
  pair <- huge.ct(
    pair.S, nlambda = 3, lambda.min.ratio = .5, verbose = FALSE
  )
  expect_equal(pair$lambda, rep(0, 3))
  expect_equal(pair$sparsity, rep(1, 3))
  .expect_ct_refittable(pair.S, pair)
})

test_that("ct default full-graph endpoint retains subnormal correlations", {
  smallest.positive = .Machine$double.xmin * .Machine$double.eps
  S = matrix(c(1, smallest.positive, smallest.positive, 1), 2, 2)

  fit = huge.ct(
    S, nlambda = 1, lambda.min.ratio = 1, verbose = FALSE
  )

  expect_identical(fit$lambda, 0)
  expect_identical(as.matrix(fit$path[[1]]), matrix(c(0, 1, 1, 0), 2, 2))
  .expect_ct_refittable(S, fit)
})

test_that("ct default intermediate thresholds retain subnormal ordering", {
  smallest.positive = .Machine$double.xmin * .Machine$double.eps
  S = diag(3)
  S[upper.tri(S)] = c(3, 2, 1) * smallest.positive
  S[lower.tri(S)] = t(S)[lower.tri(S)]

  fit = huge.ct(
    S, nlambda = 2, lambda.min.ratio = 2 / 3, verbose = FALSE
  )
  edge.count = vapply(
    fit$path,
    function(path) sum(as.matrix(path)[upper.tri(path)]),
    numeric(1)
  )

  expect_identical(
    fit$lambda, c(2, 1) * smallest.positive
  )
  expect_identical(edge.count, c(1, 2))
  .expect_ct_refittable(S, fit)
})

test_that("covariance normalization preserves symmetry across extreme scales", {
  rho = 1e-200
  variances = c(.Machine$double.xmax, 1e-300)
  covariance.value = (rho * sqrt(variances[[1]])) * sqrt(variances[[2]])
  covariance = matrix(
    c(variances[[1]], covariance.value,
      covariance.value, variances[[2]]),
    2, 2
  )

  for(S in list(covariance, covariance[2:1, 2:1])) {
    validated = huge:::.huge_validate_estimation_input(S)
    expect_identical(validated$correlation, t(validated$correlation))
    expect_equal(validated$correlation[[1, 2]], rho, tolerance = 1e-14)

    fit = huge.ct(S, lambda = rho / 2, verbose = FALSE)
    expect_identical(
      as.matrix(fit$path[[1]]), matrix(c(0, 1, 1, 0), 2, 2)
    )
  }
})

test_that("ct uses a literal strict threshold", {
  S <- diag(3)
  S[1, 2] <- S[2, 1] <- .25 + 1e-15
  S[1, 3] <- S[3, 1] <- .25

  fit <- huge.ct(S, lambda = .25, verbose = FALSE)
  path <- as.matrix(fit$path[[1]])

  expect_equal(sum(path[upper.tri(path)]), 1)
  expect_equal(path[1, 2], 1)
  expect_equal(path[1, 3], 0)
})
