test_that("mb returns valid structure", {
  set.seed(20)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "mb", verbose = FALSE)

  expect_s3_class(fit, "huge")
  expect_equal(fit$method, "mb")
  expect_equal(length(fit$path), length(fit$lambda))
  expect_equal(length(fit$sparsity), length(fit$lambda))
})

test_that("mb handles single-variable raw and covariance inputs", {
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
    do.call(huge.mb, c(args, list(verbose = FALSE)))
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
    for (index in seq_len(path_length)) {
      expect_true(inherits(fit$beta[[index]], "sparseMatrix"), info = case)
      expect_true(inherits(fit$path[[index]], "sparseMatrix"), info = case)
      expect_silent(validObject(fit$beta[[index]]))
      expect_silent(validObject(fit$path[[index]]))
      expect_equal(dim(fit$beta[[index]]), c(1L, 1L), info = case)
      expect_equal(dim(fit$path[[index]]), c(1L, 1L), info = case)
      expect_equal(sum(fit$beta[[index]]), 0, info = case)
      expect_equal(sum(fit$path[[index]]), 0, info = case)
    }
    expect_null(fit$icov, info = case)
  }

  expected_auto <- exp(seq(log(1e-3), log(2e-4), length.out = 3))
  # The generated path preserves lambda.max exactly instead of round-tripping
  # it through log()/exp(), which can move the endpoint by one ULP.
  expected_auto[[1L]] <- 1e-3
  expect_identical(fits$raw_auto$lambda, expected_auto)
  expect_identical(fits$raw_auto$lambda, fits$covariance_auto$lambda)
  expect_identical(fits$raw_explicit$lambda, c(0.5, 0.2))
  expect_identical(fits$raw_explicit$lambda, fits$covariance_explicit$lambda)
})

test_that("mb sparsity is non-decreasing", {
  set.seed(21)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "mb", verbose = FALSE)

  expect_true(all(diff(fit$lambda) < 0))
  expect_true(all(diff(fit$sparsity) >= -1e-10))
})

test_that("mb path matrices are symmetric and binary", {
  set.seed(22)
  L <- huge.generator(n = 80, d = 30, graph = "band", verbose = FALSE)
  fit <- huge(L$data, method = "mb", verbose = FALSE)

  for (k in seq_along(fit$path)) {
    p <- as.matrix(fit$path[[k]])
    expect_equal(p, t(p), info = paste("path asymmetric at k =", k))
    expect_true(all(p %in% c(0, 1)))
    expect_true(all(diag(p) == 0))
  }
})

test_that("mb with scr=TRUE works", {
  set.seed(23)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "mb", scr = TRUE, verbose = FALSE)

  expect_s3_class(fit, "huge")
  expect_true(all(diff(fit$sparsity) >= -1e-10))
})

test_that("mb sym='and' gives sparser graphs than sym='or'", {
  set.seed(24)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit_or  <- huge(L$data, method = "mb", sym = "or",  verbose = FALSE)
  fit_and <- huge(L$data, method = "mb", sym = "and", verbose = FALSE)

  # "and" should be at least as sparse as "or" at each lambda
  for (k in seq_along(fit_or$path)) {
    expect_true(sum(fit_and$path[[k]]) <= sum(fit_or$path[[k]]) + 1e-10,
                info = paste("and not sparser than or at k =", k))
  }
})

test_that("mb works across graph types", {
  set.seed(25)
  for (g in c("hub", "band", "cluster", "random")) {
    L <- huge.generator(n = 60, d = 20, graph = g, verbose = FALSE)
    fit <- huge(L$data, method = "mb", verbose = FALSE)
    expect_true(all(diff(fit$sparsity) >= -1e-10),
                info = paste("non-monotone sparsity for graph =", g))
  }
})

test_that("mb screening excludes the response under tied correlations", {
  S <- matrix(1, nrow = 4, ncol = 4)

  fit <- huge.mb(
    S, lambda = 0.5, scr = TRUE, scr.num = 1, verbose = FALSE
  )

  expect_equal(dim(fit$idx_mat), c(1L, 4L))
  for (m in seq_len(ncol(S))) {
    expect_false((m - 1L) %in% fit$idx_mat[, m])
  }
  for (path in fit$path) {
    expect_true(all(diag(as.matrix(path)) == 0))
  }

  refit <- huge.mb(
    S, lambda = 0.5, idx.mat = fit$idx_mat, verbose = FALSE
  )
  expect_equal(refit$beta, fit$beta)
  expect_equal(refit$path, fit$path)
})

test_that("mb rejects unsafe screening index matrices", {
  S <- diag(4)
  valid <- matrix(
    c(1, 2, 0, 2, 0, 1, 0, 1), nrow = 2, ncol = 4
  )

  bad_inputs <- list(
    not_a_matrix = as.vector(valid),
    wrong_columns = valid[, 1:3, drop = FALSE],
    empty = matrix(integer(), nrow = 0, ncol = 4),
    out_of_range = replace(valid, 1, 4),
    internal_sentinel = replace(valid, 1, -1),
    negative = replace(valid, 1, -2),
    duplicate = replace(valid, 2, valid[1, 1]),
    self_index = replace(valid, 1, 0),
    missing = replace(valid, 1, NA_real_),
    infinite = replace(valid, 1, Inf),
    fractional = replace(valid, 1, 1.5)
  )

  for (case in names(bad_inputs)) {
    expect_error(
      huge.mb(
        S, lambda = 0.5, idx.mat = bad_inputs[[case]], verbose = FALSE
      ),
      "idx.mat",
      info = case
    )
  }
})

test_that("mb covariance default lambda uses correlation off-diagonals", {
  S <- matrix(
    c(100, .6, 0,
      .6, 3, -.2,
      0, -.2, .5),
    nrow = 3
  )
  correlation <- cov2cor(S)
  lambda.max <- max(abs(correlation[upper.tri(correlation)]))
  expected <- exp(seq(
    log(lambda.max), log(.25 * lambda.max), length.out = 3
  ))

  fit <- huge.mb(
    S, nlambda = 3, lambda.min.ratio = .25, verbose = FALSE
  )
  expect_equal(fit$lambda, expected, tolerance = 1e-15)
})
