test_that("glasso returns valid structure", {
  set.seed(10)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)

  expect_s3_class(fit, "huge")
  expect_equal(fit$method, "glasso")
  expect_equal(length(fit$path), length(fit$lambda))
  expect_equal(length(fit$icov), length(fit$lambda))
  expect_equal(length(fit$sparsity), length(fit$lambda))
  expect_equal(length(fit$loglik), length(fit$lambda))
  expect_equal(length(fit$df), length(fit$lambda))
})

test_that("glasso lambda path is decreasing and sparsity is non-decreasing", {
  set.seed(11)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)

  expect_true(all(diff(fit$lambda) < 0))
  expect_true(all(diff(fit$sparsity) >= -1e-10))
  expect_true(all(diff(fit$df) >= 0))
})

test_that("glasso loglik values are finite", {
  set.seed(12)
  L <- huge.generator(n = 80, d = 30, graph = "band", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)
  expect_true(all(is.finite(fit$loglik)))
})

test_that("glasso loglik matches the returned precision matrices", {
  cases <- list(
    list(S = diag(c(2, 3)), lambda = 0.5),
    list(
      S = matrix(c(2, 0.6, 0,
                   0.6, 3, 0,
                   0, 0, 4), 3, 3),
      lambda = c(0.4, 0.2)
    )
  )

  for (case in cases) {
    fit <- huge(case$S, method = "glasso", lambda = case$lambda,
                verbose = FALSE)
    for (k in seq_along(fit$icov)) {
      precision <- fit$icov[[k]]
      direct <- as.numeric(determinant(precision, logarithm = TRUE)$modulus) -
        sum(case$S * precision)
      expect_equal(fit$loglik[k], direct, tolerance = 1e-10,
                   info = paste("lambda index", k))
    }
  }
})

test_that("glasso path matrices are symmetric", {
  set.seed(13)
  L <- huge.generator(n = 80, d = 30, graph = "cluster", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)

  for (k in seq_along(fit$path)) {
    p <- as.matrix(fit$path[[k]])
    expect_equal(p, t(p), info = paste("path asymmetric at k =", k))
    expect_true(all(diag(p) == 0))
  }
})

test_that("glasso icov matrices are symmetric", {
  set.seed(14)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)

  for (k in seq_along(fit$icov)) {
    ic <- fit$icov[[k]]
    expect_equal(ic, t(ic), tolerance = 1e-4,
                 info = paste("icov asymmetric at k =", k))
  }
})

test_that("glasso precision symmetrization avoids finite overflow", {
  S <- 1e-307 * matrix(c(1, .99, .99, 1), nrow = 2)
  fit <- huge.glasso(
    S, lambda = 1e-309, cov.output = TRUE, verbose = FALSE
  )

  precision <- fit$icov[[1]]
  expect_true(all(is.finite(precision)))
  expect_identical(precision, t(precision))
  expect_true(all(is.finite(fit$cov[[1]])))
  expect_true(all(is.finite(fit$loglik)))
  direct <- as.numeric(determinant(precision, logarithm = TRUE)$modulus) -
    sum(S * precision)
  expect_equal(fit$loglik[1], direct, tolerance = 1e-10)
})

test_that("glasso log-likelihood has no absolute pivot scale", {
  S <- 1e20 * matrix(c(1, .5, .5, 1), nrow = 2)
  fit <- huge.glasso(
    S, lambda = 1e19, cov.output = TRUE, verbose = FALSE
  )

  precision <- fit$icov[[1]]
  expect_true(all(is.finite(precision)))
  expect_true(all(is.finite(fit$cov[[1]])))
  direct <- as.numeric(determinant(precision, logarithm = TRUE)$modulus) -
    sum(S * precision)
  expect_true(is.finite(direct))
  expect_equal(fit$loglik[1], direct, tolerance = 1e-10)
})

test_that("glasso rejects non-finite native results", {
  for (scale in c(1e-320, 1e308)) {
    for (cov.output in c(FALSE, TRUE)) {
      expect_error(
        huge.glasso(
          diag(scale, 2), lambda = scale,
          cov.output = cov.output, verbose = FALSE
        ),
        "non-finite"
      )
    }
  }
})

test_that("glasso with scr=TRUE produces valid results", {
  set.seed(15)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", scr = TRUE, verbose = FALSE)

  expect_true(all(diff(fit$sparsity) >= -1e-10))
  expect_true(all(is.finite(fit$loglik)))
})

test_that("glasso cov.output works", {
  set.seed(16)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", cov.output = TRUE, verbose = FALSE)

  expect_true(!is.null(fit$cov))
  expect_equal(length(fit$cov), length(fit$lambda))
  for (k in seq_along(fit$cov)) {
    co <- fit$cov[[k]]
    expect_equal(dim(co), c(30, 30))
    expect_equal(co, t(co), tolerance = 1e-10)
  }
})

test_that("glasso accepts covariance matrix input", {
  set.seed(17)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  S <- cor(L$data)
  fit <- huge(S, method = "glasso", verbose = FALSE)

  expect_s3_class(fit, "huge")
  expect_true(fit$cov.input)
  expect_true(all(diff(fit$sparsity) >= -1e-10))
})

test_that("glasso default lambda ignores covariance diagonal scale", {
  expected.identity <- c(1e-3, sqrt(1e-7), 1e-4)
  diagonal <- diag(c(100, 3, .5, 8))
  diagonal.fit <- huge.glasso(
    diagonal, nlambda = 3, lambda.min.ratio = .1, verbose = FALSE,
    input.type = "covariance"
  )

  expect_equal(diagonal.fit$lambda, expected.identity, tolerance = 1e-15)
  expect_true(all(vapply(
    diagonal.fit$path, function(path) sum(path) == 0, logical(1)
  )))
  expect_true(all(vapply(
    diagonal.fit$icov, function(icov) all(is.finite(icov)), logical(1)
  )))

  S <- matrix(
    c(100, .6, 0,
      .6, 3, -.2,
      0, -.2, .5),
    nrow = 3
  )
  covariance.fit <- huge.glasso(
    S, nlambda = 3, lambda.min.ratio = .25, verbose = FALSE,
    input.type = "covariance"
  )
  expect_equal(covariance.fit$lambda, c(.6, .3, .15),
               tolerance = 1e-15)

  weak <- diag(3)
  weak[1, 2] <- weak[2, 1] <- 1e-6
  weak.fit <- huge.glasso(
    weak, nlambda = 2, lambda.min.ratio = .1, verbose = FALSE,
    input.type = "covariance"
  )
  expect_equal(weak.fit$lambda, c(1e-6, 1e-7), tolerance = 1e-20)

  explicit <- c(.55, .2)
  expect_identical(
    huge.glasso(S, lambda = explicit, verbose = FALSE)$lambda,
    explicit
  )
})

test_that("glasso auto covariance keeps its historical default lambda scale", {
  covariance <- matrix(c(.001, .002, .002, .02), nrow = 2)
  legacy.max <- max(abs(covariance - diag(2)))

  automatic <- huge.glasso(
    covariance, nlambda = 2, lambda.min.ratio = .1, verbose = FALSE
  )
  named.automatic <- huge.glasso(
    covariance, nlambda = 2, lambda.min.ratio = .1, verbose = FALSE,
    input.type = c(route = "auto")
  )
  attributed.automatic <- huge.glasso(
    covariance, nlambda = 2, lambda.min.ratio = .1, verbose = FALSE,
    input.type = structure("auto", class = "AsIs")
  )
  explicit <- huge.glasso(
    covariance, nlambda = 2, lambda.min.ratio = .1, verbose = FALSE,
    input.type = "covariance"
  )

  expect_equal(automatic$lambda, c(legacy.max, .1 * legacy.max),
               tolerance = 1e-15)
  expect_identical(named.automatic$lambda, automatic$lambda)
  expect_identical(attributed.automatic$lambda, automatic$lambda)
  expect_equal(explicit$lambda, c(.002, .0002), tolerance = 1e-15)
})

test_that("glasso raw-data default lambda remains correlation based", {
  set.seed(171)
  x <- matrix(rnorm(70 * 8), nrow = 70)
  correlation <- cor(x)
  lambda.max <- max(abs(correlation[upper.tri(correlation)]))
  expected <- exp(seq(
    log(lambda.max), log(.2 * lambda.max), length.out = 4
  ))

  fit <- huge.glasso(
    x, nlambda = 4, lambda.min.ratio = .2, verbose = FALSE
  )
  expect_equal(fit$lambda, expected, tolerance = 1e-12)
})

test_that("glasso refines symmetric precision and rejects non-SPD paths", {
  rho <- 0.9999
  covariance <- matrix(rho, 3, 3)
  diag(covariance) <- 1

  for (cov.output in c(FALSE, TRUE)) {
    expect_error(
      huge.glasso(
        covariance, lambda = c(0.1, 0.0002),
        cov.output = cov.output, verbose = FALSE
      ),
      "not positive definite"
    )
  }

  fit <- huge.glasso(
    covariance, lambda = 0.01, cov.output = TRUE, verbose = FALSE
  )
  fit.no.cov <- huge.glasso(
    covariance, lambda = 0.01, cov.output = FALSE, verbose = FALSE
  )
  precision <- fit$icov[[1]]
  expect_identical(precision, t(precision))
  expect_equal(fit.no.cov$icov[[1]], precision, tolerance = 1e-12)
  expect_gt(determinant(precision, logarithm = TRUE)$sign, 0)
  expect_true(
    max(rowSums(abs(fit$cov[[1]] %*% precision - diag(3)))) <= 1e-2
  )
})

test_that("glasso returns a coherent pair for ill-conditioned Toeplitz input", {
  d <- 20
  rho <- 0.99
  covariance <- outer(
    seq_len(d), seq_len(d), function(i, j) rho^abs(i - j)
  )
  fit <- huge.glasso(
    covariance, lambda = 0.001, cov.output = TRUE, verbose = FALSE
  )

  precision <- fit$icov[[1]]
  estimated.covariance <- fit$cov[[1]]
  expect_identical(precision, t(precision))
  expect_identical(estimated.covariance, t(estimated.covariance))
  expect_gt(determinant(precision, logarithm = TRUE)$sign, 0)
  expect_gt(determinant(estimated.covariance, logarithm = TRUE)$sign, 0)
  expect_lte(
    max(rowSums(abs(estimated.covariance %*% precision - diag(d)))),
    1e-2
  )
})

test_that("glasso refines a finite iteration-limit candidate", {
  covariance <- matrix(c(
    1.0000000000000002, .18009786022575919, -.10400558762095684,
    .76774091516199638, .18742182119843603,
    .18009786022575919, 1, .83834121674080364, .48252423066148048,
    -.066223216563695037,
    -.10400558762095684, .83834121674080364, 1, .097388083759665525,
    -.56897959252915542,
    .76774091516199638, .48252423066148048, .097388083759665525,
    1.0000000000000002, .44589188604391788,
    .18742182119843603, -.066223216563695037, -.56897959252915542,
    .44589188604391788, .99999999999999989
  ), nrow = 5, byrow = TRUE)

  expect_warning(
    fit <- huge.glasso(
      covariance, lambda = .00083834121674080362,
      cov.output = TRUE, verbose = FALSE,
      input.type = "covariance"
    ),
    "iteration limit"
  )
  expect_gt(determinant(fit$icov[[1]], logarithm = TRUE)$sign, 0)
  expect_lte(
    max(rowSums(abs(fit$cov[[1]] %*% fit$icov[[1]] - diag(5)))),
    1e-2
  )
})

test_that("glasso works across graph types", {
  set.seed(18)
  for (g in c("hub", "band", "cluster", "random")) {
    L <- huge.generator(n = 60, d = 20, graph = g, verbose = FALSE)
    fit <- huge(L$data, method = "glasso", verbose = FALSE)
    expect_true(all(diff(fit$sparsity) >= -1e-10),
                info = paste("non-monotone sparsity for graph =", g))
  }
})
