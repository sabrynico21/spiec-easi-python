test_that("huge.inference Gaussian produces valid p-values", {
  set.seed(70)
  L <- huge.generator(n = 100, d = 15, graph = "hub", g = 3, verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)
  T_hat <- tail(fit$icov, 1)[[1]]
  inf <- huge.inference(L$data, T_hat, L$theta)

  expect_true(!is.null(inf$p))
  expect_equal(dim(inf$p), c(15, 15))
  # p-values should be in [0, 1]
  expect_true(all(inf$p >= 0 & inf$p <= 1))
  # error rate should be in [0, 1]
  expect_true(inf$error >= 0 && inf$error <= 1)
})

test_that("huge.inference preserves the Gaussian one-variable boundary", {
  x <- matrix(0:3, ncol = 1)
  inf <- huge.inference(x, matrix(1, 1, 1), matrix(0, 1, 1))

  expect_equal(dim(inf$p), c(1L, 1L))
  expect_equal(inf$p[1, 1], 0.15729920705028505, tolerance = 1e-15)
  expect_identical(inf$error, 0)
})

test_that("huge.inference rejects degenerate data before arithmetic", {
  one.row <- matrix(c(1, 2, 3), nrow = 1)
  for (type in c("Gaussian", "Nonparanormal")) {
    expect_error(
      huge.inference(one.row, diag(3), matrix(0, 3, 3), type = type),
      "at least two observations",
      info = type
    )
  }

  one.variable <- matrix(0:3, ncol = 1)
  for (method in c("score", "wald")) {
    expect_error(
      huge.inference(
        one.variable, matrix(1, 1, 1), matrix(0, 1, 1),
        type = "Nonparanormal", method = method
      ),
      "at least two variables",
      info = method
    )
  }

  constant <- cbind(0:3, rep(1, 4))
  for (type in c("Gaussian", "Nonparanormal")) {
    expect_error(
      huge.inference(constant, diag(2), matrix(0, 2, 2), type = type),
      "constant column",
      info = type
    )
  }
})

test_that("huge.inference validates active options and matrix inputs", {
  set.seed(71)
  x <- matrix(rnorm(16), 8, 2)
  adj <- matrix(0, 2, 2)

  for (type in list("bad", c("Gaussian", "Nonparanormal"), 1, NA_character_)) {
    expect_error(huge.inference(x, diag(2), adj, type = type), "type")
  }
  expect_error(
    huge.inference(x, diag(2), adj, type = "Nonparanormal", method = "bad"),
    "method"
  )
  expect_no_error(huge.inference(x, diag(2), adj, method = "ignored"))

  for (alpha in list(0, -0.1, 1.1, NA_real_, Inf, TRUE, c(.05, .1))) {
    expect_error(huge.inference(x, diag(2), adj, alpha = alpha), "alpha")
  }
  expect_no_error(huge.inference(x, diag(2), adj, alpha = 1))

  expect_error(huge.inference(1:4, diag(2), adj), "data.*numeric matrix")
  bad.data <- x
  bad.data[1, 1] <- NA_real_
  expect_error(huge.inference(bad.data, diag(2), adj), "data.*finite")

  expect_error(huge.inference(x, diag(3), adj), "T.*2.*2")
  bad.t <- diag(2)
  bad.t[1, 2] <- Inf
  expect_error(huge.inference(x, bad.t, adj), "T.*finite")
  expect_error(
    huge.inference(x, diag(c(0, 1)), adj), "positive diagonal"
  )
  expect_error(
    huge.inference(x, diag(c(-1, 1)), adj), "positive diagonal"
  )

  expect_error(huge.inference(x, diag(2), matrix(0, 3, 3)), "adj.*2.*2")
  bad.adj <- adj
  bad.adj[1, 2] <- NaN
  expect_error(huge.inference(x, diag(2), bad.adj), "adj.*finite")

  sparse <- Matrix::Matrix(adj, sparse = TRUE)
  inf <- huge.inference(x, Matrix::Diagonal(2), sparse)
  expect_true(all(is.finite(inf$p)))

  sparse.data <- Matrix::Matrix(x, sparse = TRUE)
  sparse.inf <- huge.inference(sparse.data, diag(2), sparse)
  expect_identical(sparse.inf$data, sparse.data)
  expect_equal(sparse.inf$p, inf$p, tolerance = 0)
})

test_that("huge.inference rejects non-finite edge results", {
  set.seed(72)
  x <- matrix(rnorm(16), 8, 2)
  adj <- matrix(0, 2, 2)
  huge.t <- 1e308 * matrix(c(1, .5, .5, 1), 2, 2)
  tiny.t <- diag(2) * 1e-200

  for (T in list(huge.t, tiny.t)) {
    for (settings in list(
      list(type = "Gaussian", method = "score"),
      list(type = "Nonparanormal", method = "score"),
      list(type = "Nonparanormal", method = "wald")
    )) {
      expect_error(
        huge.inference(
          x, T, adj,
          type = settings$type, method = settings$method
        ),
        "finite|numerically",
        info = paste(settings$type, settings$method)
      )
    }
  }
})

test_that("Nonparanormal inference gates edges but keeps rank limits", {
  monotone <- matrix(c(0, 1, 0, 1), nrow = 2)
  score <- huge.inference(
    monotone, diag(2), matrix(0, 2, 2),
    type = "Nonparanormal", method = "score"
  )
  expect_true(all(is.nan(diag(score$p))))
  expect_true(all(is.finite(score$p[row(score$p) != col(score$p)])))

  extreme.ties <- cbind(c(1e308, -1e308, 0), c(0, 0, 1))
  for (method in c("score", "wald")) {
    inf <- huge.inference(
      extreme.ties, diag(2), matrix(0, 2, 2),
      type = "Nonparanormal", method = method
    )
    expect_true(
      all(is.finite(inf$p[row(inf$p) != col(inf$p)])),
      info = method
    )
  }
})

test_that("Gaussian inference is stable under extreme finite data scales", {
  data = cbind(
    c(-2, -1, 0, 1, 2, 3),
    c(3, -1, 2, -2, 1, 0)
  )
  reference = huge.inference(
    data, diag(2), matrix(0, 2, 2), type = "Gaussian"
  )

  for(factor in c(1e-300, 1e300)) {
    current = huge.inference(
      data * factor, diag(2), matrix(0, 2, 2), type = "Gaussian"
    )
    expect_equal(current$p, reference$p, tolerance = 1e-14)
    expect_equal(current$error, reference$error)
  }
})
