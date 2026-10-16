# Property-based invariant tests: random inputs, universal properties.
# Complements the fixed-fixture tests — these hunt for unknown bugs rather
# than guarding known fixes. Seeded loops keep runtime bounded (~seconds).

random_shape = function(seed) {
  set.seed(seed)
  list(n = sample(15:120, 1), d = sample(4:50, 1), nlambda = sample(1:8, 1))
}

test_that("all methods: outputs symmetric, sparsity in [0,1] and monotone", {
  for (seed in 1:6) {
    s = random_shape(seed)
    x = matrix(rnorm(s$n * s$d), s$n, s$d)
    for (m in c("mb", "glasso", "ct", "tiger")) {
      fit = if (m == "tiger")
        suppressWarnings(huge(x, method = m, nlambda = s$nlambda,
                              verbose = FALSE)) else
        huge(x, method = m, nlambda = s$nlambda, verbose = FALSE)
      if (m == "tiger") {
        expect_gte(length(fit$path), 1)
        expect_lte(length(fit$path), s$nlambda)
      } else {
        expect_length(fit$path, s$nlambda)
      }
      for (p in fit$path) {
        pm = as.matrix(p)
        expect_identical(max(abs(pm - t(pm))), 0)
        expect_identical(max(abs(diag(pm))), 0)
      }
      expect_true(all(fit$sparsity >= 0 & fit$sparsity <= 1))
      # lambda is decreasing along the path, so sparsity is nondecreasing
      if (s$nlambda > 1) expect_true(all(diff(fit$sparsity) >= -1e-12))
    }
  }
})

test_that("ct path is nested: larger lambda support is subset of smaller", {
  for (seed in 7:9) {
    s = random_shape(seed)
    if (s$nlambda < 2) s$nlambda = 4
    x = matrix(rnorm(s$n * s$d), s$n, s$d)
    fit = huge(x, method = "ct", nlambda = s$nlambda, verbose = FALSE)
    for (i in 1:(s$nlambda - 1)) {
      a = as.matrix(fit$path[[i]]) != 0
      b = as.matrix(fit$path[[i + 1]]) != 0
      expect_true(all(b[a]))  # every edge at larger lambda persists
    }
  }
})

test_that("glasso icov positive definite and consistent with path", {
  for (seed in 10:12) {
    s = random_shape(seed)
    x = matrix(rnorm(s$n * s$d), s$n, s$d)
    fit = huge(x, method = "glasso", nlambda = 4, verbose = FALSE)
    for (i in seq_along(fit$icov)) {
      ic = as.matrix(fit$icov[[i]])
      expect_gt(min(eigen(ic, symmetric = TRUE, only.values = TRUE)$values), 0)
      # path pattern == off-diagonal support of icov
      expect_identical(as.matrix(fit$path[[i]]) != 0,
                       (ic != 0) & (row(ic) != col(ic)))
    }
  }
})

test_that("extreme column scales do not break estimation (internal standardization)", {
  set.seed(31)
  n = 60; d = 12
  x = matrix(rnorm(n * d), n, d)
  x_scaled = sweep(x, 2, 10^seq(-6, 6, length.out = d), "*")
  f1 = huge(x, method = "mb", nlambda = 4, verbose = FALSE)
  f2 = huge(x_scaled, method = "mb", nlambda = 4, verbose = FALSE)
  # correlation-based: identical graphs regardless of column scales
  for (i in 1:4)
    expect_identical(as.matrix(f1$path[[i]]), as.matrix(f2$path[[i]]))
})

test_that("near-singular correlation input stays finite", {
  set.seed(32)
  n = 80; d = 10
  base = rnorm(n)
  x = sapply(1:d, function(j) base + rnorm(n, sd = 0.01))  # corr ~ 0.9999
  for (m in c("mb", "glasso", "tiger")) {
    if (m == "tiger") {
      expect_warning(
        fit <- huge(x, method = m, nlambda = 4, verbose = FALSE),
        "certified prefix"
      )
    } else {
      fit = huge(x, method = m, nlambda = 4, verbose = FALSE)
    }
    expect_true(all(is.finite(fit$sparsity)))
    if (!is.null(fit[["icov"]]))
      for (ic in fit$icov) expect_true(all(is.finite(as.matrix(ic))))
  }
})

test_that("nlambda = 1 works everywhere (degenerate path)", {
  set.seed(33)
  x = matrix(rnorm(50 * 8), 50, 8)
  for (m in c("mb", "glasso", "ct", "tiger")) {
    fit = huge(x, method = m, nlambda = 1, verbose = FALSE)
    expect_length(fit$path, 1)
  }
})
