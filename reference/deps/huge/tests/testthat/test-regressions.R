# Regression tests pinning bugs fixed in 2.0.0. Each test names the defect it
# guards; if one fails, the corresponding fix has been reverted or broken.

test_that("tiger icov is exactly symmetric (2.0.0: in-place symmetrization bug)", {
  set.seed(42)
  L = huge.generator(n = 60, d = 30, graph = "hub", verbose = FALSE)
  fit = huge(L$data, method = "tiger", nlambda = 5, verbose = FALSE)
  for (ic in fit$icov)
    expect_identical(max(abs(as.matrix(ic) - t(as.matrix(ic)))), 0)
})

test_that("glasso precision, path, df, sparsity, and loglik agree", {
  S = matrix(c(
    1, .517180, -.517349, .326007, .263074, -.434894, -.594540,
    .517180, 1, -.007732, .807446, -.171208, -.164722, .110568,
    -.517349, -.007732, 1, .108058, -.413103, .652356, .356858,
    .326007, .807446, .108058, 1, .131565, -.283979, -.107923,
    .263074, -.171208, -.413103, .131565, 1, -.473447, -.695847,
    -.434894, -.164722, .652356, -.283979, -.473447, 1, .390406,
    -.594540, .110568, .356858, -.107923, -.695847, .390406, 1
  ), 7, 7, byrow = TRUE)
  fit = huge(S, method = "glasso", nlambda = 8,
             lambda.min.ratio = .005, cov.output = TRUE,
             verbose = FALSE)
  denom = nrow(S) * (nrow(S) - 1)

  for (i in seq_along(fit$path)) {
    path = as.matrix(fit$path[[i]])
    precision = as.matrix(fit$icov[[i]])
    expect_equal(precision, t(precision), tolerance = 0, info = paste("lambda", i))
    expect_equal(path, t(path), tolerance = 0, info = paste("lambda", i))
    support = precision != 0
    diag(support) = FALSE
    expect_identical(path != 0, support, info = paste("lambda", i))
    edge.count = sum(path[upper.tri(path)] != 0)
    expect_identical(fit$df[i], as.integer(edge.count))
    expect_equal(fit$sparsity[i], 2 * edge.count / denom, tolerance = 0)
    determinant.info = determinant(precision, logarithm = TRUE)
    expect_identical(determinant.info$sign, 1L)
    expect_lte(max(rowSums(abs(fit$cov[[i]] %*% precision - diag(7)))),
               1e-2)
    direct = as.numeric(determinant.info$modulus) -
      sum(diag(S %*% precision))
    expect_equal(fit$loglik[i], direct, tolerance = 1e-10)
  }
})

test_that("RIC selection is scale-invariant (2.0.0: unstandardized data bug)", {
  set.seed(9)
  L = huge.generator(n = 100, d = 40, graph = "hub", verbose = FALSE)
  fit1 = huge(L$data, method = "mb", nlambda = 8, verbose = FALSE)
  fit2 = huge(L$data * 100, method = "mb", nlambda = 8, verbose = FALSE)
  set.seed(101); s1 = huge.select(fit1, criterion = "ric", verbose = FALSE)
  set.seed(101); s2 = huge.select(fit2, criterion = "ric", verbose = FALSE)
  expect_equal(s1$opt.lambda, s2$opt.lambda, tolerance = 1e-12)
  expect_equal(s1$opt.sparsity, s2$opt.sparsity)
})

test_that("R raw-data estimators preserve correlation at extreme finite scales", {
  x = cbind(c(-1, 0, 1), c(-1, 1, 0))
  minimum = .Machine$double.xmin * .Machine$double.eps
  inputs = list(
    x * 1e308,
    x * 1e-200,
    x * minimum,
    sweep(x, 2, c(1e308, 1e-200), "*")
  )

  for (method in c("ct", "mb", "glasso")) {
    explicit.reference = huge(
      x, method = method, lambda = .4, verbose = FALSE
    )
    automatic.reference = huge(
      x, method = method, nlambda = 3,
      lambda.min.ratio = .5, verbose = FALSE
    )

    for (current in inputs) {
      expect_no_warning({
        explicit = huge(
          current, method = method, lambda = .4, verbose = FALSE
        )
        automatic = huge(
          current, method = method, nlambda = 3,
          lambda.min.ratio = .5, verbose = FALSE
        )
      })
      expect_false(explicit$cov.input)
      expect_equal(explicit$lambda, explicit.reference$lambda, tolerance = 0)
      expect_equal(explicit$sparsity, explicit.reference$sparsity, tolerance = 0)
      expect_equal(automatic$lambda, automatic.reference$lambda, tolerance = 1e-14)
      expect_equal(
        automatic$sparsity, automatic.reference$sparsity, tolerance = 0
      )
      for (i in seq_along(explicit$path))
        expect_identical(
          as.matrix(explicit$path[[i]]),
          as.matrix(explicit.reference$path[[i]])
        )
      for (i in seq_along(automatic$path))
        expect_identical(
          as.matrix(automatic$path[[i]]),
          as.matrix(automatic.reference$path[[i]])
        )
    }
  }
})

test_that("R fast correlation preserves adjacent maximum-scale ULPs", {
  maximum = .Machine$double.xmax
  ulp = 2^971
  indices = cbind(
    first = 0:6,
    second = c(0, 1, 3, 6, 2, 5, 4)
  )
  x = maximum - ulp * indices

  expect_equal(huge:::.huge_fast_cor(x), cor(indices), tolerance = 1e-14)

  set.seed(3)
  n = sample(3:200, 1)
  values = rnorm(n)
  repeated = cbind(positive = values, duplicate = values, negative = -values)
  correlation = huge:::.huge_fast_cor(repeated)
  expect_identical(dimnames(correlation), list(colnames(repeated), colnames(repeated)))
  expect_lte(max(abs(correlation)), 1)
  thresholded = huge.ct(repeated, lambda = 1, verbose = FALSE)
  expect_identical(thresholded$sparsity, 0)
})

test_that("npn output columns have unit sd and names survive (2.0.0 fixes)", {
  set.seed(3)
  x = matrix(rnorm(80 * 12), 80, 12,
             dimnames = list(paste0("r", 1:80), paste0("c", 1:12)))
  for (f in c("shrinkage", "truncation")) {
    z = huge.npn(x, npn.func = f, verbose = FALSE)
    expect_equal(unname(apply(z, 2, sd)), rep(1, 12), tolerance = 1e-12)
    expect_identical(dimnames(z), dimnames(x))
  }
  # skeptic on named input errored before 2.0.0 (n-length rownames on d x d)
  k = huge.npn(x, npn.func = "skeptic", verbose = FALSE)
  expect_identical(dim(k), c(12L, 12L))
  expect_identical(rownames(k), colnames(x))
  expect_identical(colnames(k), colnames(x))
  # unknown npn.func silently returned x unchanged before 2.0.0
  expect_error(huge.npn(x, npn.func = "bogus", verbose = FALSE))
})

test_that("generator sigma is a correlation matrix with symmetric omega", {
  set.seed(7)
  L = huge.generator(n = 50, d = 25, graph = "hub", verbose = FALSE)
  expect_equal(unname(diag(L$sigma)), rep(1, 25))
  expect_identical(max(abs(L$omega - t(L$omega))), 0)
  expect_lt(max(abs(L$sigma %*% L$omega - diag(25))), 1e-10)
})

test_that("select sets opt.cov only for glasso with cov (2.0.0: partial-match bug)", {
  set.seed(3)
  x = matrix(rnorm(120 * 30), 120, 30)
  for (m in c("mb", "ct", "tiger")) {
    fit = huge(x, method = m, nlambda = 5, verbose = FALSE)
    set.seed(9); s = huge.select(fit, criterion = "ric", verbose = FALSE)
    expect_null(s[["opt.cov"]])
  }
  fit = huge(x, method = "glasso", nlambda = 5, cov.output = TRUE, verbose = FALSE)
  set.seed(9); s = huge.select(fit, criterion = "ric", verbose = FALSE)
  expect_false(is.null(s[["opt.cov"]]))
})

test_that("huge() rejects unknown method (2.0.0: returned half-empty object)", {
  x = matrix(rnorm(200), 50, 4)
  expect_error(huge(x, method = "bogus", verbose = FALSE))
})

test_that("huge.select preserves its historical positional verbose argument", {
  expect_identical(
    tail(names(formals(huge.select)), 3),
    c("rep.num", "verbose", "num.cores")
  )

  set.seed(3)
  x = matrix(rnorm(60 * 10), 60, 10)
  fit = huge(x, method = "mb", nlambda = 4, verbose = FALSE)
  set.seed(19)
  named = huge.select(fit, criterion = "ric", rep.num = 2, verbose = FALSE)
  set.seed(19)
  positional = huge.select(fit, "ric", 0.5, 0.1, NULL, 2, FALSE)
  expect_identical(positional, named)
})

test_that("stars num.cores > 1 reproduces serial results exactly", {
  skip_on_os("windows")
  with.check.core.limit = function(code) {
    previous = Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
    on.exit({
      if(is.na(previous))
        Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
      else
        do.call(
          Sys.setenv,
          setNames(list(previous), "_R_CHECK_LIMIT_CORES_")
        )
    }, add = TRUE)
    do.call(
      Sys.setenv,
      setNames(list("TRUE"), "_R_CHECK_LIMIT_CORES_")
    )
    force(code)
  }

  set.seed(3)
  x = matrix(rnorm(120 * 30), 120, 30)
  fit = huge(x, method = "mb", nlambda = 5, verbose = FALSE)
  set.seed(42); s1 = huge.select(fit, criterion = "stars", rep.num = 4,
                                 num.cores = 1, verbose = FALSE)
  set.seed(42); expect_warning(
    s2 <- huge.select(fit, criterion = "stars", rep.num = 4,
                      num.cores = 2, verbose = FALSE),
    "OpenMP.*BLAS"
  )
  s3 = with.check.core.limit({
    set.seed(42)
    expect_warning(
      huge.select(fit, criterion = "stars", rep.num = 4,
                  num.cores = 20, verbose = FALSE),
      "OpenMP.*BLAS"
    )
  })
  expect_identical(s1$opt.index, s2$opt.index)
  expect_equal(s1$variability, s2$variability, tolerance = 0)
  expect_identical(s1$opt.index, s3$opt.index)
  expect_equal(s1$variability, s3$variability, tolerance = 0)
})

test_that("stars serial accumulation matches the batch reference exactly", {
  set.seed(90)
  x = matrix(rnorm(50 * 8), 50, 8)
  fit = huge(x, method = "ct", nlambda = 4, verbose = FALSE)
  rep.num = 4
  ratio = 0.75
  subsample.size = floor(nrow(x) * ratio)

  # Reproduce the historical implementation: retain every subsample path,
  # then merge replication-by-replication within each lambda value.
  set.seed(901)
  ind.list = lapply(seq_len(rep.num), function(i)
    sample(seq_len(nrow(x)), subsample.size, replace = FALSE))
  paths = lapply(ind.list, function(ind)
    huge.ct(x[ind, ], lambda = fit$lambda, verbose = FALSE)$path)
  expected.seed = .Random.seed

  expected.merge = lapply(seq_along(fit$lambda), function(lambda.idx) {
    count = Matrix::Matrix(0, ncol(x), ncol(x))
    for(rep.idx in seq_len(rep.num))
      count = count + paths[[rep.idx]][[lambda.idx]]
    count / rep.num
  })
  expected.variability = vapply(expected.merge, function(m) {
    4 * (sum(m) - sum(m^2)) / (ncol(x) * (ncol(x) - 1))
  }, numeric(1))
  expect_true(any(expected.variability > 0))
  stars.cross = which(expected.variability >= 0.1)
  expected.index = if(length(stars.cross) == 0) {
    length(fit$lambda)
  } else {
    max(stars.cross[1] - 1, 1)
  }

  set.seed(901)
  selected = huge.select(
    fit, criterion = "stars", stars.subsample.ratio = ratio,
    rep.num = rep.num, num.cores = 1, verbose = FALSE
  )
  selected.seed = .Random.seed

  for(lambda.idx in seq_along(expected.merge))
    expect_equal(as.matrix(selected$merge[[lambda.idx]]),
                 as.matrix(expected.merge[[lambda.idx]]),
                 tolerance = 0, info = paste("lambda", lambda.idx))
  expect_true(all(vapply(selected$merge, inherits, logical(1),
                         what = "sparseMatrix")))
  expect_equal(selected$variability, expected.variability, tolerance = 0)
  expect_identical(selected$opt.index, expected.index)
  expect_equal(selected$refit, fit$path[[expected.index]], tolerance = 0)
  expect_identical(selected.seed, expected.seed)
})

test_that("mb/tiger sparse paths have sorted dgCMatrix indices (core collect_sorted)", {
  set.seed(5)
  x = matrix(rnorm(100 * 30), 100, 30)
  for (m in c("mb", "tiger")) {
    fit = huge(x, method = m, nlambda = 5, verbose = FALSE)
    for (p in fit$path) expect_true(validObject(p, complete = TRUE))
  }
})

test_that("solvers are silent normally and glasso rejects an uncertified limit case", {
  set.seed(1)
  x = matrix(rnorm(100 * 30), 100, 30)
  for (m in c("mb", "glasso", "tiger"))
    expect_no_warning(huge(x, method = m, nlambda = 5, verbose = FALSE))

  # Pathological: n < d with 5 latent factors and an absurdly small lambda.
  # The old iteration-limit result is not a usable inverse pair, so glasso
  # must reject it instead of returning it with only a warning.
  set.seed(50)
  n = 40; d = 60
  base = matrix(rnorm(n * 5), n, 5)
  xp = base[, sample(5, d, replace = TRUE)] + matrix(rnorm(n * d, sd = 0.05), n, d)
  expect_error(
    huge(scale(xp), method = "glasso", lambda = 0.001, verbose = FALSE),
    "inconsistent precision and covariance"
  )
})

test_that("ROC rejects truth matrices with only one edge class", {
  path = list(matrix(0, 3, 3))
  expect_error(huge.roc(path, matrix(0, 3, 3), verbose = FALSE), "ROC/AUC")

  complete = matrix(1, 3, 3)
  diag(complete) = 0
  expect_error(huge.roc(path, complete, verbose = FALSE), "ROC/AUC")
})

test_that("ROC AUC is invariant to equal-FPR order and duplicates", {
  graph = function(edges) {
    value = matrix(0, 4, 4)
    for(index in seq_len(nrow(edges))) {
      left = edges[index, 1]
      right = edges[index, 2]
      value[left, right] = 1
      value[right, left] = 1
    }
    value
  }

  truth = graph(matrix(c(1, 2, 1, 3), ncol = 2, byrow = TRUE))
  low = graph(matrix(c(3, 4), ncol = 2))
  high = graph(matrix(
    c(1, 2, 1, 3, 3, 4), ncol = 2, byrow = TRUE
  ))
  complete = matrix(1, 4, 4)
  diag(complete) = 0

  paths = list(
    list(low, high, complete),
    list(high, low, complete),
    list(low, high, high, complete)
  )
  auc = vapply(
    paths,
    function(path) huge.roc(path, truth, verbose = FALSE)$AUC,
    numeric(1)
  )

  expect_equal(auc, rep(.75, 3))
})

test_that("RIC certifies numerical zeros regardless of summation order (2.0.1: ATLAS)", {
  # The rotated inner products below are mathematically zero, so RIC's optimum
  # is exactly zero.  A conforming BLAS may return any value inside the dot
  # product's roundoff interval, and an optimized one (ATLAS) does: 2.0.0
  # returned a tiny positive lambda here, which bypassed the zero-lambda
  # fallback and made glasso fit singular data almost unregularized.
  # A cyclic shift relabels which rotation index sees which pairing without
  # changing the multiset of rotated inner products, so the certified optimum
  # must not move.  (An arbitrary permutation is not a valid probe here: it
  # changes which cyclic rotations exist, and therefore the true optimum.)
  x = cbind(c(-1, -1, 1, 1), c(-1, -1, 1, 1))
  shifts = list(1:4, c(2, 3, 4, 1), c(3, 4, 1, 2), c(4, 1, 2, 3))

  for(shift in shifts) {
    fit = huge(x[shift, ], method = "ct", lambda = 1, verbose = FALSE)
    selected = huge.select(
      fit, criterion = "ric", rep.num = 4, verbose = FALSE
    )
    expect_identical(selected$opt.lambda, 0, info = paste(shift, collapse = ""))
  }

  # Scaling permutes nothing but changes every rounding: still exactly zero.
  for(factor in c(1, 8, 1 / 1024, 1e6)) {
    fit = huge(x * factor, method = "ct", lambda = 1, verbose = FALSE)
    selected = huge.select(
      fit, criterion = "ric", rep.num = 4, verbose = FALSE
    )
    expect_identical(selected$opt.lambda, 0, info = format(factor))
  }
})

test_that("RIC keeps correlations that working precision can represent", {
  # Guards the other side of the certification bound: the fix must not become
  # a fixed tolerance that erases weak-but-real structure.
  #
  # This case is built to sit just above the bound, which is the only way to
  # show the bound is not over-aggressive.  For an input where some rotation
  # gives an exact zero, RIC's optimum is 0 whatever the tolerance is, so such
  # inputs cannot discriminate.  Here every cyclic rotation is nonzero and the
  # smallest is driven to a chosen multiple of the bound: `t` scales the
  # component of the second column along the rotation that would otherwise
  # vanish, so min |C[1,2]| lands at roughly 24x and 2.4x the bound.
  n = 5
  ramp = c(-1, -.5, 0, .5, 1)
  ramp = (ramp - mean(ramp)) / sd(ramp)
  rotate = function(v, r) v[((seq_along(v) - 1 + r) %% length(v)) + 1]
  set.seed(3)
  offset = rnorm(n)
  direction = rnorm(n)
  # Remove the rotation-0 component, so `t` alone sets the smallest rotation.
  base = offset - (sum(ramp * offset) / sum(ramp * direction)) * direction
  for(t in c(1e-13, 1e-14)) {
    x = cbind(ramp, base + t * direction)
    fit = huge(x, method = "ct", lambda = .5, verbose = FALSE)
    selected = huge.select(
      fit, criterion = "ric", rep.num = n, verbose = FALSE
    )
    expect_gt(selected$opt.lambda, 0)
  }

  # Weak-but-representable structure on top of an orthogonal column.
  first = c(-1, -1, 1, 1)
  orthogonal = c(-1, 1, -1, 1)
  for(strength in c(1e-3, 1e-6, 1e-9)) {
    x = cbind(first, orthogonal + strength * first)
    fit = huge(x, method = "ct", lambda = .5, verbose = FALSE)
    selected = huge.select(
      fit, criterion = "ric", rep.num = 4, verbose = FALSE
    )
    expect_equal(sum(selected$refit), 2, info = format(strength))
  }

  set.seed(404)
  L = huge.generator(n = 100, d = 30, graph = "hub", verbose = FALSE)
  for(method in c("mb", "glasso", "tiger", "ct")) {
    fit = huge(L$data, method = method, nlambda = 6, verbose = FALSE)
    set.seed(11)
    selected = huge.select(fit, criterion = "ric", verbose = FALSE)
    expect_gt(selected$opt.lambda, 0)
    expect_true(is.finite(selected$opt.sparsity), info = method)
  }
})

test_that("RIC zero routing does not depend on a bitwise-zero lambda", {
  # 2.0.0 gated the safe refit on `opt.lambda == 0`, so a residual too small
  # for glasso to certify but not bitwise zero reached the solver and raised
  # an uncaught error.  Route on the roundoff scale instead: any lambda the
  # solver cannot certify must fall back to the fitted path with a warning.
  x = cbind(c(-1, -1, 1, 1), c(-1, -1, 1, 1))
  for(method in c("glasso", "tiger")) {
    fit = huge(x, method = method, lambda = 1, verbose = FALSE)
    selected = expect_warning(
      huge.select(fit, criterion = "ric", rep.num = 4, verbose = FALSE),
      "RIC selected lambda = 0.*original fitted path",
      info = method
    )
    expect_identical(selected$opt.lambda, 0, info = method)
    expect_identical(
      as.matrix(selected$refit),
      as.matrix(fit$path[[which.min(abs(fit$lambda))]]),
      info = method
    )
  }
})
