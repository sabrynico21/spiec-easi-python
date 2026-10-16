test_that("huge.select with ebic works for glasso", {
  set.seed(50)
  L <- huge.generator(n = 80, d = 20, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)
  sel <- huge.select(fit, criterion = "ebic", verbose = FALSE)

  expect_s3_class(sel, "select")
  expect_true(!is.null(sel$refit))
  expect_true(!is.null(sel$opt.lambda))
  expect_true(!is.null(sel$opt.index))
  expect_true(sel$opt.lambda > 0)
  expect_true(sel$opt.index >= 1 && sel$opt.index <= length(sel$lambda))
  expect_equal(dim(sel$refit), c(20, 20))
})

test_that("huge.select with ric works for mb", {
  set.seed(51)
  L <- huge.generator(n = 80, d = 20, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "mb", verbose = FALSE)
  sel <- huge.select(fit, criterion = "ric", verbose = FALSE)

  expect_s3_class(sel, "select")
  expect_true(!is.null(sel$refit))
  expect_true(sel$opt.lambda > 0)
})

test_that("RIC is invariant to extreme finite uniform data scales", {
  x = cbind(
    c(-2, -1, 0, 1, 2, 3),
    c(3, -1, 2, -2, 1, 0),
    c(0, 2, -1, 3, -2, 1)
  )
  selected = lapply(c(1, 1e-300, 1e300), function(factor) {
    fit = huge(
      x * factor, method = "ct", lambda = c(.8, .2), verbose = FALSE
    )
    huge.select(
      fit, criterion = "ric", rep.num = nrow(x), verbose = FALSE
    )
  })

  expect_true(all(vapply(
    selected, function(value) is.finite(value$opt.lambda), logical(1)
  )))
  expect_equal(
    vapply(selected, function(value) value$opt.lambda, numeric(1)),
    rep(selected[[1]]$opt.lambda, 3),
    tolerance = 1e-14
  )
})

test_that("RIC handles a zero lower endpoint without weakening solvers", {
  x = cbind(
    c(-1, -1, 1, 1),
    c(-1, -1, 1, 1)
  )

  for(method in c("ct", "mb", "glasso", "tiger")) {
    fit = huge(
      x, method = method, lambda = 1, verbose = FALSE
    )
    if(method %in% c("glasso", "tiger")) {
      selected = expect_warning(
        huge.select(
          fit, criterion = "ric", rep.num = nrow(x), verbose = FALSE
        ),
        "RIC selected lambda = 0.*original fitted path",
        info = method
      )
      nearest = which.min(abs(fit$lambda))
      expect_identical(
        as.matrix(selected$refit), as.matrix(fit$path[[nearest]]),
        info = method
      )
    } else {
      selected = expect_no_warning(
        huge.select(
          fit, criterion = "ric", rep.num = nrow(x), verbose = FALSE
        )
      )
      expect_equal(sum(selected$refit), 2, info = method)
    }
    expect_identical(selected$opt.lambda, 0)
    expect_true(is.finite(selected$opt.sparsity), info = method)
  }
})

test_that("RIC zero does not erase a representable weak correlation", {
  first = c(-1, -1, 1, 1)
  orthogonal = c(-1, 1, -1, 1)
  # Below the former fixed 8-epsilon cutoff, but safely above the
  # pair-specific dot-product roundoff bound.
  x = cbind(first, orthogonal + 1.5e-15 * first)
  fit = huge(
    x, method = "ct", lambda = .5, verbose = FALSE
  )
  selected = huge.select(
    fit, criterion = "ric", rep.num = nrow(x), verbose = FALSE
  )

  expect_identical(selected$opt.lambda, 0)
  expect_equal(sum(selected$refit), 2)
})

test_that("huge.select with stars works for mb", {
  set.seed(52)
  L <- huge.generator(n = 80, d = 20, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "mb", verbose = FALSE)
  sel <- huge.select(fit, criterion = "stars", rep.num = 5, verbose = FALSE)

  expect_s3_class(sel, "select")
  expect_true(!is.null(sel$refit))
  expect_true(!is.null(sel$variability))
  expect_equal(length(sel$variability), length(sel$lambda))
  expect_true(all(sel$variability >= 0))
})

test_that("huge.select keeps TIGER RIC and rejects unsafe StARS grids", {
  set.seed(321)
  x = matrix(rnorm(120 * 20), 120, 20)
  fit = huge(x, method = "tiger", nlambda = 4,
             lambda.min.ratio = 0.4, verbose = FALSE)

  expect_error(
    huge.select(fit, criterion = "stars", rep.num = 3,
                verbose = FALSE),
    "TIGER.*StARS|common certified prefix"
  )

  set.seed(7)
  ric = huge.select(fit, criterion = "ric", rep.num = 10,
                    verbose = FALSE)
  expect_true(is.finite(ric$opt.lambda) && ric$opt.lambda > 0)
  expect_equal(dim(ric$refit), c(20, 20))
})

test_that("huge.select uses method-specific default criteria", {
  set.seed(53)
  x <- matrix(rnorm(100 * 10), nrow = 100)
  expected <- c(mb = "ric", ct = "stars", glasso = "ebic", tiger = "ric")

  for (method in names(expected)) {
    fit <- huge(x, method = method, nlambda = 3,
                lambda.min.ratio = .5, verbose = FALSE)
    set.seed(531)
    selected <- huge.select(
      fit, rep.num = 3, num.cores = 1, verbose = FALSE
    )
    expect_identical(selected$criterion, unname(expected[[method]]),
                     info = method)
    expect_false(is.null(selected$refit), info = method)
  }
})

test_that("huge.select rejects malformed estimates and criteria", {
  set.seed(54)
  x <- matrix(rnorm(80 * 8), nrow = 80)
  fit <- huge(x, method = "mb", nlambda = 3, verbose = FALSE)

  expect_error(huge.select(NULL, verbose = FALSE), "est")
  expect_error(huge.select(list(), verbose = FALSE), "est")

  incomplete <- fit
  incomplete$path <- NULL
  expect_error(huge.select(incomplete, verbose = FALSE), "est|path")

  for (criterion in list("bogus", c("ric", "stars"), 1, NA_character_)) {
    expect_error(
      huge.select(fit, criterion = criterion, verbose = FALSE),
      "criterion"
    )
  }
  expect_error(
    huge.select(fit, criterion = "ebic", verbose = FALSE),
    "ebic|glasso"
  )
})

test_that("huge.select stops cleanly for covariance input", {
  set.seed(55)
  x <- matrix(rnorm(80 * 8), nrow = 80)
  fit <- huge(cor(x), method = "glasso", nlambda = 3, verbose = FALSE)
  original.class <- class(fit)

  expect_error(
    huge.select(fit, verbose = FALSE),
    "covariance matrix"
  )
  expect_identical(class(fit), original.class)
  expect_s3_class(fit, "huge")
})

test_that("huge.select validates only active criterion parameters", {
  set.seed(56)
  x <- matrix(rnorm(80 * 8), nrow = 80)
  mb <- huge(x, method = "mb", nlambda = 3, verbose = FALSE)
  glasso <- huge(x, method = "glasso", nlambda = 3, verbose = FALSE)

  for (rep.num in list(0, -1, 1.5, NA_real_, Inf, c(2, 3), TRUE)) {
    expect_error(
      huge.select(
        mb, criterion = "ric", rep.num = rep.num, verbose = FALSE
      ),
      "rep.num"
    )
  }
  expect_error(
    huge.select(
      mb, criterion = "stars", rep.num = 2,
      num.cores = 0, verbose = FALSE
    ),
    "num.cores"
  )
  expect_error(
    huge.select(
      mb, criterion = "stars", rep.num = 2,
      stars.thresh = 0, verbose = FALSE
    ),
    "stars.thresh"
  )
  expect_error(
    huge.select(
      mb, criterion = "stars", rep.num = 2,
      stars.subsample.ratio = 0, verbose = FALSE
    ),
    "stars.subsample.ratio"
  )
  expect_error(
    huge.select(
      mb, criterion = "stars", rep.num = 2,
      stars.subsample.ratio = 1e-6, verbose = FALSE
    ),
    "stars.subsample.ratio"
  )
  expect_error(
    huge.select(
      glasso, criterion = "ebic", ebic.gamma = NA_real_,
      verbose = FALSE
    ),
    "ebic.gamma"
  )

  expect_no_error(huge.select(
    glasso, criterion = "ebic", rep.num = 0, num.cores = 0,
    stars.thresh = 0, stars.subsample.ratio = 0, verbose = FALSE
  ))
  expect_no_error(huge.select(
    mb, criterion = "ric", rep.num = 2, num.cores = 0,
    ebic.gamma = NA_real_, stars.thresh = 0,
    stars.subsample.ratio = 0, verbose = FALSE
  ))
  expect_no_error(huge.select(
    mb, criterion = "stars", rep.num = 2, num.cores = 1,
    ebic.gamma = NA_real_, verbose = FALSE
  ))
})

test_that("RIC handles an exact zero-correlation boundary", {
  x <- cbind(
    c(-1, -1, 1, 1),
    c(-1, 1, -1, 1)
  )

  for(method in c("ct", "mb", "glasso", "tiger")) {
    fit <- huge(x, method = method, nlambda = 3, verbose = FALSE)
    selected <- expect_no_warning(huge.select(
      fit, criterion = "ric", rep.num = nrow(x), verbose = FALSE
    ))

    expect_equal(selected$opt.lambda, 0, tolerance = 1e-15, info = method)
    expect_equal(sum(selected$refit), 0, info = method)
    expect_equal(selected$opt.sparsity, 0, info = method)
  }
})

test_that("RIC refits preserve explicit square-data routing", {
  x <- matrix(c(
    3, 2, 2,
    2, 4, 4,
    2, 4, 7
  ), 3, 3)
  fit <- huge(
    x, method = "tiger", lambda = 1, input.type = "data",
    verbose = FALSE
  )
  selected <- expect_no_warning(huge.select(
    fit, criterion = "ric", rep.num = nrow(x), verbose = FALSE
  ))
  data.refit <- huge.tiger(
    x, lambda = selected$opt.lambda, input.type = "data",
    verbose = FALSE
  )
  automatic.refit <- huge.tiger(
    x, lambda = selected$opt.lambda, verbose = FALSE
  )

  expect_false(selected$cov.input)
  expect_identical(
    as.matrix(selected$refit), as.matrix(data.refit$path[[1]])
  )
  expect_false(identical(
    as.matrix(selected$refit), as.matrix(automatic.refit$path[[1]])
  ))
})

test_that("RIC handles the single-variable empty-graph boundary", {
  x <- matrix(c(-1, 0, 1), ncol = 1)

  for (method in c("ct", "mb", "glasso", "tiger")) {
    fit <- huge(x, method = method, nlambda = 3, verbose = FALSE)
    selected <- huge.select(
      fit, criterion = "ric", rep.num = nrow(x), verbose = FALSE
    )

    expect_s3_class(selected, "select")
    expect_equal(selected$opt.lambda, 0, tolerance = 1e-15, info = method)
    expect_equal(dim(selected$refit), c(1L, 1L), info = method)
    expect_equal(sum(selected$refit), 0, info = method)
    expect_identical(selected$opt.sparsity, 0, info = method)

    if (method %in% c("mb", "tiger")) {
      selected_default <- huge.select(
        fit, rep.num = nrow(x), verbose = FALSE
      )
      expect_identical(selected_default$criterion, "ric", info = method)
      expect_equal(selected_default$opt.lambda, 0, tolerance = 1e-15,
                   info = method)
      expect_identical(selected_default$opt.sparsity, 0, info = method)
      expect_equal(sum(selected_default$refit), 0, info = method)
    }
  }
})

test_that("StARS rejects single-variable graph paths", {
  x <- matrix(c(-1, 0, 1), ncol = 1)

  for (method in c("ct", "mb", "glasso")) {
    fit <- huge(x, method = method, nlambda = 3, verbose = FALSE)
    expect_error(
      huge.select(
        fit, criterion = "stars", rep.num = 2, verbose = FALSE
      ),
      "StARS requires at least two variables",
      info = method
    )
  }

  tiger <- huge(x, method = "tiger", nlambda = 3, verbose = FALSE)
  expect_error(
    huge.select(tiger, criterion = "stars", rep.num = 2, verbose = FALSE),
    "TIGER.*StARS|common certified prefix"
  )
})

test_that("StARS requires a non-increasing lambda path", {
  set.seed(58)
  x <- matrix(rnorm(40 * 4), nrow = 40)

  ascending <- huge(
    x, method = "ct", lambda = c(.1, .3, .5), verbose = FALSE
  )
  expect_error(
    huge.select(ascending, rep.num = 2, verbose = FALSE),
    "non-increasing"
  )

  unordered <- huge(
    x, method = "ct", lambda = c(.5, .1, .3), verbose = FALSE
  )
  expect_error(
    huge.select(
      unordered, criterion = "stars", rep.num = 2, verbose = FALSE
    ),
    "non-increasing"
  )

  tied <- huge(
    x, method = "ct", lambda = c(.5, .5, .2), verbose = FALSE
  )
  expect_no_error(huge.select(
    tied, criterion = "stars", rep.num = 2, verbose = FALSE
  ))

  expect_no_error(huge.select(
    ascending, criterion = "ric", rep.num = 2, verbose = FALSE
  ))

  tiger <- huge(x, method = "tiger", nlambda = 3, verbose = FALSE)
  tiger$lambda <- rev(tiger$lambda)
  expect_error(
    huge.select(
      tiger, criterion = "stars", rep.num = 2, verbose = FALSE
    ),
    "TIGER.*StARS|common certified prefix"
  )
})
