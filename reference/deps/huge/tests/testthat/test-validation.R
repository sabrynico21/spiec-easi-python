test_that("generated paths reject invalid nlambda", {
  S <- diag(4)
  invalid <- list(0, -1, 1.5, NA_real_, Inf, c(2, 3), TRUE, "2")

  for (method in c("mb", "glasso", "ct", "tiger")) {
    for (value in invalid) {
      expect_error(
        huge(S, method = method, nlambda = value, verbose = FALSE),
        "nlambda",
        info = paste(method, deparse(value))
      )
    }
  }
})

test_that("generated paths reject invalid lambda.min.ratio", {
  S <- diag(4)
  invalid <- list(0, -0.1, 1.1, NA_real_, Inf, c(.1, .2), TRUE, ".1")

  for (method in c("mb", "glasso", "ct", "tiger")) {
    for (value in invalid) {
      expect_error(
        huge(
          S, method = method, nlambda = 2,
          lambda.min.ratio = value, verbose = FALSE
        ),
        "lambda.min.ratio",
        info = paste(method, deparse(value))
      )
    }
    expect_no_error(huge(
      S, method = method, nlambda = 1,
      lambda.min.ratio = 1, verbose = FALSE
    ))
  }
})

test_that("MB and glasso default lambda tails retain representable values", {
  S = matrix(c(1, 1e-200, 1e-200, 1), 2, 2)
  ratio = 1e-200
  smallest.positive = .Machine$double.xmin * .Machine$double.eps

  generated = huge:::.huge_default_lambda(
    S, d = 2, nlambda = 3, lambda.min.ratio = ratio
  )$lambda

  expect_length(generated, 3)
  expect_true(all(is.finite(generated)))
  expect_true(all(generated > 0))
  expect_true(all(diff(generated) <= 0))
  expect_equal(generated[[1]], 1e-200)
  expect_equal(log(generated[[2]]), log(1e-300), tolerance = 1e-12)
  expect_identical(generated[[3]], smallest.positive)
})

test_that("explicit lambda validates only its active branch", {
  S <- diag(4)
  valid <- list(
    mb = c(.2, .1),
    glasso = c(.2, .1),
    ct = c(.2, 0),
    tiger = c(.2, .1)
  )

  for (method in names(valid)) {
    fit <- huge(
      S, method = method, lambda = valid[[method]],
      nlambda = 0, lambda.min.ratio = 2, verbose = FALSE
    )
    expect_equal(fit$lambda, valid[[method]], info = method)
  }
})

test_that("explicit lambda rejects unsafe values", {
  S <- diag(4)
  positive.methods <- c("mb", "glasso", "tiger")
  invalid.positive <- list(
    numeric(), 0, -1, NA_real_, Inf, c(.2, NA_real_), "0.2"
  )

  for (method in positive.methods) {
    for (value in invalid.positive) {
      expect_error(
        huge(S, method = method, lambda = value, verbose = FALSE),
        "lambda",
        info = paste(method, deparse(value))
      )
    }
  }

  for (value in list(numeric(), -1, NA_real_, Inf, c(.2, NA_real_), "0.2")) {
    expect_error(
      huge(S, method = "ct", lambda = value, verbose = FALSE),
      "lambda",
      info = deparse(value)
    )
  }
  expect_no_error(huge(
    S, method = "ct", lambda = c(.2, 0), verbose = FALSE
  ))
})

test_that("explicit lambda is scalar or one-dimensional", {
  S <- matrix(c(1, .2, .2, 1), nrow = 2)
  path <- c(.9, .8, .7, .6)
  valid <- list(
    scalar = .9,
    vector = path,
    one.dim.array = array(path, dim = length(path))
  )
  invalid <- list(
    row.matrix = matrix(path, nrow = 1),
    column.matrix = matrix(path, ncol = 1),
    square.matrix = matrix(path, nrow = 2),
    byrow.matrix = matrix(path, nrow = 2, byrow = TRUE),
    three.dim.array = array(path, dim = c(2, 1, 2)),
    rank.before.domain = matrix(c(.9, .8, NaN, .6), nrow = 2)
  )

  for(method in c("ct", "mb", "glasso", "tiger")) {
    direct <- get(paste0("huge.", method))
    for(name in names(valid)) {
      value <- valid[[name]]
      expect_identical(
        direct(S, lambda = value, verbose = FALSE)$lambda,
        as.numeric(value),
        info = paste(method, name, "direct")
      )
      expect_identical(
        huge(
          S, method = method, lambda = value, verbose = FALSE
        )$lambda,
        as.numeric(value),
        info = paste(method, name, "generic")
      )
    }
    for(name in names(invalid)) {
      value <- invalid[[name]]
      expect_error(
        direct(S, lambda = value, verbose = FALSE),
        "one-dimensional",
        info = paste(method, name, "direct")
      )
      expect_error(
        huge(
          S, method = method, lambda = value, verbose = FALSE
        ),
        "one-dimensional",
        info = paste(method, name, "generic")
      )
    }
  }
})

test_that("MB and glasso lambda paths allow ties and reject increases", {
  S <- matrix(c(1, .3, .3, 1), nrow = 2)
  tied <- c(.5, .5, .2)

  for(method in c("mb", "glasso")) {
    direct <- get(paste0("huge.", method))
    expect_identical(
      direct(S, lambda = tied, verbose = FALSE)$lambda,
      tied,
      info = paste(method, "direct ties")
    )
    expect_identical(
      huge(S, method = method, lambda = tied, verbose = FALSE)$lambda,
      tied,
      info = paste(method, "generic ties")
    )

    automatic <- direct(
      S, nlambda = 3, lambda.min.ratio = 1, verbose = FALSE
    )
    replay <- direct(
      S, lambda = automatic$lambda, verbose = FALSE
    )
    expect_identical(diff(automatic$lambda), c(0, 0))
    expect_identical(replay$lambda, automatic$lambda)
    expect_identical(replay$path, automatic$path)
    expect_identical(replay$sparsity, automatic$sparsity)

    expect_error(
      direct(S, lambda = c(.1, .2), verbose = FALSE),
      "non-increasing",
      info = paste(method, "direct increase")
    )
    expect_error(
      huge(
        S, method = method, lambda = c(.5, .1, .2), verbose = FALSE
      ),
      "non-increasing",
      info = paste(method, "generic unordered")
    )
    expect_error(
      huge(
        S, method = method, lambda = c(.1, .2, NA_real_),
        verbose = FALSE
      ),
      "finite",
      info = paste(method, "finite error priority")
    )
    expect_error(
      huge(
        S, method = method, lambda = c(.1, .2, -1),
        verbose = FALSE
      ),
      "strictly positive",
      info = paste(method, "sign error priority")
    )
  }
})

test_that("native MB and glasso entries fail closed on unsafe inputs", {
  S <- matrix(c(1, .3, .3, 1), nrow = 2)
  idx <- matrix(c(1L, 0L), nrow = 1)
  native.mb <- function(lambda, nlambda = length(lambda), d = 2L,
                        matrix = S) {
    .Call(
      "_huge_SPMBgraph", matrix, lambda, as.integer(nlambda),
      as.integer(d), PACKAGE = "huge"
    )
  }
  native.mb.scr <- function(lambda, nlambda = length(lambda), d = 2L,
                            matrix = S, index = idx, nscr = 1L) {
    .Call(
      "_huge_SPMBscr", matrix, lambda, as.integer(nlambda),
      as.integer(d), index, as.integer(nscr), PACKAGE = "huge"
    )
  }
  native.glasso <- function(lambda, matrix = S) {
    .Call(
      "_huge_hugeglasso", matrix, lambda, FALSE, FALSE, FALSE,
      PACKAGE = "huge"
    )
  }
  entries <- list(
    mb = native.mb,
    mb.scr = native.mb.scr,
    glasso = native.glasso
  )

  for(name in names(entries)) {
    entry <- entries[[name]]
    expect_no_error(entry(c(.2, .2, .1)))
    expect_error(
      entry(c(.1, .2)), "non-increasing",
      info = paste(name, "increase")
    )
    for(value in list(0, -.1, NA_real_, Inf)) {
      expect_error(
        entry(value), "positive and finite",
        info = paste(name, deparse(value))
      )
    }
    expect_error(
      entry(c(.1, .2, NaN)), "positive and finite",
      info = paste(name, "finite error priority")
    )
    expect_error(
      entry(c(.1, .2, -.1)), "positive and finite",
      info = paste(name, "sign error priority")
    )
    expect_error(
      entry(numeric()), "nlambda|at least one",
      info = paste(name, "empty path")
    )
  }

  expect_error(native.mb(c(.2, .1), nlambda = 1L), "length")
  expect_error(native.mb.scr(c(.2, .1), nlambda = 1L), "length")
  extra.column <- cbind(S, c(0, 0))
  expect_error(native.mb(.1, matrix = extra.column), "d by d")
  expect_error(native.mb.scr(.1, matrix = extra.column), "d by d")
  expect_error(native.mb(.1, d = 0L), "d must be positive")
  expect_error(native.mb.scr(.1, d = 0L), "d must be positive")
  expect_error(
    native.mb.scr(.1, index = matrix(integer(), 0, 2), nscr = 0L),
    "between 1 and d - 1"
  )
  expect_error(
    native.mb.scr(.1, index = matrix(integer(), 0, 3), nscr = 0L),
    "between 1 and d - 1"
  )
  expect_error(
    native.mb.scr(.1, index = cbind(idx, 0L)),
    "nscr by d"
  )
  expect_error(native.glasso(.1, matrix = extra.column), "square")
})

test_that("registered native lambda entries reject multidimensional values", {
  S <- matrix(c(1, .2, .2, 1), nrow = 2)
  raw <- cbind(c(-1, -1, 1, 1), c(-1, 1, -1, 1))
  idx <- matrix(c(1L, 0L), nrow = 1)
  path <- c(.9, .8, .7, .6)
  entries <- list(
    mb = function(lambda) .Call(
      "_huge_SPMBgraph", S, lambda, as.integer(length(lambda)), 2L,
      PACKAGE = "huge"
    ),
    mb.scr = function(lambda) .Call(
      "_huge_SPMBscr", S, lambda, as.integer(length(lambda)), 2L,
      idx, 1L, PACKAGE = "huge"
    ),
    glasso = function(lambda) .Call(
      "_huge_hugeglasso", S, lambda, FALSE, FALSE, FALSE,
      PACKAGE = "huge"
    ),
    tiger.legacy = function(lambda) .Call(
      "_huge_SPMBgraphsqrt", raw, lambda,
      as.integer(length(lambda)), 2L, PACKAGE = "huge"
    ),
    tiger.fit = function(lambda) .Call(
      "_huge_SPMBgraphsqrtFit", S, lambda,
      as.integer(length(lambda)), 2L, TRUE, .1, PACKAGE = "huge"
    )
  )
  valid <- list(.9, path, array(path, dim = length(path)))
  invalid <- list(
    matrix(path, nrow = 1),
    matrix(path, ncol = 1),
    matrix(path, nrow = 2),
    array(path, dim = c(2, 1, 2)),
    matrix(c(.9, .8, NaN, .6), nrow = 2)
  )

  for(name in names(entries)) {
    entry <- entries[[name]]
    for(value in valid)
      expect_no_error(entry(value))
    for(value in invalid) {
      expect_error(
        entry(value), "one-dimensional", info = name
      )
    }
  }

  expect_error(
    .Call(
      "_huge_SPMBgraphsqrt", raw, matrix(.9, nrow = 1),
      0L, 2L, PACKAGE = "huge"
    ),
    "one-dimensional"
  )
})

test_that("registered native RIC entry validates dimensions and rotations", {
  X <- scale(cbind(
    c(-1, -1, 1, 1),
    c(-1, 1, -1, 1)
  ))
  native.ric <- function(
      matrix = X, d = ncol(matrix), n = nrow(matrix),
      r = seq_len(n), t = length(r)) {
    .Call(
      "_huge_RIC", matrix, as.integer(d), as.integer(n), r,
      as.integer(t), PACKAGE = "huge"
    )
  }

  expect_true(is.finite(native.ric()))
  expect_true(is.finite(native.ric(r = c(0, 1, 2, 3, 4), t = 5)))
  expect_true(is.finite(native.ric(r = c(1, 2, NaN), t = 2)))

  # Both mismatches retain more physical storage than the declared problem
  # would read if validation regressed; never probe an undersized buffer.
  expect_error(
    native.ric(
      matrix = rbind(X, c(0, 0)), d = 2L, n = 4L,
      r = 1:4, t = 4L
    ),
    "dimensions.*n by d"
  )
  expect_error(
    native.ric(
      matrix = cbind(X, 0), d = 2L, n = 4L,
      r = 1:4, t = 4L
    ),
    "dimensions.*n by d"
  )
  expect_error(native.ric(d = 0L), "n and d must be positive")
  expect_error(native.ric(n = 0L), "n and d must be positive")

  for(value in c(NA_real_, NaN, Inf, -Inf)) {
    bad <- X
    bad[1, 1] <- value
    expect_error(
      native.ric(matrix = bad), "X must contain only finite",
      info = deparse(value)
    )
  }

  for(value in list(1.5, NaN, Inf, -Inf)) {
    expect_error(
      native.ric(r = value, t = 1L), "finite integers",
      info = deparse(value)
    )
  }
  for(value in c(-1, nrow(X) + 1)) {
    expect_error(
      native.ric(r = value, t = 1L), "rotation indices",
      info = deparse(value)
    )
  }
  expect_error(
    native.ric(r = .Machine$integer.max + 1, t = 1L),
    "integer range"
  )

  expect_error(native.ric(t = 0L), "t must be positive")
  expect_error(native.ric(t = -1L), "t must be positive")
  expect_error(native.ric(r = 1:2, t = 3L), "length of r")

  singleton <- matrix(seq_len(4), ncol = 1)
  expect_identical(
    native.ric(singleton, d = 1L, n = 4L, r = 1, t = 1L),
    0
  )
  singleton[1, 1] <- NA_real_
  expect_error(
    native.ric(singleton, d = 1L, n = 4L, r = 1, t = 1L),
    "X must contain only finite"
  )
})

test_that("mb and tiger require exact sym values", {
  set.seed(571)
  x <- matrix(rnorm(60 * 6), nrow = 60)
  invalid <- list(
    "bogus", "o", NULL, character(), c("or", "and"), NA_character_, 1
  )

  for (method in c("mb", "tiger")) {
    for (value in invalid) {
      expect_error(
        huge(
          x, method = method, lambda = c(.4, .2),
          sym = value, verbose = FALSE
        ),
        "sym",
        info = paste(method, deparse(value))
      )
    }
    expect_no_error(huge(
      x, method = method, lambda = c(.4, .2),
      sym = "or", verbose = FALSE
    ))
    expect_no_error(huge(
      x, method = method, lambda = c(.4, .2),
      sym = "and", verbose = FALSE
    ))
  }
})

run_estimator_entry <- function(
  method, x, direct, input.type = NULL, lambda = .2
) {
  args <- list(x = x, lambda = lambda, verbose = FALSE)
  if(!is.null(input.type))
    args$input.type <- input.type
  if (direct) {
    do.call(get(paste0("huge.", method)), args)
  } else {
    args$method <- method
    do.call(huge, args)
  }
}

test_that("input.type disambiguates square symmetric observations", {
  x <- matrix(c(
    2, 1, 0,
    1, 2, 1,
    0, 1, 2
  ), 3, 3, byrow = TRUE)

  for(method in c("ct", "mb", "glasso", "tiger")) {
    for(direct in c(FALSE, TRUE)) {
      automatic <- run_estimator_entry(method, x, direct, lambda = 1)
      data.fit <- run_estimator_entry(
        method, x, direct, input.type = "data", lambda = 1
      )
      covariance.fit <- run_estimator_entry(
        method, x, direct, input.type = "covariance", lambda = 1
      )

      expect_true(automatic$cov.input, info = paste(method, direct))
      expect_false(data.fit$cov.input, info = paste(method, direct))
      expect_true(covariance.fit$cov.input, info = paste(method, direct))
    }
  }
})

test_that("TIGER input.type keeps correlation and lambda selection native", {
  x <- matrix(c(
    2, 1, 0,
    1, 2, 1,
    0, 1, 2
  ), 3, 3, byrow = TRUE)

  raw.native <- .Call(
    "_huge_SPMBgraphsqrtFit", x, NULL, 1L, 3L, FALSE, .1,
    PACKAGE = "huge"
  )
  covariance.native <- .Call(
    "_huge_SPMBgraphsqrtFit", x, NULL, 1L, 3L, TRUE, .1,
    PACKAGE = "huge"
  )
  raw.fit <- huge.tiger(
    x, nlambda = 1, input.type = "data", verbose = FALSE
  )
  covariance.fit <- huge.tiger(
    x, nlambda = 1, input.type = "covariance", verbose = FALSE
  )
  automatic <- huge.tiger(x, nlambda = 1, verbose = FALSE)

  expect_identical(raw.fit$lambda, raw.native$lambda)
  expect_identical(covariance.fit$lambda, covariance.native$lambda)
  expect_identical(automatic$lambda, covariance.fit$lambda)
  expect_false(identical(raw.fit$lambda, covariance.fit$lambda))

  routed <- huge:::.huge_validate_estimation_input(
    x, input.type = "covariance", prepare.covariance = FALSE
  )
  expect_true(routed$cov.input)
  expect_null(routed$covariance)
  expect_null(routed$correlation)
})

test_that("input.type rejects invalid and incompatible declarations", {
  x <- matrix(seq_len(6), 3, 2)

  for(value in list("raw", NA_character_, c("data", "auto"), 1)) {
    expect_error(
      huge.ct(x, lambda = 1, input.type = value, verbose = FALSE),
      "input.type"
    )
  }
  expect_error(
    huge.ct(
      x, lambda = 1, input.type = "covariance", verbose = FALSE
    ),
    "square"
  )
  expect_error(
    huge.tiger(
      matrix(c(1, 0, 1, 1), 2, 2),
      lambda = 1, input.type = "covariance", verbose = FALSE
    ),
    "symmetric"
  )
})

test_that("estimator entries reject malformed and non-finite inputs", {
  invalid <- list(
    vector = list(value = 1:4, message = "numeric matrix"),
    nonnumeric = list(
      value = matrix(letters[1:4], 2, 2), message = "numeric matrix"
    ),
    zero.rows = list(
      value = matrix(numeric(), 0, 2), message = "non-empty"
    ),
    zero.columns = list(
      value = matrix(numeric(), 2, 0), message = "non-empty"
    ),
    missing = list(
      value = matrix(c(1, 2, 3, NA, 5, 6), 3, 2), message = "finite"
    ),
    infinite = list(
      value = matrix(c(1, 2, 3, Inf, 5, 6), 3, 2), message = "finite"
    )
  )

  for (method in c("ct", "mb", "glasso", "tiger")) {
    for (direct in c(FALSE, TRUE)) {
      for (case in invalid) {
        expect_error(
          run_estimator_entry(method, case$value, direct),
          case$message,
          info = paste(method, if (direct) "direct" else "huge")
        )
      }
    }
  }
})

test_that("estimator entries reject undefined raw correlations", {
  invalid <- list(
    one.observation = list(
      value = matrix(c(1, 2, 3), 1, 3),
      message = "at least two observations"
    ),
    constant.column = list(
      value = cbind(seq_len(4), rep(1, 4)),
      message = "constant column"
    )
  )

  for (method in c("ct", "mb", "glasso", "tiger")) {
    for (direct in c(FALSE, TRUE)) {
      for (case in invalid) {
        expect_error(
          run_estimator_entry(method, case$value, direct),
          case$message,
          info = paste(method, if (direct) "direct" else "huge")
        )
      }
    }
  }
})

test_that("estimator entries reject invalid covariance matrices", {
  invalid <- list(
    zero.diagonal = list(
      value = diag(c(0, 1)), message = "positive.*diagonal"
    ),
    negative.diagonal = list(
      value = diag(c(-1, 1)), message = "positive.*diagonal"
    ),
    cauchy.schwarz = list(
      value = matrix(c(1, 2, 2, 1), 2, 2),
      message = "valid covariance"
    )
  )

  for (method in c("ct", "mb", "glasso", "tiger")) {
    for (direct in c(FALSE, TRUE)) {
      for (case in invalid) {
        expect_error(
          run_estimator_entry(method, case$value, direct),
          case$message,
          info = paste(method, if (direct) "direct" else "huge")
        )
      }
    }
  }
})

test_that("non-PSD covariance is accepted only when glasso regularizes it", {
  S <- matrix(c(
    1, .9, .9,
    .9, 1, 0,
    .9, 0, 1
  ), 3, 3, byrow = TRUE)

  for(method in c("ct", "mb", "tiger")) {
    for(direct in c(FALSE, TRUE)) {
      expect_error(
        run_estimator_entry(method, S, direct, lambda = .5),
        "positive semidefinite",
        info = paste(method, if (direct) "direct" else "huge")
      )
    }
  }

  for(direct in c(FALSE, TRUE)) {
    args <- list(
      x = S, lambda = .5, cov.output = TRUE, verbose = FALSE
    )
    if(direct) {
      fit <- do.call(huge.glasso, args)
    } else {
      args$method <- "glasso"
      fit <- do.call(huge, args)
    }
    precision <- fit$icov[[1]]
    covariance <- fit$cov[[1]]
    expect_gt(
      min(eigen(precision, symmetric = TRUE, only.values = TRUE)$values),
      0
    )
    expect_gt(
      min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values),
      0
    )
    expect_true(
      max(rowSums(abs(covariance %*% precision - diag(3)))) <= 1e-2
    )
  }
})

test_that("covariance validation accepts PSD and spectral roundoff", {
  roundoff <- matrix(-.5 - 1e-14, 3, 3)
  diag(roundoff) <- 1
  cases <- list(singular = matrix(1, 2, 2), roundoff = roundoff)

  for (S in cases) {
    for (method in c("ct", "mb", "glasso", "tiger")) {
      expect_no_error(huge(
        S, method = method, lambda = 1, verbose = FALSE
      ))
    }
  }
})

test_that("covariance projection preserves exact subnormal diagonal", {
  smallest <- .Machine$double.xmin * .Machine$double.eps
  S <- diag(rep(smallest, 2))
  expect_gt(smallest, 0)

  for (method in c("ct", "mb", "glasso", "tiger")) {
    for (direct in c(FALSE, TRUE)) {
      fit <- run_estimator_entry(method, S, direct)
      expect_true(fit$cov.input, info = paste(method, direct))
      expect_identical(fit$sparsity, 0, info = paste(method, direct))
      expect_identical(
        as.matrix(fit$path[[1]]), matrix(0, 2, 2),
        info = paste(method, direct)
      )
      if (!is.null(fit$icov)) {
        expect_true(
          all(is.finite(fit$icov[[1]])), info = paste(method, direct)
        )
      }
    }
  }
})

test_that("covariance validation clips only correlation roundoff", {
  S <- matrix(c(1, 1 + 5e-9, 1 + 5e-9, 1), 2, 2)
  fit <- huge.ct(S, lambda = .9, verbose = FALSE)

  expect_true(fit$cov.input)
  expect_equal(sum(fit$path[[1]]), 2)
})

test_that("accepted near-symmetric covariance is projected consistently", {
  S <- diag(3)
  S[1, 2] <- .2 + 1e-15
  S[2, 1] <- .2 - 1e-15
  expect_true(isSymmetric(S))
  projected <- S / 2 + t(S) / 2

  for (method in c("ct", "mb", "glasso", "tiger")) {
    fit <- huge(S, method = method, lambda = .2, verbose = FALSE)
    reference <- huge(
      projected, method = method, lambda = .2, verbose = FALSE
    )

    expect_identical(
      as.matrix(fit$path[[1]]), as.matrix(reference$path[[1]]),
      info = method
    )
    expect_identical(
      as.matrix(fit$path[[1]]), t(as.matrix(fit$path[[1]])),
      info = method
    )
  }
})

test_that("R covariance detection uses the native per-entry tolerance", {
  aggregate.roundoff <- matrix(c(1, 1e-14, 3e-14, 1), 2, 2)
  expect_true(isSymmetric(aggregate.roundoff))

  fit <- huge.tiger(aggregate.roundoff, lambda = 1, verbose = FALSE)
  expect_false(fit$cov.input)

  near <- matrix(c(1, .5, .5 + 1e-15, 1), 2, 2)
  expect_true(huge:::.huge_is_covariance_input(near))
  expect_true(huge.tiger(near, lambda = 1, verbose = FALSE)$cov.input)

  maximum <- .Machine$double.xmax
  large.asymmetry <- matrix(
    c(maximum, maximum, maximum * (1 - 1e-13), maximum), 2, 2
  )
  expect_true(isSymmetric(large.asymmetry))
  expect_false(huge:::.huge_is_covariance_input(large.asymmetry))

  maximum.integer <- .Machine$integer.max
  integer.asymmetry <- matrix(
    c(1L, -maximum.integer, maximum.integer, 1L), 2, 2
  )
  expect_no_warning(
    expect_false(huge:::.huge_is_covariance_input(integer.asymmetry))
  )
  expect_no_warning({
    fit <- huge.tiger(integer.asymmetry, lambda = 1, verbose = FALSE)
  })
  expect_false(fit$cov.input)
})

test_that("covariance detection is invariant to finite uniform scale", {
  raw = matrix(c(1, 4, 7, 2, 5, 8, 3, 6, 10), 3, 3)
  reference = huge.ct(raw, lambda = .95, verbose = FALSE)
  tiger.reference = huge.tiger(raw, nlambda = 1, verbose = FALSE)
  minimum = .Machine$double.xmin * .Machine$double.eps

  for (factor in c(1, 1e-15, 1e-200, minimum)) {
    current = raw * factor
    expect_false(huge:::.huge_is_covariance_input(current))

    for(method in c("mb", "glasso")) {
      expect_false(
        huge(
          current, method = method, lambda = .95, verbose = FALSE
        )$cov.input,
        info = paste(method, format(factor))
      )
    }

    fit = huge.ct(current, lambda = .95, verbose = FALSE)
    expect_false(fit$cov.input)
    expect_identical(
      as.matrix(fit$path[[1]]), as.matrix(reference$path[[1]])
    )

    tiger = huge.tiger(current, nlambda = 1, verbose = FALSE)
    expect_false(tiger$cov.input)
    expect_equal(tiger$lambda, tiger.reference$lambda, tolerance = 1e-14)

    expect_error(
      .Call(
        "_huge_SPMBgraphsqrtFit", current, 1, 1L, 3L, TRUE, .1,
        PACKAGE = "huge"
      ),
      "symmetric"
    )
  }

  near = matrix(c(1, .5, .5 + 1e-15, 1), 2, 2)
  near = near * 1e-200
  expect_true(huge:::.huge_is_covariance_input(near))
  expect_true(huge.tiger(near, lambda = 1, verbose = FALSE)$cov.input)
  expect_no_error(
    .Call(
      "_huge_SPMBgraphsqrtFit", near, 1, 1L, 2L, TRUE, .1,
      PACKAGE = "huge"
    )
  )

  extreme = matrix(
    c(.Machine$double.xmax, -1, 1, 1e-300), 2, 2
  )
  expect_false(huge:::.huge_is_covariance_input(extreme))
})

test_that("covariance symmetry tolerance uses implied-correlation scale", {
  tolerance = 100 * .Machine$double.eps
  native.tiger = function(S) {
    .Call(
      "_huge_SPMBgraphsqrtFit", S, 1, 1L, 2L, TRUE, .1,
      PACKAGE = "huge"
    )
  }

  for(variance in c(.Machine$double.xmin, 1, .Machine$double.xmax)) {
    within = matrix(c(
      variance, tolerance * variance,
      0, variance
    ), 2, 2)
    outside = matrix(c(
      variance, 101 * .Machine$double.eps * variance,
      0, variance
    ), 2, 2)

    expect_true(huge:::.huge_is_covariance_input(within))
    expect_no_error(native.tiger(within))
    expect_false(huge:::.huge_is_covariance_input(outside))
    expect_error(native.tiger(outside), "symmetric")
  }

  maximum = .Machine$double.xmax
  negligible = matrix(c(maximum, -1, 1, maximum), 2, 2)
  expect_true(huge:::.huge_is_covariance_input(negligible))
  expect_true(
    huge.tiger(negligible, lambda = 1, verbose = FALSE)$cov.input
  )
  expect_no_error(native.tiger(negligible))

  extreme = matrix(c(maximum, -1, 1, 1e-300), 2, 2)
  expect_false(huge:::.huge_is_covariance_input(extreme))
  expect_error(native.tiger(extreme), "symmetric")

  invalid.diagonal = matrix(c(0, -1, 1, 1), 2, 2)
  expect_error(native.tiger(invalid.diagonal), "positive.*diagonal")
})
