test_that("huge.generator produces valid output for all graph types", {
  set.seed(1)
  for (g in c("hub", "band", "cluster", "random", "scale-free")) {
    L <- huge.generator(n = 60, d = 20, graph = g, verbose = FALSE)
    expect_s3_class(L, "sim")
    expect_equal(dim(L$data), c(60, 20))
    expect_equal(dim(L$sigma), c(20, 20))
    expect_equal(dim(L$omega), c(20, 20))
    expect_equal(L$sigmahat, cor(L$data), tolerance = 1e-12)
    expect_equal(diag(L$sigmahat), rep(1, 20), tolerance = 0)

    # sigma should be positive definite
    eig <- eigen(L$sigma, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eig > 0))

    # theta should be symmetric binary
    theta <- as.matrix(L$theta)
    expect_equal(theta, t(theta))
    expect_true(all(theta %in% c(0, 1)))
    expect_true(all(diag(theta) == 0))

    # sparsity should be in [0, 1]
    expect_true(L$sparsity >= 0 && L$sparsity <= 1)
  }
})

test_that("huge.generator respects dimension parameters", {
  set.seed(2)
  L <- huge.generator(n = 100, d = 30, graph = "hub", verbose = FALSE)
  expect_equal(nrow(L$data), 100)
  expect_equal(ncol(L$data), 30)
  expect_equal(nrow(L$sigma), 30)
})

for (graph in c("hub", "band", "cluster", "random", "scale-free")) {
  test_that(paste("huge.generator supports d = 1 for", graph), {
    set.seed(222)
    L <- huge.generator(n = 8, d = 1, graph = graph, verbose = FALSE)

    expect_equal(dim(L$data), c(8, 1))
    expect_equal(dim(L$sigma), c(1, 1))
    expect_equal(dim(L$omega), c(1, 1))
    expect_equal(dim(L$sigmahat), c(1, 1))
    expect_equal(as.matrix(L$theta), matrix(0, 1, 1))
    expect_equal(L$sigma, matrix(1, 1, 1))
    expect_equal(L$omega, matrix(1, 1, 1))
    expect_equal(L$sigmahat, matrix(1, 1, 1))
    expect_identical(L$sparsity, 0)

    printed <- paste(capture.output(print(L)), collapse = "\n")
    expect_match(printed, "Dimension: d = 1", fixed = TRUE)
    expect_match(printed, "Sparsity level: 0", fixed = TRUE)
    expect_false(grepl("NaN", printed, fixed = TRUE))
  })
}

test_that("band graphs support the maximal meaningful bandwidth", {
  for (d in c(2, 5)) {
    set.seed(223)
    L <- huge.generator(
      n = 8, d = d, graph = "band", g = d - 1,
      verbose = FALSE
    )

    expect_equal(as.matrix(L$theta), 1 - diag(d))
    expect_identical(L$sparsity, 1)
  }
})

test_that("huge.generator rejects invalid active parameters before side effects", {
  expect_clean_error <- function(args, pattern) {
    set.seed(224)
    seed.before <- .Random.seed
    warnings <- character()
    printed <- capture.output(
      result <- withCallingHandlers(
        tryCatch(
          do.call(
            huge.generator,
            c(list(n = 8, d = 5, verbose = TRUE), args)
          ),
          error = identity
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    )

    expect_s3_class(result, "error")
    if (inherits(result, "error")) {
      expect_match(conditionMessage(result), pattern)
    }
    expect_identical(printed, character())
    expect_identical(warnings, character())
    expect_true(identical(.Random.seed, seed.before))
  }

  for (v in list(0, -1, NA_real_, Inf, TRUE, c(0.2, 0.3), "0.3")) {
    expect_clean_error(list(graph = "random", v = v), "v")
  }
  for (u in list(0, -1, NA_real_, Inf, TRUE, c(0.1, 0.2), "0.1")) {
    expect_clean_error(list(graph = "band", u = u), "u")
  }
  for (g in list(0, -1, 1.5, 3e9, NA_real_, Inf, TRUE, c(1, 2), "2")) {
    for (graph in c("hub", "cluster", "band")) {
      expect_clean_error(list(graph = graph, g = g), "g")
    }
  }
  for (prob in list(-0.1, 1.1, NA_real_, Inf, TRUE, c(0.2, 0.3), "0.3")) {
    for (graph in c("random", "cluster")) {
      expect_clean_error(list(graph = graph, prob = prob), "prob")
    }
  }
})

test_that("huge.generator rejects invalid flags before side effects", {
  expect_clean_flag_error <- function(name, value) {
    set.seed(229)
    seed.before <- .Random.seed
    args <- list(
      n = 8, d = 5, graph = "random",
      vis = FALSE, verbose = FALSE
    )
    args[name] <- list(value)
    warnings <- character()
    printed <- capture.output(
      result <- withCallingHandlers(
        tryCatch(do.call(huge.generator, args), error = identity),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    )

    expect_s3_class(result, "error")
    if (inherits(result, "error")) {
      expect_match(conditionMessage(result), name)
    }
    expect_identical(printed, character())
    expect_identical(warnings, character())
    expect_true(identical(.Random.seed, seed.before))
  }

  invalid <- list(
    NULL, NA, logical(), c(TRUE, FALSE), 0, 2, "FALSE"
  )
  for (name in c("vis", "verbose")) {
    for (value in invalid) {
      expect_clean_flag_error(name, value)
    }
  }
})

test_that("huge.generator ignores graph-inactive parameters", {
  expect_same_generation <- function(graph, name, value) {
    set.seed(225)
    expected <- huge.generator(
      n = 8, d = 5, graph = graph,
      verbose = FALSE
    )
    expected.seed <- .Random.seed

    set.seed(225)
    args <- list(n = 8, d = 5, graph = graph, verbose = FALSE)
    args[[name]] <- value
    actual <- do.call(huge.generator, args)

    expect_identical(actual, expected)
    expect_true(identical(.Random.seed, expected.seed))
  }

  for (graph in c("random", "scale-free")) {
    expect_same_generation(graph, "g", "inactive")
  }
  for (graph in c("hub", "band", "scale-free")) {
    expect_same_generation(graph, "prob", "inactive")
  }
})

test_that("huge.generator accepts parameter boundaries without changing semantics", {
  for (graph in c("random", "cluster")) {
    for (prob in c(0, 1)) {
      set.seed(226)
      L <- huge.generator(
        n = 8, d = 5, graph = graph, prob = prob,
        verbose = FALSE
      )
      expect_true(all(is.finite(L$data)))
    }
  }

  for (graph in c("hub", "cluster")) {
    set.seed(227)
    oversized <- huge.generator(
      n = 8, d = 5, graph = graph, g = 7,
      verbose = FALSE
    )
    oversized.seed <- .Random.seed
    set.seed(227)
    capped <- huge.generator(
      n = 8, d = 5, graph = graph, g = 5,
      verbose = FALSE
    )
    expect_identical(oversized, capped)
    expect_true(identical(oversized.seed, .Random.seed))
  }

  set.seed(228)
  oversized <- huge.generator(
    n = 8, d = 5, graph = "band", g = 7,
    verbose = FALSE
  )
  oversized.seed <- .Random.seed
  set.seed(228)
  capped <- huge.generator(
    n = 8, d = 5, graph = "band", g = 4,
    verbose = FALSE
  )
  expect_identical(oversized, capped)
  expect_true(identical(oversized.seed, .Random.seed))
})

test_that("huge.generator visualization restores graphics parameters", {
  plot.file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot.file, width = 6, height = 6)
  on.exit({
    grDevices::dev.off()
    unlink(plot.file)
  }, add = TRUE)
  graphics::par(
    mfrow = c(2, 3), pty = "m",
    omi = rep(0.05, 4), mai = rep(0.2, 4)
  )
  graphics::par(
    mfg = c(1, 2), cex = 1.37, mex = 0.91,
    usr = c(-2, 3, -4, 5)
  )
  graphics::par(new = FALSE)
  parameters <- c(
    "mfrow", "mfg", "pty", "omi", "mai", "cex", "mex", "usr", "new"
  )
  before <- graphics::par(parameters)

  set.seed(230)
  invisible(huge.generator(
    n = 8, d = 5, graph = "band",
    vis = TRUE, verbose = FALSE
  ))

  after <- graphics::par(parameters)
  expect_identical(after, before)
})

test_that("huge.generator rejects invalid dimensions before side effects", {
  expect_clean_error <- function(args, pattern) {
    set.seed(73)
    seed.before <- .Random.seed
    printed <- capture.output(
      result <- tryCatch(
        do.call(
          huge.generator,
          c(args, list(graph = "random", verbose = TRUE))
        ),
        error = identity
      )
    )

    expect_s3_class(result, "error")
    expect_match(conditionMessage(result), pattern)
    expect_identical(printed, character())
    expect_identical(.Random.seed, seed.before)
  }

  for (n in list(
    1, 0, -1, 1.5, 3e9, NA_real_, Inf, TRUE, c(2, 3), "2"
  )) {
    expect_clean_error(list(n = n, d = 4), "n")
  }
  for (d in list(
    0, -1, 1.5, 3e9, NA_real_, Inf, TRUE, c(2, 3), "2"
  )) {
    expect_clean_error(list(n = 8, d = d), "d")
  }
})

test_that("huge.generator rejects unknown graph types before side effects", {
  for (graph in list("unknown", c("random", "hub"), NA_character_, 1)) {
    set.seed(74)
    seed.before <- .Random.seed
    printed <- capture.output(
      result <- tryCatch(
        huge.generator(
          n = 8, d = 4, graph = graph,
          verbose = TRUE
        ),
        error = identity
      )
    )

    expect_s3_class(result, "error")
    expect_match(conditionMessage(result), "graph.*one of")
    expect_identical(printed, character())
    expect_identical(.Random.seed, seed.before)
  }
})
