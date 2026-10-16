stress_graphics_state <- function() {
  graphics::par(
    mfrow = c(2, 3), pty = "m",
    omi = rep(0.05, 4), mai = rep(0.2, 4)
  )
  graphics::par(
    mfg = c(1, 2), cex = 1.37, mex = 0.91,
    usr = c(-2, 3, -4, 5),
    xaxp = c(-2, 3, 5), yaxp = c(-4, 5, 9),
    las = 2, xaxs = "i", yaxs = "i", xpd = NA
  )
  graphics::par(new = FALSE)
}

custom_graphics_state <- function(kind) {
  graphics::par(
    mfrow = c(1, 1), pty = "m", oma = rep(0, 4),
    mar = c(5.1, 4.1, 4.1, 2.1), cex = 1, mex = 1
  )
  graphics::plot.new()
  graphics::plot.window(c(-2, 3), c(-4, 5))
  switch(
    kind,
    omi = graphics::par(omi = c(0.11, 0.22, 0.13, 0.17)),
    fig = graphics::par(fig = c(0.12, 0.83, 0.18, 0.88)),
    fin = graphics::par(fin = c(4.2, 3.8)),
    plt = graphics::par(plt = c(0.17, 0.82, 0.21, 0.79)),
    pin = graphics::par(pin = c(3.7, 3.2)),
    mai = graphics::par(mai = c(0.31, 0.44, 0.27, 0.19)),
    mar = graphics::par(mar = c(2.2, 3.1, 1.8, 0.7)),
    omd = graphics::par(omd = c(0.08, 0.91, 0.12, 0.94)),
    oma = graphics::par(
      oma = c(1.2, 2.1, 0.8, 1.3),
      mar = c(2.2, 3.1, 1.8, 0.7)
    )
  )
  graphics::plot.new()
  graphics::plot.window(c(-2, 3), c(-4, 5))
  graphics::par(xaxp = c(-2, 3, 5), yaxp = c(-4, 5, 9))
}

expect_graphics_state_unchanged <- function(
    draw, configure = stress_graphics_state) {
  plot.file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot.file, width = 6, height = 6)
  plot.device <- grDevices::dev.cur()
  on.exit({
    devices <- grDevices::dev.list()
    if (!is.null(devices) && plot.device %in% devices) {
      grDevices::dev.off(plot.device)
    }
    unlink(plot.file)
  }, add = TRUE)

  configure()
  before <- graphics::par(no.readonly = TRUE)
  set.seed(231)
  value <- force(draw)

  expect_identical(grDevices::dev.cur(), plot.device)
  expect_identical(graphics::par(no.readonly = TRUE), before)
  invisible(value)
}

test_that("plot.huge and huge.plot restore caller graphics state", {
  edge <- Matrix::Matrix(matrix(c(0, 1, 1, 0), 2, 2), sparse = TRUE)
  empty <- Matrix::Matrix(matrix(0, 2, 2), sparse = TRUE)
  fit <- structure(
    list(
      lambda = c(1, 0.5), sparsity = c(0, 1),
      path = list(empty, edge)
    ),
    class = "huge"
  )

  expect_graphics_state_unchanged(plot(fit, align = FALSE))
  expect_graphics_state_unchanged(
    invisible(capture.output(plot(fit, align = TRUE)))
  )
  expect_graphics_state_unchanged(huge.plot(edge, epsflag = FALSE))
})

test_that("generator, sim, select, and ROC plots restore graphics state", {
  edge <- Matrix::Matrix(matrix(c(0, 1, 1, 0), 2, 2), sparse = TRUE)
  sim <- structure(
    list(theta = edge, sigma = diag(2), sigmahat = diag(2)),
    class = "sim"
  )
  selected <- structure(
    list(
      cov.input = FALSE, refit = edge,
      lambda = c(1, 0.5), sparsity = c(0, 1),
      opt.lambda = 0.5, opt.sparsity = 1
    ),
    class = "select"
  )
  roc <- structure(
    list(fp = c(0, 0.25, 1), tp = c(0, 0.75, 1)),
    class = "roc"
  )
  truth <- matrix(0, 3, 3)
  truth[1, 2] <- truth[2, 1] <- 1
  complete <- 1 - diag(3)
  path <- list(matrix(0, 3, 3), complete)

  expect_graphics_state_unchanged(huge.generator(
    n = 8, d = 2, graph = "band",
    vis = TRUE, verbose = FALSE
  ))
  expect_graphics_state_unchanged(plot(sim))
  expect_graphics_state_unchanged(plot(selected))
  expect_graphics_state_unchanged(plot(roc))
  expect_graphics_state_unchanged(
    huge.roc(path, truth, verbose = FALSE)
  )
})

test_that("plotting restores custom figure, plot, margin, and outer regions", {
  roc <- structure(
    list(fp = c(0, 0.25, 1), tp = c(0, 0.75, 1)),
    class = "roc"
  )

  for (kind in c(
    "fig", "fin", "plt", "pin", "mai", "mar", "oma", "omi", "omd"
  )) {
    expect_graphics_state_unchanged(
      plot(roc),
      configure = function() custom_graphics_state(kind)
    )
  }
})

test_that("restoration does not mask an error before the first plot", {
  bad <- structure(
    list(
      theta = matrix("bad", 2, 2),
      sigma = diag(2), sigmahat = diag(2)
    ),
    class = "sim"
  )
  old.warn <- getOption("warn")
  on.exit(options(warn = old.warn), add = TRUE)
  options(warn = 2)

  expect_graphics_state_unchanged(
    expect_error(plot(bad), "numeric or logical")
  )
})

test_that("huge.plot EPS output preserves and returns to the caller device", {
  edge <- Matrix::Matrix(matrix(c(0, 1, 1, 0), 2, 2), sparse = TRUE)
  background.file <- tempfile(fileext = ".pdf")
  plot.file <- tempfile(fileext = ".pdf")
  output.dir <- tempfile(pattern = "huge-eps-")
  dir.create(output.dir)
  grDevices::pdf(background.file, width = 6, height = 6)
  background.device <- grDevices::dev.cur()
  grDevices::pdf(plot.file, width = 6, height = 6)
  caller.device <- grDevices::dev.cur()
  on.exit({
    devices <- grDevices::dev.list()
    if (!is.null(devices) && caller.device %in% devices) {
      grDevices::dev.off(caller.device)
    }
    devices <- grDevices::dev.list()
    if (!is.null(devices) && background.device %in% devices) {
      grDevices::dev.off(background.device)
    }
    unlink(background.file)
    unlink(plot.file)
    unlink(output.dir, recursive = TRUE)
  }, add = TRUE)

  stress_graphics_state()
  before <- graphics::par(no.readonly = TRUE)
  before.devices <- grDevices::dev.list()
  huge.plot(
    edge, epsflag = TRUE, graph.name = "state", cur.num = 1,
    location = output.dir
  )
  eps.file <- file.path(output.dir, "state1.eps")

  expect_identical(grDevices::dev.cur(), caller.device)
  expect_identical(grDevices::dev.list(), before.devices)
  expect_identical(graphics::par(no.readonly = TRUE), before)
  expect_true(file.exists(eps.file))
  expect_gt(file.info(eps.file)$size, 0)

  testthat::local_mocked_bindings(
    plot = function(...) stop("forced plot failure"),
    .package = "huge"
  )
  expect_error(
    huge.plot(
      edge, epsflag = TRUE, graph.name = "state", cur.num = 2,
      location = output.dir
    ),
    "forced plot failure"
  )
  expect_identical(grDevices::dev.cur(), caller.device)
  expect_identical(grDevices::dev.list(), before.devices)
  expect_identical(graphics::par(no.readonly = TRUE), before)
})
