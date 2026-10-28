#!/usr/bin/env Rscript
# Batch R-side status + discrete-output classification for the SpiecEasi stress
# runner (tests/r_reference/run_status_stress.py).
#
# For every archived case it recomputes, from the injected counts and the
# injected StARS subsample indices (no R RNG is used):
#   clr transform -> correlation -> lambda.max -> lambda path -> full-data MB
#   path -> per-replica MB paths -> StARS variability -> selected index -> refit
# and writes status (finite/nonfinite/error), opt_index, refit edge set and the
# per-edge merge frequencies at the selected lambda.
#
# Usage: Rscript stress_status_R.R <cases_dir> <status.csv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript stress_status_R.R <cases_dir> <status.csv>")
}
cases_dir <- args[[1]]
status_path <- args[[2]]

suppressWarnings(suppressMessages({
  library(jsonlite)
}))

REF_R <- Sys.getenv("SPIEC_EASI_REF_R", "reference/R")
if (!dir.exists(REF_R)) stop("reference R source not found: ", REF_R)
suppressWarnings(suppressMessages(pkgload::load_all(REF_R, quiet = TRUE)))

empty_row <- function(case) {
  data.frame(
    case = case, status = "error", error_class = "unknown",
    opt_index = NA_integer_, n_edges = NA_integer_, edges = "",
    lambda_max = NA_real_, stars_margin = NA_real_, merge_opt = "",
    seconds = NA_real_, stringsAsFactors = FALSE
  )
}

classify_case <- function(dir) {
  case <- basename(dir)
  start <- proc.time()[["elapsed"]]
  params <- fromJSON(file.path(dir, "params.json"))
  counts <- as.matrix(read.csv(file.path(dir, "counts.csv"), header = FALSE))
  idx <- as.matrix(read.csv(file.path(dir, "subsample_indices.csv"), header = FALSE))
  storage.mode(counts) <- "double"
  storage.mode(idx) <- "integer"

  p <- ncol(counts)
  reps <- nrow(idx)

  X <- t(clr(counts + 1, 1))
  S <- cor(X)
  lambda_max <- pulsar::getMaxCov(S)
  lams <- pulsar::getLamPath(
    lambda_max, lambda_max * params$lambda_min_ratio,
    as.integer(params$nlambda), log = isTRUE(params$lambda_log)
  )

  hf <- huge::huge(X, method = "mb", lambda = lams, verbose = FALSE)

  merge <- lapply(seq_along(lams), function(i)
    Matrix::Matrix(0, p, p, sparse = TRUE))
  for (k in seq_len(reps)) {
    hk <- huge::huge(X[idx[k, ], , drop = FALSE], method = "mb",
                     lambda = lams, verbose = FALSE)
    for (i in seq_along(lams)) merge[[i]] <- merge[[i]] + hk$path[[i]]
  }

  summary_raw <- vapply(seq_along(lams), function(i) {
    m <- as.matrix(merge[[i]]) / reps
    4 * sum(m * (1 - m)) / (p * (p - 1))
  }, numeric(1))
  mono <- cummax(summary_raw)
  thresh <- params$thresh
  opt_index <- if (any(mono >= thresh)) {
    max(which.max(mono >= thresh)[1] - 1, 1)
  } else {
    length(mono)
  }

  refit <- as.matrix(hf$path[[opt_index]])
  ei <- which(refit != 0 & upper.tri(refit), arr.ind = TRUE)
  edges <- if (nrow(ei)) {
    paste(sprintf("%d:%d", ei[, 1] - 1L, ei[, 2] - 1L), collapse = ";")
  } else {
    ""
  }

  mo <- as.matrix(merge[[opt_index]]) / reps
  mi <- which(upper.tri(mo), arr.ind = TRUE)
  merge_str <- if (nrow(mi)) {
    paste(sprintf("%d:%d:%.6g", mi[, 1] - 1L, mi[, 2] - 1L, mo[mi]),
          collapse = ";")
  } else {
    ""
  }

  values <- c(lambda_max, lams, refit)
  nonfinite <- any(!is.finite(values))

  data.frame(
    case = case,
    status = if (nonfinite) "nonfinite" else "finite",
    error_class = "",
    opt_index = as.integer(opt_index),
    n_edges = as.integer(nrow(ei)),
    edges = if (nonfinite) "" else edges,
    lambda_max = as.numeric(lambda_max),
    stars_margin = as.numeric(min(abs(mono - thresh))),
    merge_opt = merge_str,
    seconds = as.numeric(proc.time()[["elapsed"]] - start),
    stringsAsFactors = FALSE
  )
}

case_dirs <- sort(list.dirs(cases_dir, recursive = FALSE, full.names = TRUE))
rows <- vector("list", length(case_dirs))
for (i in seq_along(case_dirs)) {
  dir <- case_dirs[[i]]
  rows[[i]] <- tryCatch(
    classify_case(dir),
    error = function(e) {
      row <- empty_row(basename(dir))
      row$error_class <- class(e)[[1]]
      row
    }
  )
}
utils::write.csv(do.call(rbind, rows), status_path, row.names = FALSE, na = "")
