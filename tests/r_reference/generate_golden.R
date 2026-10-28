#!/usr/bin/env Rscript
# =============================================================================
# generate_golden.R - R reference golden generator for the Python SpiecEasi (MB)
#
# Runs the *reference R source* (reference/R via pkgload::load_all) on QMP and
# independent public microbiome cases and dumps, for every case:
#   - the input counts and the exact parameters,
#   - every deterministic intermediate (clr matrix, correlation matrix,
#     lambda.max, lambda path),
#   - the full-data MB solver output (beta coefficients and the symmetrised
#     adjacency per lambda, on the full data),
#   - the StARS subsample indices actually used by pulsar (so the Python port
#     never needs the R RNG),
#   - the per-replica MB adjacency paths, the per-lambda StARS variability and
#     its monotone envelope, and the selected index,
#   - the end-to-end spiec.easi() result (opt index/lambda, refit adjacency),
#   - sessionInfo().
#
# Scope: method = "mb" only (see docs/PORTING_NOTES.md, decision D2).
#
# Usage:
#   Rscript tests/r_reference/generate_golden.R
# Environment variables:
#   SPIEC_EASI_REF_R  path to the reference R source (default: reference/R)
#   GOLDEN_OUT        output directory               (default: tests/golden)
# =============================================================================

suppressWarnings(suppressMessages({
  library(jsonlite)
}))

REF_R <- Sys.getenv("SPIEC_EASI_REF_R", "reference/R")
OUT   <- Sys.getenv("GOLDEN_OUT", "tests/golden")

# Load dataset objects from the installed R data packages before pkgload::load_all
# shadows the SpiecEasi namespace with the reference source tree.
PUBLIC_HELPER <- file.path("..", "scripts", "public_microbiome_cases.R")
if (!file.exists(PUBLIC_HELPER)) stop("shared public-data selector not found: ", PUBLIC_HELPER)
source(PUBLIC_HELPER)
PUBLIC_CASES <- load_public_microbiome_cases()

if (!dir.exists(REF_R)) stop("reference R source not found: ", REF_R)
suppressWarnings(suppressMessages(pkgload::load_all(REF_R, quiet = TRUE)))
if (!requireNamespace("huge", quietly = TRUE))
  stop("the 'huge' package is required to generate goldens")

message("Loaded SpiecEasi ", as.character(utils::packageVersion("SpiecEasi")),
        " from ", REF_R)

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- small writers ---------------------------------------------------------

write_mat <- function(m, path) {
  write.table(m, path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

write_vec <- function(v, path) {
  writeLines(sprintf("%.17g", as.numeric(v)), path)
}

fmt <- function(x) sprintf("%.17g", as.numeric(x))

write_long <- function(df, path, cols) {
  if (is.null(df) || nrow(df) == 0) {
    df <- as.data.frame(stats::setNames(replicate(length(cols), character(0),
                                                  simplify = FALSE), cols))
  }
  write.table(df, path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

empty_df <- function(cols) {
  as.data.frame(stats::setNames(replicate(length(cols), integer(0),
                                          simplify = FALSE), cols))
}

# ---- deterministic synthetic counts ---------------------------------------

# Build a sparse symmetric weighted adjacency with two tiers of edge strength:
# a few STRONG edges (so lambda.max is dominated by them) plus a WEAK background
# (a nearest-neighbour chain and some extra edges). The strong/weak gap creates
# an "empty plateau" at the top of the lambda path, which is what gives StARS an
# interior optimum instead of the empty/complete extremes.
sparse_adj <- function(p, seed, weak = 0.4, strong = 0.9,
                       extra_frac = 1/3, strong_frac = 0.25) {
  set.seed(seed + 9000L)
  W <- matrix(0, p, p)
  if (p >= 2) {
    for (i in seq_len(p - 1L)) { W[i, i + 1L] <- weak; W[i + 1L, i] <- weak }
  }
  cand <- which(upper.tri(W) & W == 0, arr.ind = TRUE)
  nadd <- min(max(1L, floor(p * extra_frac)), nrow(cand))
  if (nadd > 0 && nrow(cand) > 0) {
    pick <- cand[sample(nrow(cand), nadd), , drop = FALSE]
    W[pick] <- weak
    W[pick[, 2:1, drop = FALSE]] <- weak
  }
  edges <- which(upper.tri(W) & W != 0, arr.ind = TRUE)
  ns <- min(max(1L, floor(nrow(edges) * strong_frac)), nrow(edges))
  if (ns > 0 && nrow(edges) > 0) {
    sp <- edges[sample(nrow(edges), ns), , drop = FALSE]
    W[sp] <- strong
    W[sp[, 2:1, drop = FALSE]] <- strong
  }
  W
}

# Poisson-lognormal counts from a sparse, two-tier Gaussian copula.
gen_counts <- function(n, p, seed, mu = 10, weak = 0.4, strong = 0.9,
                       zero_frac = 0) {
  W <- sparse_adj(p, seed, weak = weak, strong = strong)
  # Strictly diagonally dominant => guaranteed positive definite.
  Theta <- diag(1 + rowSums(W)) + W
  Sigma <- stats::cov2cor(solve(Theta))
  eta <- rmvnorm(n, rep(0, p), Sigma)                 # n x p Gaussian copula
  lam <- exp(log(mu) + eta)
  set.seed(seed)
  C <- matrix(rpois(n * p, lambda = as.vector(lam)), n, p)
  if (zero_frac > 0) {
    thr <- stats::quantile(C, zero_frac)
    C[C <= thr] <- 0L
  }
  storage.mode(C) <- "integer"
  dimnames(C) <- NULL
  C
}

# ---- dumpers ---------------------------------------------------------------

# beta[[i]] is d x d with beta[j, m] = coef of predictor j in the regression
# of response m (huge.mb). path[[i]] is the symmetrised 0/1 adjacency.
dump_huge_path <- function(h, d, dir, prefix) {
  bl <- list()
  for (i in seq_along(h$beta)) {
    b  <- as.matrix(h$beta[[i]])
    nz <- which(b != 0, arr.ind = TRUE)
    if (nrow(nz)) {
      bl[[length(bl) + 1]] <- data.frame(
        lambda_idx = i - 1L,
        response   = nz[, 2] - 1L,
        predictor  = nz[, 1] - 1L,
        value      = fmt(b[nz]),
        stringsAsFactors = FALSE)
    }
  }
  bl <- if (length(bl)) do.call(rbind, bl)
        else as.data.frame(stats::setNames(
          replicate(4, character(0), simplify = FALSE),
          c("lambda_idx", "response", "predictor", "value")))
  write_long(bl, file.path(dir, paste0(prefix, "_beta.tsv")),
             c("lambda_idx", "response", "predictor", "value"))

  al <- list()
  for (i in seq_along(h$path)) {
    a   <- as.matrix(h$path[[i]])
    idx <- which(a != 0 & upper.tri(a), arr.ind = TRUE)
    if (nrow(idx)) {
      al[[length(al) + 1]] <- data.frame(
        lambda_idx = i - 1L, i = idx[, 1] - 1L, j = idx[, 2] - 1L)
    }
  }
  al <- if (length(al)) do.call(rbind, al)
        else empty_df(c("lambda_idx", "i", "j"))
  write_long(al, file.path(dir, paste0(prefix, "_adj.tsv")),
             c("lambda_idx", "i", "j"))
}

dump_case <- function(name, counts, params) {
  d <- file.path(OUT, name)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  n <- nrow(counts)
  p <- ncol(counts)

  write_mat(counts, file.path(d, "counts.tsv"))
  writeLines(colnames(counts), file.path(d, "taxa.txt"))
  writeLines(rownames(counts), file.path(d, "sample_ids.txt"))
  write(toJSON(params, auto_unbox = TRUE, digits = NA),
        file.path(d, "params.json"))

  # 1. data transformation --------------------------------------------------
  X <- t(clr(counts + 1, 1))
  write_mat(X, file.path(d, "clr.tsv"))

  # 2. correlation + lambda grid -------------------------------------------
  S <- cor(X)
  write_mat(S, file.path(d, "cor.tsv"))
  lambda_max <- pulsar::getMaxCov(S)
  write_vec(lambda_max, file.path(d, "lambda_max.txt"))
  lams <- pulsar::getLamPath(lambda_max,
                             lambda_max * params$lambda_min_ratio,
                             params$nlambda, log = params$lambda_log)
  write_vec(lams, file.path(d, "lambda_path.txt"))

  # 3. full-data MB solver path --------------------------------------------
  hf <- huge::huge(X, method = "mb", lambda = lams, verbose = FALSE)
  dump_huge_path(hf, p, d, "full")

  # 4. StARS subsample indices (R RNG, captured for injection) --------------
  ratcheck <- getFromNamespace(".ratcheck", "pulsar")
  ratio    <- ratcheck(NULL, n)
  reps     <- as.integer(params$rep_num)
  set.seed(params$seed)
  idx <- replicate(reps, sample(c(1:n), floor(n * ratio), replace = FALSE),
                   simplify = FALSE)
  write_mat(do.call(rbind, idx), file.path(d, "subsample_indices.tsv"))

  # 5. per-replica MB paths + StARS variability -----------------------------
  merge <- lapply(seq_along(lams),
                  function(i) Matrix::Matrix(0, p, p, sparse = TRUE))
  rep_rows <- list()
  rep_beta_rows <- list()
  for (k in seq_len(reps)) {
    hk <- huge::huge(X[idx[[k]], , drop = FALSE], method = "mb",
                     lambda = lams, verbose = FALSE)
    for (i in seq_along(lams)) {
      merge[[i]] <- merge[[i]] + hk$path[[i]]
      a  <- as.matrix(hk$path[[i]])
      ui <- which(a != 0 & upper.tri(a), arr.ind = TRUE)
      if (nrow(ui)) {
        rep_rows[[length(rep_rows) + 1]] <- data.frame(
          rep = k - 1L, lambda_idx = i - 1L,
          i = ui[, 1] - 1L, j = ui[, 2] - 1L)
      }
      b <- as.matrix(hk$beta[[i]])
      bi <- which(b != 0, arr.ind = TRUE)
      if (nrow(bi)) {
        rep_beta_rows[[length(rep_beta_rows) + 1L]] <- data.frame(
          rep = k - 1L, lambda_idx = i - 1L,
          response = bi[, 2] - 1L, predictor = bi[, 1] - 1L,
          value = fmt(b[bi])
        )
      }
    }
  }
  rep_rows <- if (length(rep_rows)) do.call(rbind, rep_rows)
              else empty_df(c("rep", "lambda_idx", "i", "j"))
  write_long(rep_rows, file.path(d, "rep_adj_path.tsv"),
             c("rep", "lambda_idx", "i", "j"))
  rep_beta_rows <- if (length(rep_beta_rows)) do.call(rbind, rep_beta_rows)
                   else empty_df(c("rep", "lambda_idx", "response", "predictor", "value"))
  write_long(rep_beta_rows, file.path(d, "rep_beta_long.tsv"),
             c("rep", "lambda_idx", "response", "predictor", "value"))

  summary_raw <- vapply(seq_along(lams), function(i) {
    m <- as.matrix(merge[[i]]) / reps
    4 * sum(m * (1 - m)) / (p * (p - 1))
  }, numeric(1))
  mono <- cummax(summary_raw)
  opt_index <- if (any(mono >= params$thresh))
                 max(which.max(mono >= params$thresh)[1] - 1, 1)
               else length(mono)

  write_vec(summary_raw, file.path(d, "stars_summary_raw.txt"))
  write_vec(mono, file.path(d, "stars_summary_monotone.txt"))
  write_vec(opt_index, file.path(d, "opt_index.txt"))

  # edge-frequency merge at the selected lambda
  mo <- as.matrix(merge[[opt_index]]) / reps
  mi <- which(upper.tri(mo), arr.ind = TRUE)
  write_long(data.frame(i = mi[, 1] - 1L, j = mi[, 2] - 1L,
                        value = fmt(mo[mi])),
             file.path(d, "merge_opt.tsv"), c("i", "j", "value"))

  # 6. end-to-end spiec.easi (cross-check) ----------------------------------
  est <- spiec.easi(counts, method = "mb", sel.criterion = "stars",
                    verbose = FALSE, nlambda = params$nlambda,
                    lambda.min.ratio = params$lambda_min_ratio,
                    lambda.log = params$lambda_log,
                    pulsar.params = list(seed = params$seed,
                                         rep.num = reps,
                                         thresh = params$thresh))
  write_vec(getOptInd(est), file.path(d, "pipeline_opt_index.txt"))
  write_vec(getOptLambda(est), file.path(d, "pipeline_opt_lambda.txt"))
  write_mat(as.matrix(getRefit(est)), file.path(d, "pipeline_refit.tsv"))

  # consistency guard: manual stability path must match pulsar's
  if (getOptInd(est) != opt_index)
    warning(sprintf("case %s: manual opt_index (%d) != pulsar getOptInd (%d)",
                    name, opt_index, getOptInd(est)))

  writeLines(capture.output(sessionInfo()),
             file.path(d, "sessionInfo.txt"))
  message("  case ", name, ": n=", n, " p=", p, " ratio=", round(ratio, 6),
          " opt_index=", opt_index)
}

# ---- QMP real-data cases ---------------------------------------------------

qmp_rda <- Sys.getenv("QMP_RDA", "../SPRING-master/data/QMP.rda")
if (!file.exists(qmp_rda)) stop("QMP data not found: ", qmp_rda)
load(qmp_rda)
qmp <- as.matrix(QMP)
prevalence_order <- order(colSums(qmp > 0), colnames(qmp))
qmp_full <- qmp[, prevalence_order[seq_len(20)], drop = FALSE]

small_taxa <- prevalence_order[
  unique(round(seq(1, length(prevalence_order), length.out = 30)))
]
qmp_wide <- qmp[, small_taxa, drop = FALSE]
covered <- rep(FALSE, ncol(qmp_wide))
sample_rows <- integer()
while (!all(covered)) {
  candidates <- setdiff(seq_len(nrow(qmp_wide)), sample_rows)
  gains <- vapply(candidates, function(i)
    sum(qmp_wide[i, ] > 0 & !covered), integer(1))
  next_row <- candidates[which.max(gains)]
  sample_rows <- c(sample_rows, next_row)
  covered <- covered | (qmp_wide[next_row, ] > 0)
}
if (length(sample_rows) > 20L) stop("Could not construct a 20-sample QMP p>n case")
sample_rows <- c(sample_rows,
                 head(setdiff(seq_len(nrow(qmp_wide)), sample_rows),
                      20L - length(sample_rows)))
qmp_pgt_n <- qmp_wide[sample_rows, , drop = FALSE]

cases <- list(
  list(name = "qmp_n106_p20", counts = qmp_full, seed = 1101,
       dataset = "QMP", provenance = list()),
  list(name = "qmp_n20_p30", counts = qmp_pgt_n, seed = 1102,
       dataset = "QMP", provenance = list())
)

base_params <- list(method = "mb", sel_criterion = "stars",
                    nlambda = 5L, lambda_min_ratio = 0.1,
                    lambda_log = TRUE, rep_num = 4L, thresh = 0.2)

for (public_case in PUBLIC_CASES) {
  for (seed in c(3101L, 3102L, 3103L)) {
    provenance <- list(
      source = public_case$source,
      taxa_strategy = public_case$taxa_strategy,
      source_n = public_case$raw_n,
      source_p = public_case$raw_p,
      available_n = public_case$available_n,
      available_p = public_case$available_p,
      zero_fraction = public_case$zero_fraction,
      library_size_min = public_case$library_size_min,
      library_size_median = public_case$library_size_median,
      library_size_max = public_case$library_size_max
    )
    cases[[length(cases) + 1L]] <- list(
      name = sprintf("%s_seed%d", public_case$name, seed),
      counts = public_case$input,
      seed = seed,
      dataset = public_case$dataset,
      provenance = provenance
    )
  }
}

manifest <- list()
for (cs in cases) {
  counts <- cs$counts
  params <- c(base_params, cs$provenance,
              list(seed = cs$seed, n = nrow(counts), p = ncol(counts),
                   dataset = cs$dataset))
  dump_case(cs$name, counts, params)
  manifest[[cs$name]] <- params
}

write(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, digits = NA),
      file.path(OUT, "manifest.json"))
message("Wrote ", length(cases), " golden cases to ", OUT)
