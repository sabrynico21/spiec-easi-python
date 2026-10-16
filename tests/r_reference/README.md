# R reference and golden fixtures

Everything in `tests/r_reference/` and `tests/golden/` exists so that the
Python port can be validated against the **R SpiecEasi reference** without R
being installed at test time. The goldens are committed; the Python tests load
them directly.

## Layout

- `generate_golden.R` — the generator. Loads the R source in `reference/R/`
  with `pkgload::load_all()` and, for every case, writes inputs, parameters,
  all outputs and the relevant intermediates to `tests/golden/<case>/`.
- `../golden/` — the committed fixtures (regenerate deterministically, see
  below). A `MANIFEST.sha256` records the SHA-256 of every file.

## Regenerating

From the repository root (`spiec-easi-python/`):

```bash
Rscript tests/r_reference/generate_golden.R
```

Environment variables (optional):

- `SPIEC_EASI_REF_R` — path to the R source (default `reference/R`).
- `GOLDEN_OUT` — output directory (default `tests/golden`).

## Environment used to build the committed goldens

`sessionInfo()` for every case is stored next to the fixture in
`tests/golden/<case>/sessionInfo.txt`. The versions used were:

| package | version |
|---------|---------|
| R | 4.6.1 |
| SpiecEasi | 1.99.0 (from `reference/R/`) |
| huge | 2.0.1 |
| pulsar | 0.3.13 |
| Matrix | 1.7-5 |
| jsonlite | (any recent) |

Other packages loaded by SpiecEasi (`MASS`, `VGAM`, `Rcpp`, `RcppArmadillo`)
are needed to load the reference source but are not on the `mb` code path.

## What each case contains

- `counts.tsv` — input count matrix (samples × features), tab-separated, no header.
- `params.json` — the exact parameters (method, nlambda, lambda_min_ratio,
  lambda_log, rep_num, thresh, seed, n, p).
- `clr.tsv` — the transformed matrix `t(clr(counts + 1, 1))`.
- `cor.tsv` — `cor(clr)`.
- `lambda_max.txt` — `pulsar::getMaxCov(cor(clr))`.
- `lambda_path.txt` — `pulsar::getLamPath(lambda_max, lambda_max*1e-3, 20, log=TRUE)` (17 significant digits).
- `full_beta.tsv` — full-data MB coefficients, long form `lambda_idx, response, predictor, value`.
- `full_adj.tsv` — full-data MB adjacency (upper triangle), `lambda_idx, i, j`.
- `subsample_indices.tsv` — the **R RNG** subsample indices actually used by
  pulsar (row = replica, 0-based sample indices). Injected into the Python port.
- `rep_adj_path.tsv` — per-replica, per-lambda MB adjacency: `rep, lambda_idx, i, j`.
- `stars_summary_raw.txt` — per-lambda StARS variability before monotonisation.
- `stars_summary_monotone.txt` — after `cummax`.
- `opt_index.txt` — selected index from the manual reconstruction (1-based).
- `merge_opt.tsv` — edge-frequency merge at the selected lambda, `i, j, value`.
- `pipeline_opt_index.txt`, `pipeline_opt_lambda.txt`, `pipeline_refit.tsv` —
  the end-to-end `spiec.easi()` result (cross-check of the manual path).
- `sessionInfo.txt` — the R session used to produce the case.

## Determinism notes

- The synthetic counts are generated deterministically from a seed (a sparse,
  two-tier Gaussian copula → Poisson-lognormal counts), so a rerun on the same
  R version reproduces them byte-for-byte.
- The StARS subsample indices come from R's RNG, which does **not** match
  NumPy's. That is exactly why they are dumped: the Python port reuses them, so
  the R RNG is never needed on the Python side.
- All numeric dumps use `%.17g`.
