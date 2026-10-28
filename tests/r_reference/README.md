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
- `QMP_RDA` — path to the QMP `.rda` file (default `../SPRING-master/data/QMP.rda`).

## Environment used to build the committed goldens

`sessionInfo()` for every case is stored next to the fixture in
`tests/golden/<case>/sessionInfo.txt`. The versions used were:

| package | version |
|---------|---------|
| R | 4.6.1 |
| SpiecEasi | 1.99.0 (reference source plus installed data package) |
| phyloseq | 1.56.0 |
| huge | 2.0.1 |
| pulsar | 0.3.13 |
| Matrix | 1.7-5 |
| jsonlite | (any recent) |

Other packages loaded by SpiecEasi (`MASS`, `VGAM`, `Rcpp`, `RcppArmadillo`)
are needed to load the reference source but are not on the `mb` code path.

## What each case contains

- `counts.tsv` — input count matrix (samples × features), tab-separated, no header.
- `params.json` — the exact parameters (dataset, method, nlambda, lambda_min_ratio,
  lambda_log, rep_num, thresh, seed, n, p).
- `clr.tsv` — the transformed matrix `t(clr(counts + 1, 1))`.
- `cor.tsv` — `cor(clr)`.
- `lambda_max.txt` — `pulsar::getMaxCov(cor(clr))`.
- `lambda_path.txt` — `pulsar::getLamPath` using the recorded `lambda_min_ratio`
  and `nlambda` (17 significant digits).
- `full_beta.tsv` — full-data MB coefficients, long form `lambda_idx, response, predictor, value`.
- `full_adj.tsv` — full-data MB adjacency (upper triangle), `lambda_idx, i, j`.
- `subsample_indices.tsv` — the **R RNG** subsample indices actually used by
  pulsar (row = replica, 0-based sample indices). Injected into the Python port.
- `rep_adj_path.tsv` — per-replica, per-lambda MB adjacency: `rep, lambda_idx, i, j`.
- `rep_beta_long.tsv` — per-replica, per-lambda MB nonzero coefficients, used to
  diagnose whether a support flip is a true near-threshold case.
- `stars_summary_raw.txt` — per-lambda StARS variability before monotonisation.
- `stars_summary_monotone.txt` — after `cummax`.
- `opt_index.txt` — selected index from the manual reconstruction (1-based).
- `merge_opt.tsv` — edge-frequency merge at the selected lambda, `i, j, value`.
- `pipeline_opt_index.txt`, `pipeline_opt_lambda.txt`, `pipeline_refit.tsv` —
  the end-to-end `spiec.easi()` result (cross-check of the manual path).
- `sessionInfo.txt` — the R session used to produce the case.

## Determinism notes

- The inputs include QMP, American Gut Round 1, and GlobalPatterns. The shared
  preprocessing selector is `../scripts/public_microbiome_cases.R`; selected
  taxa/sample IDs and source profile are recorded per case.
- The StARS subsample indices come from R's RNG, which does **not** match
  NumPy's. That is exactly why they are dumped: the Python port reuses them, so
  the R RNG is never needed on the Python side.
- All numeric dumps use `%.17g`.

## Cases

| case | dimensions | selection |
|---|---:|---|
| `qmp_n106_p20` | 106 × 20 | all samples, 20 least-prevalent QMP taxa |
| `qmp_n20_p30` | 20 × 30 | real p > n subset, 30 taxa across prevalence ranks |
| `amgut_r1_n80_p20_seed3101..3103` | 80 × 20 | American Gut Round 1; 3 R StARS seeds |
| `amgut_r1_n20_p30_seed3101..3103` | 20 × 30 | American Gut Round 1, p > n; 3 R StARS seeds |
| `globalpatterns_n26_p15_seed3101..3103` | 26 × 15 | GlobalPatterns most-prevalent varying taxa; 3 R StARS seeds |
| `globalpatterns_n26_p40_seed3101..3103` | 26 × 40 | GlobalPatterns, p > n; 3 R StARS seeds |

All case matrices use the same preprocessing function at
`../scripts/public_microbiome_cases.R`; `taxa.txt`, `sample_ids.txt`, and
`params.json` preserve selected identifiers and source/profile metadata. The
American Gut and GlobalPatterns cases extend QMP to two independent public
study sources. `SpiecEasi` 1.99.0 and `phyloseq` 1.56.0 are needed by the R
generator for dataset extraction; they are not Python runtime dependencies.

## Seeded failure-status stress sweep

Run the registered 1,000-case SPIEC-EASI sweep (small count dimensions in six
predeclared strata: moderate sparsity, ~65% zeros, rank-deficient/constant
feature, low library sizes, `p > n`, and threshold-near weak signal) with:

```bash
.venv/bin/python tests/r_reference/run_status_stress.py \
    /tmp/spiec-stress-1000 --cases 1000 --seed 20261020
```

The runner draws one seed per case and one StARS-index draw seed per case from
a NumPy master stream, archives every case (`counts.csv`,
`subsample_indices.csv`, `params.json`), classifies each case in Python, then
runs the pinned R source (`stress_status_R.R`) on the same counts and injected
subsample indices. Both sides use `nlambda=5`, `lambda_min_ratio=0.1`,
`lambda_log=TRUE`, `rep_num=4`, `thresh=0.2`.

It writes `python_status.csv`, `r_status.csv`, `summary.json` and
`environment.json`. A case is an unexplained failure when the
`finite`/`nonfinite`/`error` status differs, or when the selected index or
refit edge set differs without meeting the predeclared boundary rule
(`|stability - thresh| <= 1e-3` for index ties; every differing edge within
0.05 of a 0.5 selection frequency). Any unexplained difference exits nonzero.
The complete case archive is retained at the output path; the summary records
per-case Python/R timings so a full sweep can be re-timed on another machine.

Recorded run: `../stress/2026-10-20-seed20261020/` (1,000 cases, seed
20261020, 135 s wall clock; 0 status mismatches, 21/1,000 `opt_index`
mismatches, 256/1,000 single-edge refit differences, all triaged to the
`lambda.max` tie documented in `docs/EQUIVALENCE.md`).
