# Equivalence with the R reference

This document states precisely what the Python port guarantees to be identical
to the R SpiecEasi `mb` pipeline, what is identical only within a numerical
tolerance, and what cannot be identical.

Reference environment: R 4.6.1, SpiecEasi 1.99.0 (source in `reference/R/`),
huge 2.0.1, pulsar 0.3.13, Matrix 1.7-5. Goldens committed in `tests/golden/`.

## Platform scope (read this first)

The equivalence results in this document were **measured on Linux x86-64
(GCC) only**. The test suite is written to run everywhere without R (the goldens
are committed), but:

- **Windows and macOS have not been tested** — not locally and not in CI. The
  CI workflow (`.github/workflows/ci.yml`) exists but has **never been run**.
- On POSIX the optional StARS process pool uses the `fork` start method; on
  Windows (and any frozen/PyInstaller build) it falls back to `spawn`, which
  changes process start-up (not the numerical result, which is aggregated in a
  fixed order). That path is covered by code inspection only.
- Different compilers, `libm` implementations and BLAS builds can change the
  last bits of floating-point results, which matters exactly at the
  `lambda.max` boundary (see "Observed divergences").

Treat every "identical" claim below as **Linux x86-64 only** until the Windows
and macOS jobs are actually executed and pass.


## Guaranteed identical (exact / bit-level on the same BLAS)

These are deterministic steps with no iterative solver and no RNG:

| step | R | tolerance in tests |
|---|---|---|
| CLR transform `t(clr(counts+1, 1))` | `normalization.R:95`, `spiec-easi.R:67` | `rtol=1e-12, atol=1e-13` |
| correlation matrix `cor(X)` | `stats::cor` (ddof = n−1) | `rtol=1e-12, atol=1e-13` |
| `lambda.max = getMaxCov(cor(X))` | `pulsar/utilities.R:130` | `rtol=1e-12` |
| lambda grid `getLamPath(...)` | `pulsar/utilities.R:111` | `rtol=1e-12` |

## Identical within solver tolerance

The MB solver and the StARS aggregation are iterative / floating-point; the
alignment follows the solver's own stop rule.

| step | tolerance | justification |
|---|---|---|
| MB coefficients (`beta`) per replica / lambda | `rtol=1e-6, atol=1e-8` | huge MB coordinate descent stops at relative change `thol=1e-4` |
| MB edge sets per replica / lambda | **exact** | adjacency is `|beta| != 0`; robust cases only (see below) |
| StARS variability and monotone envelope | `rtol=1e-6, atol=1e-8` | derived from the edge sets |
| selected index `opt_index` | **exact** | integer |
| refit adjacency (full data) | **exact** (edge set) | same solver tolerance |

## Cannot be identical (documented, not hidden)

- **RNG for StARS subsamples.** `pulsar` draws `sample(1:n, floor(n*ratio))`
  with R's RNG. NumPy's generator differs, so with the same `seed` the
  subsample index sets differ. The Python port therefore accepts
  `subsample_indices`; the golden tests always inject the R indices. Without
  injection results are statistically equivalent but not identical.
- **BLAS/LAPACK + coordinate-descent arithmetic.** Differences of order
  1e-12…1e-15 are expected; near a soft-threshold tie the two implementations
  may select an edge differently (the "degenerate" case family).
- **Eigen-decompositions** (used by the R data generators, not by the `mb`
  pipeline) are defined up to sign/basis and R's `eigen()` orders eigenvalues
  descending while `numpy.linalg.eigh` orders ascending. Only invariant
  quantities are compared.

## Robust vs degenerate golden cases

A golden case is **robust** when the selected index and the edge sets are not
sitting on a soft-threshold boundary. The generator marks the deliberately weak
case `p08_n40_weak` as degenerate (`opt_index` falls back to 1, empty graph).
Boundary cases (`opt_index == 1` with an empty selected graph) are checked for
exact emptiness rather than for interior structure.

Cases with an *interior* optimum in the committed set: `p20_n300`,
`p20_n300_zeros`, `p30_n200_ratcheck`, `p50_n200`.

## Observed divergences (measured, not assumed)

1. **Boundary flip at lambda index 0 (`lambda.max`).** In the full-data MB
   solve, the port may report one extra edge at `lambda.max` (e.g. an
   `(18, 19)` edge in `p20_n300`). Cause: R computes the solver's correlation
   matrix with `crossprod(xs)` (BLAS `dsyrk`) while NumPy uses `xs.T @ xs`
   (`dgemm`); the two differ by ~4e-16, so `max|S|` sits ~1 ulp above/below the
   `lambda.max` that the R path used. The coefficient there is at the
   machine-epsilon level, so both answers are "correct". These edges are
   tolerated by `BOUNDARY_ATOL = 1e-8` in the tests and reported, never hidden.
2. **Singular design (`p30_n20_pgtn`, p>n).** `rank(S) = 19 < 30`, so the
   nodewise Lasso solution is not unique; R and the port converge to different
   minimisers and coefficients differ by up to ~1e-3 (all lambdas). The
   **support** still matches except for the lambda-0 boundary. Coefficient
   equality is `xfail` for this case (`DEGENERATE_COEF_CASES`).

No tolerance was widened to hide either effect: divergences beyond the boundary
tolerance fail the test with the replica/lambda and the coefficient magnitudes.


## How failures are reported

`tests/test_equivalence.py` compares step by step and, on the first mismatch,
reports the step name plus, for solver steps, the replica and lambda index at
which the first divergence occurs. Tolerances are **never** widened to force a
pass; a divergence is investigated and documented here.
