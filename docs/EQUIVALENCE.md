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

A golden case is **robust** when the selected index and the edge sets do not sit
on a soft-threshold boundary. Two committed cases are deliberately degenerate:

- `p08_n40_weak` — very weak signal; no lambda is stable, `opt_index` falls back
  to 1 and the selected graph is empty.
- `p30_n20_pgtn` — `p > n`, `rank(S) = 19 < 30`: the selected graph is empty
  (`opt_index = 1`) and the pointwise coefficients are non-unique (effect 2
  above). Its invariants are checked by `test_degenerate_case_invariants`.

Cases with an *interior* optimum in the committed set: `p20_n300`,
`p20_n300_zeros`, `p30_n200_ratcheck`, `p50_n200`.

## Observed divergences (measured, not assumed)

Two effects remain. Both were investigated to root cause with measurements; no
tolerance was widened to hide either.

### 1. Boundary flip at lambda index 0 (`lambda.max`)

At the very top of the path the port may report **one extra edge** that R does
not. Measured (`|grad|/lambda` = the KKT ratio that decides the edge):

| case | `lambda[0]` | `max|S|` − `lambda[0]` (R) | `max|S|` − `lambda[0]` (port) | edge | R | port |
|---|---|---|---|---|---|---|
| `p20_n300` | 0.40311323350155553 | −10 ulp | **+5 ulp** | (18, 19) | ratio 0.9999999999999987 → inactive | 1.0000000000000007 → active |
| `p30_n20_pgtn` | 0.60681232594045142 | −4 ulp | **+1 ulp** | (27, 28) | 0.9999999999999992 → inactive | 1.0000000000000002 → active |

Cause: `spiec.easi` builds the path from R's `stats::cor(X)` but the solver runs
on huge's *own* standardised correlation (`crossprod(xs)`, BLAS `dsyrk`). The
port reproduces both, but NumPy's `xs.T @ xs` (`dgemm`) rounds ~1e-16 differently
from R's `dsyrk`, so `max|S|` lands on opposite sides of `lambda[0]`. The
decision therefore hinges on ~3–10 ulp of a quantity whose relative difference
is 7e-16. Both answers are "correct"; the edge's coefficient is ~1e-16.
Tolerated by `BOUNDARY_ATOL = 1e-8` in the tests and reported, never hidden.

### 2. Non-unique Lasso on the singular (`p > n`) design — `p30_n20_pgtn`

`n = 20`, `p = 30`, `rank(S) = 19`. The pointwise MB coefficients differ; the
*results* do not. Everything below is measured on that case (`a`–`g` of the
investigation).

**(a) Where it diverges.** 36 `(node, lambda)` pairs, **all in nodes 27 and 28**,
all at `lambda >= 1` (never at `lambda.max`). In all 36 the **support is
identical**; only the coefficient values differ (up to `1.3e-3`).

**(b) Objective.** Evaluating huge's objective
`0.5 w'Sw − S[m,:]'w + lambda||w||₁` on both solutions gives a maximum relative
difference of **7.2e-8**, and the objective is **flat along the entire segment**
`w_R + t(w_P − w_R)`, `t ∈ [0,1]` (max relative deviation 7.2e-8). Both are
minimisers — this is not a convergence gap.

**(c) KKT.** Inactive-coordinate violation is **exactly 0** for both solutions.
Active-coordinate violation is `<= 4.2e-5` for both — inside the solver's own
`thol = 1e-4`. Neither solution is "more converged" than the other.

**(d) Null space.** For the 7 pairs with `|Δβ| >= 1e-4`, the union active set has
`|A| = 20…22` while `rank(S_AA) = 19 < |A|`, and **100% of `β_R − β_Py` lies in
`null(S_AA)`** (null-space projection ratio = 1.00). Smaller differences
(`< 1e-5`, where `rank(S_AA) = |A|`) are pure solver-tolerance residue.

**(e) Solver sensitivity — decisive.** R itself, with only the **predictor order
permuted** (same data, same lambda path, same `huge` binary), moves by
`max|Δβ| = 1.3e-2` and changes the support on **5/20 lambdas**. The port is
likewise order/start/`thol`-sensitive on this case (permuting the order:
`7.4e-3`; cold start vs warm start: `4.5e-1`; `thol = 1e-8` vs `1e-4`:
`8.9e-2`). Since R is unstable *against itself*, the non-uniqueness is a property
of the problem, not of the port.

**(f) Impact on the final output.** For `p30_n20_pgtn` the StARS curve matches to
`2.2e-16` and `opt_index` is identical (1); the **refit differs by exactly the
one boundary edge of effect 1** (R: empty graph, port: `(27, 28)`), because
`opt_index = 1` selects the `lambda.max` graph. That difference is effect 1, not
effect 2. For a case with an interior optimum (`p20_n300`, `opt_index = 2`) the
refit is **identical**. The non-unique coefficients therefore do not change any
selected edge — they only make the pointwise values irreproducible.

**(g)** is effect 1 above.

**Conclusion: mathematical limit, not a bug.** The invariants that matter are
tested (`test_degenerate_case_invariants`: same support, objective agreement,
KKT within `thol`, `S_AA·Δβ` at residual level). Pointwise coefficient equality
is `xfail` for this one case, with the reason citing these measurements.



## How failures are reported

`tests/test_equivalence.py` compares step by step and, on the first mismatch,
reports the step name plus, for solver steps, the replica and lambda index at
which the first divergence occurs. Tolerances are **never** widened to force a
pass; a divergence is investigated and documented here.
