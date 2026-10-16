# SpiecEasi → Python — Porting Notes (Passo 0)

Status: **analysis only, no implementation started.** This document is the
mandatory pre-code analysis required before writing any Python.

Sources of truth (read in full before writing this):

- R package `SpiecEasi` 1.99.0 — `/home/sabrynico21/Scaricati/SpiecEasi-master`,
  copied verbatim to `reference/R/`.
- External deps, source tarballs extracted to `reference/deps/`:
  `huge` 2.0.1 and `pulsar` 0.3.13 (kept). `glmnet` 5.1 and `Matrix` 1.7-5 were
  also downloaded during the audit and then removed once `mb` was fixed as the
  only scope (neither is on the `mb` code path); their provenance and checksums
  are recorded in `reference/deps/FETCH.md`.
- Installed (verified) R stack on this machine:
  R 4.6.1, SpiecEasi 1.99.0, huge 2.0.1, pulsar 0.3.13, glmnet 5.1,
  Matrix 1.7-5, MASS 7.3-65, VGAM 1.1-14, Rcpp 1.1.1-1.1,
  RcppArmadillo 15.2.6-1, boot 1.3-32, igraph 2.3.4, phyloseq 1.56.0.

> Line numbers below are `# R: <file>:<line>` references into `reference/R/`
> (for the SpiecEasi package) and `reference/deps/<pkg>/…` for the
> dependencies. They are the anchors required in the Python comments.

---

## 0. Executive summary of the pipeline

`spiec.easi(data, method, sel.criterion, ...)` does exactly this:

1. **Transform** — `X = t(clr(data + 1, mar=1))` (per-sample CLR with a
   SpiecEasi-specific zero/centering rule). `R: spiec-easi.R:67,95`.
2. **lambda grid** — `lambda.max = getMaxCov(cor(X))` (max |off-diagonal
   Pearson correlation|), then `getLamPath(lambda.max,
   lambda.max*lambda.min.ratio, nlambda, log=TRUE)` with `lambda.min.ratio=1e-3`,
   `nlambda=20`, descending, log-spaced. `R: spiec-easi.R:188-194`,
   `pulsar/R/utilities.R:111,130`.
3. **Per-lambda estimator** — `sparseiCov` → `huge::huge(x=X, method, lambda,
   cov.output=TRUE)`. Both `mb` and `glasso` run **on the correlation matrix**
   `S = .huge_fast_cor(X)` (unit-variance standardized columns), not on `X`.
   `R: SparseICov.R:44`, `huge/R/huge-internal.R:78,228`.
   - `mb`: nodewise Lasso coordinate-descent on `S`, returns `beta` per column
     per lambda; adjacency = `sign(|beta| + |beta|ᵀ)` (`sym="or"`).
     `deps/huge/src/huge_core.cpp:679`.
   - `glasso`: FHT-style alternating coordinate descent on `S`, returns
     `icov` (precision) and `path` = `(icov != 0)` with zero diagonal.
     `deps/huge/src/huge_core.cpp:51,471`.
4. **StARS** — `pulsar(data=X, fun=sparseiCov, fargs=list(method, lambda),
   criterion="stars", thresh=0.05, subsample.ratio=NULL, rep.num=20,
   seed=NULL, ncores=1, refit=TRUE)`. `R: spiec-easi.R:216-241`,
   `pulsar/R/mcPulsarSelect.R:168`.
5. **Refit** — `pulsar::refit` re-runs the estimator on the **full** X and
   picks `path[[opt.index]]`; the returned `pulsar.refit` keeps `est` (the full
   `huge` object), `select` (the pulsar/StARS object), `lambda`, `fun`.
   `R: spiec-easi.R:243-251`, `pulsar/R/refit.R:38`.
6. **Accessors** — `getOptInd/Lambda/Merge/Stability/Net/Beta/Cov/iCov`
   (`getOptX.pulsar.refit`). `R: utilities.R:60,66`.

**Crucial correction to common assumptions:** in the *current* sources
(SpiecEasi 1.99.0 + huge 2.0.1) the `mb` method does **NOT** call `glmnet`.
`glmnet` is used only by `neighborhood.net()` for the non-default `ising` /
`poisson` exponential-family models. The `mb` path in huge 2.0.1 is a
self-contained C++ coordinate-descent Lasso (`huge::mb`). Any earlier
description that says "MB via glmnet/huge.mb" is outdated.

---

## 1. Function-by-function porting map (core)

### 1.1 `clr.default` / `clr.matrix` — `R: normalization.R:95,105`

```
clr.default(x, base=exp(1), tol=.Machine$double.eps):
    nzero = (x >= tol)
    LOG   = log(ifelse(nzero, x, 1), base)
    ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)

clr.matrix(x, mar=2):  apply(x, mar, clr.default)
```

- `mean(LOG)` is over **all** entries (zero-substituted entries contribute
  `log(1)=0`); `mean(nzero)` is the **fraction of non-zero** entries. The
  centering constant is therefore `mean(LOG)/fraction_nonzero` — this is the
  SpiecEasi-specific rule, **not** the textbook `LOG - mean(LOG)`.
- `base = e` (natural log). `tol = .Machine$double.eps ≈ 2.220446e-16`.
- Output keeps zeros in place (`ifelse(nzero, …, 0.0)`).
- Python equivalent: vectorized with `np.where`; `mean` → `np.mean` (all
  entries), no `ddof` involved. Must reproduce the exact `>= tol` test and the
  `x>=tol` short-circuit.

### 1.2 `.spiec.easi.norm` — `R: spiec-easi.R:67`

```
if matrix:  t(clr(data + 1, 1))      # per-sample clr of pseudo-counted counts
if list:    do.call(cbind, lapply(data, .spiec.easi.norm))   # multi-domain
```

- `data + 1` adds **integer 1** to every cell.
- `clr(..., 1)` ⇒ `mar=1` ⇒ clr applied **row-wise** (each row = one sample),
  returning `p × n`, then transposed to `n × p`.
- Python: `X = clr(counts + 1, axis=1)` then transpose (or apply along axis=0).

### 1.3 `getMaxCov` — `pulsar/R/utilities.R:130`

```
getMaxCov(x, cov=isSymmetric(x), abs=TRUE, diag=FALSE):
    if (!cov) x <- cov(x)
    tmp <- Matrix::triu(x, k = if (diag) 0 else 1)   # strict upper triangle
    max(if (abs) abs(tmp) else tmp)
```

- For `getMaxCov(cor(X))`: `cor(X)` is symmetric ⇒ strict upper triangle,
  `max |r_ij|`, i.e. `max` over `i<j` of `|Pearson corr|`.
- Python: `S = np.corrcoef(X.T)` (ddof = n-1, matching R `cor`), take
  `max |triu(S, k=1)|`. **Do not** use `np.cov` scaled differently.

### 1.4 `getLamPath` — `pulsar/R/utilities.R:111`

```
getLamPath(max, min, len, log=FALSE):
    if (log) { min <- log(min); max <- log(max) }
    lams <- seq(max, min, length.out=len)
    if (log) exp(lams) else lams
```

- With `log=TRUE`: descending geometric path from `lambda.max` to
  `lambda.max*lambda.min.ratio`, `nlambda` points.
- Python: `np.geomspace(lambda_max, lambda_max*ratio, n)`. (Verify
  endpoint bit-parity: R computes `seq(log(max), log(min), length.out)` then
  `exp`; keep the same operation order to stay within 1 ulp.)

### 1.5 `.ratcheck` (subsample ratio) — `pulsar/R/mcPulsarSelect.R:15`

```
if (is.null(ratio)) { if (n > 144) 10*sqrt(n)/n else 0.8 }
subsample size = floor(n * ratio)
```

### 1.6 `pulsar` — `pulsar/R/mcPulsarSelect.R:168` (core loop)

Order of operations (must be preserved exactly):

1. `n=nrow(data)`, `p=ncol(data)`; `.lamcheck(fargs$lambda)` (warn if not
   non-increasing or length < 2).
2. `subsample.ratio <- .ratcheck(subsample.ratio, n)`.
3. **Subsample indices**: `replicate(rep.num, sample(c(1:n),
   floor(n*subsample.ratio), replace=FALSE))`. `seed` set/cleared around this.
   ⇒ **R RNG; not reproducible in numpy** (Decision D1).
4. If `refit`: prepend the full index vector `1:n` (tagged `attr 'full'`).
5. `isamp <- ind.sample[1:minN]`, `minN = rep.num + refit` (or `2+refit` if
   `lb.stars`).
6. `estFun` per index: `do.call(fun, c(fargs, list(data[ind.sample,])))`; if
   `full` return whole object, else return `tmp$path`.
7. Reorder to `premerge.reord[[lambda]][[rep]]` via `.tlist`.
8. `stars.stability` → `est$stars`.
9. `opt.index <- est$stars$opt.index`.
10. If `opt.index == 1`: emit warning "Optimal lambda may be larger/smaller…".
11. `est$call`, `est$envir`, `est$est = fullmodel`.

### 1.7 `stars.stability` + `.starsind` — `pulsar/R/mcPulsarSelect.R:339,329`

```
merge[[l]]  <- Reduce("+", premerge[[l]])            # sum over reps
merge[[l]]  <- merge[[l]]/rep.num                    # edge frequency
summary[l]  <- 4*sum(merge[[l]]*(1-merge[[l]])) / (p*(p-1))
summary     <- cummax(summary)                       # monotonize
opt.index   <- .starsind(summary, thresh, offset=1)
.starsind: if any(summary >= thresh) max(which.max(summary>=thresh)-1, 1)
           else warning + length(summary)
```

- Dominant eigenvalue of the merge summand `m*(1-m)`; note the `/ (p*(p-1))`
  normaliser (full ordered pairs, includes both triangles).
- `opt.index = first index with monotone stability >= thresh, minus 1`.
- `thresh` default set to **0.05** by `spiec.easi` (`R: spiec-easi.R:230`);
  pulsar's own default is 0.1 but is overridden.

### 1.8 `huge.mb` — `huge/R/huge.mb.R:27` + `huge_core.cpp:679`

- Input `S = .huge_preprocess(x)$S = .huge_fast_cor(x)` (correlation matrix).
- `.huge_fast_cor` (`huge-internal.R:78`): standardize each column with
  `.huge_standardize` (center + unit variance, numerically careful), then
  `S = crossprod(xs)/(n-1)`, clamp to [-1,1], `diag(S)=1`.
- `.huge_default_lambda` is bypassed because `lambda` is supplied.
- Solver `huge::mb`: **nodewise Lasso** on the correlation matrix:
  - `thol = 1e-4`, `MAX_ITER = 10000` (both the outer and inner loops).
  - Sequential strong rule (`strong_thr = 2λ_i − λ_{i−1}`), active-set
    coordinate descent with `threshold_l1`, relative coefficient-change gap.
  - Warm-started across the lambda path; columns solved in parallel
    (`#pragma omp for`), **result order independent of thread count**.
  - Output per node `m`: coefficients `beta[m, ·]` at each lambda, in encoded
    `(lambda_idx, var_idx)` form.
- R layer (`huge.mb`): `beta[[i]] = G block`, `path[[i]] = abs(beta)`,
  `path[[i]] = sign(path + t(path))` for `sym="or"`
  (`sign(path * t(path))` for `"and"`), `sparsity = sum(path)/d/(d-1)`.
  **`path` is 0/1 adjacency; `beta` holds the raw coefficients.**

### 1.9 `huge.glasso` — `huge/R/huge.glasso.R:29` + `huge_core.cpp:471`

- Same correlation input `S = .huge_fast_cor(x)`.
- `huge::glasso` (FHT): alternating coordinate descent, solved per
  connected component of `{(j,k): |S_jk| > λ}` (union-find), warm-started
  across lambda.
  - `GLASSO_DEFAULT_CONVERGENCE_TOL = 1e-4`, `MAX_ITER_EXT=100`,
    `MAX_ITER_INT=MAX_ITER_ACT=10000`.
  - A refinement pass with `1e-8` tol + inverse-residual check
    (`GLASSO_INVERSE_RESIDUAL_TOL=1e-2`, trigger 0.5×).
  - Returns `icov` (precision), `cov` (if `cov.output`), `path` =
    `(icov != 0)` with **zero diagonal**, `loglik`, `sparsity`, `df`.
- `loglik[i] = log det(icov) − trace(S·icov)` (scale-equilibrated Cholesky).

### 1.10 `sparseiCov` (the pulsar `fun`) — `R: SparseICov.R:44`

```
sparseiCov(data, method, npn=FALSE, verbose=FALSE, cov.output=TRUE, ...):
    if (npn) data <- huge::huge.npn(data)
    args <- list(...); method <- switch(method, glasso="glasso", mb="mb")
    if (is.null(args$lambda.min.ratio)) args$lambda.min.ratio <- 1e-3
    do.call(huge::huge, c(args, list(x=data, method=method, verbose=verbose,
                                     cov.output=cov.output)))
```

- `npn = FALSE` in the SPIEC-EASI pipeline (never triggered).
- The `lambda` and `lambda`-path come from `fargs` (the `spiec.easi` grid).

### 1.11 `refit` + accessors — `pulsar/R/refit.R:38`, `R: utilities.R:66`

```
refit.pulsar(obj): est$est = full-data huge object;
                   est$refit[[crit]] = est$est$path[[opt.index]]
getOptX.pulsar.refit: refit→drop0(est$refit$stars); merge→drop0(select$stars$merge[[ind]]);
                      stars→select$stars$summary[ind]; lambda→est$lambda[[ind]];
                      icov/cov→est$est$icov/cov[[ind]] (glasso);
                      beta→est$est$beta[[ind]] (mb)
```

- `getOptInd` = `est$select$stars$opt.index`.

### 1.12 `symBeta` — `R: utilities.R:136`

`symBeta(beta, mode="ave"|"maxabs"|"upper"|"lower")`. **Not** used in the
default MB path (huge already symmetrizes with `or`). Kept in map for
completeness; move to a separate module, out of the default pipeline.

### 1.13 `.spiec.easi` warnings / data checks — `R: spiec-easi.R:35-52`

- Warn if rows sum to 1 (already normalized), warn if any negative, warn if a
  column is constant. These reproduce R's warnings (to stderr).

---

## 2. Methods and scope

| method  | estimator function        | deps used                              | in scope? |
|---------|---------------------------|----------------------------------------|-----------|
| `mb`    | `sparseiCov`→`huge::huge` | huge (MB CD), pulsar, Matrix           | **yes (primary)** |
| `glasso`| `sparseiCov`→`huge::huge` | huge (glasso), pulsar, Matrix          | **yes** |
| `slr`   | `sparseLowRankiCov`→ADMM  | custom C++ ADMM + `sqrtmNewt`, `softSVT`| **ask (D2)** |
| `coat`  | `coat`                    | base R + Matrix                        | **ask (D2)** |
| `ising`, `poisson` | `neighborhood.net` | **glmnet** | **ask (D2)** |
| SparCC (`sparcc`, `sparccboot`, `pval.sparccboot`) | — | VGAM, boot | **ask (D3)** |
| synthetic generators (`synth_comm_from_counts`, `rmv*`, `make_graph`, `graph2prec`, `robustPCA`) | — | VGAM, MASS, igraph-independent | **ask (D3)** |
| `multi.spiec.easi` (list input) | — | — | **ask (D3)** |
| phyloseq / otu_table input | — | Bioconductor | **propose out of scope** (D4) |

`sel.criterion`: only `stars` (and `bstars`, which maps to `stars` +
`lb.stars=ub.stars=TRUE`). `gstars` is a TODO in the R source.

---

## 3. External dependencies — exact calls, versions, options

### 3.1 `huge` 2.0.1 (GPL-2) — **critical, no Python equivalent**
Calls: `huge::huge(x=X, method="mb"|"glasso", lambda=<path>, verbose=FALSE,
cov.output=TRUE)`.
Algorithms (port candidates): `huge_core.cpp` `mb` (nodewise Lasso CD, strong
rule, thol 1e-4, max 1e4) and `glasso` (FHT + component decomposition +
refinement tol 1e-8). Both operate on the **correlation** matrix produced by
`.huge_fast_cor`/`.huge_standardize`.

### 3.2 `pulsar` 0.3.13 (GPL ≥3) — small, portable R-only
`pulsar`, `.ratcheck`, `stars.stability`, `.starsind`, `getLamPath`,
`getMaxCov`, `refit`, `opt.index`. No compiled code (`NeedsCompilation: no`).

### 3.3 `glmnet` 5.1 (GPL-2) — **only** for `ising`/`poisson` via
`neighborhood.net`. If those methods are out of scope, glmnet is not needed.

### 3.4 `Matrix` 1.7-5 (GPL-2/≥3) — sparse containers, `triu`, `drop0`,
`forceSymmetric`, `summary`, `which`. Portable to `scipy.sparse` / dense numpy,
but `triu/tril/drop0` semantics must be matched.

### 3.5 `MASS` (GPL-2/≥3) `ginv`; `VGAM` (GPL-2/≥3) `rdiric/qzipois/...`;
`Rcpp`/`RcppArmadillo` (ADMM, slr only). `boot` (SparCC bootstrap).

### 3.6 Licensing implications (report only — user chooses, D5)
- SpiecEasi: **GPL (≥3)**. huge: **GPL-2**. pulsar: **GPL (≥3)**.
  glmnet: **GPL-2**. Porting this code makes the Python port a **derivative
  work ⇒ must be GPL-compatible** (likely GPL-3+). MIT is **incompatible**.
- The old Python repo shipped an MIT `LICENSE` — that must be replaced.

---

## 4. Dependency-strategy options (the (a)/(b)/(c) question)

For each externally-delegated numeric routine:

**(a) Re-implement the exact algorithm in NumPy/Numba/C++.**
- +True algorithmic parity; same tolerances, same grid, same stop rules.
- +No GPL/BSD dialect problem beyond the port's own licence.
- −Largest effort; must match CD sweep order & floating-point behaviour.

**(b) Use a Python library only if numerically equivalent.**
- `sklearn.linear_model.Lasso` / `graphical_lasso`:
  **NOT equivalent.** Different standardisation (sklearn centers but does not
  rescale to unit variance by default; `Lasso` objective `1/(2n)‖y−Xw‖²` with
  intercept, huge solves `½wᵀSw − bᵀw + λ‖w‖₁` with no intercept on the
  correlation scale), different lambda parameterisation, different stop rule
  (sklearn `tol` on dual gap vs huge relative coefficient change), different
  tie/active-set handling. ⇒ **rule out**.
- No maintained Python package reproduces huge's `mb` + `glasso` exactly.

**(c) Declare non-equivalence explicitly for that piece.**
- +Honest; smaller effort.
- −Fails the user's requirement that `mb`/`glasso` results equal R.

**Recommendation:** (a) for `huge::mb` and `huge::glasso`, and (a) for the
pulsar StARS logic (it is small pure R). This is the only route that satisfies
"strictly equivalent".

---

## 5. Exact equivalence hazards / uncertainty register

| # | Hazard | Handle |
|---|--------|--------|
| U1 | R `sample()` RNG ≠ numpy RNG ⇒ different StARS subsamples | inject indices (D1) |
| U2 | BLAS/LAPACK + C++ solver FMA/order ⇒ ~1e-12…1e-15 diffs; iter count ±1 near thresholds; stability/edges near threshold can flip | document; golden "robust" classification |
| U3 | `eigen()` R returns **descending** eigenvalues, `numpy.eigh` **ascending**; eigenvectors sign/basis ambiguous | compare only sign/basis-invariant quantities |
| U4 | `Matrix::triu(x,k=1)` on a base matrix returns sparse; `max()` over implicit zeros unaffected for `abs` | re-check for `abs=FALSE` |
| U5 | huge MB boundary: at `λ = λ.max` the model can be empty or keep the max-correlation edge depending on summation order | golden case classification |
| U6 | `.huge_standardize` uses a custom shift to avoid overflow; identical only up to fp | treat as 1e-12-level; document |
| U7 | `clr.default` centering `mean(LOG)/mean(nzero)` is non-standard | implement literally (no "fix") |
| U8 | `cov2cor`/`cor` ddof = n-1 | `np.corrcoef`/`np.cov(..., ddof=1)` |
| U9 | Coefficient symmetrisation `or` vs `and`; `path` is 0/1, `beta` raw | keep both outputs |
| U10 | parallel `#pragma omp` in huge → results must be thread-count invariant (verified by huge's own tests) | assert in our port |

---

## 6. Decisions log (RESOLVED with user, 16 Oct 2026)

- **D1 — random subsampling.** Replicate pulsar's `sample(c(1:n), floor(n·r),
  replace=FALSE)` is impossible to reproduce across R/numpy RNG. We will
  (i) allow `subsample_indices=None` injection, (ii) otherwise generate with
  `numpy.random.Generator(seed)`. Golden tests **always inject** the indices
  captured from R. Documented as a known non-equivalence when no injection.
- **D2 — method scope. DECIDED: `mb` only.** `glasso`, `slr`, `coat`,
  `ising`, `poisson` are **out of scope**. This removes the `glmnet` dependency
  entirely and removes the ADMM C++ code. `sparseiCov(method="glasso")` is not
  ported.
- **D3 — accessory scope. DECIDED: only what `spiec.easi(method="mb")`
  needs.** No SparCC (`sparcc`/`sparccboot`/`pval.sparccboot`), no synthetic
  generators (`rmv*`, `synth_comm_from_counts`, `make_graph`, `graph2prec`,
  `robustPCA`), no `multi.spiec.easi`, no `symBeta`/`triu2diag` unless needed
  internally. `symBeta` is not used by the `mb` path.
- **D4 — input types. DECIDED: numeric matrices only**, same as the other
  estimators available to the app (SparCC/CCLasso): raw count matrices.
  No phyloseq/otu_table S3 dispatch.
- **D5 — licence. DECIDED: GPL-3.0-or-later** (most appropriate given the
  derivative-work obligation from GPL≥3 / GPL-2 sources).
- **D6 — name mapping.** R→Python snake_case mapping (to be listed in
  README): `spiec.easi`→`spiec_easi`, `sel.criterion`→`sel_criterion`,
  `lambda.min.ratio`→`lambda_min_ratio`, `pulsar.params`→`pulsar_params`,
  `subsample.ratio`→`subsample_ratio`, `rep.num`→`rep_num`,
  `lambda.log`→`lambda_log`, `pulsar.select`→`pulsar_select`.
- **D7 — repo. DECIDED:** nested git repository inside
  `spiec-easi-python/`; also update the app `README.md` / `CHANGELOG.md` and
  restore the `python/analysis.py` imports once the new API exists.

---

## 7. Proposed repository layout (to confirm)

```
spiec-easi-python/
  pyproject.toml            # runtime: numpy (+ scipy only if needed); dev: pytest, ruff, pytest-cov
  README.md  CHANGELOG.md  LICENSE
  src/spieceasi/
    __init__.py  __main__.py            # CLI: python -m spieceasi
    transforms.py                        # clr, spiec_easi_norm      R: normalization.R, spiec-easi.R:67
    lambda_path.py                       # get_max_cov, get_lam_path R: pulsar/utilities.R
    solvers/mb.py  solvers/glasso.py     # ports of huge::mb / huge::glasso
    stars.py                             # pulsar/stars.stability
    pipeline.py                          # spiec_easi + refit + accessors
    _validation.py  _types.py
  reference/R/                           # vendored SpiecEasi source (read-only)
  reference/deps/{huge,pulsar}/          # vendored dep sources (read-only) + FETCH.md
  tests/r_reference/generate_golden.R    # R golden generator
  tests/golden/                          # committed goldens (R not needed in CI)
  tests/test_equivalence.py              # pytest
  docs/{PORTING_NOTES,EQUIVALENCE,PACKAGING,DIFFERENCES}.md
```

---

## 8. What I still need from you (see the AskUserQuestion that follows)

1. **(a)/(b)/(c)** for `huge::mb` and `huge::glasso` (recommend **a**).
2. Method scope: `mb`+`glasso` only, or also `slr` / `coat` / `ising` /
   `poisson`?
3. Accessory scope: SparCC, synthetic generators, `multi.spiec.easi`?
4. phyloseq/otu_table input — exclude?
5. Licence choice (GPL-3+ recommended given the sources) — user must pick.
6. Repo location: keep the standalone library as a nested folder inside
   `electron-analysis-app` (its own git repo) or a separate directory?
