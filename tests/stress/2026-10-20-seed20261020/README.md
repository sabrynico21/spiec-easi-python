# SPIEC-EASI MB stress sweep — 20 Oct 2026, seed 20261020

Predeclared seeded stress suite required by
`docs/EQUIVALENCE_ACCEPTANCE_CONTRACT.md` ("predeclared seeded stress suite of
at least 1,000 small cases per method family, weighted toward sparse, singular,
and threshold-near inputs").

## Command

```bash
cd spiec-easi-python
.venv/bin/python tests/r_reference/run_status_stress.py \
    /tmp/spiec-stress-1000 --cases 1000 --seed 20261020
```

Runner: `tests/r_reference/run_status_stress.py`; R side:
`tests/r_reference/stress_status_R.R` (pinned source in `reference/R/`, loaded
with `pkgload::load_all`).

## Environment

| item | value |
|---|---|
| date | 2026-10-20 |
| R | 4.6.1 |
| Python | 3.14 (venv `spiec-easi-python/.venv`) |
| NumPy | see `environment.json` |
| platform | Linux x86-64, 8 logical CPUs, serial run |
| settings | `nlambda=5`, `lambda_min_ratio=0.1`, `lambda_log=TRUE`, `rep_num=4`, `thresh=0.2` |

## Predeclared design

Six equally represented strata (`case_index % 6`): moderate sparsity (~15%
zeros), sparse (~65% zeros), rank-deficient (constant feature) + ~35% zeros,
low library sizes + ~45% zeros, `p > n`, and threshold-near weak signal + ~30%
zeros. Per-case seeds and StARS subsample indices are drawn from one NumPy
master stream; the indices are injected into both languages, so the comparison
isolates arithmetic rather than resampling. Case tree: 1,000 directories, 3,000
files, checksummed in `cases.sha256` (the cases themselves stay in the run
output directory and are regenerable from the seed).

## Result

| metric | value |
|---|---:|
| cases | 1,000 |
| status (`finite`/`nonfinite`/`error`), Python | 1000 / 0 / 0 |
| status (`finite`/`nonfinite`/`error`), R | 1000 / 0 / 0 |
| status mismatches | **0** |
| `opt_index` mismatches | **21** (all at p = 2) |
| `opt_index` mismatches passing the predeclared boundary rule | 0 |
| refit edge-set mismatches | **256** |
| edge-set mismatches passing the predeclared boundary rule | 87 |
| unexplained by the predeclared rule (exit code 1) | **190** |

Runner exit code: **1** (unexplained mismatches present). The rule was not
relaxed after seeing the results.

## Post-run triage (not a rule change)

1. **Every one of the 256 edge-set differences is a single edge**, and in all
   256 that edge is exactly the maximum-`|cor|` feature pair, i.e. the only
   edge reachable at `lambda = lambda.max = max|cor(X)|`. 254 of the 256 are
   cases where both languages selected index 1 (the refit *is* the lambda.max
   model); the other 2 are Python `opt_index=5` / R `opt_index=1`, where R's
   refit at lambda.max still carries that pair.
2. **All 21 `opt_index` mismatches are p = 2 designs** and swap between the two
   extreme indices (14 × Python 5 / R 1, 7 × Python 1 / R 5). For p = 2 the
   StARS variability is quantised (0, 0.75 or 1 for `rep_num=4`), so one
   replica's lambda.max edge flips the selection between "smallest lambda" and
   "lambda.max" with no interior optimum. Dumping R's per-replica supports
   shows the tie directly, e.g. `case-00055`: R replicas hold the single edge
   at lambda.max in 3/4 draws, Python in 0/4.
3. **Matched-grid control.** Re-running Python with R's exact
   `lambda_max`/grid removes 7 of the 21 `opt_index` mismatches (and their edge
   differences) and removes the lambda.max edge in 5 of 10 sampled edge
   mismatches. The rest persist on R's grid, so the remaining divergence is
   inside `huge`'s MB coordinate descent at the `lambda.max` boundary — the
   deviation already documented in `docs/PORTING_NOTES.md` and
   `docs/EQUIVALENCE.md` (BLAS summation vs R's scalar loop) and in the
   committed goldens (step 9: exactly 1 divergent edge in 8/20 cases when R's
   `opt_index == 1`).

So the 190 cases flagged "unexplained" by the predeclared edge rule all belong
to the single, previously documented lambda.max tie class; none exposes a new
failure mode. The predeclared rule keyed on StARS selection frequency, which
does not describe an MB support tie at the soft-threshold boundary; that is a
rule design gap, recorded here rather than patched afterwards. Any replacement
criterion (e.g. "differing edge is the max-`|cor|` pair and `opt_index == 1` in
both languages", or a coefficient-boundary test at lambda.max) must be
registered *before* the next run and the sweep re-executed under it.

## Timing (this machine, serial)

| phase | seconds | per case |
|---|---:|---:|
| Python classification (1,000 cases) | 22.9 | mean 23 ms, median 19 ms, p95 56 ms |
| R classification (1,000 cases, incl. ~10 s `pkgload` startup) | 110.4 | mean ~110 ms |
| **total wall clock** | **135.2 (≈ 2 min 16 s)** | |

Projection: cost scales roughly linearly in `nlambda × rep_num`. At the
`spiec.easi` defaults (`nlambda=20`, `rep_num=20`) the same 1,000-case sweep
would take roughly 20 × 135 s ≈ **45 min** (≈ 40 min Python + R) on this
machine; a smaller case count remains a valid smoke test.

## Files

- `summary.json` — machine-readable run result (counts, timings, previews).
- `python_status.csv`, `r_status.csv` — per-case status, selected index, refit
  edges, lambda.max, StARS margin, merge frequencies at the selected lambda.
- `environment.json` — Python/NumPy/platform/R versions, CPU count.
- `cases.sha256` — SHA-256 of every archived case file (3,000 entries).
