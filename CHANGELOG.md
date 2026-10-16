# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Initial repository scaffold for a from-scratch Python port of the R
  SpiecEasi `mb` pipeline.
- **Core implementation** (`src/spieceasi/`): `transforms` (`clr`,
  `spiec_easi_norm`), `lambda_path` (`get_max_cov`, `get_lam_path`),
  `preprocess` (`huge_standardize`, `huge_fast_cor`), `mb` (exact port of
  `huge::mb`'s coordinate descent), `stars` (`pulsar` StARS), `pipeline`
  (`spiec_easi` + accessors), `cli` (`python -m spieceasi`, JSON in/out) and
  `_validation`.
- **Equivalence tests** (`tests/test_equivalence.py`) comparing against the R
  goldens step by step, plus `tests/test_validation.py` (validation + CLI) and
  `tests/test_properties.py` (structure + synthetic recovery).
- **Optional multiprocessing** over StARS subsamples (`n_jobs`); results are
  identical for any worker count.
- `reference/R/` — vendored copy of the SpiecEasi R source (the port's source of
  truth) and `reference/deps/` — vendored sources of `huge` 2.0.1,
  `pulsar` 0.3.13, `glmnet` 5.1 and `Matrix` 1.7-5 used to document the exact
  delegated algorithms, tolerances and lambda grids.
- `docs/PORTING_NOTES.md` — the Passo 0 function-by-function analysis, the
  external-dependency audit, the equivalence-hazard register and the decision
  log (scope = `mb` only; exact-algorithm port; GPL-3.0-or-later).
- `tests/r_reference/generate_golden.R` — R golden generator (loads
  `reference/R/` via `pkgload::load_all`), with 8 deterministic cases committed
  under `tests/golden/`.
- `tests/r_reference/README.md` — regeneration instructions and the pinned R /
  package versions.
- GPL-3.0 full text as `LICENSE`.

### Removed

- The previous standalone `spiec-easi-python` implementation (source, tests,
  golden fixtures, docs, benchmarks, examples, packaging) was deleted; the port
  is being rebuilt from the R source.
