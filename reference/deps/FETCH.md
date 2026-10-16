# Vendored reference sources

This directory contains the **source of the upstream projects this port is
based on**, kept only so the exact algorithms, defaults, tolerances and lambda
grids can be re-checked offline. Nothing here is imported, built or shipped at
runtime.

The port targets `method = "mb"` only (see `docs/PORTING_NOTES.md`, decision D2).

## Kept: sources actually ported

### SpiecEasi 1.99.0 — the pipeline itself

- Location in this repo: `reference/R/` (R sources, C++ sources, tests, vignettes).
- Upstream: <https://github.com/zdk123/SpiecEasi>
- Version: `1.99.0` (see `reference/R/DESCRIPTION`; it is the development version
  that was installed on the reference machine).
- Licence: **GPL (>= 3)**.
- Ported from it: `spiec.easi.default`, `.spiec.easi.norm`, `clr.default` /
  `clr.matrix`, `.data.checks`, `getOptX.pulsar.refit`, and the pulse of the
  pipeline (transform → lambda grid → StARS → refit).

### huge 2.0.1 — Meinshausen–Bühlmann neighbourhood selection

- Location in this repo: `reference/deps/huge/`
- Upstream: CRAN <https://cran.r-project.org/package=huge>
- Tarball: <https://cran.r-project.org/src/contrib/huge_2.0.1.tar.gz>
- **SHA-256 (tarball):** `ca525f3d306a98520962d36d64d2dc3c08f7c1b387d44a2938a88c80ec76f133`
- Licence: **GPL-2** (see `reference/deps/huge/DESCRIPTION`).
- Ported from it: `.huge_standardize`, `.huge_fast_cor`, `.huge_preprocess`,
  `huge.mb` and the C++ `huge::mb` coordinate descent
  (`src/huge_core.cpp`: `mb`, `mb_partial_residual`, `mb_refine_active`,
  `threshold_l1`, `collect_sorted`).

### pulsar 0.3.13 — StARS model selection

- Location in this repo: `reference/deps/pulsar/`
- Upstream: CRAN archive (the package is archived) and
  <https://github.com/zdk123/pulsar>
- Tarball: <https://cran.r-project.org/src/contrib/Archive/pulsar/pulsar_0.3.13.tar.gz>
- **SHA-256 (tarball):** `3703c0930ffdc46a3472c0dc26b23c08a92d31e2ad18ac5f1604dd256509a7f8`
- Licence: **GPL (>= 3)** (see `reference/deps/pulsar/DESCRIPTION`).
- Ported from it: `getLamPath`, `getMaxCov`, `.ratcheck`, `stars.stability`,
  `.starsind`, the `pulsar` subsampling loop and `refit.pulsar`.

## Removed: fetched for the audit, not needed for `mb`

These were downloaded while deciding the dependency strategy (options a/b/c in
`docs/PORTING_NOTES.md`) and then deleted once `method = "mb"` was fixed as the
only scope, because neither package is on the `mb` code path.

| package | version | tarball URL | SHA-256 (tarball) | licence |
|---|---|---|---|---|
| glmnet | 5.1 | <https://cran.r-project.org/src/contrib/glmnet_5.1.tar.gz> | `a6c23c906806f89a451de30b202e5bf088577f4da260d0b7f8887a7df3f22778` | GPL-2 |
| Matrix | 1.7-5 | <https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.7-5.tar.gz> | `308b0edd7bcb023a8bb50cc5cec98e66b3658c94a5badfb1ce6024e595fa8a83` | GPL-2 \| GPL-3 |

`glmnet` is only used by SpiecEasi's `neighborhood.net()` for the non-default
`ising`/`poisson` models (`reference/R/R/SparseICov.R:113`). `Matrix` provides
sparse containers that this port replaces with NumPy/SciPy-free code; its
*algorithms* are not used.

## Licence and attribution

- This project is distributed under **GPL-3.0-or-later** (see `/LICENSE`). It is
  a derivative work of the sources above, so their notices apply.
- `SpiecEasi` (GPL≥3) and `pulsar` (GPL≥3) are compatible with GPL-3.
- **`huge` is declared `GPL-2`** in its `DESCRIPTION`. Read strictly that is
  *GPL-2-only*, which is **not** compatible with GPL-3 for a combined work. This
  is the same arrangement the upstream `SpiecEasi` package itself relies on (it
  is `GPL (>= 3)` and imports `huge`), but it is a real obligation to be aware
  of. See `docs/PORTING_NOTES.md` §3.6 and §D5.
- If a different licensing outcome is required, the options are (i) keep GPL-3
  and treat the risk as accepted/documented, (ii) contact the `huge`
  maintainers for a relicense, or (iii) re-derive the MB solver from the
  published algorithm rather than from the `huge` source (a much larger
  clean-room effort). This is a decision for the project owner.

## How to refresh

```bash
curl -sSLO https://cran.r-project.org/src/contrib/huge_2.0.1.tar.gz
curl -sSLO https://cran.r-project.org/src/contrib/Archive/pulsar/pulsar_0.3.13.tar.gz
sha256sum huge_2.0.1.tar.gz pulsar_0.3.13.tar.gz   # must match the hashes above
tar xzf huge_2.0.1.tar.gz
tar xzf pulsar_0.3.13.tar.gz
```
