# spieceasi (Python)

A faithful, standalone Python port of the **SPIEC-EASI** `mb` pipeline from the
R package [`SpiecEasi`](https://github.com/zdk123/SpiecEasi) (Kurtz *et al.*,
2015, *PLoS Comput Biol*), plus the exact numeric routines it delegates to
(`huge` for Meinshausen–Bühlmann neighbourhood selection, `pulsar` for StARS).

The goal is **strict numerical equivalence with the R reference**, not a
re-imagining: the same data transformations, the same lambda grid, the same
solver tolerances and iteration limits, the same StARS sub-sampling and
selection rules. Where something cannot be identical (notably the RNG used for
the StARS subsamples), it is documented and made injectable.

> **Scope.** Only `method = "mb"` is implemented. `glasso`, `slr`, `coat`,
> `ising`/`poisson`, SparCC and the synthetic-data generators of the R package
> are intentionally out of scope. See `docs/PORTING_NOTES.md`.

> **Equivalence scope.** The step-by-step equivalence against R
> (`docs/EQUIVALENCE.md`) was measured on **Linux x86-64 (GCC) only**. The test
> suite runs anywhere without R (the goldens are committed), but **Windows and
> macOS have not been tested** and the CI workflow, although present, has never
> been run. Treat the equivalence claims as Linux x86-64 only.

> **Status.** Under construction. See `docs/PORTING_NOTES.md`,
> `docs/EQUIVALENCE.md` and `CHANGELOG.md`.

## Install

```bash
pip install -e ".[dev]"
```

Runtime dependency: `numpy` only.

## Usage (Python)

```python
import numpy as np
from spieceasi import spiec_easi

# counts: (n_samples, n_features) raw count matrix
result = spiec_easi(counts, method="mb", sel_criterion="stars")
adj = result["refit"]          # selected adjacency
opt = result["opt_index"]      # selected lambda index
```

Parameter names follow the R source, converted to snake_case:

| R (SpiecEasi / pulsar) | Python |
|---|---|
| `data` | `data` |
| `method` | `method` |
| `sel.criterion` | `sel_criterion` |
| `verbose` | `verbose` |
| `pulsar.select` | `pulsar_select` |
| `pulsar.params` | `pulsar_params` |
| `lambda.log` | `lambda_log` |
| `lambda.min.ratio` | `lambda_min_ratio` (inside `pulsar_params`/`**kwargs`) |
| `nlambda` | `nlambda` |
| `subsample.ratio` | `subsample_ratio` |
| `rep.num` | `rep_num` |
| `thresh` | `thresh` |
| `seed` | `seed` |

The same defaults as R are used: `method="mb"`, `sel_criterion="stars"`,
`nlambda=20`, `lambda_min_ratio=1e-3`, `lambda_log=True`, `rep_num=20`,
`thresh=0.05`.

### Determinism

The StARS subsample indices are generated from R's RNG in the reference. They
cannot be reproduced with NumPy, so they can be injected:

```python
result = spiec_easi(counts, subsample_indices=indices)  # shape (rep_num, size)
```

When not injected, they are drawn with `numpy.random.default_rng(seed)`.
Results that depend on the subsampling (stability, selected lambda, final edges)
are then **statistically** equivalent to R but not identical.

### Parallelism

StARS runs `rep_num` independent subsample fits. Pass `n_jobs > 1` to run them
with a process pool; results are identical for any worker count (the subsample
results are aggregated in a fixed order):

```python
result = spiec_easi(counts, pulsar_params={"rep_num": 20}, n_jobs=4)
```

`n_jobs=1` (the default) is fully serial.

## Usage (CLI / Electron)

The CLI reads a JSON request from a file or stdin and writes a JSON response to
stdout (errors go to stderr with a non-zero exit code; progress, if any, is
reported on stderr only, never stdout):

```bash
python -m spieceasi --input request.json > response.json
spieceasi --input request.json
```

Request (JSON):

```json
{
  "counts": [[12, 0, 5], [7, 3, 0]],
  "params": {
    "method": "mb",
    "sel_criterion": "stars",
    "nlambda": 20,
    "lambda_min_ratio": 0.001,
    "rep_num": 20,
    "thresh": 0.05,
    "seed": 17
  }
}
```

Response (JSON):

```json
{
  "lambda_path": [0.4031, 0.2801, 0.1947, "..."],
  "opt_index": 3,
  "opt_lambda": 0.1947,
  "stability": [0.013, 0.031, 0.035, "..."],
  "edges": [[0, 4, 0.95], [1, 7, 0.90]],
  "refit": {"format": "edge_list", "edges": [[0, 4], [1, 7]]},
  "diagnostics": {"n": 30, "p": 4, "method": "mb", "converged": true}
}
```

- Indices are 0-based; `edges` are `[i, j, weight]` where `weight` is the StARS
  edge frequency at the selected lambda; `refit` is the selected 0/1 graph.
- `diagnostics.converged` is `false` when the MB solver hit its iteration limit
  (huge would warn "mb solver reached its iteration limit").

Calling from Node (`child_process.spawn`):

```javascript
const { spawn } = require("child_process");
const child = spawn("spieceasi", ["--input", "request.json"]);
let out = "", err = "";
child.stdout.on("data", (d) => (out += d));
child.stderr.on("data", (d) => (err += d));
child.on("close", (code) => {
  if (code !== 0) throw new Error(`spiec_easi failed: ${err}`);
  const result = JSON.parse(out);
});
```

Sparse outputs (adjacency, merge/stability matrices) are serialised as explicit
edge lists `[i, j, weight]` with 0-based indices.

## Reference

- Article: Kurtz ZD, Müller CL, Miraldi ER, Littman DR, Blaser MJ, Bonneau RA.
  *Sparse and Compositionally Robust Inference of Microbial Ecological
  Networks.* PLoS Comput Biol 11(5): e1004226, 2015.
- Original R package: <https://github.com/zdk123/SpiecEasi>
- `huge`: <https://github.com/Gatech-Flash/huge>
- `pulsar`: <https://github.com/zdk123/pulsar>

## Licence

GPL-3.0-or-later. The port is a derivative work of GPL-licensed code
(`SpiecEasi` GPL≥3, `huge` GPL-2, `pulsar` GPL≥3); see `docs/PORTING_NOTES.md`.
