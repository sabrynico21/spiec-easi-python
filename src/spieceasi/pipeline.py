"""The SPIEC-EASI pipeline (``method="mb"``).

Faithful port of:

- ``spiec.easi.default`` (R: spiec-easi.R:149)
- ``.spiec.easi.norm``   (R: spiec-easi.R:67)
- ``sparseiCov``         (R: SparseICov.R:44)  - the pulsar ``fun``
- ``pulsar``             (pulsar/R/mcPulsarSelect.R:168)
- ``refit.pulsar``       (pulsar/R/refit.R:38)

Only ``method="mb"`` is ported; ``glasso``/``slr``/``coat``/``ising``/
``poisson`` are out of scope (docs/PORTING_NOTES.md, D2).
"""

from __future__ import annotations

import math
import warnings
from typing import Any

import numpy as np

from ._validation import validate_counts
from .lambda_path import get_lam_path, get_max_cov
from .mb import mb
from .preprocess import huge_fast_cor
from .stars import stars_stability
from .transforms import spiec_easi_norm

_DEFAULT_LAMBDA_MIN_RATIO = 1e-3
_DEFAULT_NLAMBDA = 20
_DEFAULT_REP_NUM = 20
_DEFAULT_THRESH = 0.05


def _ratcheck(subsample_ratio: float | None, n: int) -> float:
    """``.ratcheck`` (pulsar/R/mcPulsarSelect.R:15)."""
    if subsample_ratio is None:
        if n > 144:
            return 10.0 * math.sqrt(n) / n
        return 0.8
    return float(subsample_ratio)


def r_cor(x: np.ndarray) -> np.ndarray:
    """``stats::cor`` (Pearson, ddof = n-1)."""
    return np.corrcoef(x, rowvar=False)


def _solve_subsample(
    job: tuple[np.ndarray, np.ndarray, str],
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Solve the MB path for one subsample (module-level for multiprocessing)."""
    Xk, lambdas, sym = job
    betas_k, paths_k, _ = mb(huge_fast_cor(Xk), lambdas, sym=sym)
    return paths_k, betas_k


def _pool_context() -> Any:
    """Multiprocessing context for the subsample pool.

    ``fork`` is used on Linux only: it avoids Python 3.14's ``forkserver``
    default (which re-imports the entry point and fails when the caller runs
    from stdin) and needs no importable ``__main__``. Everywhere else
    (Windows, macOS) the platform default is used, i.e. ``spawn``; that path
    requires an entry point guarded by ``if __name__ == "__main__"`` and
    ``multiprocessing.freeze_support()`` (both are present in
    :mod:`spieceasi.__main__` and :func:`spieceasi.cli.main`).

    The ``spawn`` path has **not been executed** in this repository (the
    equivalence work was Linux x86-64 only); see ``docs/EQUIVALENCE.md``.
    """
    import multiprocessing as mp
    import sys

    if sys.platform.startswith("linux"):
        return mp.get_context("fork")
    return mp.get_context("spawn")


def spiec_easi(
    data: Any,
    method: str = "mb",
    sel_criterion: str = "stars",
    verbose: bool = False,
    pulsar_select: bool = True,
    pulsar_params: dict[str, Any] | None = None,
    lambda_log: bool = True,
    nlambda: int = _DEFAULT_NLAMBDA,
    lambda_min_ratio: float = _DEFAULT_LAMBDA_MIN_RATIO,
    lambda_max: float | None = None,
    sym: str = "or",
    subsample_indices: Any = None,
    keep_replica_betas: bool = False,
    n_jobs: int = 1,
    progress_callback: Any = None,
) -> dict[str, Any]:
    """Run the SPIEC-EASI ``mb`` pipeline.

    Parameters mirror the R function (snake_case). ``pulsar_params`` accepts the
    StARS settings ``thresh``, ``subsample_ratio``, ``rep_num``, ``seed`` and
    ``ncores`` (``ncores`` is accepted for API parity but the current
    implementation is serial).

    ``subsample_indices`` (optional, 0-based, shape ``(rep_num, size)``) injects
    the StARS subsample indices. When omitted they are drawn with
    ``numpy.random.default_rng(seed)``; because R and NumPy RNGs differ, results
    are then statistically - not numerically - equivalent to R.
    """
    if method != "mb":
        raise NotImplementedError(
            f'method={method!r} is not ported; only "mb" is supported'
        )
    if sel_criterion not in ("stars", "bstars"):
        raise ValueError(f"Unknown selection criterion: {sel_criterion!r}")

    pulsar_params = dict(pulsar_params or {})
    thresh = pulsar_params.pop("thresh", _DEFAULT_THRESH)
    subsample_ratio = pulsar_params.pop("subsample_ratio", None)
    rep_num = int(pulsar_params.pop("rep_num", _DEFAULT_REP_NUM))
    seed = pulsar_params.pop("seed", None)
    n_jobs = int(pulsar_params.pop("ncores", n_jobs))
    if pulsar_params:
        raise ValueError(f"Unrecognized pulsar_params: {sorted(pulsar_params)}")

    counts = validate_counts(data)

    # 1. transform
    X = spiec_easi_norm(counts)
    n, p = X.shape

    # 2. lambda grid (max|cor(X)|, R's cor - not huge's standardised one)
    if lambda_max is None:
        lambda_max = get_max_cov(r_cor(X))
    lambdas = get_lam_path(
        lambda_max, lambda_max * lambda_min_ratio, nlambda, log=lambda_log
    )

    # 3. StARS subsampling + per-replica MB paths
    ratio = _ratcheck(subsample_ratio, n)
    size = int(math.floor(n * ratio))
    if subsample_indices is None:
        rng = np.random.default_rng(seed)
        indices = [rng.choice(n, size=size, replace=False) for _ in range(rep_num)]
    else:
        indices = [np.asarray(a, dtype=np.intp) for a in subsample_indices]
        rep_num = len(indices)

    # full-data fit (the refit / est$est model)
    S_full = huge_fast_cor(X)
    beta_full, path_full, hit_full = mb(S_full, lambdas, sym=sym)

    premerge: list[list[np.ndarray]] = [[] for _ in range(nlambda)]
    replica_betas: list[list[np.ndarray]] = []

    jobs = [(X[indices[k]], lambdas, sym) for k in range(rep_num)]
    if n_jobs > 1 and rep_num > 1:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(
            max_workers=min(n_jobs, rep_num), mp_context=_pool_context()
        ) as ex:
            solved = list(ex.map(_solve_subsample, jobs))
    else:
        solved = [_solve_subsample(job) for job in jobs]

    for k, (paths_k, betas_k) in enumerate(solved):
        if progress_callback is not None:
            progress_callback(k + 1, rep_num)
        if keep_replica_betas:
            replica_betas.append(betas_k)
        for i in range(nlambda):
            premerge[i].append(paths_k[i])

    stars = stars_stability(premerge, thresh, rep_num, p)
    opt_index = int(stars["opt_index"])

    if opt_index == 1:
        direction = "larger" if bool(np.any(stars["summary"] >= 0.1)) else "smaller"
        warnings.warn(
            f"Optimal lambda may be {direction} than the supplied values",
            stacklevel=2,
        )

    # 4. refit on the full data at the selected index (1-based -> 0-based)
    refit = path_full[opt_index - 1]

    select = {
        "stars": stars,
        "est": {"beta": beta_full, "path": path_full, "method": "mb"},
        "call": {
            "fun": "sparseiCov",
            "fargs": {"method": "mb", "lambda": lambdas},
        },
    }
    result = {
        "lambda": lambdas,
        "lambda_max": float(lambda_max),
        "opt_index": opt_index,
        "opt_lambda": float(lambdas[opt_index - 1]),
        "stability": stars["summary"],
        "merge": stars["merge"],
        "refit": refit,
        "beta": beta_full,
        "path": path_full,
        "select": select,
        "est": select["est"],
        "fun": "sparseiCov",
        "method": "mb",
        "hit_max_iter": bool(hit_full),
        "data": X,
    }
    if keep_replica_betas:
        result["_replica_betas"] = replica_betas
    return result


# ---------------------------------------------------------------------------
# Accessors (R: utilities.R:60 getOptX.pulsar.refit)
# ---------------------------------------------------------------------------


def get_opt_ind(result: dict[str, Any]) -> int:
    """``getOptInd`` -> ``est$select$stars$opt.index``."""
    return int(result["opt_index"])


def get_opt_lambda(result: dict[str, Any]) -> float:
    """``getOptLambda``."""
    return float(result["lambda"][get_opt_ind(result) - 1])


def get_opt_merge(result: dict[str, Any]) -> np.ndarray:
    """``getOptMerge`` -> edge-frequency merge at the selected lambda."""
    return result["merge"][get_opt_ind(result) - 1]


def get_stability(result: dict[str, Any]) -> float:
    """``getStability`` -> StARS summary at the selected lambda."""
    return float(result["stability"][get_opt_ind(result) - 1])


def get_opt_net(result: dict[str, Any]) -> np.ndarray:
    """``getOptNet`` / ``getRefit`` -> the selected adjacency."""
    return result["refit"]


def get_opt_beta(result: dict[str, Any]) -> np.ndarray:
    """``getOptBeta`` -> full-data MB coefficients at the selected lambda."""
    return result["beta"][get_opt_ind(result) - 1]
