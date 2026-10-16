"""Lambda-path utilities.

Faithful ports from the R ``pulsar`` package:

- ``getMaxCov`` (pulsar/R/utilities.R:130)
- ``getLamPath`` (pulsar/R/utilities.R:111)
"""

from __future__ import annotations

from typing import Any

import numpy as np


def _is_symmetric(x: np.ndarray) -> bool:
    """Mirror of R ``isSymmetric`` (numeric tolerance)."""
    if x.ndim != 2 or x.shape[0] != x.shape[1]:
        return False
    return bool(np.allclose(x, x.T))


def get_max_cov(
    x: Any,
    cov: bool | None = None,
    abs: bool = True,  # noqa: A002 - name mirrors the R argument
    diag: bool = False,
) -> float:
    """Maximum (absolute) value of a covariance/correlation matrix.

    R: ``getMaxCov`` (pulsar/R/utilities.R:130)::

        if (!cov) x <- cov(x)
        tmp <- Matrix::triu(x, k = if (diag) 0 else 1)
        max(if (abs) abs(tmp) else tmp)

    ``cov`` defaults to ``isSymmetric(x)``. The strict upper triangle (``k=1``)
    is used unless ``diag`` is set.
    """
    x = np.asarray(x, dtype=float)
    if cov is None:
        cov = _is_symmetric(x)
    if not cov:
        x = np.cov(x, rowvar=False, ddof=1)
    k = 0 if diag else 1
    upper = np.triu(x, k=k)
    if abs:
        upper = np.abs(upper)
    return float(np.max(upper))


def get_lam_path(max: float, min: float, len: int, log: bool = False) -> np.ndarray:  # noqa: A002
    """Descending lambda path, equally or log-spaced.

    R: ``getLamPath`` (pulsar/R/utilities.R:111)::

        if (log) { min <- log(min); max <- log(max) }
        lams <- seq(max, min, length.out=len)
        if (log) exp(lams) else lams
    """
    if max < min:
        raise ValueError("Did you flip min and max?")
    if log:
        min = float(np.log(min))
        max = float(np.log(max))
    lams = np.linspace(max, min, num=len)
    if log:
        return np.exp(lams)
    return lams
