"""Data transformations.

Faithful ports of the SpiecEasi R functions:

- ``clr.default`` / ``clr.matrix``  (R: normalization.R:95,105)
- ``.spiec.easi.norm``              (R: spiec-easi.R:67)

The SpiecEasi CLR is *not* the textbook centred log-ratio: it centres by
``mean(LOG) / mean(x >= tol)`` (the mean of the log-values divided by the
fraction of non-zero entries). This is reproduced literally - see
``docs/DIFFERENCES.md``.
"""

from __future__ import annotations

from typing import Any

import numpy as np

# R: .Machine$double.eps
_DOUBLE_EPS = float(np.finfo(np.float64).eps)


def clr(
    x: Any,
    base: float = float(np.e),
    tol: float = _DOUBLE_EPS,
) -> np.ndarray:
    """Centred log-ratio transform of a 1-D vector.

    R: ``clr.default`` (normalization.R:95)

    ``nzero <- (x.f >= tol); LOG <- log(ifelse(nzero, x.f, 1), base);
    ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)``
    """
    x = np.asarray(x, dtype=float)
    nzero = x >= tol
    log_x = np.where(nzero, x, 1.0)
    log_x = np.log(log_x) / np.log(base)
    frac = float(np.mean(nzero))
    if frac == 0.0:
        return np.zeros_like(x)
    centre = float(np.mean(log_x)) / frac
    return np.where(nzero, log_x - centre, 0.0)


def clr_rows(
    x: Any,
    base: float = float(np.e),
    tol: float = _DOUBLE_EPS,
) -> np.ndarray:
    """Apply :func:`clr` independently to each row of a 2-D array.

    Equivalent to ``clr.matrix(x, mar=1)`` (R: normalization.R:105).
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 2:
        raise ValueError("clr_rows expects a 2-D array")
    nzero = x >= tol
    log_x = np.where(nzero, x, 1.0)
    log_x = np.log(log_x) / np.log(base)
    frac = nzero.mean(axis=1, keepdims=True)
    centre = np.where(frac > 0.0, log_x.mean(axis=1, keepdims=True) / frac, 0.0)
    out = np.where(nzero, log_x - centre, 0.0)
    return np.where(frac > 0.0, out, 0.0)


def clr_columns(
    x: Any,
    base: float = float(np.e),
    tol: float = _DOUBLE_EPS,
) -> np.ndarray:
    """Apply :func:`clr` independently to each column.

    Equivalent to ``clr.matrix(x, mar=2)`` (R: normalization.R:105), the R
    default margin.
    """
    return clr_rows(np.asarray(x, dtype=float).T, base=base, tol=tol).T


def spiec_easi_norm(data: Any) -> np.ndarray:
    """SpiecEasi data transformation: ``t(clr(data + 1, 1))``.

    R: ``.spiec.easi.norm`` (spiec-easi.R:67). Input ``data`` is an
    ``n_samples x n_features`` count matrix; the output has the same shape.

    R computes ``t(clr(data + 1, 1))``. ``clr.matrix(..., mar=1)`` applies
    :func:`clr` to each row and (because R's ``apply`` over rows transposes the
    result) returns ``features x samples``; the outer ``t`` restores
    ``samples x features``. :func:`clr_rows` already keeps the row layout, so
    no transpose is needed here.
    """
    data = np.asarray(data, dtype=float)
    return clr_rows(data + 1.0)
