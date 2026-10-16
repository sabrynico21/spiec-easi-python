"""huge data preprocessing.

Exact ports of the huge 2.0.1 internals the ``mb`` path relies on:

- ``.huge_is_covariance_input`` (huge/R/huge-internal.R:90)
- ``.huge_standardize``        (huge/R/huge-internal.R:47)
- ``.huge_fast_cor``           (huge/R/huge-internal.R:78)
- ``.huge_preprocess``         (huge/R/huge-internal.R:228)

``huge.mb`` does **not** run on R's ``cor()``: it standardises each column with
``.huge_standardize`` (a numerically careful centre/scale) and forms the
correlation matrix as ``crossprod(xs) / (n - 1)``, clamped to [-1, 1]. This
differs from ``stats::cor`` at the 1e-16 level, which matters exactly at
``lambda.max`` (the boundary of the MB path), so it must be reproduced.
"""

from __future__ import annotations

from typing import Any

import numpy as np

_DOUBLE_EPS = float(np.finfo(np.float64).eps)


def huge_is_covariance_input(x: np.ndarray) -> bool:
    """``.huge_is_covariance_input`` (huge-internal.R:90)."""
    if x.ndim != 2 or x.shape[0] != x.shape[1]:
        return False
    values = x * 1.0
    if not np.all(np.isfinite(values)):
        return False
    tolerance = 100.0 * _DOUBLE_EPS
    dimension = values.shape[0]
    if dimension <= 1:
        return True
    diagonal = np.abs(np.diag(values))
    diagonal_root = np.sqrt(diagonal)
    for column in range(1, dimension):
        rows = np.arange(column)
        left = values[rows, column]
        right = values[column, rows]
        different = left != right
        if not np.any(different):
            continue
        rows = rows[different]
        left = left[different]
        right = right[different]
        covariance_scale = diagonal_root[rows] * diagonal_root[column]
        equal_diagonal = diagonal[rows] == diagonal[column]
        covariance_scale = np.where(equal_diagonal, diagonal[column], covariance_scale)
        reference = np.maximum(
            np.maximum(np.abs(left), np.abs(right)), covariance_scale
        )
        left = left / reference
        right = right / reference
        scale = np.maximum(np.abs(left), np.abs(right))
        threshold = np.where(scale <= tolerance, tolerance, tolerance * scale)
        if np.any(np.abs(left - right) > threshold):
            return False
    return True


def huge_standardize(x: Any) -> np.ndarray:
    """``.huge_standardize`` (huge-internal.R:47): unit-variance centred columns."""
    x = np.asarray(x, dtype=float)
    n, d = x.shape
    xs = np.zeros((n, d))
    for j in range(d):
        source = x[:, j] * 1.0
        column_scale = float(np.max(np.abs(source)))
        if not np.isfinite(column_scale) or column_scale == 0.0:
            raise ValueError(
                "Data must have finite, positive sample standard deviations."
            )
        reference = float(source[0])
        direct = (
            (source == 0.0) | (reference == 0.0) | ((source < 0.0) == (reference < 0.0))
        )
        delta = source / column_scale - reference / column_scale
        delta[direct] = (source[direct] - reference) / column_scale
        delta = delta - np.mean(delta)
        centered_scale = float(np.max(np.abs(delta)))
        if not np.isfinite(centered_scale) or centered_scale == 0.0:
            raise ValueError(
                "Data must have finite, positive sample standard deviations."
            )
        delta = delta / centered_scale
        column_sd = float(np.std(delta, ddof=1))
        if not np.isfinite(column_sd) or column_sd <= 0.0:
            raise ValueError(
                "Data must have finite, positive sample standard deviations."
            )
        xs[:, j] = delta / column_sd
    return xs


def huge_fast_cor(x: Any) -> np.ndarray:
    """``.huge_fast_cor`` (huge-internal.R:78).

    Correlation via standardised crossprod.
    """
    x = np.asarray(x, dtype=float)
    try:
        xs = huge_standardize(x)
    except ValueError:
        n = x.shape[0]
        return np.corrcoef(x, rowvar=False) * (1.0 if n > 1 else 1.0)
    n = xs.shape[0]
    S = (xs.T @ xs) / (n - 1)
    S[S > 1.0] = 1.0
    S[S < -1.0] = -1.0
    np.fill_diagonal(S, 1.0)
    return S


def huge_preprocess(x: Any) -> dict[str, Any]:
    """``.huge_preprocess`` (huge-internal.R:228).

    Returns the standardized correlation matrix ``S`` together with ``n``,
    ``d`` and the ``cov.input`` flag. Only the data-matrix route is ported
    (covariance input carries the same code path but is not used by the ``mb``
    pipeline, which always passes the raw transformed data).
    """
    x = np.asarray(x, dtype=float)
    n, d = x.shape
    cov_input = huge_is_covariance_input(x)
    if cov_input:
        raise NotImplementedError(
            "covariance-matrix input is not part of the ported mb pipeline"
        )
    S = huge_fast_cor(x)
    if not np.all(np.isfinite(S)):
        raise ValueError("Raw data x cannot produce a finite correlation matrix.")
    return {"S": S, "n": n, "d": d, "cov_input": cov_input}
