"""Input validation (mirrors the R ``.data.checks`` and huge's guards)."""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np


def validate_counts(data: Any) -> np.ndarray:
    """Validate an ``n_samples x n_features`` count matrix.

    Reproduces the R warnings from ``.data.checks`` (spiec-easi.R:35) and the
    hard errors from huge's ``.huge_validate_estimation_input``.
    """
    arr = np.asarray(data)
    if arr.ndim != 2:
        raise ValueError("data must be a 2-D (n_samples x n_features) matrix")
    if not np.issubdtype(arr.dtype, np.number):
        raise ValueError("data must be numeric")
    arr = arr.astype(float, copy=False)
    if arr.shape[0] < 2:
        raise ValueError("at least 2 samples are required")
    if arr.shape[1] < 2:
        raise ValueError("at least 2 features are required")
    if not np.all(np.isfinite(arr)):
        raise ValueError("data must contain only finite values (no NaN/Inf)")

    if bool(np.any(arr < 0)):
        warnings.warn(
            "Negative values detected, but raw counts are expected", stacklevel=2
        )
    row_sums = arr.sum(axis=1)
    if np.allclose(row_sums, np.ones_like(row_sums)):
        warnings.warn("Data is normalized, but raw counts are expected", stacklevel=2)

    constant = np.all(arr == arr[0:1, :], axis=0)
    if bool(np.any(constant)):
        raise ValueError("data contains a constant (zero-variance) column")
    return arr
