"""Meinshausen-Buhlmann neighbourhood selection.

Exact port of ``huge::mb`` (huge 2.0.1), i.e. the C++ routine ``huge::mb`` in
``huge/src/huge_core.cpp:679`` and its helpers ``mb_partial_residual``
(line 640) and ``mb_refine_active`` (line 651), together with the R wrapper
``huge.mb`` (huge/R/huge.mb.R:27).

The estimator solves, for every node ``m``, a nodewise Lasso over the
*correlation* matrix ``S`` (``huge.mb`` calls ``.huge_preprocess`` which returns
the correlation matrix via ``.huge_fast_cor``)::

    min_w  0.5 w' S w - S[m, :]' w + lambda ||w||_1      (w_m fixed to 0)

via active-set coordinate descent with the sequential strong rule
(``strong_thr = 2*lambda_i - lambda_{i-1}``), a warm start along the lambda
path, ``thol = 1e-4`` and ``max_iter = 10000``. The update

    r_j = S[m, j] + w_j - sum_k S[j, a_k] w_{a_k}
    w_j = threshold_l1(r_j, lambda)

is the source's exact form; the ``+ w_j`` cancels the diagonal term
``S[j, j] w_j`` because ``diag(S) = 1``.

The inner loops run on plain Python floats (``S.tolist()`` rows): this is the
same operation order as the C++ and is ~1.6x faster than the equivalent NumPy
element expressions because it avoids per-coordinate call overhead.
"""

from __future__ import annotations

from typing import Any

import numpy as np

# R: huge_core.cpp:683-684 (mb)
_THOL = 1e-4
_MAX_ITER = 10000


def _threshold_l1(x: float, thr: float) -> float:
    """``threshold_l1`` (huge_core.cpp:18)."""
    if x > thr:
        return x - thr
    if x < -thr:
        return x + thr
    return 0.0


def _mb_single_node(
    Srows: list[list[float]], m: int, lambdas: np.ndarray
) -> tuple[list[np.ndarray], list[np.ndarray], bool]:
    """Solve the nodewise Lasso for node ``m`` along the whole lambda path.

    Returns ``(predictors, coefficients, hit_max_iter)`` where the two lists are
    indexed by lambda and each holds the non-zero predictors and their values.
    """
    d = len(Srows)
    Sm = Srows[m]
    idx_i = [True] * d
    idx_i[m] = False
    w0 = [0.0] * d
    w1 = [0.0] * d
    grad_abs = [abs(v) for v in Sm]
    idx_a = [0] * d
    size_a = 0

    preds: list[np.ndarray] = []
    coefs: list[np.ndarray] = []
    hit_any = False

    for i in range(len(lambdas)):
        lam = float(lambdas[i])
        strong_thr = 2.0 * lam - float(lambdas[i - 1]) if i > 0 else 2.0 * lam
        strong = [idx_i[j] and grad_abs[j] > strong_thr for j in range(d)]

        iter_ext = 0
        certified = False
        while (not certified) and iter_ext < _MAX_ITER:
            gap_ext = 1
            while gap_ext != 0 and iter_ext < _MAX_ITER:
                size_a_prev = size_a
                for j in range(d):
                    if idx_i[j] and strong[j]:
                        # mb_partial_residual: S[m,j] - sum_k S[j,a_k] w[a_k]
                        r = Sm[j]
                        row = Srows[j]
                        for kk in range(size_a):
                            r -= row[idx_a[kk]] * w0[idx_a[kk]]
                        grad_abs[j] = abs(r)
                        wj = _threshold_l1(r, lam)
                        w1[j] = wj
                        if wj != 0.0:
                            idx_a[size_a] = j
                            size_a += 1
                            idx_i[j] = False
                        w0[j] = wj
                gap_ext = size_a - size_a_prev

                # mb_refine_active: cyclic CD over the active set
                gap_int = 1.0
                iter_int = 0
                while gap_int > _THOL and iter_int < _MAX_ITER:
                    change_sum = 0.0
                    coef_sum = 0.0
                    for jj in range(size_a):
                        widx = idx_a[jj]
                        row = Srows[widx]
                        r = Sm[widx] + w0[widx]
                        for kk in range(size_a):
                            r -= row[idx_a[kk]] * w0[idx_a[kk]]
                        wnew = _threshold_l1(r, lam)
                        w1[widx] = wnew
                        coef_sum += abs(wnew)
                        change_sum += abs(wnew - w0[widx])
                        w0[widx] = wnew
                    gap_int = change_sum / coef_sum if coef_sum > 0.0 else 0.0
                    iter_int += 1

                junk = 0
                for jj in range(size_a):
                    widx = idx_a[jj]
                    if w1[widx] == 0.0:
                        junk += 1
                        idx_i[widx] = True
                        strong[widx] = True
                    else:
                        idx_a[jj - junk] = widx
                size_a -= junk
                iter_ext += 1

            # Certification: full pass over every inactive coordinate.
            violations = 0
            for j in range(d):
                if idx_i[j]:
                    r = Sm[j]
                    row = Srows[j]
                    for kk in range(size_a):
                        r -= row[idx_a[kk]] * w0[idx_a[kk]]
                    grad_abs[j] = abs(r)
                    wj = _threshold_l1(r, lam)
                    w1[j] = wj
                    if wj != 0.0:
                        idx_a[size_a] = j
                        size_a += 1
                        idx_i[j] = False
                        strong[j] = True
                        violations += 1
                    w0[j] = wj
            certified = violations == 0
        if not certified:
            hit_any = True

        active = idx_a[:size_a]
        preds.append(np.array(active, dtype=np.intp))
        coefs.append(np.array([w1[j] for j in active]))

    return preds, coefs, hit_any


def mb(
    S: Any,
    lambdas: Any,
    sym: str = "or",
) -> tuple[list[np.ndarray], list[np.ndarray], bool]:
    """Meinshausen-Buhlmann graph estimation on a correlation matrix.

    Returns ``(betas, paths, hit_max_iter)``. ``betas[i][j, m]`` is the
    coefficient of predictor ``j`` in the regression of response ``m`` at lambda
    index ``i`` (matching ``huge.mb``'s ``beta``). ``paths[i]`` is the
    symmetrised 0/1 adjacency (``sign(|beta| + |beta|')`` for ``sym="or"``,
    ``sign(|beta| * |beta|')`` for ``"and"``).
    """
    S = np.asarray(S, dtype=float)
    if S.ndim != 2 or S.shape[0] != S.shape[1]:
        raise ValueError("S must be a square matrix")
    d = S.shape[0]
    lambdas = np.asarray(lambdas, dtype=float)
    if lambdas.ndim != 1 or lambdas.size == 0:
        raise ValueError("lambdas must be a non-empty 1-D sequence")
    if sym not in ("or", "and"):
        raise ValueError('sym must be exactly one of "or" or "and"')

    Srows = S.tolist()
    betas = [np.zeros((d, d)) for _ in range(lambdas.size)]
    hit = False
    for m in range(d):
        preds, coefs, node_hit = _mb_single_node(Srows, m, lambdas)
        hit = hit or node_hit
        for i in range(lambdas.size):
            pj = preds[i]
            if pj.size:
                betas[i][pj, m] = coefs[i]

    paths: list[np.ndarray] = []
    for B in betas:
        A = np.abs(B)
        A = A + A.T if sym == "or" else A * A.T
        paths.append(np.sign(A))
    return betas, paths, hit
