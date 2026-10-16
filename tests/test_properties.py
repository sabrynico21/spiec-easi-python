"""Structural properties and synthetic-recovery tests."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

from spieceasi.pipeline import spiec_easi
from spieceasi.transforms import clr

GOLDEN = Path(__file__).parent / "golden"


@pytest.mark.slow
@pytest.mark.parametrize("case", ["p20_n300", "p50_n200"])
def test_refit_adjacency_structure(case: str) -> None:
    counts = np.loadtxt(GOLDEN / case / "counts.tsv")
    params = json.loads((GOLDEN / case / "params.json").read_text())
    idx = np.loadtxt(GOLDEN / case / "subsample_indices.tsv", dtype=int) - 1
    res = spiec_easi(
        counts,
        nlambda=params["nlambda"],
        lambda_min_ratio=params["lambda_min_ratio"],
        pulsar_params={"rep_num": params["rep_num"], "thresh": params["thresh"]},
        subsample_indices=idx,
        n_jobs=int(os.environ.get("SPIEC_EASI_TEST_JOBS", "4")),
    )
    A = res["refit"]
    assert np.allclose(A, A.T)  # symmetric
    assert np.all(np.diag(A) == 0)  # zero diagonal
    assert np.all((A == 0) | (A == 1))  # 0/1 adjacency
    # merge (edge frequencies) in [0, 1] and symmetric
    M = res["merge"][res["opt_index"] - 1]
    assert np.allclose(M, M.T)
    assert np.all(M >= -1e-12) and np.all(M <= 1 + 1e-12)


def test_clr_handles_zeros_and_round_trips_positive() -> None:
    x = np.array([0.0, 1.0, 4.0, 0.0, 9.0])
    out = clr(x)
    assert out[0] == 0.0 and out[3] == 0.0  # zeros stay zero
    assert np.all(np.isfinite(out))
    # clr of a constant positive vector is zero
    assert np.allclose(clr(np.full(5, 3.0)), 0.0)


def _simulate_sparse(
    n: int, p: int, seed: int, mu: float = 12.0, weak: float = 0.4, strong: float = 0.8
) -> tuple[np.ndarray, set]:
    """Sparse two-tier Gaussian copula -> Poisson-lognormal counts.

    Mirrors the golden generator: a nearest-neighbour chain plus extra edges
    (``weak``), a subset promoted to ``strong``, and an *empirical* Gaussian draw
    (unit-variance whitening, like ``rmvnorm(empirical = TRUE)``). Returns the
    counts and the true edge set.
    """
    rng = np.random.default_rng(seed)
    A = np.zeros((p, p))
    for i in range(p - 1):
        A[i, i + 1] = A[i + 1, i] = weak
    cand = [(a, b) for a in range(p) for b in range(a + 1, p) if A[a, b] == 0]
    for idx in rng.choice(len(cand), size=max(1, int(p / 3)), replace=False):
        a, b = cand[int(idx)]
        A[a, b] = A[b, a] = weak
    edges = [(a, b) for a in range(p) for b in range(a + 1, p) if A[a, b] != 0]
    for idx in rng.choice(
        len(edges), size=max(1, int(len(edges) * 0.3)), replace=False
    ):
        a, b = edges[int(idx)]
        A[a, b] = A[b, a] = strong

    sigma = np.eye(p) + A
    eig = np.linalg.eigvalsh(sigma)
    if eig.min() <= 1e-6:
        sigma = sigma + (abs(eig.min()) + 1e-3) * np.eye(p)
    sd = np.sqrt(np.diag(sigma))
    sigma = sigma / np.outer(sd, sd)

    z = rng.standard_normal((n, p))
    zc = z - z.mean(0)
    u, _, vt = np.linalg.svd(zc, full_matrices=False)
    zw = u @ vt
    zw = zw / zw.std(axis=0, ddof=1)  # unit variance (rmvnorm empirical=TRUE)
    eta = zw @ np.linalg.cholesky(sigma + 1e-10 * np.eye(p)).T
    counts = rng.poisson(np.exp(np.log(mu) + eta)).astype(float) + 1.0
    return counts, {(min(a, b), max(a, b)) for a, b in edges}


@pytest.mark.slow
def test_recovery_on_known_graph() -> None:
    counts, true_edges = _simulate_sparse(n=300, p=20, seed=11)
    res = spiec_easi(
        counts,
        nlambda=20,
        pulsar_params={"rep_num": 10, "seed": 7},
        n_jobs=int(os.environ.get("SPIEC_EASI_TEST_JOBS", "4")),
    )
    A = res["refit"]
    p = counts.shape[1]
    got = {(a, b) for a in range(p) for b in range(a + 1, p) if A[a, b] != 0}
    assert got, "StARS selected an empty graph on a strong sparse signal"
    tp = len(got & true_edges)
    precision = tp / len(got)
    recall = tp / len(true_edges)
    # StARS selects the sparsest stable graph: high precision, lower recall.
    assert precision >= 0.7, f"precision {precision:.2f} too low"
    assert recall >= 0.2, f"recall {recall:.2f} too low"
