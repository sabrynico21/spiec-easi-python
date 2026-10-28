"""Equivalence tests against the committed R goldens.

The goldens in ``tests/golden/`` are produced by the R reference
(``tests/r_reference/generate_golden.R``). These tests load the input, inject the
R subsample indices, run the Python port and compare step by step. R is not
required.

Tolerances (documented in ``docs/EQUIVALENCE.md``; never widened to force a
pass):

- deterministic steps (transform, correlation, lambda grid): rtol 1e-12
- MB coefficients, StARS variability: rtol 1e-6 / atol 1e-8 (the MB coordinate
  descent stops at relative change 1e-4)
- edge sets and opt index: exact, except an edge whose coefficient magnitude is
  within ``BOUNDARY_ATOL`` in *both* implementations (a soft-threshold boundary
  tie at the machine-epsilon level).
"""

from __future__ import annotations

import json
import os
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from spieceasi.lambda_path import get_lam_path, get_max_cov
from spieceasi.mb import mb
from spieceasi.pipeline import spiec_easi
from spieceasi.preprocess import huge_fast_cor
from spieceasi.transforms import spiec_easi_norm

GOLDEN = Path(os.environ.get("SPIEC_EASI_GOLDEN_DIR", Path(__file__).parent / "golden"))
CASES = sorted(p.name for p in GOLDEN.iterdir() if p.is_dir())

TOL_DET = {"rtol": 1e-12, "atol": 1e-13}
TOL_COEF = {"rtol": 1e-6, "atol": 1e-8}
BOUNDARY_ATOL = 1e-8

# The QMP p>n case has rank(S)=19 < p=30. Its R/Python coefficients currently
# agree to 1.2e-14 absolute, but the design remains singular. Keep objective/KKT
# checks as mathematical invariants in addition to the pointwise comparison.
SINGULAR_DESIGN_CASES = {"qmp_n20_p30"}


def _lasso_objective(S: np.ndarray, b: np.ndarray, w: np.ndarray, lam: float) -> float:
    """The nodewise-Lasso objective huge minimises: 0.5 w'Sw - b'w + lam|w|_1."""
    return float(0.5 * w @ S @ w - b @ w + lam * np.abs(w).sum())


def _lasso_kkt(
    S: np.ndarray, b: np.ndarray, w: np.ndarray, m: int, lam: float
) -> tuple[float, float]:
    """Max KKT violation (active, inactive) for node ``m`` (excluding m itself)."""
    keep = np.arange(S.shape[0]) != m
    r = (b - S @ w)[keep]
    ww = w[keep]
    act = np.nonzero(ww)[0]
    v_act = (
        0.0 if act.size == 0 else float(np.max(np.abs(r[act] - lam * np.sign(ww[act]))))
    )
    inact = np.setdiff1d(np.arange(r.size), act)
    v_inact = (
        0.0
        if inact.size == 0
        else float(np.max(np.maximum(0.0, np.abs(r[inact]) - lam)))
    )
    return v_act, v_inact


_CACHE: dict[str, dict[str, Any]] = {}


def _dir(case: str) -> Path:
    return GOLDEN / case


def _counts(case: str) -> np.ndarray:
    return np.loadtxt(_dir(case) / "counts.tsv")


def _params(case: str) -> dict[str, Any]:
    return json.loads((_dir(case) / "params.json").read_text())


def _beta_long(path: Path) -> dict[tuple[int, int, int], float]:
    out: dict[tuple[int, int, int], float] = {}
    for line in path.read_text().splitlines():
        if line.strip():
            li, m, j, v = line.split("\t")
            out[(int(li), int(m), int(j))] = float(v)
    return out


def _adj_long(path: Path) -> dict[int, set[tuple[int, int]]]:
    out: dict[int, set[tuple[int, int]]] = {}
    for line in path.read_text().splitlines():
        if line.strip():
            li, i, j = line.split("\t")
            out.setdefault(int(li), set()).add((int(i), int(j)))
    return out


def _rep_adj_long(path: Path) -> dict[tuple[int, int], set[tuple[int, int]]]:
    out: dict[tuple[int, int], set[tuple[int, int]]] = {}
    for line in path.read_text().splitlines():
        if line.strip():
            k, li, i, j = line.split("\t")
            out.setdefault((int(k), int(li)), set()).add((int(i), int(j)))
    return out


def _edges_from_beta(beta: np.ndarray) -> set[tuple[int, int]]:
    d = beta.shape[0]
    A = np.abs(beta) + np.abs(beta).T
    return {(a, b) for a in range(d) for b in range(a + 1, d) if A[a, b] != 0}


def _run(case: str) -> dict[str, Any]:
    """Run the pipeline for a case, injecting the R subsample indices (cached)."""
    if case in _CACHE:
        return _CACHE[case]
    params = _params(case)
    idx = np.loadtxt(_dir(case) / "subsample_indices.tsv", dtype=int) - 1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = spiec_easi(
            _counts(case),
            method="mb",
            sel_criterion="stars",
            nlambda=params["nlambda"],
            lambda_min_ratio=params["lambda_min_ratio"],
            lambda_log=params["lambda_log"],
            pulsar_params={
                "rep_num": params["rep_num"],
                "thresh": params["thresh"],
                "seed": params["seed"],
            },
            subsample_indices=idx,
            keep_replica_betas=True,
            n_jobs=int(os.environ.get("SPIEC_EASI_TEST_JOBS", "4")),
        )
    _CACHE[case] = res
    return res


def _check_edges(
    label: str,
    got: set[tuple[int, int]],
    ref: set[tuple[int, int]],
    got_beta: np.ndarray | None,
    ref_beta: dict[tuple[int, int], float] | None,
) -> None:
    """Exact edge comparison, allowing only machine-epsilon boundary ties."""
    diff = got ^ ref
    for a, b in sorted(diff):
        gm = 0.0
        pm = 0.0
        if ref_beta is not None:
            gm = max(abs(ref_beta.get((a, b), 0.0)), abs(ref_beta.get((b, a), 0.0)))
        if got_beta is not None:
            pm = max(abs(got_beta[a, b]), abs(got_beta[b, a]))
        if max(gm, pm) > BOUNDARY_ATOL:
            raise AssertionError(
                f"{label}: edge {(a, b)} diverges (gold|coef|={gm:.3e}, "
                f"port|coef|={pm:.3e}) - not a boundary tie"
            )


# ---------------------------------------------------------------------------
# Deterministic steps (fast)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("case", CASES)
def test_transform(case: str) -> None:
    assert _params(case)["dataset"] in {
        "QMP",
        "American Gut Project Round 1",
        "GlobalPatterns / Caporaso et al. 2011",
    }
    np.testing.assert_allclose(
        spiec_easi_norm(_counts(case)),
        np.loadtxt(_dir(case) / "clr.tsv"),
        **TOL_DET,  # type: ignore[arg-type]
    )


@pytest.mark.parametrize("case", CASES)
def test_correlation(case: str) -> None:
    X = spiec_easi_norm(_counts(case))
    np.testing.assert_allclose(
        np.corrcoef(X, rowvar=False),
        np.loadtxt(_dir(case) / "cor.tsv"),
        **TOL_DET,  # type: ignore[arg-type]
    )


@pytest.mark.parametrize("case", CASES)
def test_lambda_max_and_path(case: str) -> None:
    params = _params(case)
    X = spiec_easi_norm(_counts(case))
    lmax = get_max_cov(np.corrcoef(X, rowvar=False))
    np.testing.assert_allclose(
        lmax,
        float(np.loadtxt(_dir(case) / "lambda_max.txt")),
        **TOL_DET,  # type: ignore[arg-type]
    )
    np.testing.assert_allclose(
        get_lam_path(
            lmax,
            lmax * params["lambda_min_ratio"],
            params["nlambda"],
            log=params["lambda_log"],
        ),
        np.loadtxt(_dir(case) / "lambda_path.txt"),
        **TOL_DET,  # type: ignore[arg-type]
    )


# ---------------------------------------------------------------------------
# Full-data MB solver (fast)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("case", CASES)
def test_full_beta(case: str) -> None:
    X = spiec_easi_norm(_counts(case))
    lams = np.loadtxt(_dir(case) / "lambda_path.txt")
    betas, _, _ = mb(huge_fast_cor(X), lams, sym="or")
    got = np.stack([betas[i].T for i in range(len(lams))])  # [lambda, resp, pred]
    ref = np.zeros_like(got)
    for (li, m, j), v in _beta_long(_dir(case) / "full_beta.tsv").items():
        ref[li, m, j] = v
    np.testing.assert_allclose(got, ref, **TOL_COEF)  # type: ignore[arg-type]


@pytest.mark.parametrize("case", CASES)
def test_full_adjacency(case: str) -> None:
    X = spiec_easi_norm(_counts(case))
    lams = np.loadtxt(_dir(case) / "lambda_path.txt")
    betas, paths, _ = mb(huge_fast_cor(X), lams, sym="or")
    gold_adj = _adj_long(_dir(case) / "full_adj.tsv")
    gold_beta_rows = _beta_long(_dir(case) / "full_beta.tsv")
    p = X.shape[1]
    for i in range(len(lams)):
        got = {(a, b) for a in range(p) for b in range(a + 1, p) if paths[i][a, b] != 0}
        gold_edges = gold_adj.get(i, set())
        ref_beta = {
            (a, b): max(
                abs(gold_beta_rows.get((i, a, b), 0.0)),
                abs(gold_beta_rows.get((i, b, a), 0.0)),
            )
            for (a, b) in (got | gold_edges)
        }
        _check_edges(f"{case} full lambda {i}", got, gold_edges, betas[i], ref_beta)


# ---------------------------------------------------------------------------
# Singular (p > n) case: coefficient and optimization invariants
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("case", sorted(SINGULAR_DESIGN_CASES))
def test_singular_case_invariants(case: str) -> None:
    """Check KKT/objective invariants for the real rank-deficient QMP design.

    Checks, for every node and every lambda >= 1 (lambda 0 is the machine-ulp
    boundary, covered by ``test_full_adjacency``):

    * identical support;
    * the objective (as huge minimises it) agrees between R and the port;
    * both solutions satisfy KKT within the solver's own ``thol = 1e-4``;
    * if coefficients differ, ``S_AA @ (beta_R - beta_Py)`` is at the solver's
      residual level (a null-space difference rather than a different optimum).
    """
    X = spiec_easi_norm(_counts(case))
    lams = np.loadtxt(_dir(case) / "lambda_path.txt")
    S = huge_fast_cor(X)
    d = S.shape[0]
    assert np.linalg.matrix_rank(S) < d
    betas_py, _, _ = mb(S, lams, sym="or")

    rows = _beta_long(_dir(case) / "full_beta.tsv")
    beta_R = np.zeros((len(lams), d, d))
    for (li, m, j), v in rows.items():
        beta_R[li, m, j] = v
    beta_Py = np.stack([b.T for b in betas_py])

    worst_obj = 0.0
    worst_kkt = 0.0
    worst_nullres = 0.0
    checked = 0
    for li in range(1, len(lams)):
        for m in range(d):
            b = S[m].copy()
            wR = beta_R[li, m].copy()
            wP = beta_Py[li, m].copy()
            wR[m] = wP[m] = 0.0
            lam = float(lams[li])
            supR = set(np.nonzero(wR)[0].tolist())
            supP = set(np.nonzero(wP)[0].tolist())
            assert supR == supP, f"{case} lambda {li} node {m}: support differs"
            checked += 1
            fR = _lasso_objective(S, b, wR, lam)
            fP = _lasso_objective(S, b, wP, lam)
            worst_obj = max(worst_obj, abs(fP - fR) / max(abs(fR), 1e-30))
            for w in (wR, wP):
                va, vi = _lasso_kkt(S, b, w, m, lam)
                assert vi == 0.0, f"{case} lambda {li} node {m}: KKT inactive violated"
                worst_kkt = max(worst_kkt, va)
            if np.max(np.abs(wR - wP)) < 1e-12:
                continue
            A = sorted(supR | supP)
            SAA = S[np.ix_(A, A)]
            worst_nullres = max(
                worst_nullres, float(np.max(np.abs(SAA @ (wR - wP)[A])))
            )

    assert checked > 0
    assert worst_obj <= 1e-6, f"objective mismatch {worst_obj:.2e}"
    assert worst_kkt <= 1e-4, f"KKT violation {worst_kkt:.2e} exceeds thol"
    assert worst_nullres <= 1e-4, f"S_AA dbeta {worst_nullres:.2e} too large"


# ---------------------------------------------------------------------------
# Pipeline (slow): opt index, stability, refit, per-replica adjacency
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("case", CASES)
def test_pipeline_opt_index(case: str) -> None:
    res = _run(case)
    assert res["opt_index"] == int(np.loadtxt(_dir(case) / "opt_index.txt"))


@pytest.mark.slow
@pytest.mark.parametrize("case", CASES)
def test_pipeline_stability(case: str) -> None:
    res = _run(case)
    np.testing.assert_allclose(
        res["stability"],
        np.loadtxt(_dir(case) / "stars_summary_monotone.txt"),
        **TOL_COEF,  # type: ignore[arg-type]
    )


@pytest.mark.slow
@pytest.mark.parametrize("case", CASES)
def test_pipeline_refit(case: str) -> None:
    res = _run(case)
    gold = np.loadtxt(_dir(case) / "pipeline_refit.tsv")
    p = _counts(case).shape[1]
    got = _edges_from_beta(res["refit"])
    ref = {(a, b) for a in range(p) for b in range(a + 1, p) if gold[a, b] != 0}
    opt = res["opt_index"] - 1
    gold_beta_rows = _beta_long(_dir(case) / "full_beta.tsv")
    ref_beta = {
        (a, b): max(
            abs(gold_beta_rows.get((opt, a, b), 0.0)),
            abs(gold_beta_rows.get((opt, b, a), 0.0)),
        )
        for (a, b) in (got | ref)
    }
    _check_edges(f"{case} refit", got, ref, res["beta"][opt], ref_beta)


@pytest.mark.slow
@pytest.mark.parametrize("case", CASES)
def test_per_replica_adjacency(case: str) -> None:
    res = _run(case)
    gold = _rep_adj_long(_dir(case) / "rep_adj_path.tsv")
    for k, betas_k in enumerate(res["_replica_betas"]):
        for i, beta in enumerate(betas_k):
            got = _edges_from_beta(beta)
            ref = gold.get((k, i), set())
            _check_edges(f"{case} replica {k} lambda {i}", got, ref, beta, None)


@pytest.mark.slow
@pytest.mark.parametrize("case", ["qmp_n106_p20", "qmp_n20_p30"])
def test_jobs_determinism(case: str) -> None:
    """Serial and pooled runs must be identical (multiprocessing is exact)."""
    params = _params(case)
    idx = np.loadtxt(_dir(case) / "subsample_indices.tsv", dtype=int) - 1

    def run(n_jobs: int) -> dict[str, Any]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return spiec_easi(
                _counts(case),
                nlambda=params["nlambda"],
                lambda_min_ratio=params["lambda_min_ratio"],
                pulsar_params={
                    "rep_num": params["rep_num"],
                    "thresh": params["thresh"],
                },
                subsample_indices=idx,
                n_jobs=n_jobs,
            )

    a = run(1)
    b = run(4)
    assert a["opt_index"] == b["opt_index"]
    np.testing.assert_array_equal(a["refit"], b["refit"])
    np.testing.assert_array_equal(a["stability"], b["stability"])
    for merge_a, merge_b in zip(a["merge"], b["merge"], strict=True):
        np.testing.assert_array_equal(merge_a, merge_b)
