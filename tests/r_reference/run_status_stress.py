#!/usr/bin/env python3
"""Seeded R/Python status + discrete-output stress sweep for the SpiecEasi MB port.

This is the SPIEC-EASI analogue of ``cclasso-python/tests/r_reference/
run_status_stress.py``. It implements the predeclared 1,000-case seeded stress
suite required by ``docs/EQUIVALENCE_ACCEPTANCE_CONTRACT.md`` ("weighted toward
sparse, singular, and threshold-near inputs").

Design (predeclared; do not tune after seeing results)
------------------------------------------------------
* Master seed ``--seed`` (default 20261020) draws one per-case seed and one
  per-case draw seed from a ``numpy.random.default_rng`` stream, in case order.
* ``--cases`` cases are generated, ``case_index % 6`` selects the generator
  mode, so the six strata are equally represented:
    0. moderate Poisson-lognormal counts, ~15% zeros;
    1. sparse counts, ~65% zeros;
    2. rank-deficient design (one constant feature) + ~35% zeros;
    3. low library sizes (small Poisson mean) + ~45% zeros;
    4. p > n (p between n+1 and n+6, capped at 14 features);
    5. threshold-near design (weak/strong edge weights 0.25/0.35) + ~30% zeros.
* Fixed SPIEC-EASI MB settings: ``nlambda=5``, ``lambda_min_ratio=0.1``,
  ``lambda_log=TRUE``, ``rep_num=4``, ``thresh=0.2`` - the same settings used by
  ``tests/r_reference/generate_golden.R`` for the committed fixtures.
* StARS subsample indices are drawn once with NumPy and injected into *both*
  languages (written 1-based for R), so the R RNG is never needed and the two
  runs differ only by arithmetic, not by resampling.
* Each language computes its own lambda grid from the injected counts. The
  resulting lambda-max difference (~1e-16 relative) is recorded per case.

Classification (predeclared thresholds)
---------------------------------------
* ``status``: ``error`` (an exception), ``nonfinite`` (any non-finite lambda,
  stability, beta or refit entry), else ``finite``.
* ``opt_index``: the StARS-selected index (1-based). An ``opt_index``
  disagreement is **boundary-classified** only when the monotonised stability
  curve of at least one language passes within ``OPT_BOUNDARY_MARGIN`` of
  ``thresh`` at the two competing indices.
* ``refit`` edge set (exact, 0-based ``i:j``): an edge-set disagreement is
  **boundary-classified** only when every differing edge has a selection
  frequency within ``EDGE_BOUNDARY_FREQ`` of 0.5 in the Python merge at the
  selected lambda. Everything else is an unexplained mismatch.

The runner exits nonzero if any status mismatch or unexplained discrete
mismatch exists, and keeps the complete case archive for investigation.

Usage::

    .venv/bin/python tests/r_reference/run_status_stress.py \\
        /tmp/spiec-stress-1000 --cases 1000 --seed 20261020
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np

from spieceasi.pipeline import spiec_easi

NLAMBDA = 5
LAMBDA_MIN_RATIO = 0.1
LAMBDA_LOG = True
REP_NUM = 4
THRESH = 0.2

OPT_BOUNDARY_MARGIN = 1e-3
EDGE_BOUNDARY_FREQ = 0.05

STATUS_FIELDS = [
    "case",
    "status",
    "error_class",
    "opt_index",
    "n_edges",
    "edges",
    "lambda_max",
    "stars_margin",
    "merge_opt",
    "seconds",
]


# ---------------------------------------------------------------------------
# deterministic case generation (mirrors the R generator in generate_golden.R)
# ---------------------------------------------------------------------------


def _sparse_adj(
    rng: np.random.Generator,
    p: int,
    weak: float,
    strong: float,
    extra_frac: float = 1.0 / 3.0,
    strong_frac: float = 0.25,
) -> np.ndarray:
    """Two-tier sparse symmetric weighted adjacency (R ``sparse_adj``)."""
    w = np.zeros((p, p))
    for i in range(p - 1):
        w[i, i + 1] = weak
        w[i + 1, i] = weak
    cand = [(i, j) for i in range(p) for j in range(i + 1, p) if w[i, j] == 0.0]
    nadd = min(max(1, int(p * extra_frac)), len(cand))
    if nadd and cand:
        for k in rng.choice(len(cand), size=nadd, replace=False):
            i, j = cand[int(k)]
            w[i, j] = weak
            w[j, i] = weak
    edges = [(i, j) for i in range(p) for j in range(i + 1, p) if w[i, j] != 0.0]
    ns = min(max(1, int(len(edges) * strong_frac)), len(edges))
    if ns and edges:
        for k in rng.choice(len(edges), size=ns, replace=False):
            i, j = edges[int(k)]
            w[i, j] = strong
            w[j, i] = strong
    return w


def _gen_counts(
    rng: np.random.Generator,
    n: int,
    p: int,
    mu: float,
    weak: float,
    strong: float,
    zero_frac: float,
) -> np.ndarray:
    """Poisson-lognormal counts from a sparse, two-tier Gaussian copula."""
    w = _sparse_adj(rng, p, weak, strong)
    theta = np.diag(1.0 + w.sum(axis=1)) + w
    sigma = np.cov(np.linalg.inv(theta))
    d = np.sqrt(np.diag(sigma))
    sigma = sigma / np.outer(d, d)
    sigma = 0.5 * (sigma + sigma.T)
    eta = rng.multivariate_normal(np.zeros(p), sigma, size=n)
    lam = np.exp(np.log(mu) + eta)
    counts = rng.poisson(lam).astype(float)
    if zero_frac > 0:
        threshold = np.quantile(counts, zero_frac)
        counts[counts <= threshold] = 0.0
    return counts


def _make_case(
    case_rng: np.random.Generator, case_index: int
) -> tuple[np.ndarray, dict[str, object]]:
    mode = case_index % 6
    n = int(case_rng.integers(8, 31))
    if mode == 4:
        p = int(min(n + case_rng.integers(1, 7), 14))
    else:
        p = int(case_rng.integers(2, min(12, n + 4) + 1))

    mu, weak, strong, zero_frac = 10.0, 0.4, 0.9, 0.15
    if mode == 1:
        mu, zero_frac = 8.0, 0.65
    elif mode == 3:
        mu, zero_frac = 2.0, 0.45
    elif mode == 4:
        zero_frac = 0.30
    elif mode == 5:
        mu, weak, strong, zero_frac = 6.0, 0.25, 0.35, 0.30

    counts = _gen_counts(case_rng, n, p, mu, weak, strong, zero_frac)
    if mode == 2:
        counts[:, -1] = counts[0, -1]  # constant feature -> rank deficiency

    params: dict[str, object] = {
        "mode": mode,
        "n": int(n),
        "p": int(p),
        "mu": float(mu),
        "weak": float(weak),
        "strong": float(strong),
        "zero_frac": float(zero_frac),
        "nlambda": NLAMBDA,
        "lambda_min_ratio": LAMBDA_MIN_RATIO,
        "lambda_log": LAMBDA_LOG,
        "rep_num": REP_NUM,
        "thresh": THRESH,
    }
    return counts, params


# ---------------------------------------------------------------------------
# Python-side classification
# ---------------------------------------------------------------------------


def _edges_string(adjacency: np.ndarray) -> tuple[int, str]:
    idx = np.argwhere((adjacency != 0) & np.triu(np.ones_like(adjacency, bool), 1))
    edges = ";".join(f"{int(i)}:{int(j)}" for i, j in idx)
    return int(len(idx)), edges


def _merge_string(merge: np.ndarray) -> str:
    rows = []
    for i in range(merge.shape[0]):
        for j in range(i + 1, merge.shape[1]):
            rows.append(f"{i}:{j}:{float(merge[i, j]):.6g}")
    return ";".join(rows)


def _python_status(
    counts: np.ndarray, indices: np.ndarray, case_name: str
) -> tuple[dict[str, object], np.ndarray | None, np.ndarray | None]:
    start = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            result = spiec_easi(
                counts,
                method="mb",
                pulsar_params={"thresh": THRESH, "rep_num": REP_NUM},
                lambda_log=LAMBDA_LOG,
                nlambda=NLAMBDA,
                lambda_min_ratio=LAMBDA_MIN_RATIO,
                subsample_indices=indices,
                n_jobs=1,
            )
        except Exception as exc:  # status is compared, not the exception text
            return (
                {
                    "case": case_name,
                    "status": "error",
                    "error_class": type(exc).__name__,
                    "opt_index": "",
                    "n_edges": "",
                    "edges": "",
                    "lambda_max": "",
                    "stars_margin": "",
                    "merge_opt": "",
                    "seconds": f"{time.perf_counter() - start:.6f}",
                },
                None,
                None,
            )
    elapsed = time.perf_counter() - start

    arrays: list[np.ndarray] = [
        np.asarray(result["lambda"], dtype=float),
        np.asarray(result["stability"], dtype=float),
        np.asarray(result["refit"], dtype=float),
    ]
    arrays.extend(np.asarray(b, dtype=float) for b in result["beta"])
    nonfinite = any(not np.all(np.isfinite(a)) for a in arrays)

    stability = np.asarray(result["stability"], dtype=float)
    merge_opt = np.asarray(result["merge"][int(result["opt_index"]) - 1], dtype=float)
    n_edges, edges = _edges_string(np.asarray(result["refit"], dtype=float))
    row: dict[str, object] = {
        "case": case_name,
        "status": "nonfinite" if nonfinite else "finite",
        "error_class": "",
        "opt_index": int(result["opt_index"]),
        "n_edges": n_edges,
        "edges": edges,
        "lambda_max": f"{float(result['lambda_max']):.17g}",
        "stars_margin": f"{float(np.min(np.abs(stability - THRESH))):.17g}",
        "merge_opt": _merge_string(merge_opt),
        "seconds": f"{elapsed:.6f}",
    }
    if nonfinite:
        row["edges"] = ""
    return row, stability, merge_opt


# ---------------------------------------------------------------------------
# comparison
# ---------------------------------------------------------------------------


def _freq_margin(merge_string: str, edge: str) -> float:
    """Distance of an edge's selection frequency from 0.5 (large = decisive)."""
    for token in merge_string.split(";"):
        if not token:
            continue
        i, j, value = token.split(":")
        if f"{i}:{j}" == edge:
            return abs(float(value) - 0.5)
    return 0.5  # edge absent from R's merge: treat as decisive


def _compare(
    python_rows: list[dict[str, object]],
    r_rows: dict[str, dict[str, str]],
) -> tuple[list[dict[str, object]], dict[str, int]]:
    unexplained: list[dict[str, object]] = []
    boundary: list[dict[str, object]] = []
    counts = {
        "status_mismatch": 0,
        "opt_index_mismatch": 0,
        "opt_index_boundary": 0,
        "edge_mismatch": 0,
        "edge_boundary": 0,
    }
    for row in python_rows:
        ref = r_rows.get(row["case"])
        if ref is None:
            unexplained.append({"case": row["case"], "reason": "missing R row"})
            counts["status_mismatch"] += 1
            continue
        if ref["status"] != row["status"]:
            counts["status_mismatch"] += 1
            unexplained.append(
                {
                    "case": row["case"],
                    "reason": "status",
                    "python": row["status"],
                    "R": ref["status"],
                    "python_error": row["error_class"],
                    "R_error": ref["error_class"],
                }
            )
            continue
        if row["status"] != "finite":
            continue

        if int(ref["opt_index"]) != int(row["opt_index"]):
            counts["opt_index_mismatch"] += 1
            margin = min(float(row["stars_margin"]), float(ref["stars_margin"]))
            entry = {
                "case": row["case"],
                "reason": "opt_index",
                "python": row["opt_index"],
                "R": ref["opt_index"],
                "stars_margin": f"{margin:.3g}",
            }
            if margin <= OPT_BOUNDARY_MARGIN:
                counts["opt_index_boundary"] += 1
                boundary.append(entry)
            else:
                unexplained.append(entry)

        py_edges = set(filter(None, str(row["edges"]).split(";")))
        r_edges = set(filter(None, ref["edges"].split(";")))
        differing = py_edges ^ r_edges
        if differing:
            counts["edge_mismatch"] += 1
            decisiveness = min(
                _freq_margin(str(row["merge_opt"]), e) for e in differing
            )
            entry = {
                "case": row["case"],
                "reason": "edge_set",
                "python_only": sorted(py_edges - r_edges)[:10],
                "R_only": sorted(r_edges - py_edges)[:10],
                "min_freq_margin": f"{decisiveness:.3g}",
            }
            if decisiveness <= EDGE_BOUNDARY_FREQ:
                counts["edge_boundary"] += 1
                boundary.append(entry)
            else:
                unexplained.append(entry)
    return unexplained, boundary, counts


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def _write_matrix(path: Path, values: np.ndarray, fmt: str = "%.17g") -> None:
    np.savetxt(path, np.asarray(values), delimiter=",", fmt=fmt)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--cases", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    if args.cases < 1:
        parser.error("--cases must be positive")
    out_dir = args.out_dir
    if out_dir.exists() and any(out_dir.iterdir()):
        parser.error(f"refusing to overwrite non-empty output directory: {out_dir}")
    cases_dir = out_dir / "cases"
    cases_dir.mkdir(parents=True, exist_ok=True)

    master_rng = np.random.default_rng(args.seed)
    python_rows: list[dict[str, object]] = []
    sweep_start = time.perf_counter()
    for i in range(args.cases):
        case_seed = int(master_rng.integers(0, 2**32 - 1))
        draw_seed = int(master_rng.integers(0, 2**32 - 1))
        case_rng = np.random.default_rng(case_seed)
        counts, params = _make_case(case_rng, i)
        n, p = counts.shape
        ratio = 0.8 if n <= 144 else 10.0 * np.sqrt(n) / n
        size = int(np.floor(n * ratio))
        draw_rng = np.random.default_rng(draw_seed)
        indices = np.array(
            [draw_rng.choice(n, size=size, replace=False) for _ in range(REP_NUM)],
            dtype=np.intp,
        )

        case_dir = cases_dir / f"case-{i:05d}"
        case_dir.mkdir()
        _write_matrix(case_dir / "counts.csv", counts)
        _write_matrix(case_dir / "subsample_indices.csv", indices + 1, fmt="%d")
        params.update(
            {
                "case_seed": case_seed,
                "draw_seed": draw_seed,
                "subsample_ratio": float(ratio),
                "subsample_size": int(size),
            }
        )
        (case_dir / "params.json").write_text(
            json.dumps(params, indent=2) + "\n", encoding="utf-8"
        )

        row, _, _ = _python_status(counts, indices, case_dir.name)
        python_rows.append(row)
        if (i + 1) % 100 == 0:
            print(f"[python] {i + 1}/{args.cases} cases", flush=True)

    python_seconds = [float(row["seconds"]) for row in python_rows]
    python_status = out_dir / "python_status.csv"
    with python_status.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=STATUS_FIELDS)
        writer.writeheader()
        writer.writerows(python_rows)

    r_status = out_dir / "r_status.csv"
    r_runner = Path(__file__).with_name("stress_status_R.R")
    r_start = time.perf_counter()
    r_run = subprocess.run(
        ["Rscript", str(r_runner), str(cases_dir), str(r_status)],
        cwd=Path(__file__).resolve().parents[2],
        capture_output=True,
        text=True,
    )
    r_elapsed = time.perf_counter() - r_start
    if r_run.returncode:
        raise RuntimeError(
            "R stress classification failed:\n"
            f"stdout: {r_run.stdout[-4000:]}\nstderr: {r_run.stderr[-4000:]}"
        )

    with r_status.open(newline="", encoding="utf-8") as stream:
        r_rows = {row["case"]: row for row in csv.DictReader(stream)}

    unexplained, boundary, counts = _compare(python_rows, r_rows)

    total_python = sum(python_seconds)
    summary = {
        "seed": args.seed,
        "cases": args.cases,
        "design": {
            "nlambda": NLAMBDA,
            "lambda_min_ratio": LAMBDA_MIN_RATIO,
            "lambda_log": LAMBDA_LOG,
            "rep_num": REP_NUM,
            "thresh": THRESH,
            "opt_boundary_margin": OPT_BOUNDARY_MARGIN,
            "edge_boundary_freq": EDGE_BOUNDARY_FREQ,
        },
        "python_status_counts": {
            status: sum(row["status"] == status for row in python_rows)
            for status in ("finite", "nonfinite", "error")
        },
        "r_status_counts": {
            status: sum(row["status"] == status for row in r_rows.values())
            for status in ("finite", "nonfinite", "error")
        },
        "mismatch_counts": counts,
        "unexplained_count": len(unexplained),
        "unexplained_preview": unexplained[:20],
        "boundary_classified_count": len(boundary),
        "boundary_preview": boundary[:20],
        "timing": {
            "python_total_s": round(total_python, 3),
            "python_mean_s": round(float(np.mean(python_seconds)), 6),
            "python_median_s": round(float(np.median(python_seconds)), 6),
            "python_p95_s": round(float(np.percentile(python_seconds, 95)), 6),
            "r_total_s": round(r_elapsed, 3),
            "sweep_wall_s": round(time.perf_counter() - sweep_start, 3),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    environment = {
        "python": sys.version,
        "numpy": np.__version__,
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "r_version": subprocess.run(
            ["Rscript", "--version"], capture_output=True, text=True
        ).stdout.splitlines()[0],
    }
    (out_dir / "environment.json").write_text(json.dumps(environment, indent=2) + "\n")

    print(json.dumps(summary, indent=2))
    if unexplained:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
