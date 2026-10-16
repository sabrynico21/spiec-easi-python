"""Command-line interface for the SPIEC-EASI ``mb`` pipeline.

Reads a JSON request from ``--input FILE`` (or stdin), writes a JSON response to
stdout. Errors go to stderr with a non-zero exit code; nothing else is ever
written to stdout, so the output can be parsed directly by the Electron
backend. See ``README.md`` for the schema.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any

import numpy as np

from .pipeline import spiec_easi


def _edges_from_adjacency(
    adj: np.ndarray, weights: np.ndarray | None
) -> list[list[float]]:
    d = adj.shape[0]
    out: list[list[float]] = []
    for a in range(d):
        for b in range(a + 1, d):
            if adj[a, b] != 0:
                w = float(weights[a, b]) if weights is not None else 1.0
                out.append([a, b, w])
    return out


def run_request(request: dict[str, Any]) -> dict[str, Any]:
    """Execute one JSON request and return the response dict."""
    if "counts" not in request:
        raise ValueError("request must contain a 'counts' field")
    counts = np.asarray(request["counts"], dtype=float)
    params = dict(request.get("params", {}))

    method = params.pop("method", "mb")
    sel_criterion = params.pop("sel_criterion", "stars")
    lambda_log = params.pop("lambda_log", True)
    nlambda = int(params.pop("nlambda", 20))
    lambda_min_ratio = float(params.pop("lambda_min_ratio", 1e-3))
    lambda_max = params.pop("lambda_max", None)
    sym = params.pop("sym", "or")
    pulsar_params = params.pop("pulsar_params", None)
    subsample_indices = params.pop("subsample_indices", None)
    if params:
        raise ValueError(f"Unrecognized params: {sorted(params)}")

    result = spiec_easi(
        counts,
        method=method,
        sel_criterion=sel_criterion,
        verbose=False,
        pulsar_params=pulsar_params,
        lambda_log=lambda_log,
        nlambda=nlambda,
        lambda_min_ratio=lambda_min_ratio,
        lambda_max=lambda_max,
        sym=sym,
        subsample_indices=subsample_indices,
    )

    opt = int(result["opt_index"]) - 1
    merge_opt = result["merge"][opt]
    refit = result["refit"]
    n, p = result["data"].shape
    return {
        "lambda_path": [float(v) for v in result["lambda"]],
        "opt_index": int(result["opt_index"]),
        "opt_lambda": float(result["opt_lambda"]),
        "stability": [float(v) for v in result["stability"]],
        "edges": _edges_from_adjacency(refit, merge_opt),
        "refit": {"format": "edge_list", "edges": _edges_from_adjacency(refit, None)},
        "diagnostics": {
            "n": int(n),
            "p": int(p),
            "method": "mb",
            "converged": not bool(result["hit_max_iter"]),
        },
    }


def main(argv: list[str] | None = None) -> int:
    # Required when the optional StARS process pool uses the "spawn" start
    # method (Windows, and frozen builds): without it a re-imported/packaged
    # main module would re-run the CLI in every worker. Harmless otherwise.
    import multiprocessing

    multiprocessing.freeze_support()

    parser = argparse.ArgumentParser(
        prog="spieceasi",
        description="SPIEC-EASI (mb) network inference: JSON in, JSON out.",
    )
    parser.add_argument(
        "--input", "-i", default=None, help="path to the JSON request (default: stdin)"
    )
    args = parser.parse_args(argv)

    try:
        if args.input is None:
            raw = sys.stdin.read()
        else:
            with open(args.input, encoding="utf-8") as fh:
                raw = fh.read()
        request = json.loads(raw)
        response = run_request(request)
    except Exception as exc:  # noqa: BLE001 - report any failure on stderr
        print(f"spieceasi: error: {exc}", file=sys.stderr)
        return 2

    json.dump(response, sys.stdout)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
