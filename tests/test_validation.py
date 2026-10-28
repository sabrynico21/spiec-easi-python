"""Input-validation and CLI tests (fast)."""

from __future__ import annotations

import json

import numpy as np
import pytest

from spieceasi.cli import main, run_request
from spieceasi.pipeline import spiec_easi


def _counts(n: int = 30, p: int = 4, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.poisson(5.0, size=(n, p)).astype(float) + 1.0


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_rejects_non_2d() -> None:
    with pytest.raises(ValueError, match="2-D"):
        spiec_easi(np.arange(10.0))


def test_rejects_nan() -> None:
    c = _counts()
    c[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        spiec_easi(c, pulsar_params={"rep_num": 2})


def test_rejects_too_few_samples() -> None:
    with pytest.raises(ValueError, match="2 samples"):
        spiec_easi(np.ones((1, 3)), pulsar_params={"rep_num": 2})


def test_raw_constant_count_column_is_not_rejected_before_clr() -> None:
    c = _counts()
    c[:, 2] = 7.0
    result = spiec_easi(
        c,
        nlambda=3,
        pulsar_params={"rep_num": 2, "seed": 4},
    )
    assert np.isfinite(result["lambda"]).all()
    assert np.isfinite(result["stability"]).all()


def test_warns_on_negative() -> None:
    c = _counts()
    c[0, 0] = -1.0
    with pytest.warns(UserWarning, match="Negative"):
        spiec_easi(c, pulsar_params={"rep_num": 2})


def test_rejects_unknown_method() -> None:
    with pytest.raises(NotImplementedError):
        spiec_easi(_counts(), method="glasso")


def test_rejects_unknown_pulsar_param() -> None:
    with pytest.raises(ValueError, match="Unrecognized"):
        spiec_easi(_counts(), pulsar_params={"bogus": 1})


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_run_request_small() -> None:
    req = {
        "counts": _counts().tolist(),
        "params": {"nlambda": 5, "pulsar_params": {"rep_num": 3}},
    }
    resp = run_request(req)
    assert set(resp) >= {
        "lambda_path",
        "opt_index",
        "opt_lambda",
        "stability",
        "edges",
        "refit",
        "diagnostics",
    }
    assert len(resp["lambda_path"]) == 5
    assert resp["diagnostics"]["method"] == "mb"
    assert resp["diagnostics"]["n"] == 30 and resp["diagnostics"]["p"] == 4
    # edges are 0-based, i < j, weight in [0, 1]
    for i, j, w in resp["edges"]:
        assert 0 <= i < j < 4
        assert 0.0 - 1e-12 <= w <= 1.0 + 1e-12
    # refit is the same support as edges
    assert {tuple(e[:2]) for e in resp["refit"]["edges"]} == {
        (i, j) for i, j, _ in resp["edges"]
    }


def test_main_reads_stdin_writes_stdout(monkeypatch, capsys) -> None:
    req = {
        "counts": _counts().tolist(),
        "params": {"nlambda": 4, "pulsar_params": {"rep_num": 2}},
    }
    monkeypatch.setattr("sys.stdin", __import__("io").StringIO(json.dumps(req)))
    rc = main([])
    assert rc == 0
    out = capsys.readouterr().out
    parsed = json.loads(out)
    assert parsed["opt_index"] >= 1


def test_main_error_exit_code(monkeypatch, capsys) -> None:
    monkeypatch.setattr("sys.stdin", __import__("io").StringIO("{not json"))
    rc = main([])
    assert rc == 2
    err = capsys.readouterr().err
    assert "error" in err.lower()
