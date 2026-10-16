"""StARS model selection.

Faithful port of the R ``pulsar`` functions used by SpiecEasi:

- ``stars.stability`` (pulsar/R/mcPulsarSelect.R:339)
- ``.starsind``       (pulsar/R/mcPulsarSelect.R:329)

``spiec.easi`` sets ``thresh = 0.05`` (pulsar's own default of 0.1 is
overridden).
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np


def starsind(summary: np.ndarray, thresh: float, offset: int = 1) -> int:
    """``.starsind`` (pulsar/R/mcPulsarSelect.R:329).

    First index (1-based) at which the (monotonised) summary reaches ``thresh``,
    minus ``offset``, floored at 1. If the threshold is never reached, warns and
    returns ``length(summary)``.
    """
    ge = np.asarray(summary) >= thresh
    if bool(np.any(ge)):
        first = int(np.argmax(ge)) + 1  # which.max(...)=first TRUE, 1-based
        return max(first - offset, 1)
    warnings.warn("Optimal lambda may be outside the supplied values", stacklevel=2)
    return int(len(summary))


def stars_stability(
    premerge: list[list[np.ndarray]],
    stars_thresh: float | None,
    rep_num: int,
    p: int,
) -> dict[str, Any]:
    """``stars.stability`` (pulsar/R/mcPulsarSelect.R:339).

    ``premerge`` is indexed ``[lambda][replica]`` and holds the 0/1 adjacency
    matrices of each replica. Returns a dict with ``merge`` (edge frequencies),
    ``summary`` (monotonised variability), ``opt_index``, ``criterion`` and
    ``thresh``.
    """
    if stars_thresh is None:
        stars_thresh = 0.1

    merge: list[np.ndarray] = []
    for lam_reps in premerge:
        acc = np.zeros_like(lam_reps[0], dtype=float)
        for a in lam_reps:
            acc = acc + a
        merge.append(acc)

    summary = np.zeros(len(merge))
    for i, m in enumerate(merge):
        m = m / rep_num
        merge[i] = m
        summary[i] = 4.0 * np.sum(m * (1.0 - m)) / (p * (p - 1))

    summary = np.maximum.accumulate(summary)  # monotonise
    opt_index = starsind(summary, stars_thresh)
    return {
        "merge": merge,
        "summary": summary,
        "opt_index": opt_index,
        "criterion": "stars.stability",
        "thresh": stars_thresh,
    }
