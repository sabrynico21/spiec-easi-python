"""spieceasi — faithful Python port of the R SpiecEasi `mb` pipeline.

The implementation mirrors the R source in ``reference/R/`` (SpiecEasi) and the
numeric routines it delegates to (``huge`` for the Meinshausen–Bühlmann
neighbourhood-selection Lasso, ``pulsar`` for the StARS model selection).

See ``docs/PORTING_NOTES.md`` (porting map), ``docs/EQUIVALENCE.md`` (what is
guaranteed identical) and ``docs/DIFFERENCES.md`` (deviations).
"""

from __future__ import annotations

from .lambda_path import get_lam_path, get_max_cov
from .pipeline import (
    get_opt_beta,
    get_opt_ind,
    get_opt_lambda,
    get_opt_merge,
    get_opt_net,
    get_stability,
    spiec_easi,
)
from .transforms import clr, spiec_easi_norm

__all__ = [
    "__version__",
    "spiec_easi",
    "clr",
    "spiec_easi_norm",
    "get_max_cov",
    "get_lam_path",
    "get_opt_ind",
    "get_opt_lambda",
    "get_opt_merge",
    "get_opt_net",
    "get_opt_beta",
    "get_stability",
]

__version__ = "0.1.0"
