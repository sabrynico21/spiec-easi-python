# Differences from the R source and from the paper

This file records every place where the implementation deliberately follows the
R **code** over the article, or where the Python port differs from R for a
documented reason. It is a living document.

## Code vs. article

- *TBD* — to be filled while implementing. The R code is the source of truth; any
  divergence between the code and Kurtz *et al.* (2015) is noted here.

## R behaviour intentionally reproduced (not "fixed")

- The SpiecEasi `clr` does **not** use the textbook centred-log-ratio
  transformation. `clr.default` (`normalization.R:95`) centres by
  `mean(LOG) / mean(x >= tol)`, i.e. divides the mean of the log-values by the
  fraction of non-zero entries. This is reproduced literally.
- `getMaxCov` uses the strict upper triangle of the correlation matrix
  (`k = 1`), i.e. the diagonal is excluded.

## Python-port deviations

- **Subsample RNG** — see `EQUIVALENCE.md`. Indices are injectable.
- *TBD* — any deviation introduced while porting (e.g. dtype handling, memory
  layout) is recorded here with its impact on results.
