# huge 2.0.1

* Fixed a BLAS-dependent rounding issue in RIC selection. The rotated inner
  products are computed with matrix products, and a conforming BLAS may return
  any value inside a dot product's roundoff interval, so a mathematically zero
  regularization parameter could come back as a tiny positive value. RIC now
  certifies such values to exact zero using a pair-specific forward-error
  bound, which is scale-aware and keeps correlations the working precision can
  represent. Selected lambda values are therefore reproducible across BLAS
  implementations. Reported by CRAN's ATLAS additional check.
* `huge.select()` (R) and `huge_select()` (Python) no longer gate the
  zero-lambda refit on a bitwise-zero lambda. A residual too small for a
  solver to certify but not exactly zero previously bypassed the safety
  fallback, which made graphical lasso fit rank-deficient data almost
  unregularized and raise a non-finite estimate error. Routing now uses the
  same roundoff scale in both packages.

# huge 2.0.0

* Added explicit input routing for square symmetric observation matrices:
  R accepts `input.type = "data"` or `"covariance"`, and Python accepts
  `input_type="data"` or `"covariance"`; the compatible default `"auto"`
  retains symmetry-based detection. RIC and StARS refits preserve resolved
  raw-data routing. TIGER now sends the original matrix and resolved flag
  directly to C++, where covariance validation, correlation construction,
  automatic lambda selection, and fitting use one native correlation matrix.
  For compatibility with existing glasso callers, auto-detected covariance
  input retains its historical diagonal-sensitive default lambda scale;
  explicit covariance routing uses the corrected off-diagonal scale.
* Python MB and TIGER now consume the core's compressed column support
  directly instead of materializing and discarding an
  `nlambda * d * d` dense coefficient cube. The private native entries retain
  their default dense output for compatibility. RIC also limits its OpenMP
  team to the number of requested rotations, preventing idle workers from
  each allocating a `d * d` scratch matrix. R StARS uses at most two forked
  workers, limits huge's native OpenMP code to one thread in each child,
  disables recursive `mclapply`, and retains the historical positional
  `verbose` argument.
* Repaired the standalone CMake package: it now finds and exports BLAS,
  installs the real public headers and package configuration, respects the
  OpenMP option at link time, and has static/shared installed-consumer smoke
  coverage. The R package now explicitly requests C++17. Its macOS OpenMP
  probe locates R, compiles C++ rather than C, and loads and executes a
  representative parallel loop before enabling OpenMP.
* Python release checks now install the generated source distribution outside
  the checkout, smoke-test native MB and TIGER, exercise compressed binding
  output on every published CPython ABI from 3.9 through 3.14, and test the
  installed standalone core on Linux and macOS. Unsupported Python source
  platforms fail with an explicit message instead of receiving Linux flags.
* Fixed `huge.tiger()` covariance input being interpreted by the native solver
  as `d` raw observations. TIGER now builds one correlation matrix in C++,
  derives its default lambda path there, and solves the correlation-domain
  square-root Lasso with full KKT certification. Raw data and `cov(data)` now
  produce the same path for the same lambda values. Generated paths stop with
  a warning at the longest certified prefix instead of returning unstable
  low-lambda iterates; uncertifiable user-supplied values raise an error.
* TIGER's native automatic lambda path now reserves its `1e-3` fallback for
  identity and one-variable inputs. Weak nonzero correlations use their true
  maximum magnitude, preserving the established path scale in both R and
  Python while keeping lambda selection entirely in C++.
* TIGER's C++ automatic lambda grid now handles an unrepresentable
  `lambda_max * lambda_min_ratio` in log space and saturates only values below
  the smallest positive double. Subnormal correlations or ratios no longer
  produce zero/`NaN` lambda paths or spurious certification failures.
* Explicit TIGER lambda paths now have one cross-language contract: values
  must be positive, finite, and non-increasing, with ties allowed. Automatic
  tied paths can therefore be replayed through R and Python, while increasing
  paths are rejected by both frontends and both native TIGER entries.
* Explicit path validation now matches method semantics across the public R
  and Python interfaces. MB and glasso require positive, finite,
  non-increasing paths with ties allowed. CT retains its independent-threshold
  contract: finite non-negative values, including zero, run in supplied order;
  StARS alone requires that path to be non-increasing.
* Explicit lambda inputs now accept only a numeric scalar (including a NumPy
  0-D value) or a non-empty one-dimensional sequence. Multidimensional inputs
  fail closed instead of being flattened in R's column-major order or Python's
  row-major order; valid paths retain their method-specific rules.
* Direct native MB and glasso entries now validate non-empty, positive,
  finite, non-increasing lambda paths before allocating or solving. R native
  wrappers also reject mismatched matrix dimensions and path lengths, closing
  out-of-bounds reads; Python's native CT threshold entry now rejects
  non-square matrices and invalid thresholds while preserving supplied order.
* Fixed `huge()` tiger method returning asymmetric `icov` matrices: the
  in-place symmetrization averaged against already-averaged entries.
* Fixed the shared glasso core returning numerically asymmetric precision and
  graph matrices at finite solver tolerance. Each component is now projected
  onto the symmetric matrices before it is returned or reused as a warm start;
  path, `df`, sparsity, and log-likelihood are derived from that same projected
  precision matrix in both R and Python. Half-scaled addition avoids overflow
  when averaging large finite precision entries. Scale-equilibrated Cholesky
  now computes the log-determinant and rejects non-positive-definite precision
  estimates without an absolute pivot cutoff. The shared core also rejects
  non-finite results and covariance/precision pairs whose inverse residual is
  too large before they can be reused as a warm start or returned through R
  or Python. Components whose symmetric projection loses inverse consistency
  receive one tighter solve before final certification.
* Fixed Python raw one-variable inputs being collapsed to a scalar correlation
  by NumPy. CT, MB, glasso, and TIGER now return well-formed `1 x 1` paths for
  valid single-variable observations. Raw inputs with only one observation or
  a constant column are rejected before an undefined sample correlation can
  reach a solver.
* Aligned Python covariance auto-detection with R's machine-precision symmetry
  tolerance. Nearly symmetric square raw observations are no longer
  misrouted as covariance matrices. TIGER still constructs its correlation
  matrix and selects automatic lambda values in C++ after receiving the
  corrected raw/covariance flag.
* Fixed Python RIC selection changing at very small data scales. Standardizing
  for the rotation criterion now preserves every positive finite sample
  standard deviation instead of replacing values below `1e-12` with one.
* Fixed Python covariance-to-correlation conversion changing CT graphs under
  very small uniform rescalings. The conversion now uses the true positive
  diagonal, rejects non-positive diagonals and pairwise Cauchy--Schwarz
  violations, and clips only floating-point roundoff at the correlation bound.
  Covariance symmetrization also avoids overflowing on large finite entries.
* Fixed Python MB/glasso automatic lambda grids applying the `1e-3` zero-signal
  fallback to weak but nonzero correlations or covariances. TIGER automatic
  lambda selection remains entirely in C++.
* Fixed Python MB/glasso rejecting `lambda_min_ratio = 1`, despite the public
  `(0, 1]` contract. These automatic paths now retain the requested length as
  exact ties; explicit non-increasing MB/glasso paths may also contain ties so
  StARS can refit them, while increasing paths are still rejected.
* Fixed R `huge.mb()` and `huge.tiger()` dropping single-variable sparse paths
  to vectors. Their raw and covariance interfaces now return valid `1 x 1`
  paths with zero, finite sparsity; MB also constructs zero-edge sparse paths
  without invalid indices.
* Added shared R estimator input validation. CT, MB, glasso, and TIGER now
  reject non-matrix, empty, or non-finite input, and reject raw data with fewer
  than two observations or a constant column before correlation or solver
  code runs. Covariance input must also have a positive finite diagonal and
  satisfy pairwise Cauchy--Schwarz bounds, apart from `1e-8` correlation
  roundoff. Accepted near-symmetric covariance matrices are stably projected
  before CT, MB, or glasso uses them, preventing one-sided threshold graphs.
  TIGER still passes the original validated matrix to C++, where symmetric
  correlation construction and automatic lambda selection remain.
* Covariance validation for CT, MB, and TIGER now also requires the whole
  matrix to be positive semidefinite instead of checking only pairwise
  Cauchy--Schwarz bounds. Singular covariance matrices remain valid within a
  scale-aware roundoff tolerance. TIGER repeats this check in C++ before
  native correlation, automatic lambda selection, and fitting. Glasso still
  accepts indefinite pairwise covariance estimates, but returns them only
  when regularization produces a finite, positive-definite precision estimate
  and a certified covariance/precision pair.
* TIGER's C++ covariance symmetrization no longer adds two large finite
  entries before scaling. Valid covariance matrices near the floating-point
  maximum now reach native correlation construction and lambda selection
  without intermediate overflow.
* TIGER covariance normalization now applies the larger inverse standard
  deviation first. Representable weak correlations between variables with
  extremely different variances no longer underflow to zero or trigger the
  native identity-matrix lambda fallback.
* Covariance projection in the R and Python frontends now preserves entries
  that are already exactly symmetric. In particular, the smallest positive
  subnormal diagonal is no longer halved to zero before validation.
* R covariance auto-detection and direct native TIGER covariance input now use
  the same per-entry, scale-aware symmetry tolerance as Python. Materially
  asymmetric matrices can no longer be routed as covariance input and then
  averaged into a different, apparently valid problem.
* Covariance symmetry tolerance in R, Python, and native TIGER is now measured
  on the implied-correlation scale. Raw/covariance routing therefore remains
  stable across representable finite uniform rescalings, so tiny square raw
  samples are no longer misclassified as covariance matrices. TIGER still
  constructs correlation matrices and selects automatic lambda values in C++.
* TIGER now power-of-two rescales and reference-centers raw columns in C++,
  using compensated sums before constructing their correlation matrix. Finite
  data from the smallest subnormal through the floating-point maximum no
  longer look constant or distort native automatic lambda selection.
* R and Python raw-correlation construction for correlation thresholding, MB,
  and graphical lasso now reference-centers and normalizes columns before
  forming correlations. It no longer overflows or underflows at extreme
  finite scales; R also clips correlation roundoff to its mathematical range.
* TIGER with StARS now fails fast with a clear message. Different subsamples
  can certify different prefixes of the supplied lambda grid; selection is
  unsafe until a common certified-prefix protocol is implemented. TIGER with
  RIC remains supported. R StARS also rejects single-variable paths, for which
  its variability denominator has no edge pairs.
* Python `huge_generator().sigmahat` now matches R's empirical correlation
  matrix, and rejects `n = 1` where sample correlation is undefined.
* `huge.roc()` and Python `huge_roc()` now reject truth matrices containing
  only edges or only non-edges, for which ROC/AUC is undefined.
* Python source and wheel distributions now include the GPL-2.0-only license;
  missing native build dependencies cause a build failure rather than a
  silently unusable pure-Python wheel.
* Python release publishing is now gated on the full parity suite, artifact
  content checks, a clean wheel install, and the declared minimum dependency
  versions.
* Python path sparsity and EBIC edge counting no longer densify every graph
  along the regularization path.
* Python StARS now streams completed subsample paths into its frequency counts
  and accumulates sparse coordinates without dense per-graph temporaries.
  Counts use the smallest safe unsigned integer type for `rep_num` (normally
  one byte instead of eight), then convert one lambda layer at a time for the
  unchanged variability calculation.
  CT and MB additionally store only the condensed upper triangle; glasso
  retains both directions for compatibility with previously created or
  caller-constructed paths. StARS now rejects inputs with fewer than two
  variables, for which its variability denominator is zero. Python StARS also
  rejects a subsample ratio whose floored size is below two, rather than
  silently increasing the requested subsample.
* R StARS now merges each serial subsample path immediately instead of retaining
  every replication until all fits finish; parallel execution is unchanged.
* R StARS now rejects increasing or unordered lambda paths before subsampling;
  its first-threshold-crossing rule assumes a non-increasing regularization
  path. Tied CT thresholds remain supported, and RIC/EBIC are unchanged.
* Python numeric inputs now use strict finite-scalar validation; explicit
  lambda paths ignore inactive grid controls, CT accepts a zero threshold,
  and small-dimensional generator defaults cap edge probabilities at one.
* Python ROC, inference, and plotting now reject non-finite secondary matrices
  and report non-numeric conversion failures consistently.
* Python inference now rejects fewer than two observations, constant columns,
  non-positive precision diagonals, and numerically non-finite edge p-values.
  Gaussian inference supports the valid one-variable boundary with the same
  finite result as R; nonparanormal inference requires at least two variables.
* R inference now applies the same fail-fast data, dimension, precision
  diagonal, and finite edge-p-value contract. It preserves Gaussian
  one-variable inference and legal nonparanormal score-test diagonal `NaN`s,
  and handles finite rank differences near floating-point limits safely.
* R `huge.generator()` now rejects invalid `n`, `d`, and graph names before
  printing or consuming random numbers. Empirical correlation requires at
  least two generated observations.
* R `huge.generator()` now supports the one-variable boundary for every graph
  type with zero sparsity, and band graphs support bandwidth `d - 1` without
  dimension-dropping errors.
* R `huge.generator()` now validates active `v`, `u`, `g`, and `prob` values
  before printing or consuming random numbers. Graph-inactive `g` and `prob`
  arguments are ignored, and oversized group counts no longer allocate empty
  groups.
* R `huge.generator()` now requires scalar, non-missing logical values for
  `vis` and `verbose`, rejecting invalid flags before output, plotting, or RNG
  consumption.
* R `huge.generator(vis = TRUE)` now restores the caller's graphics parameters
  after visualization, including when plotting exits with an error.
* R `plot.huge()` and `huge.plot()` now restore dependent graphics parameters,
  including custom figure, plot, and outer regions, in a stable order. EPS
  output closes its own device on exit and reliably returns to the original
  caller when multiple graphics devices are open.
* Remaining R plotting entry points (`plot.sim()`, `plot.select()`,
  `huge.roc()`, and `plot.roc()`) now use the same dependency-aware
  graphics-state restoration.
* Python Gaussian inference now ignores its nonparanormal-only method option,
  and non-EPS plotting ignores inactive output-file arguments.
* Fixed `huge.select()` RIC criterion scale-dependence: the rotation
  criterion saw unstandardized data while the lambda path is defined on the
  correlation scale, so rescaling the input silently changed the selected
  graph. The Python interface now also preserves R's empty-graph boundary,
  including single-variable and exactly uncorrelated inputs, instead of
  passing a zero lambda back to solvers that require a positive value. Also
  fixed the R shortcut, which compared against `max(cor(data))` (always 1),
  and made single-variable RIC report zero rather than `NaN` sparsity.
* Fixed `huge.npn(npn.func = "skeptic")` erroring on input with row names
  (n-length row names were assigned to the d x d output); the skeptic
  correlation matrix now carries variable names on both dimensions.
* Fixed Python nonparanormal skeptic estimation for one- and two-variable
  input. It now uses a standardized rank crossproduct, preserves the valid
  nonconstant `1 x 1` boundary, and rejects fewer than two observations or
  undefined correlations caused by a constant column.
* Solvers now warn when an iteration budget is exhausted instead of
  silently returning partially converged estimates (new `hit_max_iter`
  channel from the shared C++ core).
* Added `huge.select(..., num.cores = k)`: StARS subsamplings can now fit in
  parallel via `parallel::mclapply`. The R interface starts at most two forked
  workers, limits huge's package-owned OpenMP regions to one thread in each
  child, and returns results identical to the serial path. External BLAS
  thread settings remain controlled by the R installation. New Imports:
  parallel.
* Speedups: RIC selection now uses BLAS matrix products (about 10-20x);
  mb adds a sequential strong screening rule with full KKT certification
  (about 1.2x, and estimates now satisfy exact KKT);
  glasso updates covariance columns over the active set only (about 1.15x)
  and avoids strided memory access in its log-likelihood computation.
* Speedups in the R layer: correlation matrices are computed via BLAS
  `crossprod` (about 40x for the correlation step; `huge(method = "ct")`
  about 2x end to end), the nonparanormal skeptic uses rank crossproducts
  (about 5x), `huge.roc` uses precomputed edge masks (about 3x), and
  `huge.generator` uses symmetric eigenvalue computation and Cholesky-based
  inversion (about 1.6x at d = 2000).
* Reduced memory: mb/tiger graph symmetrization and StARS variability now
  stay in sparse matrices instead of densifying d x d intermediates. R and
  Python glasso adapters derive graph support from precision matrices instead
  of retaining a duplicate dense native path; Python CT converts native
  thresholds one lambda at a time.
* Fixed graphics parameter leaks in `plot.huge`, `plot.select`, `plot.roc`,
  `huge.roc`, and `huge.plot` (`par` is now restored on exit).
* Removed the unused `maxdf` argument from internal C++ wrappers.
* Documentation fixes: corrected dimension notations, the ct `nlambda`
  default (20, not 30), the `scr` description in `huge.glasso`, and several
  typos.
* Added an automated check that the C++ core copies shared with the pyhuge
  Python package remain identical (`tools/check_core_mirrors.sh`, CI).

# huge 1.6

* Fixed `huge.roc()` F1 score: precision was computed from rates instead of
  counts, inflating F1 in sparse graphs (TP rate and FP rate were correct;
  only F1 was affected).
* Fixed `huge.select()` missing support for tiger method (no default criterion,
  RIC/StARS refit, or opt.icov extraction).
* Fixed C++ performance regression for glasso (3-4x) and tiger (7x) vs 1.4:
  replaced scalar loops with BLAS ddot calls, precomputed column norms, and
  reverted glasso screening to proven basic method.
* Fixed `plot.sim()` graphics parameter leak (`par()` shadowing and missing
  `on.exit()` restore).
* Fixed `ebic.score` field naming mismatch in documentation.
* Fixed `scr` parameter documentation (incorrectly stated MB not supported).
* Fixed `align` parameter documentation in `plot.huge` (logically inverted).
* Removed dead `fit$rss` matrix from `huge.mb()` and `huge.tiger()`.
* Added CI workflow for R CMD check on Linux and Windows.

# huge 1.5

* BLAS acceleration for graphical lasso and TIGER solvers (2-3x speedup at d >= 200).
* Added strong screening rule for graphical lasso.
* Removed RcppEigen dependency; the package now uses lightweight C++ headers with direct BLAS calls.
* Added testthat test suite.
* Added TIGER to DESCRIPTION text.
* Fixed Makevars.win BLAS linking.
* Fixed nonparanormal transform (`huge.npn`) rownames bug.
* Fixed documentation typos and parameter reference errors.
