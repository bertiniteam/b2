
# bertini2 Changelog

All notable changes to this project will be documented in this file.

The format is based on [CHANGELOG.md][CHANGELOG.md]
and this project adheres to [Semantic Versioning][Semantic Versioning].

<!--
_______________________________________________________________________________

## [1.0.0] - 2026-04-02

Preparation for pypi release with github workflow

### Added

- github workflow for pypi and github release

### Changed

- `publish-to-test-pypi.yml` for handling the comments

### Changed

* merged the pull request for github ci release by @hkmoon in https://github.com/hkmoon/b2/pull/1
* windows release preparation
  * `size_t` is translated into `unsigned long` in linux, mac while `unsigned long long` in windows 10: `core/include/bertini2/eigen_extensions.hpp` and `core/test/classes/start_system_test.cpp` are modified
  * use `clang` of LLVM in Windows since MSVC has different compiling way for `template`
  * use `--no-isolation` for `scikit-build` in Windows
* For linux wheel naming convention, we cannot use x86_64, x86_i386 anymore for pypi repository. https://peps.python.org/pep-0600/
  * use `auditwheel` for it

### New Contributors
* @hkmoon made their first contribution in https://github.com/hkmoon/b2/pull/1

_______________________________________________________________________________

-->

<!--
_______________________________________________________________________________
TEMPLATE

## [major.minor.patch] - yyyy-mm-dd

A message that notes the main changes in the update.

### Added

### Changed

### Deprecated

### Fixed

### Removed

### Security

### New Contributors

_______________________________________________________________________________

-->


_______________________________________________________________________________

## [3.5.0] - unreleased

A correctness fix to the `MakeMovingHomotopy` guards: they decided function identity on a
*presentation* rendering, which silently refused valid homotopies.

### Removed

- **A System no longer carries a precision.**  `System::precision(unsigned)`,
  `System::precision()` and the `precision_` member are all gone, in C++ and in Python.
  Evaluation happens at the precision of the point it is handed, so there is nothing for a
  caller to set and nothing to keep in sync -- that hand-alignment
  (`target = max(point, system, ambient); system.precision(target); ...`) was the complaint
  in #377, and refusing to evaluate on a mismatch could wedge a System outright, with no
  escape through either setter.
  The "materialized at" tag now lives with each holder of multiprecision values -- every
  block's working coefficients, the patch, the SLP's memory -- each of which returns
  immediately when already there, so `SetVariables` fans out on every evaluation for the
  cost of a few integer compares.  The elision is deliberately per-holder and not at the
  System: holders can legitimately disagree (that disagreement IS #377), and a System-level
  cache would skip the very repair such a case needs.
  There is also nothing to prepare and nothing to fan out.  Every evaluable type -- each of
  the four blocks, the patch, and the SLP -- self-aligns to the precision of the point it is
  handed, under one uniform `SyncPrecision`.  Two of them already did this independently;
  the change makes the pattern and the name uniform, and the patch had a commented-out
  assert demanding callers match its precision, which it now honours by aligning itself.
  See ADR-0057.
  **Archive format changed**: `precision_` is no longer serialized (transient evaluation
  state, which the same `serialize` already excludes elsewhere), so a System archived by an
  older build will not load into a newer one.  Boost archives carry Systems between MPI ranks
  of one run; durable storage is the records/JSON path and is unaffected.

### Added

- **A showpiece, "The Seam"**, on what path tracking does at a branch cut.  Every pixel of a
  parameter grid is one tracked path of `log(x) - y = 0`, `y^2 - x + c = 0`, and every step of
  every path is exposed onto one image: arrivals sweep away, and the paths that fail crawl onto
  the negative real axis and stall there, drawing a line that is in no equation.  The line is
  where the numeric logarithm jumps, which is a convention rather than a theorem, and the only
  reason it is visible is that the library hands back the complete trajectory of a path that
  failed -- position by position, with how far along the path time it got.  Measured and stated
  on the page: it is not a precision artifact, since double precision and eighty digits agree on
  every one of 799 sampled parameters, arrivals and deaths alike.
- **A tutorial on tracking an analytic homotopy** ("Tracking an analytic homotopy", under
  "Doing things manually").  Sine has no degree, so no start system exists and no count bounds
  its roots -- but a homotopy needs neither, only start points you choose.  The page tracks
  `sin(x) = 1/2` from five zeros of sine and checks the endpoints against `arcsin`, shows why
  the precision model has to be chosen rather than defaulted, and walks into a branch point at
  `sin(x) = 1` where the power series endgame reports cycle number 2 -- correct, and explained
  by Weierstrass preparation, which says an analytic family branches like an algebraic one near
  a finite endpoint.  Doctested, with a figure of the paths in the complex plane drawn from real
  tracked data, and plain about what this does not give you: no solution count, no completeness,
  and no theory for paths that run to infinity.  The neighbouring tutorial has promised this
  since it was written, and then changed the subject.
- The endgames announce every sample point they compute along a path (`ComputedSamplePoint`,
  carrying the point and its time), from both the power series and the Cauchy endgame, and a
  `SampleSequenceCollector` observer (C++ and Python, `bertini.endgame.*.SampleSequenceCollector`)
  gathers them: the sequence of points on the path at shrinking times, kept apart from the
  Cauchy loop's circle points and from the successive approximations of the root, with run
  boundaries so one collector attached to a solver's endgame can tell one path from the next.
  This is what lets a caller watch a quantity such as a Jacobian's singular values as a
  function of distance to the root instead of judging it at one point.  Every endgame run now
  announces `Initializing`; an adaptive endgame that needs a higher precision before its first
  approximation recomputes its sample window at the new precision and announces
  `SamplesRecomputedAtHigherPrecision`, so the collector drops the superseded samples as the
  endgame does, and counts every precision increase (`num_precision_increases`).  (#361, in part)
- The records archive is reloadable.  Every system and homotopy a solve records is stored as its
  exact canonical encoding (the text its content digest is the hash of); there is now a reader for
  that text: `System.from_canonical(text)` rebuilds a system from it, `System.canonical_encoding()`
  produces it, and `bertini.records.load_system(digest, directory)` fetches and rebuilds an
  archived one.  The rebuilt system's content digest equals the archived one, and that equality is
  checked on load.  Classic (Bertini 1) text plays no part in persistence.  (C++:
  `System::FromCanonicalEncoding`, `records::LoadSystem`, `node::DecodeCanonical`.)
- `SolutionMetaData.singular_values`: the singular values of the target system's Jacobian at the
  endpoint, largest first, at the endpoint's own precision -- the spectrum `condition_number` was
  already computed from and threw away.  Published raw, archived with each path record and
  restored on recall, so a caller can read a numerical rank at a tolerance of its own choosing;
  no rank verdict is baked in, because no single tolerance suits every endpoint.  (#409)
- `SolutionMetaData.crossing_unresolved`: true for a path that was flagged as crossing another at
  the endgame boundary and whose crossing re-tracking did not resolve.  Its success codes may
  still read `Success`; this flag is the per-path form of the solver's "the affected solutions
  may be wrong" warning, and it is recomputed on recall from the archived boundary data.  Also a
  column of `to_dataframe()`.  (#365)
- Python bindings for two endgame accessors that already existed in C++ but were unreachable:
  `previous_approximation()` and `approximate_error()`, alongside the already-bound
  `final_approximation()`.  Together they give the pair of successive root approximations the
  endgame's own convergence test compares, plus the infinity norm between them -- a second
  sample of the root at a known, coarser accuracy, which is what lets a caller judge how a
  derived quantity (the singular values of a Jacobian, say) behaves as the approximation
  improves, rather than thresholding it at one point.

### Changed

- **The adaptive-precision error bounds are named after what they bound.**  On
  `AMPConfig`, `Phi` is now `jacobian_eval_error_bound`, `Psi` is now
  `function_eval_error_bound`, and `epsilon` is now `linear_solve_error_bound` -- in C++ and
  in Python (`phi`, `psi` and `epsilon` as attribute names are gone).  `SetPhiPsiFromBounds`
  is `SetErrorBoundsFromDegreeAndCoefficient` (`set_error_bounds_from_degree_and_coefficient`
  in Python) and `SetBoundsAndEpsilonFrom` is `SetBoundsFrom` (`set_bounds_from`).  A Greek
  letter told a caller nothing, which matters more now that a refusal invites them to set
  these values by hand.  **No digest moved**: the canonical encoding keeps the old spellings,
  because that text is a preimage rather than presentation, so `b2cfgenc` is untouched and
  every existing record still recalls.
- **Adaptive precision refuses a system that is not polynomial, instead of miscalibrating.**
  Its error bounds are derived from the system's degree, which such a system does not have,
  and the derivation used to produce a plausible-looking Jacobian bound and a **negative**
  function bound; every path then died before the endgame with
  `FailedToSelectPrecisionAndStepsize`, after precision had climbed to 290 bits.  Tracking
  an analytic homotopy at fixed precision was, and remains, fine -- the refusal names that
  way out, and the other one: supply the two evaluation error bounds yourself, which
  `bertini.HomotopySolver` now accepts as `amp_config=`.  The open question of what those
  bounds should be for an analytic system is issue #439.
- **Printing is a family of dialects, one per audience** (ADR-0059).  `str(node)` and
  `str(system)` are the Python spelling: powers as `**`.  `repr(node)` is the exact Python
  spelling: `eval` of it in the `bertini` namespace rebuilds an equal node at full precision
  (`real_mp('digits', precision)`, `Complex('re', 'im', precision)`, `Rational('p/q')`; integers
  and values a Python literal holds exactly stay bare).  `node.to_classic()`,
  `System.to_classic_input()` and the CLI write Bertini 1's spelling: `^`, and a complex
  constant as `(re+im*I)` -- the only complex form Bertini 1 reads.  The `(re,im)` pair form is
  gone from every output, and so are the two Python shims that patched around it (the `^` to
  `**` character replace in `repr`, the `(re,im)` regex in `parse.system`).  `real_mp(text,
  precision)` and `Complex(re, im, precision)` construct a value at a chosen precision.
- **The system encoding version is `b2sysenc/2`.**  Canonicalization breaks multidegree ties
  between operands on the canonical encoding instead of on printed text, and the encoder writes
  operands in that order.  Records written by earlier versions carry `b2sysenc/1` and are not
  recalled; the "records read forever" promise is withdrawn until the identity of the algorithm
  itself is part of a record's ask (#420).
- **The default endgame is the power series endgame**, matching Bertini 1's documented default
  (`EndgameNum: 1`), at every choice point: the configuration default (so the CLI and a classic
  input without `endgamenum`), and the Python factories `ZeroDimSolver`, `HomotopySolver`,
  `user_homotopy`, `parameter_sweep` and `bertini.solve`.  Measured on a regeneration workload,
  the Cauchy default cost 190 s where power series takes 2.5 s, because slowly diverging paths
  are a Cauchy-specific pathology.  Cauchy remains available by explicit request
  (`endgamenum: 2`, `endgame='cauchy'`).  A solve that relied on the default is now a different
  computation, so its records are new asks rather than recalls.  See ADR-0058.

### Fixed

- **A degree that does not exist is refused rather than laundered into a number.**  A function
  that is not a polynomial reports degree -1, which is the right answer and a catastrophic
  operand, and four places treated it as one.  `System::DegreeBound()` took a maximum, so the
  sentinel vanished behind any function of positive degree -- a system of one polynomial and
  one sine reported a degree bound of 2 and said nothing; it now refuses.  `System::Randomize`
  compensates differing degrees with powers of a homogenizing variable, and its target degree,
  starting at zero, clamped the sentinel to a definite claim of degree zero, after which a
  transcendental function was quietly multiplied by one such power; it now refuses, and
  whether squaring up means anything for an analytic system is issue #440.  `BlendBlock` did
  the same clamping, so a homotopy onto a sine reported itself as a constant; it now
  propagates the sentinel.
- **`System::IsPolynomial()` looked only at the declared variable groups**, never at the
  ungrouped variables, and a variable outside every group has degree 0 with respect to every
  group.  So `sin(t)*x^2` with `t` ungrouped reported degree 2 and claimed to be polynomial,
  and a system with no variable groups at all was vacuously polynomial whatever it contained.
  Every guard built on the check inherited the hole, including the zero-dimensional solver's
  refusal of non-polynomial systems.  Removing the ungrouped-variable concept entirely is
  issue #442, targeted at 3.6.
- A power with a non-integer exponent and no variable in it -- `5^(1/2)`, which is how Bertini 1
  input spells a square root -- answered "not homogeneous" while reporting degree 0.  A system
  with such a constant among its coefficients (the Barth sextic with the golden ratio written
  as `(5^(1/2)+1)/2`) therefore homogenized to nothing and `AutoPatch()` refused the result
  with "requesting to AutoPatch a system which is not homogenized".  A variable-free power is a
  constant and is now homogeneous whatever its exponent, in agreement with its degree.  A
  table-driven test pins the variable-free form of every operator (`sqrt`, `exp`, `log`, the
  trigonometric functions, rational powers, and compositions of them) as degree 0, polynomial
  and homogeneous, both on its own and as a coefficient of a system that must homogenize and
  patch.  (#419)
- Negating a sum inflated its degree with respect to a variable group: `NegateOperator`
  inherited `UnaryOperator`'s group degree, which summed the per-variable degrees -- correct
  only for a single monomial -- so `-(x^2+y^2)` reported degree 4, and `System::DegreeBound()`,
  which sizes adaptive precision, followed it.  This is the defect #397 found in
  `PowerOperator`, in the other class that inherited the sum.  Negation now passes its
  operand's degree and multidegree through, and the non-polynomial rule (a variable-free
  operand makes a constant; anything else is not a polynomial) lives once in `UnaryOperator`
  instead of in four identical copies in `sqrt`, `exp`, `log` and the trigonometric operators.
- The tracker's path-truncation check measured the 2-norm of the point while the endgame's
  `Security::max_norm` check and the post-processing `endpoint_finite_threshold` measure the
  infinity norm (the largest coordinate), so the three thresholds did not measure the same
  quantity: a path whose coordinates all stayed under `path_truncation_threshold` was truncated
  once `sqrt(n)` carried its 2-norm over the line.  The tracker now uses the infinity norm too,
  as Bertini 1 does.  The never-incremented `num_total_steps_taken_` member is gone
  (`NumTotalStepsTaken()` already computed the sum).  (#404)
- `RandomConjugateOrthonormalMatrix(rows, cols)` factored a square matrix of the *larger*
  dimension and truncated it, so an `8 x 4908` randomization matrix for an isosingular deflation
  became a `4908 x 4908` multiprecision factorization that did not finish in 900 s.  It now
  factors a matrix sized to the request (the longer side by the shorter, transposed when the
  shape is wide): the same distribution, at O(max * min^2) instead of O(max^3) -- milliseconds.
  Square requests draw and factor exactly as before; non-square ones consume fewer random draws
  and so differ from the old recipe for the same seed.  (#401)
- `to_classic_input()` left out the `pathvariable t;` declaration, so the classic text of a
  homotopy did not parse back: the functions used a variable the file never declared.  A
  homotopy is a system too; the declaration is now emitted with the variable groups, and a
  homotopy round-trips through the classic writer and parser.  (#366)
- An affine products-of-linears block (`add_products_of_linears`, the shape of a regeneration
  deformation) did not homogenize: its `Homogenize` was a no-op written for the m-homogeneous
  start system's already-homogeneous form, and its `IsHomogeneous` always answered true, so a
  homogenized system kept an affine block with the wrong variable count and could neither be
  expanded to nodes ("variable count mismatch") nor tracked projectively.  The block now folds
  each factor's constant column onto the homogenizing variable, exactly as the linear-forms
  block does, and reports homogeneity from its shape.  (#376)
- Four Python-boundary defects where a caller error killed the interpreter or came back as a raw
  Boost.Python error, now Python exceptions that name the problem: a start point with the wrong
  number of coordinates, or a homotopy with more functions than variables, handed to
  `HomotopySolver` aborted the process from inside the tracker (`ValueError` at construction;
  #369, #383, #386); a non-square matrix, or a right-hand side of the wrong length, handed to
  `bertini.linalg.solve` corrupted the heap inside Eigen's LU (`numpy.linalg.LinAlgError`, as the
  double path already raised, and the native layer refuses too; #390); `System.eval` did not
  accept a plain list as the point (#367); a `Slice` was not accepted by `System.add` or by
  `moving_homotopy` where a system of linear forms was meant (#372, #381).
- Two solver-configuration traps (#392): `get_config()` refused the very names `config_names()`
  lists (it takes them now, alongside the class); and `final_tolerance` set directly on the
  endgame was silently reverted at every `solve()` by the solver's own copy.  The solver's value
  now flows into the endgame at setup and whenever the solver's value changes, so whichever was
  set last wins and nothing is reverted without a word.
- Building a large expression no longer costs its expansion (the second half of #417, closing
  it): canonicalization broke multidegree ties by printing the operands, and the printer walks
  the tree -- a 46-node graph printed as 196 KB at every level.  Ties are now broken on the
  canonical encoding, which is linear in the graph, computed only when two operands tie.
- Classic input with Bertini 1 comments (`%` to the end of the line) is accepted by every
  parse entry point (`System(text)`, `parse.system`, the CLI), including a comment that
  mentions `INPUT` or `END;`: comments are removed in C++ before the file wrappers are
  unwrapped and the grammar runs.  The Python text-scanning shim that looked for the INPUT
  section is gone.  (#407)
- Asking for a random slice with more linear forms than variables (`Slice.random_complex(vars,
  dim)` with `dim > len(vars)`, and the real and through-point forms alike) silently produced
  dependent rows; it is now a `ValueError` (C++ `std::invalid_argument`) that says so.  (#380)
- `MakeMovingHomotopy` no longer rejects a valid deformation whose two moving endpoints merely
  *print* the same.  Both of its guards (a fixed function duplicated in the moving rows, and a
  moving row identical at both endpoints) compared functions via `operator<<`, whose default
  stream precision is **6 significant digits** -- so two genuinely different rows agreeing to 6
  digits compared equal and the homotopy was refused with a message asserting the endpoints were
  the same row.  The collision is ~1e-6 *relative*, so it bit at every coordinate scale.  Identity
  is now decided on `node::CanonicalEncoding` -- the exact-value encoding the content digests are
  built on (ADR-0042) -- while `operator<<` is still used for the human-readable message text.
  Found by a surface cell decomposition whose slice values were 1.9e-5 apart at a magnitude of
  2409, and reproduced at order one (0.162749 vs 0.1627491).  (bertiniteam/b2#391)
- The power series endgame left `previous_approximation_` holding a COPY of
  `final_approximation_` after every successful run.  It assigned the two at the bottom of its
  convergence loop while testing the loop condition at the top, so the assignment ran one final
  time on the way out.  The Cauchy endgame never had this -- it returns from its acceptance gate
  before the corresponding assignment -- so the two endgames disagreed about their own post-run
  state.  Power series now matches Cauchy.  Consequences: `PreviousApproximation()` is now a
  genuine predecessor for both endgames, and `ZeroDimSolver`'s reported
  `accuracy_estimate_user_coords` -- computed as the distance between the final approximation and
  the previous one -- is no longer identically zero for power-series solves, which had it
  reporting an exactly-perfect accuracy for every such path.
- `EndgameBase::approximate_error_` was left uninitialized, so `ApproximateError()` read an
  indeterminate value before any run.  Now initialized to infinity, which is the only safe
  sentinel: the convergence gates compare it in both directions, and NaN -- which loses every
  relational comparison -- would make the power series loop's `error > tolerance` test false and
  skip the loop entirely, reporting instant success.

_______________________________________________________________________________

## [3.4.0] - 2026-07-16

A one-call way to build a linear slice that passes through a chosen point, an endgame-hardening
sweep (power series and Cauchy both), and friendlier start-point handling in the Python layer.

### Added

- **`Slice.through_point(variables, point, dim=1, coefficients=None, real=False, orthogonal=True, homogeneous=False)`**
  — the single place to make a slice through a given point (#343).  With `coefficients=None` (the default) it
  draws a random block of `dim` linear forms (complex, or real with `real=True`; orthonormalized when
  `orthogonal`) and sets the constant column so every form vanishes at `point`; pass `coefficients` (a
  bare, non-augmented block) to use exactly those directional coefficients.  Backed by the new C++
  primitives `Slice::ThroughPoint` and a `through_point` option on `Slice::RandomReal` /
  `Slice::RandomComplex`.  `homogeneous=True` builds a projective slice through the point (rows
  orthogonal to it); it is random-only.
- **`HomotopySolver` accepts start points in any faithful numeric representation** (#350, fixing #347):
  multiprecision scalars, Python/numpy numbers, and *constant* symbolic nodes (an ndarray of
  `symbolics.Complex` previously crashed with a raw eigenpy converter error).  Start points are
  transported values — the tracker refines them — so lossy doubles are welcome here; expressions still
  containing variables are refused as the math errors they are, with an error that names the variables.
  `complex_mp` now constructs explicitly from a Python `complex` (still no implicit conversion).

### Fixed

- **The power-series endgame's Hermite interpolation now evaluates the actual Hermite
  interpolant** (#353).  The Horner walk over the doubled node list advanced at half speed, evaluating
  a different (lower-order) interpolant — the limit was still correct, but convergence order was
  degraded.  This changes computed results at agreeing inputs: approximations land measurably closer
  to the truth (the old test oracle values were themselves off and have been re-pinned exactly).
- **`max_cycle_number` is enforced as a ceiling, not a floor** (#353).  The bound was applied with
  `max()`, and a near-unity sample ratio could push the estimate through an unsigned conversion of
  infinity (UB).  Clamped before conversion; cycle-number candidates now default sanely on
  degenerate samples.
- **NaN is a failure, never `Converged`** (#354).  IEEE comparison semantics made every NaN
  comparison false, so a NaN correction step exited the convergence loop as success, and the
  security valve (`norm > max_norm`) was blind to NaN norms.  Both endgames now fail fast on NaN
  approximations (new `bertini::ContainsNaN` in `eigen_extensions.hpp` — component-wise, because
  multiprecision complex NaN compares *equal* to itself) and the valve is NaN-aware.
- **Slow divergers truncate honestly in the Cauchy endgame** (#355).  Paths diverging to infinity
  slower than the shrinking time zones could grind precision escalation for minutes before dying.
  The security valve now also arms below B1's `cycle_cutoff_time` and watches the *loop floor* — the
  minimum dehomogenized norm over the Cauchy loop samples — which exceeds `max_norm` only when the
  entire loop is beyond it (single-sample spikes on legitimate paths cannot trip it).
- **Singularity classification uses the endpoint's spectral-norm condition number** (#344, the B1
  `CondNumThreshold` spec), instead of a mixed-norm estimate that mislabeled borderline endpoints.
- **`frequency_of_CN_estimation` was inert** (#345): the tracker's condition-number refresh counter
  was passed by value, so the estimate never refreshed at the configured cadence.
- **`max_precision_used` is harvested from failed endgames too** (#349); previously a path that
  failed after escalating precision reported as if it had never left double.

### Changed

- **Binomial start points are computed at working precision** (#346) — about 3.2× faster on small
  total-degree solves, with start-point accuracy unchanged (start points are transported values).

_______________________________________________________________________________

## [3.3.2] - 2026-07-13

A `metadata_for` robustness patch (feed a solver result straight back in and it resolves), plus a
documentation showcase and a friendlier docs site.

### Fixed

- **`metadata_for(point)` no longer reports a spurious "more than one distinct solution cluster"
  for a point taken from the solver's own results** (#338).  It now matches against the multiplicity
  REPRESENTATIVES (one per distinct solution) rather than the complete set of coincident copies, and
  defaults the match tolerance to `final_tolerance` (the accuracy each endpoint is computed to) instead
  of the looser same-point clustering tolerance.  Because representatives are at least the same-point
  tolerance apart, a `final_tolerance` window holds at most one — so a solution fed straight back
  resolves to its representative (carrying its `.multiplicity`) and can never be called ambiguous, even
  at a tolerance far tighter than the clustering scale.  This bit hardest on singular points.  See
  ADR-0056.

### Added

- **`solver.same_point_tolerance()`** — the clustering tolerance
  (`final_tolerance × same_point_tolerance_multiplier`), the scale the solver uses to group coincident
  endpoints into multiplicities.  `default_point_match_tolerance()` now returns `final_tolerance`.
- **`metadata_for(..., representatives_only=False)`** — a debugging view that matches against every
  endpoint, including non-representative multiplicity copies.

### Documentation

- **New "Showpieces" gallery** (#341): beautiful renders generated entirely from real tracked data.
  *The Monodromy Loom* — a 3-D braid of solution paths as a parameter loops the discriminant, showing
  a family's monodromy (its Galois action). *The Flight Recorder* — one hard path to a
  multiplicity-35 singular point with the adaptive-precision tracker's full telemetry (endgame spiral,
  precision staircase into mpfr, condition blow-up, step-size sawtooth), plus a system-level view of
  all 35 paths converging.
- **`bertini2.org` lands on the current release** (#340): the docs-site root now redirects straight to
  the latest version instead of a version chooser (which moves to `/versions.html`).

## [3.3.1] - 2026-07-13

A bindings-robustness patch.

### Fixed

- **Storing a multiprecision complex into a real numpy array no longer crashes the
  interpreter** (#336). Assigning a `complex_mp` that carries a nonzero imaginary part into a
  `float64` array — e.g. `M = np.zeros(...); M[i, j] = solver.real_solutions()[k][0]`, whose
  solution coordinates carry ~1e-13 imaginary noise — routed through a `complex_mp → double`
  cast that threw a C++ exception (`"Could not convert imaginary number to scalar."`) *out of*
  numpy's C cast loop, calling `std::terminate()` → SIGABRT. The cast now mirrors numpy's
  builtin `complex128 → float64` behavior (keep the real part, discard the imaginary part).
  See ADR-0055.

### Internal

- The `eigenpy_numpy` numpy-interop test suite was never collected by pytest (its filename
  matched neither `test_*.py` nor `*_test.py`); renamed to `eigenpy_numpy_test.py` so it runs,
  and added regression coverage for the mp → narrower-dtype casts.

## [3.3.0] - 2026-07-12

A large pass on the symbolic and solve/result ergonomics: symbolic substitution and
richer differentiation on the function tree, a friendlier solve/settings surface, and a
typed result taxonomy — plus internal tolerance typing and a documentation refresh.

### Added

- **Symbolic `subs`, richer differentiation, eval ergonomics** (#327): leaf-level
  `node.subs({var: expr})` symbolic substitution; `differentiate(x, 2)` /
  `differentiate([x, x, y])` repeated / sequence differentiation; positional eval with an
  explicit coordinate ordering and a `strict=False` option; constant-power folding (so
  `I**2 → -1`); `node.simplify()`; and a Python `__repr__` that renders `**`.
- **`solver.solve()` returns a `SolveResult`** (#330), and `solver.result()` re-derives it:
  the bare `ZeroDimSolver` / `HomotopySolver` record themselves, so they hand back the same
  records-aware result that `bertini.solve` returns.
- **`ZeroDimResult`** (#331): the typed *answer* of a zero-dimensional solve — the distinct
  finite solutions plus `real` / `singular` / `nonsingular` / `at_infinity` / `nonsolutions`
  views. `SolveResult` is now a records decorator around `.answer`; both `ZeroDimSolver` and
  `HomotopySolver` produce one.
- **Config settings at solve / construction** (#330): `solve(**settings)` and
  `ZeroDimSolver(..., **settings)` (keywords or a `settings={...}` dict), plus
  `get_settings(as_dict=True)` for a flat `{field: value}` view.
- **Random symbolic constants** (#330): `random_real(symbolic=True)` /
  `random_complex(symbolic=True)` / `random_vector(..., symbolic=True)` and
  `symbolics.random_real()` / `random_complex()` produce random constant *nodes*.
- **`metadata_for(pt)` needs no explicit tolerance** (#330): it defaults to the solver's own
  same-point tolerance (`default_point_match_tolerance()`).

### Changed

- **`precision=` is now an integer number of digits; `mptype=` is the precision model**
  (`'double'` / `'multiple'` / `'adaptive'`) everywhere — `ZeroDimSolver`, `bertini.solve`,
  `HomotopySolver`, `user_homotopy` (#330). A *string* `precision=` is still honored as the
  old model alias for one release, with a `DeprecationWarning`.
- **Clearer errors**: a SymPy expression passed to `add_function` now points at
  `bertini.sympy_bridge.from_sympy`; `bertini.solve` on a positive-dimensional
  (under-determined) system raises a clear "needs numerical irreducible decomposition, not
  yet implemented" instead of a raw solver error (#330, #331).
- **`NumErrorT` for tolerance typing** (#329): the error / tolerance type now lives in
  `num_traits.hpp` and types the tolerance parameters across the trackers, endgames, and
  point-comparison helpers (readability only — `NumErrorT` is `double`).
- An informative `repr(solver)` replaces the default object line (#330); the CLI splash drops
  its stale primary-authors block (#332).

### Fixed

- **A `Complex` (constant) node is accepted as a coefficient** (#326, #328), matching the
  existing `complex_mp` behavior.
- **Publishing downloads only the wheel artifacts** (#325): the Doxygen `docs-cpp` artifact
  no longer leaks into the release upload.

### Documentation

- The classic continuation-cartoon figure now auto-fits the finite paths, uses a
  total-degree linear-product start so start points do not overlap, and marks divergences
  with ∞ (inside the axes) and the start system with a green play triangle (#329).

_______________________________________________________________________________

## [3.2.0] - 2026-07-10

Multiprecision linear algebra and a more flexible randomization on the library
side, plus a substantial CI / build-infrastructure pass that markedly shortens
the release cycle.

### Added

- **`bertini.linalg`** (#317): multiprecision linear algebra — LU, QR, and SVD —
  exposed by instantiating eigenpy's own decomposition visitors on the `real_mp`
  / `complex_mp` matrix types (reuse, not re-export). One dtype-agnostic surface
  that also handles `float64` / `complex128` via NumPy. See ADR-0054.
- **`bertini.precision(A, n)`** (#317): set the working precision of an entire
  vector or matrix in a single call.
- **Randomization to any codimension** (#315): `System.randomize(codimension=…)`
  and the `bertini.randomize(system, codimension)` free function generalize the
  square case; all prior `randomize()` calls are unchanged. See the ADR-0025
  amendment.

### Changed

- **Windows drops the `/WHOLEARCHIVE` link workaround** (#287): the
  `ExplicitRKPredictor` Butcher tables are now C++17 inline members in the
  header, so the Windows test executables link normally. See ADR-0052.
- **Faster CI**: wheel builds now run concurrently with the C++ tests — they are
  independent full compiles that shared no artifacts — and the Doxygen build runs
  in parallel with the wheel in the docs workflow (#316, ADR-0053). The Windows
  C++ test build now runs through ccache (#322).
- A **wheel-free Python docstring lint** (`tools/py_doclint.py`) now runs beside
  the C++ Doxygen lint as part of the cheap doc-lint gate (#316).

### Fixed

- Docstrings and tutorials: escaped absolute-value bars that reStructuredText
  misread as substitution references, which had broken the documentation build
  (#313).

_______________________________________________________________________________

## [3.1.0] - 2026-07-09

Quality-of-life and correctness release on top of 3.0.0: a thorough NumPy
interoperability pass for the multiprecision dtypes, a batch of Python UI
ergonomics discovered while writing a real-cellular-decomposition notebook, a
records-recall escape hatch, and CI / documentation-infrastructure work.

### Added

- **NumPy interoperability for `real_mp` / `complex_mp`** (#306): full ufunc
  coverage, sorting slots, safe reductions (`sum`, `prod`, …), and
  tolerance-based comparisons against `float64`; element access returns owned
  copies. See ADR-0051.
- **Python UI ergonomics** (#293–#304, #305):
  - list form `x, y, z = bertini.variables(['x', 'y', 'z'])` and variadic
    `System.add_variable_group(x, y, z)`;
  - random factories `random_real()`, `random_complex()`,
    `random_vector(n, real=…)`, all visible under `bertini.random`;
  - `bertini.is_distinct_up_to(p, q, tol)` — infinity-norm point comparison,
    accepting multiprecision and double vectors alike;
  - `merge_multiplicities=` on every solution accessor (**default True**);
  - `solver.metadata_for(point)` with a call-shape-determined return type;
  - `System.functions()`, `System.copy_functions()`, `System.clone()`;
  - a `group=` projection kwarg on every solution getter and `to_dataframe()`,
    plus the `System.coordinates_of(point, group)` primitive behind it;
  - a sympy auto-sympify bridge, so `sympy.Matrix([...nodes...]).det()` works;
  - `node.eval(point | dict | array)` (previously keyword-only);
  - `variable_group @ coefficients` dot-product sugar;
  - NumPy-native helpers `bertini.real` / `imag` / `abs` / `conj` / `round` /
    `sum` / `norm` / `is_real`.
- **`ZeroDimConfig.recall`** (#308, default `True`): set `False` to force a
  fresh re-track even when an identical ask is already recorded — the escape
  hatch for path observers, benchmarking, and re-verification. Transient: it
  does not affect the run's identity digest.

### Changed

- Solution accessors now **merge multiplicities by default**; pass
  `merge_multiplicities=False` for the raw per-path endpoints (#299).
- CI builds against **prebuilt dependency artifacts** — a custom manylinux
  image plus macOS Boost / eigenpy tarballs (ADR-0049, #282) — cutting build
  time and flakiness.
- Documentation is served from a branch-source `docs-store` with per-version
  snapshots (ADR-0050, #291, #292).

### Fixed

- Random seeding: `set_random_seed` now governs every draw, and real
  projection directions come out actually real (#294).
- NumPy 2.5 reduction use-after-free / uninitialized-slot hazards in the
  eigenpy bindings (#306).
- Windows wheel builds no longer spuriously fail on a precompiled-header
  mtime race (`-fno-pch-timestamp`, #306).

_______________________________________________________________________________

## [3.0.0] - 2026-07-07

The first stable release of the modernized Bertini 2: a rebuilt C++17 core and
a much friendlier Python interface (`import bertini`), with the `bertini2` CLI
shipped inside the wheel. Wheels for CPython 3.10 – 3.14 on Linux
(manylinux_2_34), macOS (arm64), and Windows. This release consolidates ~70
internal PRs; the full themed index, per-PR links, and the upstream issues it
closes are in #238.

### Added

- **Durable, resumable output with provenance.** Every solve writes a
  plain-text structured output directory — content-addressed inputs, one
  results file per run, walkable provenance chains. Solves consult it first, so
  a killed run *finishes* on rerun instead of restarting. Content-identity
  digests mean `seed=42` reproduces the exact same homotopy forever, across
  machines and versions.
- **Content-addressed function trees and systems** — hash-consing / interning,
  symbolic Jacobian, `Seal()` — so equal objects share an identity.
- **Block-structured systems**, with the multihomogeneous start system exposed
  to Python; **user-defined homotopies** with start points; automatic square-up
  and filtering of overdetermined inputs; first-class randomization.
- **sympy bridge** — exact two-way conversion and round-trip solving.
- `to_dataframe()` for pandas; Unicode / emoji identifiers and string-valued
  config fields.
- **The `bertini2` CLI ships inside the wheel** (`bertini2` on your PATH after
  `pip install`), emitting Bertini 1.7-compatible solution files.
- Executable, **doctest-verified** tutorials; a C++ documentation-lint gate;
  and ADRs for load-bearing design decisions. Docs at https://bertini2.org.

### Changed

- **Flattened, discoverable public Python API**; build systems from a list of
  functions (`System.add_functions`) and get solutions back in your own
  coordinates by default; a revamped configuration model.
- **Parallel-by-default solving** on shared memory — **no MPI required** — with
  MPI serialization and reproducibility fixes for multihomogeneous solves at
  scale.
- **~10× faster multiprecision evaluation** (tiered SLP arithmetic,
  allocation-free eval, common-subexpression elimination, a stateful in-place
  multiprecision LU solver) and **5 – 7× faster well-conditioned
  adaptive-precision solves** by staying in double precision where it is
  provably safe (cyclic-5: **11s → 3.5s**).
- Rewritten predict / correct (per-track condition probe, pure kernel);
  adaptive-numeric-type AMP endgames.

### Fixed

- Cauchy endgame **security and pole-zone guards** — never reports success at a
  non-root.

_______________________________________________________________________________

## [2.0.2] - 2026-05-22

Packaging and CI maintenance following the 2.0.1 PyPI debut.

### Changed

- Factored out the documentation workflow and added a `ref` input for tag
  rebuilds (#223, #224).
- Extended the Python support matrix in CI and updated cibuildwheel (#229).
- Minor 2.0.1 follow-up fixes and MPI sync force-push handling (#226, #227).
- Version bump to 2.0.2 and README Python-version refresh (#225, #230).

Full changelog: <https://github.com/bertiniteam/b2/compare/v2.0.1...v2.0.2>

_______________________________________________________________________________

## [2.0.1] - 2026-05-16

First release under the `bertini2` PyPI name. This is the consolidation of
several months of cross-platform packaging work, dependency-compatibility
fixes, documentation, and a final round of precision-handling fixes in the
system / SLP path. Intel macOS is dropped from the supported platform list
for this release (see Removed below).

### Added

- Wheel distributions on PyPI for Linux, macOS (Apple Silicon), and Windows
  across Python 3.9 – 3.13. Linux wheels use the `manylinux_2_28` image; macOS
  and Linux wheels are produced via `cibuildwheel`; Windows wheels bundle
  required DLLs via `delvewheel` and a `windows_dll_manager.py` helper that
  registers DLL search paths before importing `_pybertini`.
- GitHub Pages documentation site: C++ API via Doxygen, Python API via Sphinx.
- Version numbers (b2, GMP, MPFR, Eigen, Boost) exposed from the Python
  package so users can query them at runtime.
- `CHANGELOG.md` itself, following the *Keep a Changelog* format.

### Changed

- **Package renamed: `pybertini` → `bertini` → `bertini2`.** The current PyPI
  name is `bertini2`. Update your `pip install` and import statements
  accordingly; the Python module still imports as `import bertini`.
- Version is now read from `pyproject.toml` by CMake, removing the two-place
  manual sync that previously drifted.
- Build system: `scikit-build-core` configures CMake for the wheel build;
  `pyproject.toml` is the single source of truth for build inputs. The Linux
  wheel build rebuilds Boost.Python and eigenpy per target Python version so
  each wheel ships a matching `libboost_python3X`.
- CI matrix expanded to cover Ubuntu, macOS (Apple Silicon), and Windows for
  Python 3.9 – 3.13. PRs targeting `develop` and pushes to `develop` / `main`
  run the full matrix; other branches run a fast Ubuntu + macos-14 + Python
  3.11 matrix.
- CI build now uses Boost 1.90 (was 1.87). Boost 1.90 pre-seeds
  `thread_default_precision` from the global default and guards against 0,
  removing a class of MPFR abort hazards.
- `boost_system` is now conditionally linked for Boost < 1.89 only (Boost 1.89
  dropped the separate library).
- Eigen requirement updated to `3.3...3.4` with the macOS install switched to
  Homebrew's `eigen@3` formula.
- Tests across platforms now use a unified precision-handling pattern,
  eliminating per-platform skips and conditional precision rituals.
- `MACOSX_DEPLOYMENT_TARGET` co-varies with the runner version so wheels
  produced on macos-14 are loadable on the same OS family.
- TestPyPI publishes on every push to `develop`; PyPI publishes on tagged
  releases (`v*.*.*`) with Sigstore signing and a GitHub Release.

### Fixed

- `System::operator+=` and `System::operator*=` now invalidate the cached
  derivatives flag (`is_differentiated_`), so a subsequent `eval` rebuilds the
  SLP instead of evaluating against stale derivatives. The companion mutators
  `Reorder`, `Simplify`, and `ClearVariables` also invalidate the cache.
- `SLPCompiler::Compile` now seeds the mpfr default precision from the SLP's
  own `precision_` immediately before growing the `mpfr_complex` memory block.
  This prevents a `mpfr_init2(x, 0)` abort when `thread_default_precision()`
  is left at 0 on a fresh thread under Boost ≥ 1.87.
- Python module init seeds `thread_default_precision` so the first MPFR
  construction on the interpreter thread is always valid.
- `eval` overload registration order in the Python bindings changed so the
  double-precision overload is tried before the multi-precision overload for
  ambiguous inputs (e.g. NumPy int64 arrays); this avoids the eigenpy
  `Vec<mpfr>` extractor probing MPFR construction with precision 0.
- `EIGEN_MAKE_ALIGNED_OPERATOR_NEW` moved to the `public` section of the
  `System` class (was incorrectly placed in a non-public section).
- A second `find_package(Boost)` no longer clobbers `Boost_LIBRARIES`.
- The `_pybertini` target is always created (previously it was only created
  when bertini was the top-level CMake project, breaking out-of-tree builds).
- Various MSVC-specific build fixes: `/bigobj`, `/EHsc`, explicit-type
  workarounds for template instantiation issues, Release-config library
  linking.

### Removed

- **Intel macOS (`macos-15-intel`) is not built or tested in this release.**
  A SIGABRT in pytest involving MPFR/Boost surfaces only on Intel runners and
  cannot be reproduced on the maintainer's development hardware. Intel wheels
  will return once a reproducer or upstream fix is in hand. Users on Intel
  Macs should pin to a 1.0.x release or build from source.
- Stale build artifacts and outdated Python-binding documentation removed
  from the repo.

_______________________________________________________________________________

## [1.0.3] - 2025-05-16

Preparation for pypi release with github workflow

### Changed

* make it compatible for Windows

_______________________________________________________________________________

## [1.0.2] - 2025-05-07

Preparation for pypi release with github workflow

### Added

- github workflow for pypi and github release

### Changed

- `publish-to-test-pypi.yml` for handling the comments correctly

### Changed

* merged the pull request for github ci release by @hkmoon in https://github.com/hkmoon/b2/pull/1
* windows release preparation
    * `size_t` is translated into `unsigned long` in linux, mac while `unsigned long long` in windows 10: `core/include/bertini2/eigen_extensions.hpp` and `core/test/classes/start_system_test.cpp` are modified
    * use `clang` of LLVM in Windows since MSVC has different compiling way for `template`
    * use `--no-isolation` for `scikit-build` in Windows
* For linux wheel naming convention, we cannot use x86_64, x86_i386 anymore for pypi repository. https://peps.python.org/pep-0600/
    * use `auditwheel` for it

### New Contributors
* @hkmoon made their first contribution in https://github.com/hkmoon/b2/pull/1

_______________________________________________________________________________

## [1.0.1] - 2025-05-06

This is the initial version of the project.

### Added

- The base project

[CHANGELOG.md]: https://keepachangelog.com/en/1.1.0/
[Semantic Versioning]: http://semver.org/

<!-- markdownlint-configure-file {
    "MD022": false,
    "MD024": false,
    "MD030": false,
    "MD032": false
} -->
<!--
    MD022: Blanks around headings
    MD024: No duplicate headings
    MD030: Spaces after list markers
    MD032: Blanks around lists
-->
