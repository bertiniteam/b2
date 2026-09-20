
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

## [4.0.0] - unreleased

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

- **A system can say which of its coordinates are auxiliary** (#403).  Bertini judges a point by
  its largest coordinate, three times over: the tracker truncates a path past
  `path_truncation_threshold`, the endgame abandons one past `Security.max_norm`, and the solver
  calls an endpoint past `endpoint_finite_threshold` infinite.  Several standard constructions
  carry coordinates that exist for the construction rather than the answer -- a critical-point
  system in null-vector form needs a patch on its null vector, and that patch fixes a
  normalization nobody chose, so the block's direction is meaningful and its magnitude is an
  artifact.  Judged on such a coordinate, genuine solutions were declared at infinity, and there
  is no threshold that would have been right, because the quantity is not governed.
  `System.set_auxiliary_variable_groups` and `System.set_auxiliary_coordinates` name the
  coordinates the question is not about; `System.is_finite` and `System.is_real` are the
  judgements, and the system owns them, so the tracker, the endgame and the classifier stop each
  deriving their own.  Realness had the identical defect one line from finiteness and is fixed
  with it.  Nothing else changes: an auxiliary coordinate is tracked, recorded and returned
  exactly as before.  Which coordinates are auxiliary is part of the system's content identity,
  since the tracker truncates on the ones that are not -- so declaring one changes the digest and
  makes a new ask.  See ADR-0062 and the new "Auxiliary coordinates" tutorial.
- **The path-crossing verdict reaches the records** (#365).  The midpath check compares every path
  against every other at the endgame boundary and re-tracks the ones that appear to have jumped
  onto a neighbour; a path it finally gives up on has carried `crossing_unresolved` in the solver's
  metadata since this release's earlier work, but nothing persisted it.  A records directory could not be
  asked which paths were flagged, and recall brought such an endpoint back with the flag cleared --
  the one point the library could not vouch for, arriving vouched for.  Now the flag is part of the
  path record and comes back with it, and each recording solve also writes one `midpath` record
  carrying the check itself: whether it passed, how many crossings it found, how many re-tracks it
  spent.  Read them with the new `bertini.crossing_checks()` and the `crossing_unresolved` column
  of `bertini.tracks()`.  A run with no `midpath` record never reached the check -- a solve cut
  short skips it -- which is deliberately distinguishable from a run that checked and found
  nothing.  Records written before this simply lack the field, and read as unflagged, which is all
  they were ever able to say.
- **A configuration reference page: every setting there is, with its default** (#406).  There was
  no single place that answered "what settings exist, and what do they default to?".  The class
  listings name fields but no defaults, the config classes are spread over three modules with no
  index, and a handful of tracker settings are set by method and belong to no config struct at
  all, so every "here are the configs" listing misses them by construction.  The new page, in the top-level reference navigation, collects all of it and is
  **read out of the library when the docs are built**: adding a config field, changing a default or
  rewriting a field's docstring updates the page with no documentation change.  A hand-maintained
  table drifts from the code silently, which is what the issue was filed about.  The settings with
  no config struct are hand-written, in their own section, since no amount of introspecting config
  classes will find a method.
- **Ctrl-C stops a solve.**  A solve releases the interpreter lock and runs on the calling
  thread, so a keyboard interrupt was noted by CPython and ignored until the solve finished;
  in a notebook the kernel was simply trapped, and killing the process was the only exit.
  Now the solve installs a signal handler for its duration, every path being tracked notices
  the request between steps, the solve unwinds cooperatively, and `KeyboardInterrupt` is
  raised -- with the solver left exactly as the stop found it.  Finished paths stay finished
  and readable; paths in flight are abandoned and read `SuccessCode.ExternallyTerminated`, a
  value that had been in the enum, unproduced, for years; paths not yet begun stay
  `NeverStarted` and are not recorded.  `solver.was_stopped_early()` and
  `solver.num_paths_never_started()` say what happened, and the midpath check is skipped on a
  partial set rather than comparing paths that never ran.  Because the records hold every
  finished path, re-running the same solve recalls them and tracks only the rest: an
  interrupt is resumable.  Also programmatic: `bertini.request_stop()` from another thread
  stops the running solve without raising.  In C++, `bertini::RequestStop()`,
  `StopRequested()`, `ClearStopRequest()` and the RAII `ScopedStopRequest`; a bare tracker
  honours the request too, since it is the tracker that checks.  Not applied to the MPI
  solve.
- **Wall-clock budgets.**  `max_path_wall_clock_duration` on the solver (seconds, 0 = none)
  gives every path a budget; a path that has not finished when it runs out is abandoned
  between steps with `SuccessCode.WallClockLimitReached`, stamped with where it got to, and
  recorded with the budget that stopped it.  `max_solve_wall_clock_duration` budgets the whole
  `solve()` call instead: once it runs out no further path starts and any in flight is
  abandoned, and the solve reads exactly as if Ctrl-C had been pressed (`was_stopped_early()`,
  `ExternallyTerminated` / `NeverStarted`, re-tracked on recall).  The overrun is at most one
  step.  Both are plain numbers of seconds and both are off by default.  A bare tracker
  can be limited on its own: `set_max_wall_clock_duration(seconds)` /
  `clear_max_wall_clock_time()` in Python, `SetMaxWallClockTime(time_point)` /
  `SetMaxWallClockDuration(duration)` / `ClearMaxWallClockTime()` in C++.  The primitive is a
  deadline, not a duration, because an endgame issues hundreds of tracking calls for one path
  and a per-call budget would restart with each.  The budget is deliberately not part of the
  run's identity: a 60 s and a 61 s run are the same ask.  Wall-clock time is machine
  dependent; the recorded stamp is the machine-independent account.
- **A tutorial on stopping a solve and budgeting paths** ("Stopping a solve, and giving paths a
  budget", under "Solver settings & multiprecision"): Ctrl-C and `request_stop`, the stamp an
  abandoned path leaves, the per-path and whole-solve budgets, and how the recall policy reads
  an abandonment.  Doctested with budgets of a nanosecond and an hour, so it is the same on
  every machine.
- **A recall policy, in its own config.**  What already-recorded work counts as done is now
  `RecordsConfig.recall`, a `RecallPolicy`: `Nothing` tracks every path fresh (the old
  `recall = False`), `Completed` (the default) reuses completed paths and re-tracks abandoned
  ones -- except a path a wall-clock budget cut off, which is reused while the current budget
  asks no more patience than the one that abandoned it, and re-tracked once the budget goes up
  -- and `Everything` reuses every recorded outcome as recorded.  An interrupted path is always
  re-tracked.  A bool is still accepted (`True` is the default policy).  It moved out of
  `ZeroDimConfig` because the question is the same for every algorithm that records paths.
- **A path that did not succeed says where it got to.**  Its metadata now carries
  `latest_path_point`, the point the tracker was at when it gave up (in the solver's internal
  coordinates), with `final_time_used` holding the matching time -- previously a failed path
  reported a time of zero and no point.  Every path also carries `num_successful_steps` and
  `num_failed_steps`, the predictor-corrector steps over the whole path, pre-endgame and
  endgame together, so the cost of a path is visible without timing it; and
  `max_precision_used` is honest for fixed-precision solves too, where it read zero.  The
  stamp is recorded with the path and comes back on recall.  Two corrections came with it: a
  path whose endgame failed used to be handed the previous path's approximation as its
  "solution" (the accessor was stale); its solution slot is now empty, and the stamp is where
  its point lives.  In C++ the tracker exposes the tally that made this possible:
  `ResetCumulativeStepCounts()`, `CumulativeSuccessfulSteps()`, `CumulativeFailedSteps()`,
  counting across `TrackPath` calls (an endgame issues hundreds per path) until told to start
  over.
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

- **The tracker's own settings are a config** (#457).  The predictor, the tracking tolerance and
  the path truncation threshold were bare members of the tracker, each with its own setter and
  getter and no config struct.  They were therefore absent from `get_settings`/`set_settings`, from
  the settings digest, and from the generated configuration reference; `path_truncation_threshold`
  had no route through a solver at all, so the only way to change it was to reach for the tracker
  and call a method.  They are `tracking::TrackerConfig` now, so `solver.update(predictor=...,
  tracking_tolerance=..., path_truncation_threshold=...)` works, a settings bundle carries them,
  and the reference page lists them.  The methods still work and read and write the config.  The
  tracker derives what it had been caching -- the digits implied by the tolerance, the predictor
  object -- when it starts a path rather than when a setter is called, so a value set by any route
  takes effect.  In the same move `final_tolerance` lost its duplicate: it was a field of both the
  solver's `TolerancesConfig` and the endgame's `EndgameConfig`, with the solver's pushed onto the
  endgame at solve time, and it now lives only on the endgame, which is what achieves it.  No field
  name is shared by two configs any more.  Setting `final_tolerance` on a solver still works and
  reaches the endgame's config.  **Identity:** both moves change the settings text a solve folds
  in, so `b2cfgenc/4`'s registry line was edited in place -- 4.0.0 is unreleased, which is exactly
  the case ADR-0061 covers.
- **Encoding versions now count encodings, not commits: at most one per released version.**
  `b2sysenc/<n>` and `b2cfgenc/<n>` are part of the digest preimage, so moving one makes every
  existing record a different ask -- and the old "append a line, never edit one" rule made the
  number climb with development rather than with releases.  The evidence was already in the
  registry: `b2cfgenc/1` and `/2` were both minted while 3.0.0 was being developed and neither
  ever shipped, so the number claimed three encodings existed when one had ever been written.
  Now each registry line records the release it first shipped in, no two lines may claim the
  same release, and while a release is unreleased its line is edited in place rather than
  superseded; once it ships it is history and is never touched.  The consequence for a user is
  a promise the version can keep: `b2sysenc/2` means "the second encoding b2 has ever written",
  and a record names the release that wrote it.  The registry also gained a column for the
  *composition* of the settings text -- which configs a solve folds in -- because the existing
  hash could not see it: `b2cfgenc/3` and `/4` are byte-identical there, `/4` having been
  exactly a composition change.  (ADR-0061)
- **Path-crossing detection joined the ask: config encoding `b2cfgenc/4`.**  `MidPathConfig`
  decides when two paths at the endgame boundary count as the same point, and so decides which
  paths get re-tracked and what the solve returns -- but it was not in the settings text, so two
  solves that disagreed about `same_point_tolerance` were the same ask and would recall each
  other's results.  It is in the text now, appended last, since the order is extended by
  appending and never by reordering.  **Every settings digest therefore changes**: runs recorded
  by an earlier version are a different ask now and will be recomputed rather than recalled.
  The gap that let this sit unnoticed is closed too.  The golden fixture and the version registry
  both watch the encoders -- what one config struct turns into -- and neither watched the
  composition, which configs a solver folds into its settings text.  A new test pins that list in
  order, so adding, removing or reordering a config is a named failure rather than a silent
  change of identity.
- **One homotopy builder, named for the mathematics: `straight_line_homotopy`** (#371).  It
  replaces both `blend_homotopy` and `coefficient_parameter_homotopy`, which built the same
  object and differed only in gamma.  "Blend" named the implementation, a blend block, rather
  than the mathematics; and a coefficient-parameter homotopy is not a different construction, it
  is this one with gamma fixed at 1, so the second name claimed a specialization it did not have.
  `straight_line_homotopy(target, start, *, path_variable='t', gamma=None)` is
  `(1-t)*target + gamma*t*start`, with the gamma trick by default.  **Both old names are gone**:
  read `blend_homotopy(a, b)` as `straight_line_homotopy(a, b)`, and
  `coefficient_parameter_homotopy(a, b)` as `straight_line_homotopy(a, b, gamma=1)`, where the
  `gamma=1` is now visible rather than hidden in the choice of function.  `gamma` also accepts a
  plain int, so `gamma=1` needs no `Integer` wrapper.
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
  operands in that order; the encoding also carries which coordinates a system has declared
  auxiliary (#403).  Both changes land under the one version, which is the point of ADR-0061:
  a release gets one encoding, not one per change.  Records written by earlier versions carry
  `b2sysenc/1` and are not recalled; the "records read forever" promise is withdrawn until the
  identity of the algorithm itself is part of a record's ask (#420).
- **The default endgame is the power series endgame**, matching Bertini 1's documented default
  (`EndgameNum: 1`), at every choice point: the configuration default (so the CLI and a classic
  input without `endgamenum`), and the Python factories `ZeroDimSolver`, `HomotopySolver`,
  `user_homotopy`, `parameter_sweep` and `bertini.solve`.  Measured on a regeneration workload,
  the Cauchy default cost 190 s where power series takes 2.5 s, because slowly diverging paths
  are a Cauchy-specific pathology.  Cauchy remains available by explicit request
  (`endgamenum: 2`, `endgame='cauchy'`).  A solve that relied on the default is now a different
  computation, so its records are new asks rather than recalls.  See ADR-0058.

### Fixed

- **`InEGOperatingZone` means the same thing in both endgames** (#402).  The event says the path
  has reached the asymptotic regime, where the Puiseux model the endgame is built on dominates --
  which is what makes the samples usable quantitatively, for anyone fitting a power law against
  `|t|` to see what vanishes at the root.  The Cauchy endgame emitted it once, on a real test: its
  c/k estimates settling.  The power series endgame emitted it at the end of every successful
  advance, with no test at all, so there it meant "advanced", and a consumer that trusted it got a
  silently wrong answer under one flavor and a right one under the other.  Both run the same test
  now.  It reads only the geometrically spaced approach samples, which every flavor already keeps
  -- the power series endgame's `samples_` and the Cauchy endgame's `pseg_samples_` are the same
  window under two names -- so it lives on the base endgame, along with the fixed probe direction
  that keeps consecutive estimates comparable.  `minimum_for_c_over_k_stabilization` and
  `num_needed_for_stabilization` move from `CauchyConfig` to `EndgameConfig` with it, which is
  where a setting every endgame reads belongs; scripts that set them on the Cauchy config need the
  one-line change.  Every config digest moves, under the existing `b2cfgenc/4` rather than a new
  version, per ADR-0061.
- **Four settings had no default at all.**  `SharpeningConfig.sharpendigits` and
  `RegenerationConfig`'s three slice tolerances were declared without initializers, so a
  default-constructed config held indeterminate values and reading one was undefined behaviour --
  and there was no default to document, which is how the omission surfaced while generating the
  configuration reference.  They are initialized now: sharpening off (`0`, Bertini 1's default),
  and the slice-moving tolerances matching the main tracking tolerances they shadow, which is
  where Bertini 1 takes its own `SliceTol*` defaults from.  Nothing computes differently: the
  classic parser was the only thing that ever set them, and regeneration is scaffolding.
- The classic parser reports two mistakes it used to let through (#441).  A name declared and
  never defined -- `function f, g;` with only `f` defined -- reached a `std::map::at` lookup and
  surfaced as `map::at` (an `IndexError` in Python), naming neither the forgotten name nor the
  problem; it now says which name was declared and never defined, and what it was declared as.
  A second `pathvariable` declaration was accepted, quietly overwriting the first; a system has
  exactly one path variable, so the second declaration is now refused.  `System::AddPathVariable`
  stays permissive, so a caller building a system programmatically may still change its mind --
  the enforcement is on the declaration, in the input file.
- Every setting a solve uses is settable through the solver's own settings surface (#364).  The
  tracker's (`max_step_size`, the Newton counts, the precision configs) and the endgame's
  (`num_sample_points`, `sample_factor`, the security and flavour settings) used to be reachable
  only by fetching the sub-object and round-tripping its config -- `solver.get_endgame()
  .get_endgame_settings()`, modify, set back -- which is easy to get wrong and unreachable from
  any code path that only carries a settings dict.  Now `solver.update(num_sample_points=6,
  max_step_size="0.05")`, `solve(**settings)`, `configure()`, `get_settings()` and
  `set_settings()` all reach them, routed to whichever object holds each field.  One field name
  is shared by two configs -- `final_tolerance`, on the solver's tolerances and on the endgame --
  and the solver's own keeps winning, with the endgame's reachable as
  `configure(endgame={'final_tolerance': ...})`.  Endgames also gained the config surface
  trackers and solvers already had (`get_config`, `set_config`, `config_types`, and with them
  `update` / `configure` / `get_settings`).
- A misspelled setting suggests the one you meant.  A settings call names its fields as keywords,
  so a typo is silent where it is written and the error message is the only place it can be
  caught; listing the valid names says what exists, not what was meant, and those lists run to
  dozens of entries.  `solver.update(final_tolerence=...)` now answers `Did you mean
  'final_tolerance'?`, for fields, for config names, and whether the field belongs to the solver,
  its tracker or its endgame.  Nothing close by means no guess, and the valid names are still
  listed either way.
- The last settings that no caller could name are nameable (#364).  `same_point_tolerance`, which
  governs path-crossing detection, was held only by the midpath checker and so appeared in no
  config list the solver publishes; the solver holds it now and pushes it into the checker at the
  start of each solve, which makes it settable like anything else and gives the setting one home
  rather than two.  `EndgameConfig.refine_when_increasing_precision`, and `ZeroDimConfig`'s
  `initial_ambient_precision` and `path_variable_name`, existed in C++ but had never been exported
  to Python at all.  The path variable's name is read when a solver builds its homotopy, so it
  belongs on a config you hand to a constructor; its docstring says so, since setting it on a
  solver that already built one would leave the two disagreeing.
- An exact-rational setting takes every exact spelling (#364).  `sample_factor` is a rational,
  and accepted only a `rational_mp`: a string went to the float parser and came back with
  `Unable to parse string "1/10" as a valid floating point number`, which is a confusing thing to
  be told about a field that is not a float.  It now takes `'1/10'` and `'0.1'` -- the same
  number, both exact, since a decimal is read as the rational it denotes and never as a rounded
  binary float -- and an `int`, a `fractions.Fraction`, or a `rational_mp`.  A Python float is
  still refused, now saying why.  Rationals are also picklable now, so a settings bundle
  containing one still travels.
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
