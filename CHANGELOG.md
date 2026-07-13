
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

## [3.3.2] - 2026-07-13

A `metadata_for` robustness patch: feed a solver result straight back in and it resolves.

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
