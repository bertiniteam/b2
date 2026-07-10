# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Bertini 2 (b2) is a C++17 numerical algebraic geometry library with Python bindings. It implements homotopy continuation for solving polynomial systems, including path tracking, endgames (power series, Cauchy), and start systems (total degree, multihomogeneous). The Python package is published to PyPI as `bertini`.

## Build Commands

### C++ Core Library (CMake + Ninja)

```bash
# Configure (from repo root)
cmake -DENABLE_UNIT_TESTING=ON -G Ninja -B build -S .

# Build everything (core library, bindings, exe, tests)
cmake --build build --target all --config Release

# For iterative C++ work, build just the core library first -- it's much faster, and
# the heavy Boost.Python/eigenpy bindings (`_pybertini`) only need rebuilding when the
# bindings themselves change:
cmake --build build --target bertini2 --config Release

# Run all C++ tests
ctest --test-dir build/core
```

### Python Wheel (scikit-build-core)

```bash
python3 -m build --wheel
```

The build system uses `scikit-build-core` (configured in `pyproject.toml`). The wheel build invokes CMake internally.

### macOS Dependencies (Homebrew)

```bash
brew install gmp mpfr libmpc eigen@3 eigenpy boost boost-python3
```

### Linux/CI Dependencies (conda)

Use `environment.yml` (Ubuntu) or `environment-win.yml` (Windows) with conda/mamba/micromamba.

## Running Tests

### C++ Tests (Boost.Test)

Tests are defined in `core/CMakeLists.txt`. Each test suite is a separate executable:

```bash
# All tests
ctest --test-dir build/core

# Individual test executables (after build)
./build/core/test_classes
./build/core/test_blackbox
./build/core/test_classic
./build/core/test_endgames
./build/core/test_generating
./build/core/test_nag_algorithms
./build/core/test_nag_datatypes
./build/core/test_pool
./build/core/test_tracking_basics
./build/core/test_settings
```

### Python Tests

`pytest` is the single way to run the Python tests (the suites are plain pytest
functions + fixtures; the old `unittest` `TextTestRunner` aggregator scripts are gone):

```bash
pytest python/test/
```

The multiprecision default precision is **global mutable state**
(`bertini.default_precision(n)`). An **autouse fixture in `python/test/conftest.py`**
(`_reset_precision`) resets it to a known baseline (`DEFAULT_TEST_PRECISION = 30`) before
every test and restores it afterward, so no test can inherit a neighbor's precision — do
**not** re-introduce per-test `default_precision(...)` setup. To override the precision for
a specific test, use the `precision` fixture (parametrize it indirectly, e.g.
`@pytest.mark.parametrize("precision", [30, 50, 80], indirect=True)` with a
precision-derived tolerance). When adding or debugging precision-sensitive tests, run the
file on its own (`pytest python/test/classes/<file>.py`) to confirm it does not depend on
cross-test state.

### Tutorial doctests and doc figures

Tutorial code blocks are executable and CI-tested: `python -m sphinx -b doctest source
build/doctest` from `python/docs/` (the output goes under `python/docs/build/`, which is
gitignored).  Behavior changes that alter tutorial output (path counts, verdicts, printed
tables) will fail there -- update the prose *and* the doctest expectations together.

Committed tutorial figures ship as **both** `.png` and `.svg` and are regenerated ONLY
through `tools/refresh_doc_artifacts.py`.  Mind the flags: the default (no flags) runs the
**timing benchmark tables**, not plots -- to redraw one figure use
`python tools/refresh_doc_artifacts.py --plots --only <name>`.

## Architecture

The project has three layers, built in order:

1. **`core/`** -- C++ shared library (`libbertini2`). Header-only-heavy design under `core/include/bertini2/`. Key subsystems:
   - `function_tree/` -- Expression tree (nodes, operators, symbols) for *building and representing* polynomial systems. Nodes no longer evaluate: node-level recursive evaluation was removed -- the **SLP (`straight_line_program`) is the sole evaluator** (compile a system once, evaluate the compiled program). `Function`/`Handle` are gone; `NamedExpression` is the sole root node. See ADR-0027 (SLP: immutable Program + per-thread Memory) and ADR-0028 (named-node taxonomy).
   - `system/` -- Polynomial system construction, start systems (`start/total_degree_linear_product.hpp`, `start/mhom.hpp`), patches, slices
   - `trackers/` -- Path tracking (fixed-precision and adaptive-precision trackers, predictors, Newton correctors)
   - `endgames/` -- Power series and Cauchy endgames for singular endpoint handling
   - `nag_algorithms/` -- Higher-level algorithms (zero-dim solve; numerical irreducible decomposition is *framework scaffolding* -- not yet implemented, its `Solve()` throws)
   - `records/` -- The structured output directory (record schema `b2rec/1`, spec at `docs/records/b2rec-1.md`): durable, self-describing run records with full provenance and resume-by-recall.  `OutputDirectory` = append-only JSONL `history/` + content-addressed `definitions/` under kind folders (`systems/`, `configs/`, `givens/`, ...); `README.txt`/`INDEX.txt`/`results.json` are derived, rebuildable views, never truth.  Solvers write through the emission seam in `nag_algorithms` (`RecordTo(...)` or the ambient `BERTINI_RECORDS_DIR`); `bertini.solve` is *ensure-answered* -- paths already recorded for an identical ask (system digest + settings digest + seed) are **recalled**, not recomputed.  See ADR-0042..0047 and "Persistent digests" below.
   - `io/parsing/` -- Boost.Spirit Qi parsers for classic Bertini input format.  `io/json_writer.hpp` renders a System's parts (variable groups, functions, blocks, patches) as JSON for the records archive -- classic syntax is an INPUT/compat format only (it cannot express block structure) and never appears inside records.
   - `blackbox/` -- CLI executable entry point (CMake target `bertini2_exe`, binary named `bertini2`)

2. **`python_bindings/`** -- Boost.Python + eigenpy bindings producing `_pybertini` native module. Each `*_export.cpp` wraps the corresponding C++ subsystem. Depends on `eigenpy` for NumPy/Eigen interop.

3. **`python/bertini/`** -- Pure Python package that wraps `_pybertini` into a user-friendly API. Submodules mirror the C++ structure: `function_tree`, `system`, `tracking`, `endgame`, `parse`, `nag_algorithm`, `multiprec`, etc.

## Persistent Digests -- the forever contract

Systems and configurations have **stable cross-session identities**: SHA-256 digests over
versioned canonical text encodings (`b2sysenc/<n>` for Systems, ADR-0042; `b2cfgenc/<n>`
for configs, ADR-0043; seeds are rooted per `b2rand/1`, ADR-0044).  Records reference
objects by digest, so equal objects must digest equally *forever* -- across machines,
compilers, and versions.  Rules that follow:

- **The canonical texts are digest preimages, never presentation.**  Exact values only
  (doubles as IEEE-754 bit patterns `d64:<16 hex>`, rationals via exact `.str()`, enums by
  fixed name tables).  Human-readable JSON views are *derived* from them and free to
  change only if the transform is invertible.  Never "improve" an encoding for
  readability.
- **Adding a field to a config struct REQUIRES extending its encoder** (they are
  hand-maintained mirrors, like `serialize`), and an identity-affecting encoder change
  REQUIRES, in the same commit: bump the version token, regenerate the golden digest
  fixture (`core/test/classes/data/{config,system}_digest_fixture.txt` -- run the test,
  it prints the new digests), and **append** a line to the version registry
  (`core/test/classes/data/{config,system}_encoding_versions.txt` -- the failing test
  prints the keyspace hash to append).  Registries are append-only: never edit an
  existing line; line *k* carries version suffix *k*.  Tests enforce all of this.
- Deliberately excluded from identity: the RNG seed (its own slot in the ask, beside the
  config digest), `ZeroDimConfig::num_threads` (thread count must not change what was
  computed), and all transient eval state.

## Key Dependencies

- **GMP/MPFR/MPC** -- Arbitrary-precision arithmetic (found via custom CMake modules in `cmake/`)
- **Eigen 3** -- Linear algebra. **Not** pinned in cmake (`find_package(Eigen3)`, no version floor). In practice the version is coupled to the eigenpy build: the wheel CI builds **eigen 3.4.0** and then builds eigenpy against it (a dev env may use newer, e.g. `eigen=5.0.1`). Newer Eigen is welcome -- we *want* upstream improvements -- but it must be matched by an eigenpy built against the same Eigen (they share Eigen types across the binding ABI).
- **Boost** (serialization, filesystem, log, graph, regex, timer, chrono, thread, unit_test_framework, python) -- no minimum version pinned in cmake; `boost_system` is conditionally linked for Boost < 1.89 (header-only from 1.89). Boost.Python is ABI-locked to one CPython version, so it is built **per target Python** -- but this now happens **once, up front, in the prebuilt CI deps** (Linux image / macOS tarballs, ADR-0049), *not* recompiled in every wheel run.
- **eigenpy** -- Eigen/NumPy bridge for Python bindings. **Prebuilt** into the CI deps (image/tarballs), not built-from-source per run, against the chosen Eigen -- eigenpy and bertini must use the *same* Eigen. The version is single-sourced in `.github/ci-deps-versions.env` (`EIGENPY_VERSION`, currently `3.13.0`); eigenpy >= 3.13 sets the Python floor (>= 3.10).
- **jrl-cmakemodules** -- CMake helper macros (auto-fetched via FetchContent if not found)

## Build System Notes

- The root `CMakeLists.txt` uses `jrl-cmakemodules` (fetched automatically). It currently only adds `core/` as a subdirectory; `python_bindings/` and `python/` subdirectory calls are commented out (the wheel build via scikit-build-core handles them).
- `pyproject.toml` configures scikit-build-core: wheel packages from `python/bertini/`, build dir is `bld/`.
- Cross-platform: Linux builds in a **custom prebuilt manylinux image** (`ghcr.io/bertiniteam/b2-manylinux-deps`, ADR-0049) + `auditwheel`; macOS uses Homebrew + **prebuilt Boost/eigenpy tarballs** (the `ci-deps` release); Windows uses conda-forge (which already ships prebuilt Boost/eigenpy) + clang-cl (MSVC has template compilation issues).
- `-Werror` is disabled globally. `-pedantic` is stripped from flags.

## CI/CD

- `.github/workflows/build_and_test.yml` -- Builds wheels on Ubuntu/macOS/Windows and runs tests. Triggered by pull requests and pushes to `develop`/`main`.  Docs-only changes (`**/*.md`, `python/docs/**`, `Doxyfile`, `VERSION`) are in `paths-ignore` — but beware: for `pull_request` events GitHub evaluates the filter against the PR's **entire** diff, so a docs-only *push to a PR that also carries code* still re-runs CI (and, via the concurrency group, cancels the PR's in-flight run — batch docs commits with code pushes). Only pushes to `develop`/`main` and PRs whose whole diff is docs-only are skipped.  MPI is verified per-platform: the CLI smoke test requires `mpirun -n 2` to produce the *same solution count* as serial, for both start-system families.  **Wheel builds run CONCURRENTLY with the C++ test jobs, not behind them** — they are independent full compiles (prebuilt deps + cibuildwheel vs. boost-base-from-source + ctest; different ccache namespaces), so the only gate they share is the cheap `doc_lint`. A C++-test failure still fails the run at the workflow level (and blocks publish); it just no longer serializes the wheels behind ~7 min (unix) / ~22 min (windows) of tests. Do **not** re-add `test_cpp_*` to the wheel jobs' `needs`.
- `.github/workflows/doc_lint.yml` -- A cheap wheel-free doc-correctness gate that **the build matrix depends on** (it runs first; if it fails, nothing compiles). Two parallel jobs, one per language — **run both locally before pushing docs/source changes**, or CI will bounce the whole build:
  - **C++ (`bash tools/doclint.sh`)** — needs `doxygen` on PATH (`brew install doxygen`). Two passes: (1) *correctness* — `@param` names must match signatures, no doc blocks on removed signatures, no unresolved `\ref`/`\cite` (zero tolerance); (2) *undocumented ratchet* — the count in `tools/doc_undocumented_baseline.txt` may only **decrease** (currently `0`, so **every new public C++ entity — including each struct data member — needs a Doxygen comment**, e.g. `///< ...`). `bash tools/doclint.sh --update-baseline` locks in a reduction.
  - **Python (`python tools/py_doclint.py`)** — pure-stdlib `ast`, no build/wheel/pip. Same two-pass shape over `python/bertini/**/*.py`: (1) *correctness* — a numpydoc `Parameters` entry naming a parameter absent from the signature fails (zero tolerance); (2) *undocumented ratchet* — undocumented public entities (module / public class / public function / public method; dunders excluded) counted against `tools/py_doc_undocumented_baseline.txt`, may only **decrease**. `python tools/py_doclint.py --update-baseline` locks in a reduction. Only lints the pure-Python layer; native `_pybertini` docstrings (Boost.Python strings in `python_bindings/src/*_export.cpp`) are invisible to a wheel-free tool.
- `.github/workflows/publish.yml` -- Publishes to TestPyPI on `develop` push, PyPI on version tags (`v*.*.*`) with Sigstore signing and GitHub Releases.

### Linux wheel test coverage

Linux wheels are built inside the **custom prebuilt deps image** (`ghcr.io/bertiniteam/b2-manylinux-deps`, an `manylinux_2_34`/AlmaLinux 9 base with Boost+eigenpy baked in; set via `CIBW_MANYLINUX_X86_64_IMAGE`; ADR-0049). The **full pytest suite runs on all three platforms** — on Linux it runs *inside* that container via `CIBW_TEST_COMMAND_LINUX`, and on macOS/Windows via the host-runner test jobs. `mpi4py` is installed in every test env (the image ships OpenMPI), so the MPI test modules **run at 1 rank rather than skip** — the suites are 0-skip on all three platforms.

This was not always so: for a while Linux ran an import smoke test only, because the suite was SIGABRT/SIGSEGV-crashing — a crash *misattributed* to the older `manylinux_2_28` container's MPFR 3.1.6. The real cause is a **version-independent** bug (uninitialized `mpfr`/`mpc` numpy slots), now fixed in the bindings. Do **not** try to fix Linux test crashes by bumping MPFR or the manylinux image (that was tried and does not work) or by building MPFR from source (specifically out of bounds). See `docs/adr/0006-eigenpy-uninitialized-numpy-slot-guards.md` for the fix and `docs/adr/0003-manylinux-no-full-pytest.md` for the (now reversed) smoke-test stopgap and its history.

## Python Bindings — Known Pitfalls

### eigenpy writable Ref + adjacent scalar (ADR-0001)

Never place a writable `Eigen::Ref<Vec<mpc_complex>>` argument **adjacent** to a `mpc_complex const&` scalar argument in a Boost.Python binding. eigenpy's from-Python converter for the writable Ref writes into a static rvalue-converter slot that overlaps with the storage for adjacent `const&` scalars, corrupting them. The corrupted `mpc_complex` then triggers `MPFR_ASSERTN` → SIGABRT.

**Rule:** If a binding takes a writable `Eigen::Ref<Vec<ComplexT>>` and also needs scalar `ComplexT` args, pass the scalars **by value**:

```cpp
// WRONG — start_time/end_time get corrupted
SuccessCode wrap(Eigen::Ref<Vec<ComplexT>> result,
                 ComplexT const& start_time,   // ← adjacent const& scalar
                 ComplexT const& end_time);

// CORRECT — by-value copy is taken before the Ref converter runs
SuccessCode wrap(Eigen::Ref<Vec<ComplexT>> result,
                 ComplexT start_time,           // ← by value
                 ComplexT end_time);
```

Single-argument bindings and read-only `Vec<T> const&` bindings are unaffected. See `docs/adr/0001-eigenpy-writable-ref-scalar-by-value.md`.

## Architecture Decision Records

`docs/adr/` contains ADRs for load-bearing design decisions — where the *why* would not be obvious from reading the code. Check there before undoing anything that looks strange.

**Assess ADR need as part of every PR — both when opening it and before closing/merging it.** Ask: does this change *make* a load-bearing decision (a non-obvious tradeoff someone might later "helpfully" revert), or does it *supersede/invalidate* an existing ADR? If yes, the ADR belongs in the **same PR** as the change:
- **New decision** → add `docs/adr/NNNN-title.md` (next free number; template in `docs/adr/README.md` — Status/Date/Context/Decision/Consequences) **and** a row in the README index table. Write the *Consequences* as a guardrail ("do not re-add X", "do not move Y back to a .cpp").
- **Supersedes an existing ADR** → set the old one's `**Status:**` to `Superseded by ADR-NNNN` (and annotate its index row `_(superseded by NNNN)_`), then write the new one.
- Bar for a *new* ADR: a reviewer reading only the diff couldn't reconstruct *why*, or the obvious-looking "cleanup" reintroduces a bug the change fixed. Pure mechanics, refactors, and bug fixes with a regression test usually don't need one.

## Conventions

- C++ standard: C++17. Headers use `.hpp` extension.
- License: GPL v3 with additional terms (see `licenses/`, `core/ADDITIONAL_GPL_TERMS`).
- Version single source of truth is the top-level **`VERSION`** file; `pyproject.toml` reads it *dynamically* (scikit-build-core), and `publish.yml`'s `check_version` asserts the release tag matches it. Read at runtime via `importlib.metadata.version("bertini2")`. Bumping requires a PR (develop is ruleset-guarded), but `VERSION` is in `paths-ignore` so a version-only PR runs no CI.
