# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Bertini 2 (b2) is a C++17 numerical algebraic geometry library with Python bindings. It implements homotopy continuation for solving polynomial systems, including path tracking, endgames (power series, Cauchy), and start systems (total degree, multihomogeneous). The Python package is published to PyPI as `bertini`.

## Build Commands

### C++ Core Library (CMake + Ninja)

```bash
# Configure (from repo root)
cmake -DENABLE_UNIT_TESTING=ON -G Ninja -B build -S .

# Build
cmake --build build --target all --config Release

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

## Architecture

The project has three layers, built in order:

1. **`core/`** -- C++ shared library (`libbertini2`). Header-only-heavy design under `core/include/bertini2/`. Key subsystems:
   - `function_tree/` -- Expression tree (nodes, operators, symbols) for representing polynomial systems
   - `system/` -- Polynomial system construction, start systems (`start/total_degree.hpp`, `start/mhom.hpp`), patches, slices
   - `trackers/` -- Path tracking (fixed-precision and adaptive-precision trackers, predictors, Newton correctors)
   - `endgames/` -- Power series and Cauchy endgames for singular endpoint handling
   - `nag_algorithms/` -- Higher-level algorithms (zero-dim solve, numerical irreducible decomposition)
   - `io/parsing/` -- Boost.Spirit Qi parsers for classic Bertini input format
   - `blackbox/` -- CLI executable entry point (`bertini2_exe`)

2. **`python_bindings/`** -- Boost.Python + eigenpy bindings producing `_pybertini` native module. Each `*_export.cpp` wraps the corresponding C++ subsystem. Depends on `eigenpy` for NumPy/Eigen interop.

3. **`python/bertini/`** -- Pure Python package that wraps `_pybertini` into a user-friendly API. Submodules mirror the C++ structure: `function_tree`, `system`, `tracking`, `endgame`, `parse`, `nag_algorithm`, `multiprec`, etc.

## Key Dependencies

- **GMP/MPFR/MPC** -- Arbitrary-precision arithmetic (found via custom CMake modules in `cmake/`)
- **Eigen 3.3** -- Linear algebra (pinned to v3.3)
- **Boost** (serialization, filesystem, log, graph, regex, timer, chrono, thread, unit_test_framework, python) -- Boost >= 1.82 required; `boost_system` is conditionally linked for Boost < 1.89
- **eigenpy** -- Eigen/NumPy bridge for Python bindings
- **jrl-cmakemodules** -- CMake helper macros (auto-fetched via FetchContent if not found)

## Build System Notes

- The root `CMakeLists.txt` uses `jrl-cmakemodules` (fetched automatically). It currently only adds `core/` as a subdirectory; `python_bindings/` and `python/` subdirectory calls are commented out (the wheel build via scikit-build-core handles them).
- `pyproject.toml` configures scikit-build-core: wheel packages from `python/bertini/`, build dir is `bld/`.
- Cross-platform: Linux uses manylinux Docker + `auditwheel`; macOS uses Homebrew; Windows uses conda + clang-cl (MSVC has template compilation issues).
- `-Werror` is disabled globally. `-pedantic` is stripped from flags.

## CI/CD

- `.github/workflows/build_and_test.yml` -- Builds wheels on Ubuntu/macOS/Windows and runs tests. Triggered by pull requests and pushes to `develop`/`main`.
- `.github/workflows/publish.yml` -- Publishes to TestPyPI on `develop` push, PyPI on version tags (`v*.*.*`) with Sigstore signing and GitHub Releases.

### Linux wheel test coverage

Linux wheels are built inside a `manylinux_2_34` container (AlmaLinux 9, MPFR 4.1; set via `CIBW_MANYLINUX_X86_64_IMAGE`). The **full pytest suite runs on all three platforms** — on Linux it runs *inside* that container via `CIBW_TEST_COMMAND_LINUX`, and on macOS/Windows via the host-runner test jobs.

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

## Conventions

- C++ standard: C++17. Headers use `.hpp` extension.
- License: GPL v3 with additional terms (see `licenses/`, `core/ADDITIONAL_GPL_TERMS`).
- Version is tracked in `pyproject.toml` (the `version = "..."` line), read at runtime via `importlib.metadata.version("bertini2")`.
