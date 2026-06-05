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

- `.github/workflows/build-and-publish-to-pypi.yml` -- Builds wheels on Ubuntu/macOS/Windows, publishes to TestPyPI on `develop` push, PyPI on version tags (`v*.*.*`).
- Pushes to `develop` trigger TestPyPI publish; tagged releases go to PyPI with Sigstore signing and GitHub Releases.

## Conventions

- C++ standard: C++17. Headers use `.hpp` extension.
- License: GPL v3 with additional terms (see `licenses/`, `core/ADDITIONAL_GPL_TERMS`).
- Version is tracked in `python/bertini/_version.py` and `pyproject.toml`.
