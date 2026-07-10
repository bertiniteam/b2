# ADR-0053: Wheel builds run concurrently with the C++ tests, not behind them

**Status:** Accepted
**Date:** 2026-07-10

## Context

`build_and_test.yml` builds the release wheels (`build_macos_ubuntu_wheels`,
`build_windows_wheels`) and runs the C++ unit-test suites (`test_cpp_unix`,
`test_cpp_windows`). The wheel jobs originally declared `needs: [..., test_cpp_*]`, so a wheel
build could not start until its platform's C++ tests finished.

That gate was a **false dependency**. The two jobs share no artifacts and are entirely
independent full compiles:

- `test_cpp_*` builds Boost-base from source and runs `ctest`; its compiler cache lives in the
  `ccache-cpp-*` namespace.
- the wheel jobs build via `cibuildwheel` against the **prebuilt** deps image/tarball (ADR-0049);
  their cache lives in the `ccache-wheels-*` namespace.

Nothing produced by the test job feeds the wheel job. The gate only serialized them — worst on
Windows, where ~22 min of tests ran *before* ~22 min of wheel build (~44 min) that could have run
alongside it (~22 min).

## Decision

Drop the C++-test jobs from the wheel jobs' `needs`. Keep only the genuine inputs plus the cheap
`doc_lint` fail-fast gate:

- `build_macos_ubuntu_wheels: needs [doc_lint, set_matrix, set_versions]`
- `build_windows_wheels: needs [doc_lint]` (it hardcodes its own Python matrix and reads no
  `set_matrix`/`set_versions` outputs)

The wheels now build **concurrently** with the C++ tests, gated only by the ~1-minute doc-lint.

Correctness of the run as a whole is unchanged: `test_cpp_*` are still jobs in the workflow, so a
C++-test failure still fails the run at the workflow level — and therefore still blocks
`publish.yml` (which invokes `build_and_test.yml` via `workflow_call`) from reaching PyPI. The
only thing lost is the runner-minutes saved when tests fail *and* wheels would have been skipped;
we trade those for roughly-halved wall-clock on the common (green) path.

The docs workflow (`build_docs.yml`) got the same treatment: the Doxygen C++ build, which needs
no wheel, was split into a parallel `cpp_docs` job.

The publish path is deliberately **not** loosened: `publish.yml` still waits on the entire matrix
before deploying — a real release should verify the whole build first.

## Consequences

- **Do not re-add `test_cpp_unix` / `test_cpp_windows` to the wheel jobs' `needs`.** It reads as
  "don't build wheels if tests fail," but it only re-serializes two independent compiles; the
  workflow-level result already enforces test success before publish.
- Wheels are now built even on a run whose C++ unit tests fail. That is intended: the wheel's own
  in-container / host pytest still runs and still catches build breakage, and total wall-clock is
  the priority on the green path.
- The `doc_lint` gate is retained in front of every heavy job (including the wheels) as the one
  cheap, genuinely-worth-it fail-fast — a doc typo still fails the run in ~1 minute without
  burning the matrix.
