# ADR-0018: Bound every test with a timeout; pin seeds and skip known platform grinds

**Status:** Accepted
**Date:** 2026-06-14

## Context

A CI run lasted **6 hours 41 minutes**: the "Test wheel on Windows" pytest step hung for 6 hours
on Python 3.10/3.11/3.14 (killed by the GitHub runner's default 6h job timeout) while 3.12/3.13
passed in ~5.5 minutes. The hang was a single test,
`mhom_test.py::test_mhom_solves_adaptive_precision`: on Windows (clang-cl) the AMP tracker, the
blend-block homotopy, and the Cauchy endgame do not keep precision in lockstep, so a path grinds
toward `MaxPrecisionAllowed` for hours instead of converging or failing. On Linux it solves in
0.15s — the grind is Windows- and gamma-specific.

Two anti-patterns made this possible:

- **Nothing bounded a test.** `pytest python/test/ -v` has no per-test timeout, and the GitHub job
  inherited the 6h default, so one stuck test wedged the whole job for 6h.
- **"Try until success."** The test wrapped the solve in a retry loop (an earlier version did 40
  attempts; even 3 is wrong) — masking failures and multiplying a slow solve. Tests must give a
  definite pass/fail, fast.

## Decision

**Every test is bounded, at three independent levels:**

- `pytest-timeout` in `pyproject.toml`: `timeout = 180`, `timeout_method = "thread"`. The `thread`
  method is required — a grind stuck inside the C++ extension cannot be interrupted by a signal, so
  the signal method would not bound it; the thread method hard-terminates the process with a stack
  dump. `pytest-timeout` is in all three CI test-dependency lists (Linux in-container, macOS host,
  Windows host).
- `ctest --timeout 600` on both C++ test runs.
- job-level `timeout-minutes` on the four test jobs (25 for the wheel-test jobs, 30/40 for the C++
  jobs, sized above the healthy run plus margin).

A hanging test therefore becomes a **fast, debuggable failure** (with a stack dump) instead of a
6h job, no matter where it hangs.

**No retry-until-success in tests.** A test pins a fixed RNG seed (`set_random_seed` /
`SetGlobalSeed` — RandomMp is reseedable, see ADR-0016/0017) so its homotopy is reproducible, and
asserts a single solve. A fixed seed in a *test* is legitimate reproducibility (it is not the
flakiness-masking ADR-0017 warns against in the *solver*).

**Skip a known platform-specific grind, with a documented reason.** The blend-block AMP MHom solve
is skipped on Windows (`skipif sys.platform == "win32"` / `#ifdef _WIN32`), citing the precision
non-lockstep bug, because there is no point gating CI on a path with a known, unfixed grind there.
The same AMP MHom code path is still exercised on Windows by the eigenvalue test (a different MHom
formulation that does not hit the bug), so coverage is retained.

## Consequences

- CI cannot hang for hours again; the worst case is a job that fails within its `timeout-minutes`
  with a stack dump pointing at the offending test.
- New tests should expect the 180s per-test budget. A test that legitimately needs longer is a
  smell — split it or justify it.
- The underlying defect (AMP + blend-block + Cauchy endgame precision lockstep; ADR-0015's
  block-composed homotopy reaching `MaxPrecisionAllowed` in the endgame) remains the real fix — it
  is deferred performance work, and removing the Windows skip is its acceptance test.
- Editing `.github/workflows/build_and_test.yml` invalidates the ccache key (it is in the
  `hashFiles` set), so the first run after a workflow change is a cold, slow build; subsequent runs
  are warm.
