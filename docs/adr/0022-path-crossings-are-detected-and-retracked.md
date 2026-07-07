# ADR-0022: Path crossings at the endgame boundary are detected and re-tracked, not tolerated

**Status:** Accepted
**Date:** 2026-06-16

## Context

The zero-dim solve has always included a *midpath check* (`MidpathChecker`,
`EGBoundaryAction`): at the endgame boundary it compares every pair of path
endpoints and, if two land on the same point (within a relaxed tolerance), it is
supposed to re-track the offending paths with tightened settings. Two distinct
paths reaching the *same* point at the boundary is a probability-0 event, so when
it happens it is a **path crossing** — a symptom of under-resolved tracking (most
often a too-low-order predictor) — not a benign coincidence.

The machinery did not work, for two reasons:

- **Inverted `same_start`.** `MidpathChecker::Check` computed
  `same_start = (||start_i - start_j|| > tol)`, i.e. true when the start points
  *differ*. Since `CrossedPath` sets `rerun = !same_start`, a genuine crossing
  (distinct starts, coincident endpoints) was flagged `rerun == false` and **never
  re-tracked**. The comparison must be `<`.
- **`Check` never reset its state.** `passed_` was only ever set `false`, and
  `crossed_paths_` was never cleared. `EGBoundaryAction` calls `Check` in a loop
  (once before any resolve, again after each `MidpathResolve`); on the second call
  it returned a stale `false` and re-tracked stale indices, so the loop always
  burned all attempts and never observed a clean pass.

Separately, ADR-0017 ("flakiness is tolerances, not seeding") had framed boundary
duplicates as a tolerance-flakiness nuisance to *measure around* by counting
distinct points (dividing by multiplicity). That framing — and the "harmless
duplicate" language in the `solving_at_scale` tutorial — is misleading: a boundary
duplicate is a crossing the algorithm should *fix*, not merely tolerate.

## Decision

- **Fix the two bugs.** `same_start` uses `<`; `Check` resets `passed_ = true` and
  clears `crossed_paths_` at entry.
- **Make the resolve a bounded, escalating action — never an open-ended loop.**
  `EGBoundaryAction` re-tracks crossed paths at most
  `ZeroDimConfig::max_num_crossed_path_resolve_attempts` times (default 2), so it
  always terminates. `max_num_crossed_path_resolve_attempts == 0` means
  *detect-and-report only* (the conservative, Bertini-1-style choice).
- **Escalate the remedy, don't just shrink tolerance.** Each attempt tightens the
  tracking tolerance **and** raises the ODE predictor to the default (RKF45) if a
  low-order predictor is in use — a too-low-order predictor (notably Euler) is the
  usual root cause, and tightening tolerance on a first-order predictor is far less
  effective than moving to a higher order.
- **Unresolved crossings stay observable, and we keep tracking them.** After the
  attempts are exhausted, the still-crossed paths still go through the endgame (we
  never drop more answers than Bertini 1 would), but the outcome is recorded in a
  `MidpathCheckReport` (`passed`, `num_crossings_detected`, `num_resolve_attempts`,
  `crossed_path_indices`), retrievable via `ZeroDim::EndgameBoundaryMetadata()`
  (Python: `endgame_boundary_metadata`), and a warning is emitted. We do **not**
  add a `SuccessCode::MidpathCrossing`: those paths tracked fine, they just
  collided, so the crossing fact belongs in the report, not in per-path
  `success_code`.
- **Refine ADR-0017, don't contradict it.** Tolerances/precision are still the
  right fix; this work simply *automates* that fix inside `EGBoundaryAction` and
  stops calling duplicates harmless. The distinct-solution count remains a valid
  *verification*, but it is a check, not an excuse. "Harmless duplicate" language is
  removed from the tutorials, examples, and comments.

## Consequences

- The accessor `EndgameBoundaryData()` is renamed `EndgameBoundarySolutions()`
  (Python `endgame_boundary_data` → `endgame_boundary_solutions`); a new
  `EndgameBoundaryMetadata()` exposes the report. This is a breaking rename, but the
  old name had only internal callers.
- A solve that hits a crossing now self-corrects (with the default predictor it
  generally does not hit one at all); a caller can detect an unresolved crossing
  programmatically rather than discovering a wrong count later.
- See ADR-0017 (tolerances vs. seeding) and ADR-0015 (same-point classification),
  which this ADR refines and complements.

## Future work (deferred — intentionally not built here)

The per-attempt remedy wants to become a **configurable, ordered remedy list**
(tighten tolerance, more Newton iterations, bump predictor, shrink min step size,
…), applied in sequence until exhausted, so a user can compose a resolution
strategy. For now the escalation is a fixed, minimal two-remedy step; getting the
strategy abstraction right is out of scope. Termination is guaranteed regardless,
because the attempt count is bounded.
