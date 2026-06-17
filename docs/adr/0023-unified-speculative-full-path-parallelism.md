# ADR-0023: ZeroDim uses one per-path primitive and a speculative-full-path parallel model

**Status:** Accepted
**Date:** 2026-06-17

## Context

The distributed ZeroDim solve was a **two-phase, barrier-separated** manager-worker model:

1. Phase 1 — workers track every path start→boundary; each boundary result is serialized back to
   the manager.
2. A barrier, then the manager alone runs the midpath (path-crossing) check and any re-tracks while
   the workers sit idle.
3. Phase 2 — the manager serializes each boundary point back out to a worker, which runs the endgame.

Three problems came out of the crossed-paths work (ADR-0022):

- **A precision-transfer bug by construction of the handoff.** The Phase-2 task serialized the
  boundary point and stepsize but dropped the boundary *precision*. An adaptive-precision worker
  therefore resumed the endgame at double precision, and a distributed AMP solve diverged from the
  serial one — cyclic-5 returned a precision-degraded subset instead of all 70 finite solutions. The
  bug existed *only because* boundary state crossed the worker→manager→worker handoff.
- **A mid-solve barrier + idle-manager window.** The slowest pre-endgame path stalled every worker
  before any endgame began, and all workers idled during the manager-only midpath check / re-track
  (flagged in a code comment as a known bottleneck).
- **Two separate per-path code paths** (a serial flow and a distributed worker flow) that set up
  precision/tolerance/stepsize differently and could drift.

## Decision

**One per-path primitive.** The four near-duplicate per-path routines collapse onto two shared
bodies, `ExecuteBeforeEG` / `ExecuteDuringEG`, driven through a reference-context
(`BeforeEGContext` / `DuringEGContext`). The serial flow passes its member tracker/endgame/observers;
a worker thread passes its cloned `PathThreadState`. They standardize on `SetThreadPrecision`
(thread-local), correct on the main thread and on worker threads alike. `ExecuteOnePath` composes
the two into a whole path (set pre-endgame tolerance → track to boundary → set endgame tolerance →
run endgame).

**Speculative-full-path model.** A worker (or the serial loop) executes one *whole* path as the unit
of work. Because the boundary state never leaves the executing context, the precision-transfer bug
*cannot recur* — it is designed out, not patched.

**Crossing resolution is bounded and parallel.** After a round of whole-path execution, the manager
runs the midpath check on the collected boundary points; crossed paths are re-dispatched as ordinary
whole-path tasks to the same worker pool (re-tracking is path-independent — no idle-manager
bottleneck), with escalated settings, bounded by `max_num_crossed_path_resolve_attempts`
(`0` = detect-and-report only). A one-int `MPI_Bcast` between rounds tells workers whether to run
again and to apply the same `EscalateRetrackSettings` the serial flow uses. Serial and distributed
share the per-path primitive and the same resolution logic.

**One MPI message type.** `PathBeforeEGResult` / `Phase2Task` / `PathDuringEGResult` collapse to a
single `FullPathResult` (boundary data + final solution + endgame metadata); the task is just the
path index.

## Consequences

- A distributed solve now produces the same correct answer as serial on generic problems: cyclic-5
  recovers all 70 distinct finite solutions in **both** double and adaptive precision, serial and
  distributed. Guarded by new acceptance tests in `python/test/parallel/test_mpi_zerodim.py`.
- In the speculative model a crossed path's (now-wasted) endgame runs before the crossing is
  detected, so a re-track redoes the *whole* path including a second endgame. Crossings are
  probability-0 events, so the expected extra cost is ~0 and the barrier is removed on *every* solve;
  a deliberately crossing-heavy config (e.g. the `crossed_paths` tutorial) pays more, bounded by the
  attempt count.
- `ReinitializeInitialStepSize(true)` is re-asserted at the start of each pre-endgame track: a path's
  endgame now runs immediately before the next path's pre-endgame on the same tracker and leaves a
  tiny step, which would otherwise make the next path crawl from the start time.

## Known limitation (separate from this work)

Serial and distributed are **not** bit-for-bit identical on pathological systems (e.g. Griewank-Osborn
at a tight endgame boundary): per-path `max_precision_used` and the classification of *infinite/singular*
paths can differ. The cause is that the distributed solve re-forms the homotopy from the broadcast
seed (ADR-0016) and that re-form is not bit-identical to the serially-constructed homotopy, so the
two solve slightly different (but equally valid) random homotopies. The *finite* solutions — the
gamma-independent answer — agree. Making the re-form reproduce construction bit-for-bit is a separate
reproducibility fix tracked against ADR-0016, not the parallelism model.

See ADR-0016 (broadcast homotopy seed), ADR-0017 (tolerances, not seeding), ADR-0022 (detect and
re-track crossings).
