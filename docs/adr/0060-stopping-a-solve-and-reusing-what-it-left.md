# 0060 — Stopping a solve, and reusing what it left behind

**Status:** Accepted
**Date:** 2026-09-19 (designed with the maintainer 2026-09-18; landed with 3.5.0, PR #449)

## Context

Until 3.5 a solve could not be stopped.  In a notebook, Ctrl-C was noted by CPython and
ignored until the solve finished, because the solve runs on the calling thread with the
interpreter lock released and `KeyboardInterrupt` is raised from an evaluation loop that is
not running.  There was also no way to give up on a path that was taking too long, and a
path that failed for any reason reported a time of zero and no point: it was simply missing.

Adding these touched three things that are easy to get wrong later, because each has an
obvious-looking alternative that is worse.

## Decision

### 1. One cooperative stop, checked between steps, never an exception across a worker

A process-wide lock-free flag (`bertini::RequestStop()` and friends) and a per-tracker
wall-clock deadline (`SetMaxWallClockTime`) are both checked at the top of every tracking
iteration, before the tracker's own budgets.  A path in flight returns
`ExternallyTerminated` or `WallClockLimitReached` from there; a path not yet begun is never
begun and stays `NeverStarted`.  Nothing is thrown, no thread is preempted, and the solve does
not move off the calling thread.

Why not a watcher thread or a solve on a worker thread: a fresh thread's multiprecision
working precision can be zero and aborts inside mpfr on first use; the session guard restores
thread-local random streams belonging to whichever thread runs it; and MPI is initialised
funneled.  Why not an exception: it would leave the per-path result slots half written, which
is exactly the state the user then inspects.  The cost of the cooperative check is that a
stop overruns by at most the step in flight.

The deadline is a point in time, not a duration, because it must hold across `TrackPath`
calls: an endgame issues hundreds of them for one path, and a per-call budget would start
over with each.

### 2. Every path says where it got to, and an abandoned path is recorded

Every path's metadata carries `final_time_used` (written after every stage), the step tally
(`num_successful_steps`, `num_failed_steps`, over the whole path), the precision reached, and,
for a path that did not succeed, `latest_path_point`.  No observer is attached to obtain this: the
tracker already holds its position when it returns, and it already counts its steps; the base
tracker gained a second pair of step counters that only the caller resets.  An abandoned path
IS recorded, with its code, its stamp, and the budget that stopped it.  A path never begun is
NOT recorded: its absence is explained once by the solve having been cut short.

The maintainer's words, which settle the first half: "an abandoned path with no record is
NOT more honest.  it does not record what happened, why that result is missing.  it's just
silently missing."

### 3. Budgets are out of the ask; what to reuse is the user's policy

The wall-clock budgets (`max_path_wall_clock_duration`, `max_solve_wall_clock_duration`)
are excluded from the configuration digest, beside `num_threads`.  A budget says how long to
wait for an answer, not what the answer is: a 60 s and a 61 s run are the same ask.  Encoding
them would have forced a config-encoding version bump, regenerated every fixture digest, and
invalidated every existing record, for nothing.

That leaves the question of whether a recorded abandonment may stand in for tracking the
path again.  The rule: **a path record is reusable when the run that would reuse it asks no
more of that path than the run that produced it.**  A completed path (success, divergence, a
failure under settings that are part of the ask) always qualifies.  An abandonment was ended
for a reason deliberately kept out of the ask, so whether to reuse it is the user's call, made
through `RecordsConfig::recall`, a `RecallPolicy`:

- `Nothing`: track every path fresh (the old `recall = false`).
- `Completed` (default): reuse completed paths; re-track abandoned ones -- except a path a
  per-path budget cut off, which is reused while the current budget is no larger than the
  recorded one, and re-tracked once it goes up, which is exactly when paying again is the
  point.  An interrupted path is always re-tracked: an interrupt has no budget to compare.
- `Everything`: reuse every recorded outcome as recorded.

The policy lives in its own `RecordsConfig`, not on `ZeroDimConfig`, because the question is
the same for every algorithm that records paths.  The decision is made inside the recall
routine, so every topology (serial, threaded, MPI manager) applies it.

The whole-solve budget is treated exactly as an interrupt: a path it cuts off reads
`ExternallyTerminated`, not `WallClockLimitReached`, because the path was abandoned for a
reason that has nothing to do with the path.  The tracker cannot tell the two deadlines apart,
so the solver rewrites the code when the solve deadline was the one that bound.

Step and precision budgets need none of this machinery: they are inside the ask already, so a
different budget is a different run by identity and the old abandonment is never consulted.

## Consequences

- Do not add a watcher thread, move the solve off the calling thread, or throw across a
  worker to make stopping "faster".  The one-step overrun is the accepted price.
- Do not make the wall-clock budgets part of the configuration digest, and do not add a
  machine identity to records to make the wall-clock comparison exact across machines (the
  maintainer declined: it likely needs user consent).  The recorded stamp -- steps, precision,
  time reached -- is the machine-independent account; the wall-clock comparison is exact on
  one machine and a heuristic across machines, and the docs say so.
- Do not record paths that were never started, and do not stop recording abandoned ones.
- Do not fold the recall policy back into a per-algorithm config, and do not hard-code a
  reuse rule for abandonments: it is the user's choice.  Keep the bool spelling working
  (`True` = `Completed`, `False` = `Nothing`).
- New `SuccessCode` values are appended at the END of the enum: the integers are part of the
  `b2rec/1` record contract.  `WallClockLimitReached` is currently last.
- A path whose endgame fails now has an EMPTY solution slot.  Before this the slot held the
  previous path's approximation, because the endgame writes its approximation only on
  convergence.  Do not "restore" the unconditional copy.
- Deferred, deliberately: a retry workflow (re-run only the failed or abandoned paths under
  different settings; a records-identity question, since changing a setting changes the ask),
  the deterministic work budget of #400, and whether any budget should ever come out of the
  ask on purpose.  Not applied to the MPI solve: a signal reaches one rank.
