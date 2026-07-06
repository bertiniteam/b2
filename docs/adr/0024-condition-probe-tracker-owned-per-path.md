# ADR-0024: The condition-number random direction is refreshed once per path, by the tracker

**Status:** Accepted
**Date:** 2026-06-17

## Context

The adaptive-precision machinery estimates a path's condition number by solving `J x = r` for a
random probe direction `r` and forming `||J|| * ||J^{-1}||`.  That probe, `rand_temp_`, was drawn in
`NewtonCorrector::ChangeSystem` -- i.e. **once per `SetSystem`** -- via `RandomOfUnits`, which reads
this thread's RNG engine (`bertini::ThreadEngine()`).

That draw was at an **uncontrolled RNG point**.  The zero-dim algorithm reseeds the thread RNG per
path (`ReseedThisThread(idx)`) so tracking is reproducible, but the probe was drawn before any such
reseed, at whatever engine state `SetSystem` happened to leave.  That state differs between a serial
solve and a distributed worker (which has just called `SetGlobalSeed` of the broadcast seed), so the
two got **different probes**.  For well-conditioned paths the probe doesn't change the outcome; for
ill-conditioned/near-singular paths it changes adaptive-precision decisions, so serial and
distributed diverged there.  This surfaced while making distributed solves reproducible (ADR-0023,
authoritative pi (#156), broadcasting rank 0's systems, rank 0 shipping start points).

The probe's correct lifecycle is **per path**: one direction, held fixed for the *entire* track of a
point -- all Newton steps and the endgame's sample-circle sub-tracks -- which also keeps condition
estimates comparable across steps.  A first attempt that refreshed the probe inside
`Tracker::TrackPath` was **wrong and reverted**: the endgame issues many `TrackPath` calls per point,
so it redrew mid-track, churned the RNG, and regressed bit-identity.

## Decision

Refresh the probe **once per path, triggered by the algorithm, through the tracker** -- not in
`ChangeSystem`, not in `TrackPath`, not per Newton step:

- `NewtonCorrector::RefreshRandomDirection()` (re)draws the probe; `ChangeSystem` calls it once for
  initialization.
- `Tracker::RefreshConditionDirection()` forwards to the corrector.
- `ZeroDim::ExecuteBeforeEG` calls it once, right after the per-path `ReseedThisThread(idx)` and
  `SetThreadPrecision`, before tracking -- and **not** again for the endgame.  So the direction is
  drawn from an RNG state determined solely by the path index, identical whether the path runs in a
  serial loop or on a distributed worker, and held for the whole track (pre-endgame + endgame).

## Consequences

- The probe is now a reproducible, per-path quantity.  Bit-identity that already held on
  well-conditioned finite solutions is preserved (cyclic-5 finite hash unchanged), and more
  ill-conditioned cases now agree serial vs distributed.  All suites pass (C++, 53 serial pytest, 6
  MPI under `mpirun -n 3`).
- **The probe was not the whole story.** Controlling it did *not* fully close the Griewank-Osborn
  seed-1 case, so there is at least **one more uncontrolled non-deterministic input** on deeply
  singular paths still to find (low priority: a measure-zero case whose finite answer is already
  correct in both modes).
- **Ownership is still split:** the corrector owns the probe storage and the tracker merely triggers
  the refresh.  The cleaner end state -- the random direction fully owned by the tracker/track and
  the corrector a **pure kernel** (`(J, f, probe) -> (step, condition estimate)`, no hidden RNG or
  scratch) -- is left to the **planned predict/correct rewrite**, whose interface is meant to express
  exactly a direction that persists over a whole track.  The seed-1 residual is best hunted there too.

See ADR-0023 (unified speculative-full-path parallelism), ADR-0016 (broadcast homotopy seed),
ADR-0017 (tolerances, not seeding).
