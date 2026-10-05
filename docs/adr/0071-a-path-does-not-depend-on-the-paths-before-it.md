# ADR-0071: A path does not depend on the paths tracked before it

**Status:** Accepted
**Date:** 2026-10-05

## Context

A tracker and an endgame are reused path after path: by the zero-dim solver, on its own
objects in a serial solve and on one clone per thread in a threaded one, and by anyone driving
them by hand.  Each path's random draws were already deterministic, because the solver reseeds
the thread's generator from the global seed and the path index before each stage
(`ReseedThisThread`, ADR-0044).  Results were still irreproducible under threads (#378), because
state that is not random flowed from one path into the next:

- The adaptive tracker built a track's first step at `current_precision_`, which still held the
  precision the previous track ended in, and set the step's precision before the assignment
  that replaced it.
- The downward conversion to double re-precisioned only the end time.  The step size, the time
  and its increment kept the multiple precision the track came down from.
- The predictor and corrector rounded the condition-number probe in place at every change of
  precision, so a drop destroyed digits that a later rise could not restore.
- The endgame's c/k probe was drawn the first time an endgame needed it and kept for the
  object's life.  An escalation then moved it out of the double lane, so the next run in double
  drew a new one mid-path.
- An adaptive endgame cleared only the lane a run used, so a run in double left an earlier run's
  multiprecision times in place, and `LatestTime` read them.

With one thread, a path's predecessors are fixed by the path order, so a serial solve was
reproducible and still depended on that order.  With a pool, the scheduler picks each thread's
predecessors.

## Decision

A path's result is a function of the systems, the settings, the seed and the path's index, and
of nothing a previous path left behind.  The rule holds whatever runs the paths: one thread or
many, one process or MPI ranks.  It is also met by state hygiene in the tracker and endgame,
not by building fresh objects for each path.  A fresh object would hide the problem from the
solver and leave it in place for anyone reusing a tracker or an endgame by hand.

- Whatever a track or run sets up, it sets from its own inputs: the start precision, the start
  point, the stored configuration.  Never from the value a previous track left in a member.
- A random direction is stored as drawn.  Every working copy is derived from that copy at the
  precision it is needed in, never rounded in place.  The endgame's c/k probe is drawn in double,
  so its value is exact in every lane and at every precision.
- A random direction that should be fresh per path is drawn at the start of the run, after the
  per-path reseed, never lazily mid-path.
- Containers kept per lane are cleared in every lane when a run begins.

## Consequences

- **Do not store a value that the next path will read without resetting it.**  When adding a
  member to a tracker or an endgame, decide where each track or run gives it its value.
- **Do not round a stored random direction in place.**  Derive the working copy from the stored
  copy.
- **Do not draw random numbers lazily on first use.**  The draw then lands in whichever path
  happens to be first.
- **Do not replace this with a fresh tracker per path.**  That would pass the solver's tests and
  leave the trap for hand users.
- Results changed from 4.0.0.dev1 in the last bits and, occasionally, in step counts.  The old
  results depended on path order, so no earlier result is the reference.
- Pinned by tests at every level, each of which fails against the code without its fix:
  - `a_track_does_not_depend_on_the_track_before_it` and
    `the_step_size_is_at_the_working_precision_after_a_track` (`amp_tracker_test.cpp`);
  - `the_probe_direction_survives_a_drop_in_precision` (`newton_correct_test.cpp`);
  - `cauchy_run_does_not_depend_on_the_run_before_it` and
    `power_series_run_does_not_depend_on_the_run_before_it`, which run for every tracker type;
  - `threaded_solve_is_bit_identical_to_serial_*` (`threaded_solve.cpp`), which compares
    endpoints and every metadata field across 1, 2, 3 and 8 threads.
  Python counterparts are in `python/test/zero_dim/reproducibility_test.py`, and the
  distributed-equals-serial test is in `python/test/parallel/test_mpi_zerodim.py`.
- The guarantee is bit-identity for one build on one machine.  Different compilers, libraries
  or instruction sets may round differently.
