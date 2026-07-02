# ADR-0046: The solver records seam — solve() is ensure-answered

**Status:** Accepted
**Date:** 2026-07-03

## Context

The structured-output-directory arc's rung 4: solvers should record every path as it
completes and consult the records before computing, so that a fresh run, a crashed-run
resume, and a completed no-op rerun are the same code path (resume is memoization, not
a mode), for both faces of bertini2 (CLI and Python) and all three topologies (serial,
threaded, MPI).  Rungs 1–3 supplied the pieces: config digests (ADR-0043) make the ask
identity computable; seed-rooted randomness (ADR-0044) makes the rebuilt homotopy
identical so recorded endpoints are valid on rerun; the OutputDirectory (ADR-0045) is
the storage.

## Decision

The seam lives in `HomotopySolver` (the shared continuation engine, ADR-0040), keyed to
one architectural observation: **`StoreFullPathResult` is the single manager-side
installation point for completed paths in every topology** — the MPI manager loop, the
thread pool's main-thread install (the pack/store split), and re-track rounds all pass
through it, and it runs only on the main/manager thread.  Therefore:

- **A track record IS a serialized `FullPathResult`**
  (`records/solver_recording.hpp`): boundary data + endpoint + all endgame metadata,
  scalars exact and human-readable (doubles as `%.17g`, bit-exact round-trip;
  multiprecision as full-precision decimal + precision).  Emission is a tail call in
  `StoreFullPathResult` (plus explicit emits on the direct-install serial path and the
  re-track lambda) — single-writer by construction, no locks.
- **Hydration replays records through the same installer**: recorded paths are decoded
  to `FullPathResult`s and passed to `StoreFullPathResult`, so hydrated state is
  identical to computed state BY CONSTRUCTION — boundary data included, so the midpath
  crossing check works on a resumed run.  `Solve()`/`RunParallel()` then dispatch only
  the missing indices.  A re-track appends a fresh record; last-per-index wins on read.
- **The ask** (`RecordsAsk`): op + tracker/endgame kind (stable `kRecordName` strings
  on TrackerTraits/AlgoTraits — never typeid) + target `ContentDigest` + a
  `SettingsDigest` over (ZeroDimConf, Tolerances, AutoRetrack, PostProcessing — fixed
  documented order; extend by appending) + the global seed.  The run id is a hash
  prefix of the ask.
- **Attachment**: explicit `RecordTo(directory)`, or ambient via the
  `BERTINI_RECORDS_DIR` environment variable — attached on the manager rank only;
  workers never touch the records.  When neither is present, nothing records: the
  test suites stay clean, and surface-level default-on is rung 5's product decision.

## Consequences

- **The CLI, killed and rerun, finishes — with zero new flags** (demonstrated: SIGKILL
  mid-run, rerun completes, one run header, every path recorded once).  Under
  `mpirun`, records are written by the manager alone (one history file) and an MPI
  rerun fully hydrates — identical `main_data`, zero recomputation.  Note: OpenMPI
  does not forward arbitrary environment variables to ranks; use
  `mpirun -x BERTINI_RECORDS_DIR ...` (only the manager needs it, but forwarding is
  the reliable habit).
- Tested (test_nag_algorithms/zero_dim_records): record-then-hydrate equivalence to
  1e-14; PARTIAL directory resumes computing only the missing paths and matches the
  uninterrupted solve; a different seed is a different ask (two run headers, one
  directory); ambient env attachment.
- Hydration trusts the records (coordinates are a cache of the recorded computation);
  classification (PostEGAction) reruns over hydrated metadata identically.
- The `hydrating_` flag suppresses re-emission during replay; emission requires an
  attached directory and a run id, so mis-ordered calls are inert rather than wrong.
- Known scope edges (deliberate): `SharpeningConfig`/`MidPathConfig` are not yet in
  the settings digest (appending them later changes future asks only); the
  solutions-facing results handle and auto-declared `result` records are rung 5.
