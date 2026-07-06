# ADR-0016: Distributed ZeroDim broadcasts the homotopy seed so all ranks agree

**Status:** Accepted
**Date:** 2026-06-14

## Context

`solve(communicator=...)` runs the solver constructor on **every** rank, and the
constructor draws random data: the patch, the total-degree start-system
coefficients, and the homotopy `gamma`. Each rank therefore built a **different**
homotopy, so "path index `i`" denoted a different path on each worker. The manager
collected one endpoint per index from whichever worker computed it and stitched
together a scrambled, mostly-wrong result — cyclic-5 returned ~17 distinct
solutions instead of 70. The solve completed and reported a plausible path count;
it was **silently wrong**, with no crash or error.

Per-path tracking randomness was already deterministic in the path index
(`ReseedThisThread`), so the *only* source of divergence was the per-rank homotopy
construction.

## Decision

At the start of `RunParallel` (before `PreSolveChecks`), broadcast the manager's
RNG seed to all ranks, reseed, and re-run the system-management policy's
`SystemSetup` so every rank forms the **identical** start system and homotopy
(and re-points the tracker at it):

```cpp
unsigned long seed = parallel::IsManager() ? GetGlobalSeed() : 0ul;
MPI_Bcast(&seed, 1, MPI_UNSIGNED_LONG, 0, comm);
SetGlobalSeed(seed);
SystemManagementPolicy::SystemSetup(this->template Get<ZeroDimConf>().path_variable_name);
num_start_points_ = StartSystem().NumStartPoints();
GetTracker().SetSystem(Homotopy());
```

Combined with the deterministic per-path RNG, this makes the distributed solve
byte-for-byte consistent with the serial solve.

## Consequences

- Any future per-rank randomness that affects the start system or homotopy must be
  **derived from the broadcast seed** (or itself broadcast). Drawing random data
  independently per rank reintroduces this bug.
- A distributed correctness regression manifests as a **wrong distinct-solution
  count**, not a crash. Keep a count-against-known-value correctness check in the
  MPI tests and in the distributed examples (see ADR-0017 on counting *distinct*
  solutions) — a path tally alone will not catch it.
- This is why the cyclic scaling example verifies the recovered distinct finite
  count equals the mathematically known cyclic-n value.
