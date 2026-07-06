# ADR-0017: Solve flakiness is a tolerance problem, not a seeding problem

**Status:** Accepted
**Date:** 2026-06-14

## Context

With the default tracking tolerances (`newton_before_endgame` 1e-5,
`newton_during_endgame` 1e-6) and a fresh random `gamma` each run, a
double-precision total-degree solve occasionally **missed** a solution (a
borderline path fails and a root goes unfound) or returned a **duplicate** (two
paths converge to the same solution). The finite-solution count jittered run to
run — cyclic-5 gave 69 / 70 / 71 / 72 — serial *and* distributed.

The tempting "fix" is to pin a random seed so runs reproduce. That is the **wrong
fix**: it merely selects a `gamma` that happens to work and hides the fragility
from future changes (e.g. performance work that perturbs the numerics). The
flakiness is real and should be removed, not masked.

## Decision

- **Fix flakiness with tolerances/settings, not seeding.** Tightening
  `newton_before_endgame` / `newton_during_endgame` (e.g. to 1e-7 / 1e-8) makes the
  solve reliably find every solution across random `gamma` draws. The right knob is
  the tracker tolerance (and, where relevant, precision), not the RNG.
- **Seeding is for reproducibility only** — identical runs for timing comparisons
  or debugging — never as the remedy for a flaky solve.
- **Measure correctness by counting *distinct* solutions**, the deterministic and
  mathematically meaningful quantity, not the raw path tally (which an occasional
  duplicate makes jitter even when the solve is correct). Two endpoints are the
  same per the same-point rule in ADR-0015.

## Consequences

- When a solve is flaky, reach for tolerances/precision first. If a fixed seed is
  the *only* thing that makes it pass, the real bug is still there — do not ship
  that as the fix.
- Examples and benchmarks must stay correct under a random `gamma`; a tutorial that
  only works at one seed is a red flag.
- Relevant for the upcoming performance work: changing step-size/precision logic
  must not be "validated" by a lucky seed — re-check correctness across several
  random draws (distinct-solution count against the known value).
- ADR-0029 later improved the *conditioning* of `gamma` (unit modulus, generated at
  maximum precision). That is complementary numerical hygiene, not a contradiction:
  a better-conditioned `gamma` is still no substitute for fixing a flaky solve with
  tolerances/precision, and correctness is still measured by the distinct-solution
  count across random draws.
