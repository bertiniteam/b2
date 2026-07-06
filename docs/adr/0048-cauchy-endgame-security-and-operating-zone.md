# ADR-0048: The Cauchy endgame's divergence handling — watch the endpoint, gate on the operating zone

**Status:** Accepted
**Date:** 2026-07-06

## Context

The Cauchy endgame extrapolates a path's endpoint from a circle of tracked samples
around the target time (the Cauchy integral / roots-of-unity mean).  Two mechanisms
decide a path's fate near infinity, and getting the interaction right is subtle enough
that it broke a textbook system:

- **Acceptance**: when do we *trust* a converged approximation as a root?
- **Security truncation**: when do we cut a path off as diverging (`SecurityMaxNormReached`)?

Three prior decisions collided:

1. The **junk-success bug** (found 2026-07-03 via the structured output directory's own
   `function_residual`; present in Bertini 1): on an unpatched affine homotopy toward a
   *deficient* target, a path is a Laurent pole `x(t) = a/t + b`.  The Cauchy mean
   annihilates the pole term *exactly* (roots-of-unity identity), so consecutive
   approximations agree at the finite constant `b` — and the endgame minted **Success at
   a non-root** (residual ~1).  Fixed in PR #70 with two changes: a pole-component
   *operating-zone* acceptance gate (refuse convergence while the loop's negative-mode
   "pole mass" is significant), **and** a pole-growth *truncation* (truncate after N
   rounds of geometrically growing pole mass), **and** switching the security max-norm
   check from the extrapolated endpoint to the loop **samples**.

2. That combination **regressed cyclic-6** (a textbook system; Bertini 1 and pre-#70 b2
   both find its full **156** nonsingular solutions).  Post-#70 b2 lost 3–6 of them,
   seed-dependently.  The records' per-path store made the loss auditable against B1's
   ground truth; an endgame bisect pinned the causes.

## Decision

Three rules, replacing #70's samples-watch + pole-growth-truncation:

- **The security check watches the extrapolated ENDPOINT, not the loop samples.**  A
  genuinely-infinite endpoint's approximation grows to infinity *faster* than any finite
  loop sample along the way (the endpoint *is* the infinity; the samples are all finite),
  so the endpoint is the earliest, strongest divergence signal.  Watching the samples
  over-truncated finite paths whose samples transiently spiked.

- **Security truncation fires ONLY in the operating zone.**  Outside the zone the loop
  may encircle a pole or a branch point of the cover, and the Cauchy mean is *garbage* —
  its norm means nothing.  Truncating on that garbage killed clean convergent paths whose
  endgame-boundary loop transiently enclosed a branch point (mass grows a round or two,
  then dies as the radius shrinks past it) — the cyclic-6 residual.  Two consecutive
  *in-zone* rounds above `max_norm` truncate honestly; an out-of-zone round restarts the
  count.

- **No active pole-growth truncation.**  The operating-zone *acceptance gate* alone cures
  junk-success: a Laurent pole's mass never vanishes, so the gate never accepts it — the
  path instead runs to the endgame's natural terminal condition and returns non-Success.
  The pole-growth truncation was both unnecessary (the gate suffices) and the dominant
  cause of the cyclic-6 loss, so it — and its `num_pole_growth_rounds_before_truncation`
  config field — are removed (`b2cfgenc/2 → b2cfgenc/3`).

## Consequences

- **cyclic-6 recovers its full 156 nonsingular solutions** on every tested seed, matching
  Bertini 1; junk-success stays cured (the named endgame regressions pass with the
  pole-growth truncation gone).
- **At security level 0, an at-infinity path may EITHER truncate (`diverged`) OR converge
  to an infinite endpoint (`success`, classified infinite)** — both are valid, and which
  happens depends on the endgame's exact convergence.  Records and tests assert the robust
  invariant: every path is an answer (never `failed`), `success + diverged` covers all,
  and `CoarsePathStatus` maps the SuccessCodes to the vocabulary directly.  Level 1
  computes the infinite endpoints without truncating.
- This is the records arc paying a debt forward: the durable per-path store is what
  *caught* the regression (`function_residual`, coordinate audit) and let Bertini 1
  adjudicate ground truth.  See the arc doc (`arcs/structured-output-directory.md`, round
  10) and PR #70 for the prior design this refines.
