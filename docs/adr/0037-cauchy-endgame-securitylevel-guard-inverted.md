# ADR-0037: The Cauchy endgame's SecurityLevel guard was inverted (diverging paths ground at escalating precision)

**Status:** Accepted

## Context

Benchmarking the b2 CLI against Bertini 1.7 on zero-dim solves (cyclic5, adaptive precision)
showed b2 ~100× slower. An instrumented investigation found:

- The over-escalation to multiprecision is **entirely in the endgame** (≈89% of corrector
  criterion checks, ≈97% of the multiprecision arithmetic cost); main-path tracking runs in double
  and already matches b1.
- The escalation is driven solely by AMP **Criterion C**, responding to a Jacobian that is
  **genuinely near-singular**: an independent full-inverse norm matched the random-probe estimate of
  `‖J⁻¹‖` to within ~6%, reaching `~10²⁰`. (So the condition probe is faithful; the matrix really is
  near-singular.)
- The near-singular points belong to paths whose **homogenizing coordinate `x[0] → 0`** — i.e.
  paths **diverging to infinity**. As `x[0]` collapses, the homogenized Jacobian degenerates, so
  Criterion C correctly demands ever more precision — but it is *wasted* work on doomed paths.

Bertini 1 truncates such paths cheaply (it reported "Truncated infinite paths: 50 — set
SecurityLevel to 1"). b2 was *not* truncating them. Root cause: in `endgames/cauchy.hpp` the
SecurityMaxNorm divergence check was guarded by `if (SecuritySettings().level)` — i.e. it only ran
when `level != 0`. But the **default is `level = 0`**, and the PowerSeries endgame guards the same
check with `if (level <= 0)`. An in-code NOTE already flagged the two as "appear inconsistent."
Bertini 1 semantics (authoritative): **SecurityLevel ≤ 0 ⇒ truncate paths going to infinity**
(the safe default); **level ≥ 1 ⇒ keep tracking them** so the at-infinity endpoints are computed.
So b2's Cauchy guard was inverted: at the default level it never truncated, and ground every
diverging path through the full Cauchy endgame at escalating precision.

## Decision

1. Fix the guard in `endgames/cauchy.hpp` to `if (this->SecuritySettings().level <= 0)`, matching
   the PowerSeries endgame and Bertini 1. At the default SecurityLevel 0, diverging paths are now
   truncated (returning `SuccessCode::SecurityMaxNormReached`).
2. Keep `SecurityMaxNormReached` classified as a **path failure**, not an infinite endpoint — a
   truncation is "we gave up on this path", not a computed at-infinity solution (this mirrors
   Bertini 1, which lists them as truncated and tells you to raise SecurityLevel to compute them).
   `ZeroDim::Report()` therefore continues to count only `GoingToInfinity` as a divergence.
3. To *compute* infinite endpoints, set `SecurityConfig::level >= 1` so paths are tracked to their
   (infinite) endpoint and classified `GoingToInfinity`. The `infinite_solutions_at_infinity` test
   now sets `level = 1` to exercise that path.

## Consequences

- cyclic5 (serial, mptype 2): **25.3s → 11.1s** (~2.3×), with all 70 finite solutions unchanged
  and matching Bertini 1 to <1e-8. The remaining gap vs b1 (~0.2s) is the per-step double-precision
  machinery overhead, tracked separately.
- Default behavior is now b1-faithful: diverging paths are truncated cheaply instead of ground
  through the MP endgame. Users who want at-infinity endpoints raise SecurityLevel (as in b1).
- Note left for follow-up: `SecurityConfig::max_norm`'s default (`1e4`) carries a "wrong default
  value" comment; and `InfiniteSolutions()` still lists `SecurityMaxNormReached` points — both are
  pre-existing and out of scope for this fix.
