# ADR-0067: The endgames differ in root estimation only; security is asked before acceptance

**Status:** Accepted
**Date:** 2026-09-29

## Context

The security ceiling is one rule: at security level 0, two consecutive endpoint approximations
whose dehomogenized infinity norm is above `Security.max_norm` truncate the path with
`SecurityMaxNormReached`.  It is a statement about where the path is going, and it is made
before anyone asks whether the path has converged.

The power series endgame applies it that way: after every extrapolation, before its loop
re-tests convergence.  The Cauchy endgame did not.  Its loop tested acceptance first and
returned `Success`, and consulted the ceiling only for rounds that had not been accepted.  So a
path that converged to a finite point above the ceiling had two answers:

| `x^2 = 1e8, y = 2x`, default `max_norm` 1e4 | finite solutions | codes |
|---|---|---|
| power series | 0 of 2 | `SecurityMaxNormReached` |
| Cauchy | 2 of 2 | `Success` |

Since 4.0 the power series endgame is the default (ADR-0058), so the difference became what a
user sees when they switch endgames: solutions appear or vanish.  It was found by the cellular
decomposition port, whose contract test for the accuracy estimates solves exactly such a system.

Moving the security block above the acceptance test was not enough.  Cauchy's count of
consecutive approximations runs over rounds in the operating zone only, because outside it the
Cauchy mean is not an approximation of anything (truncating on it lost genuine cyclic-6
solutions, 156 -> 155).  Acceptance, though, needs a single in-zone round.  A path whose first
in-zone round was already above the ceiling had a count of one and was accepted.

## Decision

The two endgames differ in how they estimate the root, and in nothing else.  For security:

1. **Security is asked before acceptance**, in both Cauchy loops (`RunImpl` and the adaptive
   one), as in the power series endgame.
2. **An approximation above the ceiling is never accepted** at security level 0.  It is the
   first of the two that truncate, or the second.
3. Cauchy's count still runs over **rounds in the operating zone** (or below
   `cycle_cutoff_time`), and still watches the larger of the endpoint norm and the loop floor.
   Those are properties of how Cauchy estimates: an out-of-zone mean is not an endpoint
   approximation, so it neither counts toward truncation nor can be accepted.

The user-visible invariant: **at security level 0, no endgame returns `Success` for an endpoint
above `max_norm`.**  At level 1 neither endgame truncates.

## Consequences

- A solve whose solutions lie above `max_norm` needs the ceiling raised or security level 1,
  under either endgame.  That was already true of the default endgame; it is now true of both.
- Do not move Cauchy's acceptance test back above the security block, and do not drop the
  `!above_ceiling` clause from it: either change brings back `Success` above the ceiling.
- Do not remove the operating-zone condition from Cauchy's count in the name of making the two
  endgames textually alike.  It is what keeps cyclic-6 at 156 solutions.
- A new endgame applies the same rule in the same place.
- Pinned by `security_ceiling_is_the_same_rule_in_both_endgames`
  (`core/test/nag_algorithms/zero_dim.cpp`), over both trackers and both endgames, and by
  `python/test/zero_dim/security_ceiling_test.py`.  The C++ test lowers the ceiling instead of
  raising the roots, because the double precision tracker is not reliable on a badly scaled
  system and the test is not about that.
