# ADR-0056: metadata_for matches against representatives, at final_tolerance

**Status:** Accepted
**Date:** 2026-07-13

## Context

`metadata_for(point)` (issue #302) looks a solution's metadata up by its point, so a user can take a
result the solver returned and ask for its diagnostics. The first implementation matched `point`
against the **complete** solution list (every endpoint, including the non-representative copies of a
multiplicity cluster) and used the solver's **same-point clustering tolerance**
(`final_tolerance × same_point_tolerance_multiplier`) as the default match window. It then threw
`"the point matches more than one distinct solution cluster; reduce tol"` if two representatives fell
inside the window.

This broke the most basic round-trip: **feed a point the solver just handed you back in, and it can
report ambiguity** — even though the solver had already clustered those endpoints as the *same*
solution. It bit on singular points, whose coincident copies scatter, and it got worse (not better)
at a *tighter* user tolerance, because the same-point window is exactly the one scale wide enough to
span from one cluster into a neighbor. A point the solver classified as not-distinct being called
ambiguous on lookup is self-contradictory.

## Decision

Two coupled changes to `MetadataForPoint` / `CoincidentMetadataForPoint`:

1. **Match against the multiplicity representatives only** (one candidate per distinct solution), not
   the complete set of copies. A `representatives_only=false` argument keeps the match-everything view
   for debugging, but it is not the default.
2. **Default the match tolerance to `final_tolerance`** (the accuracy each endpoint is computed to),
   *not* the same-point clustering tolerance. `DefaultPointMatchTolerance()` now returns
   `final_tolerance`; the clustering tolerance is exposed separately as `SamePointTolerance()`
   (`final_tolerance × same_point_tolerance_multiplier`), and `ComputeMultiplicities` and
   `CoincidentMetadataForPoint` use it.

The invariant that makes this correct: `ComputeMultiplicities` marks a point a distinct representative
only when it is at least `SamePointTolerance` from every earlier one, so **any two representatives are
≥ `SamePointTolerance` apart**. `final_tolerance` is that divided by `same_point_tolerance_multiplier`
(≥ 1), so a match window of `final_tolerance` around any point contains **at most one representative**.
A point taken from the solver's own solution lists (which return representatives by default,
`merge_multiplicities=true`) therefore matches its own representative and can *never* be reported
ambiguous. `MetadataForPoint` returns that representative (carrying `.multiplicity`);
`CoincidentMetadataForPoint` finds the representative, then gathers its whole cluster (everything
within `SamePointTolerance` of it) so the returned count is the multiplicity, independent of `tol`.

## Consequences

- Feeding any solver result straight back into `metadata_for` resolves to its representative — the
  contract users expect. The "more than one distinct cluster" throw is now only reachable with a
  deliberately coarse user-supplied `tol`, where it is the honest answer.
- **Do not restore the same-point tolerance as the `metadata_for` default, and do not match against
  the complete solution set by default.** Either one reintroduces the spurious ambiguity: they are the
  two halves of the same bug. The representative-only match at `final_tolerance` is load-bearing.
- `DefaultPointMatchTolerance()` changed value (same-point → final_tolerance). It is not an identity
  input (not digested), so no digest/version bump is required; the C++ and Python tests that pinned
  the old formula were updated, and `same_point_tolerance()` is newly exposed for callers who want the
  clustering scale.
- The multiplicity classification itself is unchanged — this is purely how a point is *looked up*
  after the fact. A wrong multiplicity is a `ComputeMultiplicities` tolerance question (see ADR-0017),
  separate from this lookup fix.
