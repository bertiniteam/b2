# ADR-0062: the system judges its own points, and auxiliary coordinates are part of its identity

**Status:** Accepted
**Date:** 2026-09-20

## Context

Bertini asks "how big is this point" three times, with three thresholds: the tracker truncates a
path exceeding `path_truncation_threshold`, the endgame abandons one exceeding
`Security::max_norm`, and the solver classifies an endpoint past `endpoint_finite_threshold` as at
infinity.  Realness is a fourth question of the same shape, over the imaginary parts.

All four measured the whole dehomogenized point.  Several standard constructions carry coordinates
that exist for the construction rather than the answer — a critical-point system in null-vector
form, `f; v^T M; patch`, where `v` must carry its own patch or the system admits the zero null
vector.  That patch fixes a normalization nobody chose on geometric grounds, so `v`'s direction is
the meaningful object and its magnitude is an artifact: not intentionally divergent, simply not
governed.  Judging a point on the larger of two independently uncontrolled scales is a kind error
rather than a tuning problem, because no threshold can be correct for an ungoverned quantity
(#403).

A caller can always replace the classification.  A caller cannot un-abandon a path the tracker
truncated partway.

## Decision

**A system says which of its coordinates are auxiliary, and answers the questions itself.**

- `System::SetAuxiliaryVariableGroups` and `System::SetAuxiliaryCoordinates` declare them, by FIFO
  group index and by index into a point in user coordinates.  The two are unioned; empty means
  none, which is the previous behaviour exactly.  Indices are into USER coordinates, so they
  survive `Homogenize()` and `AutoPatch()`.
- `System::IsFinite(point, threshold)` and `System::IsReal(point, threshold)` are the judgements,
  and `InfinityNormOfDehomogenized` — which already documented itself as the single canonical
  measurement — honours the declaration.  The tracker and the classifier had each hand-rolled the
  measurement beside that chokepoint; they call the system now.  Thresholds stay where they are, in
  their own configs: what was missing was never a setting, it was that nobody owned the judgement.
- **It is part of the system's content identity** (`b2sysenc`).  The tracker truncates on the
  coordinates that are not auxiliary, so two systems differing only here do not compute the same
  thing and must not recall each other's records.  Landing in unreleased 3.5, this edits
  `b2sysenc/2` in place per ADR-0061.
- Declaring every coordinate auxiliary is refused: a point with nothing to judge would be
  unconditionally finite and real.

Rejected: **a caller-supplied predicate.**  A function cannot join the ask, so the same question
would recall a differently-classified answer with nothing in the record able to explain it — the
records-identity flaw of #420, deliberately widened.

Rejected: **judging on the first variable group by default.**  In a multihomogeneous system every
group is geometric, so "the first group is the real one" is false for most multi-group systems.
The grouping says where the patches are, not which coordinates carry meaning.

Rejected: **a config read by all three layers.**  It was the first design, and the System is the
better owner: it already knows its coordinates, and a config would have had to be plumbed to the
tracker, the endgame and the solver separately while meaning the same thing in all three.

## Consequences

- Do not reintroduce a direct norm-and-compare at a call site.  Four call sites, one
  implementation; a fifth caller asks the system.
- b2 attaches no meaning to "auxiliary" beyond the exclusion.  The coordinate is tracked, recorded
  and returned as before.  Do not let it grow into a general notion of variable roles without
  deciding what else it would mean.
- Adding or removing a declaration changes the digest, so it is a new ask and is recomputed.  This
  is correct and should not be "optimized" away.
- `System`'s copy constructor is hand-written member by member — a third hand-maintained mirror
  beside `serialize` and the canonical encoder — and forgetting a member there means every `Clone`
  silently drops it.  That is how this feature first appeared to work on the caller's own object
  and do nothing inside the solver.  A member added to `System` must be added in four places: the
  class, `serialize`, the copy constructor, and (if identity-affecting) the canonical encoder and
  its decoder.
