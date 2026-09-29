# ADR-0069: An accuracy estimate names its coordinates, and the bare name is retired

**Status:** Accepted
**Date:** 2026-09-29

## Context

A solution's metadata carried two accuracy estimates.  Both are the distance, in the infinity
norm, between the endgame's last two approximations of the root.

| field | coordinates |
|---|---|
| `accuracy_estimate` | the solver's internal ones: homogenized, on the patch |
| `accuracy_estimate_user_coords` | the user's: dehomogenized |

The plain name went to the quantity fewer people want.  A user reads `accuracy_estimate` as
the accuracy of their solution in their own coordinates, because those are the only
coordinates they think in.  It was not that.  The internal estimate is what the endgame
compares with `final_tolerance`, so it behaves like a number of correct digits; the absolute
error in the user's variables is that, times roughly the scale of the solution.  A solution of
order 1000 with six correct digits has an absolute error near `1e-3`.  The two differ by the
scale, and nothing in the name said which was which.

The cellular decomposition port found the consequence.  It reads the internal estimate as a
relative error, and once multiplied a tolerance by the scale a second time.

## Decision

Both estimates are named for their coordinates, and there is no field named plain
`accuracy_estimate`.

- `accuracy_estimate_internal_coords`: the former `accuracy_estimate`.
- `accuracy_estimate_user_coords`: unchanged.
- `accuracy_estimate` is **retired**.  In Python, reading it raises a `RuntimeError` that names
  both successors.  It is not an `AttributeError`, because `getattr(md, 'accuracy_estimate',
  None)` swallows that and the caller goes on with `None`.

The bare name is retired, not given to the user-coordinates estimate.  The same name
returning a different number is the failure a rename exists to prevent: every reader would
keep running, on a value off by the scale of its solution, with no error anywhere.

In the records each estimate is written under its new key.  The reader accepts the earlier
key `accuracy_estimate` for the internal one, so a record written before the rename loads.

## Consequences

- Code that read `accuracy_estimate` fails loudly and is told what to read.
- **Do not re-add `accuracy_estimate` as an alias** for either estimate in the 4.x line.  If the
  plain name is ever to mean the user's coordinates, that is for a later major release, once
  nothing is left reading it with the old meaning.
- **Do not change the retirement to an `AttributeError`**, however conventional that is for a
  missing attribute.
- A new quantity that exists in both coordinate systems is named for the one it is in.
- `accuracy_digits` is derived from the internal estimate and keeps its name: a count of
  digits has no coordinates.
- The estimates are distances between successive approximations.  Whether they bound the true
  error is a separate question and is not settled here: the triangle inequality makes the
  estimate the error of the PREVIOUS approximation, give or take the error of the last, so
  it bounds the true error from above only when the last round at least halved the error.
- Pinned by `python/test/zero_dim/accuracy_estimate_names_test.py` and by
  `the_accuracy_estimates_are_recorded_under_names_that_say_their_coordinates` in
  `core/test/nag_algorithms/zero_dim_records.cpp`.
