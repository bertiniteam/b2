# ADR-0013: Solutions are reported in user coordinates by default; internal coordinates are an explicit choice

**Status:** Accepted
**Date:** 2026-06-11

## Context

The zero-dim algorithm clones the supplied system (`CloneGiven`, policies.hpp),
homogenizes and auto-patches the clone, and tracks in those coordinates. Solutions
were stored AND returned in that internal representation: a 2-variable affine system
yields 3-vectors `[h, x, y]`. Nothing in the python surface dehomogenized, the
zero-dim python test asserted only solution *counts*, and the representation leaked
solver implementation into every consumer — during the sympy-interop investigation
it cost real time to recognize the "wrong-looking" values as patched projective
coordinates.

Both representations are legitimately needed: user coordinates for reading off
answers, internal coordinates for continuing work (start points for further
tracking, re-using the same patch, refinement).

## Decision

- **The unqualified accessor means user coordinates.** A user probably expects
  solutions in their own variables; that is the default mode.
- **Language-appropriate idiom for the toggle, same semantics:**
  - python: pandas-style bool kwarg — `solutions(user_coords=True)`; passing
    `False` is *declining* the usual behavior to get internal coordinates.
  - C++: explicit method names — `SolutionsUserCoords()` /
    `SolutionsInternalCoords()` (no mysterious bools at C++ call sites). The old
    `FinalSolutions()` is gone ("Final" was filler; solutions are solutions).
- **Reference semantics, no recompute:** the user-coordinates container is
  dehomogenized once, lazily, into a solver-owned cache (invalidated when a solve
  starts), returned by const reference; the python binding wraps it
  (`return_internal_reference`) so repeated `solutions()[i]` neither recomputes nor
  copies the container.
- **The lift completes the round trip:** `System::HomogenizePoint` /
  `homogenize_point` takes a user-coordinates point to internal coordinates —
  insert the homogenizing coordinate (value 1) per affine group in FIFO layout,
  then rescale onto the system's patch. It is the inverse of `DehomogenizePoint`.
- **Coordinate labels come from the system that owns the representation:**
  `variable_ordering` on your original system labels user coordinates;
  `target_system()` (the prepared clone, now exposed) labels internal ones. No
  separate ordering API.

## Consequences

- **Breaking change in python:** `solutions()` previously returned internal
  coordinates. Code that wants the old behavior writes
  `solutions(user_coords=False)`. The new value-asserting zero-dim tests pin both
  representations and the round trip.
- House style for future representation/mode toggles (e.g. an eventual NID
  surface): bool kwarg in python with the default being what a user expects;
  named methods in C++.
- The lazy cache assumes post-solve serial access, like the algorithm's other
  accessors; both `Solve()` and `RunParallel()` invalidate it on entry.
