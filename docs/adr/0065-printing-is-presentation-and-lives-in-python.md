# ADR-0065: printing an object is presentation, and it lives in the python layer

**Status:** Accepted
**Date:** 2026-09-21

## Context

106 of the 188 public classes printed the default `<bertini._pybertini.tracking.AMPTracker
object at 0x7f3c...>`: the type you already knew, and an address you cannot use.  That covered
every tracker, every endgame, all 38 observer event types, the observers and collectors, the
linear algebra decompositions and the decomposition scaffolding (#99).

The count overstates the work.  Most of those are one template exposed at three precisions, so
roughly fifteen functions cover all 106.

Two places to put them.  In the bindings, `.def("__repr__", ...)` once per exported template,
ideally backed by a C++ `operator<<` so the same text serves both languages.  Or in the python
layer, attaching `__repr__` to the bound classes at import -- which works, because a
Boost.Python class accepts the assignment and `str()` follows `repr()` for any type that does
not define its own.

## Decision

**The text lives in `bertini/_repr.py`, attached at import.**  A repr is presentation: what a
user should see when they print an object, which is a different question from what the object
is, and it changes on different grounds.  Keeping it in one python file puts every such
decision in one place, where it can be read as a whole and changed without a rebuild.  A C++
`operator<<` would serve a different audience with a different need; when one is wanted it can
be written then, without this file having pre-empted its shape.

**Machinery is described in angle brackets; values keep the call form.**  A tracker reads
`<AMPTracker: homotopy in 2 variables, RKF45, tolerance 1e-05, adaptive precision at 30 digits,
at t=0.31 after 47 steps>`.  None of these objects could be rebuilt from any text -- most hold a
C++ object no python expression can name -- and the brackets say so, where `AMPTracker(...)`
would imply a round trip that does not exist.  Configs, metadata and `StraightLineHomotopy`
already read `SteppingConfig(initial_step_size=..., ...)`, which is close to something you could
type back, and they keep it.

**A tracker's description is live.**  It reports where the path is now, not only how it was
configured, because the moment you most want to print a tracker is inside an observer callback
mid-path.

**No repr may raise.**  Every accessor goes through a helper that swallows failure and omits
that clause, since they all reach into C++ state that may not exist yet.  A repr that throws
turns an ordinary `print` in a debugging session into a traceback about printing.

## Consequences

- **Do not move this text into the bindings as a "proper" home** without deciding separately
  what C++ callers need from `operator<<`.  The python spelling is the python audience's, and
  the two need not match.
- **Do not make any of these reprs round-trippable**, or name them as if they were.  The
  objects hold C++ state that no text reconstructs; the angle brackets are the honest signal.
- A test walks every public class and fails if any of them inherits the default repr
  (`python/test/classes/repr_test.py`), so a new binding arrives with a description or not at
  all.  A type that genuinely has nothing to report still names itself: `<GoryDetailLogger>`.
- Several of these types expose no state at all, so their descriptions are bare names.  That is
  a gap in what the bindings expose, not in this file -- when such a type gains accessors, its
  description should gain the corresponding clause.
