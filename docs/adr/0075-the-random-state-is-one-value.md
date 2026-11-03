# ADR-0075: The random state is one value, and every draw comes from one stream

**Status:** Accepted
**Date:** 2026-11-03

## Context

A long computation that is written out between its stages must continue, after a restart, with the numbers an uninterrupted run would have drawn (#518).  The cellular decomposition of the Barth decic on the 3-sphere runs for many hours in eight stages, and its stages draw from b2's global random state: `random_vector`, `AdvanceSolveSeed` (formerly `DeriveSolveSeed`), the gammas and patches of every solve.  b2 could set a seed and read it back, but could not get the state the draws had reached, nor set it.

The state was split across interfaces.  ADR-0044 moved every draw onto the pinned `b2rand/1` draw stream but kept the thread's `std::mt19937` engine "for any non-identity consumer", seeded alongside the stream; nothing drew from it any more, yet `SetGlobalSeed`, `ReseedThisThread` and the solve's `ScopedSessionState` still seeded, saved and restored it.  A caller saving "the random state" had to know to take both, and the global seed and the solve-seed ordinal besides.

## Decision

- **One stream.**  The `mt19937` engine and `ThreadEngine()` are removed.  Every draw comes from the thread's `b2rand/1` draw stream.  No draw changes: none used the engine.
- **One value.**  `bertini::RandomState` (`random.hpp`) is everything the calling thread's next draws depend on: `global_seed`, `seed_chosen`, `solve_ordinal`, and `draw_stream`.  `GetRandomState()` returns the calling thread's; `SetRandomState(state)` sets it, so the draws after setting a state are the draws that followed getting it.  States compare with `==`.
- **The calling thread only.**  Per-path and per-worker streams are not part of it: each is reseeded from the global seed and its index when it is used (`ReseedThisThread`, `DerivedWorkerSeed`, ADR-0071).
- **Public fields.**  A caller may build a state field by field, including one no run reached.
- **Text is serialization, at the edge.**  `ToText()` and `FromText(text)` write and read versioned plaintext, `b2randstate/1`, one labelled line per field (the draw stream's key and current block in hex), for plaintext records.  `FromText` refuses any other text with `std::invalid_argument`.  The state version goes with the draw derivation: a new `b2rand` version means a new `b2randstate` version.
- **Archives carry the text.**  Boost.Serialization (`save`/`load`) writes the `ToText` text, so every archive -- text, binary, an MPI message -- carries the one versioned form, and a load refuses a state of another version.
- **Python.**  `bertini.random.RandomState` with the same fields, `to_text()`, `from_text()`, `==`, and pickling as its text; `get_random_state()` and `set_random_state(state)`, named beside `set_random_seed`.

## Consequences

- Do not add a second source of randomness beside the draw stream.  A new kind of draw is a new function of the stream; a new piece of per-thread random state is a new field of `RandomState`, with its line in the text and a new `b2randstate` version.
- Do not save or restore part of the state by hand.  Code that brackets a computation's randomness takes a `RandomState`.  `ScopedSessionState` (a solve leaves the session as it found it) still restores only the draw stream and the precision, deliberately not the seed or the ordinal: a seedless solve rekeys to a derived seed and an MPI solve sets rank 0's, and whether a solve should leave those as it found them is an open question, not decided here.
- "The calling thread only" may prove wrong: a program that draws setup randomness on several threads of its own has a state per thread, and a run continued from one of them draws differently on the others.  If such a program appears, the state grows to cover it; until then the restriction is documented on `RandomState`.
- A state written by one version of b2 is read only by a version with the same `b2randstate` version; a record carrying one is otherwise asked again, as records are across encoding versions.
