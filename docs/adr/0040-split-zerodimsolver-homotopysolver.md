# ADR-0040: Split ZeroDim into HomotopySolver (engine) + ZeroDimSolver (algorithm)

**Status:** Accepted

## Context

`algorithm::ZeroDim<Tracker, Endgame, System, Policy>` conflated two different jobs, selected by a
system-management *policy* template parameter:

- **The continuation engine** — given a homotopy and a source of start points, track each path
  through the tracker and endgame, resolve crossings, classify the endpoints, report. This is ~90%
  of the class. The `policy::RefToGiven` configuration *was* this: the user supplies
  target/start/homotopy by reference and `SystemSetup` is a no-op (the "user homotopy" path).
- **The zero-dim algorithm** — given a polynomial system, form a start system and a homotopy, then
  solve. This was `policy::CloneGiven`: clone + homogenize + `FormStart` (via an injected
  `StartSystemFactory`) + `FormHomotopy`, plus `ConsistencyCheck`.

PR #48 had already de-templated the start-system *type* (held polymorphically via a factory) but
left the policy parameter in place. The two roles also want different *names*: downstream algorithms
(slice-moving, regeneration, monodromy, NID) are all "track a set of points through a homotopy" and
should call the engine directly, not masquerade as a zero-dim solve — exactly the wart issue #261
flagged about `user_homotopy`.

## Decision

Replace the one policy-parameterized class with **two concrete classes**, retiring the policy:

- **`HomotopySolver<Tracker, Endgame, System>`** is the engine. It holds the homotopy, start
  system, and target by *reference* (the old `RefToGiven` storage, now its only storage) and owns
  all the machinery: pre-endgame tracking, endgame, crossing resolution, classification, reporting,
  threaded + MPI solve. The MPI system-broadcast became a virtual `DistributeSystems()` hook
  (a no-op here: a user-supplied homotopy is the user's cross-rank responsibility).

- **`ZeroDimSolver<Tracker, Endgame, System>`** is the algorithm. It clones the user's target,
  homogenizes/patches it, builds a start system (via the injected factory) and a homotopy, then
  drives the engine. It **is-a** `HomotopySolver`, with a first private base `OwnedHomotopy<System>`
  that builds and owns the systems *before* the engine base references them (base-initialization
  order = declaration order). It overrides `DistributeSystems()` to broadcast its owned systems.

`HomotopySolver` is a first-class name in C++ and Python; `ZeroDimSolver` replaces the old `ZeroDim`.
The old `ZeroDim` / `...UserHomotopy` names are deleted outright (no alias) — a deliberate loud break.

### Why inheritance, not a composed member

The natural reading of "ZeroDimSolver calls HomotopySolver" is a member + delegation. We chose
inheritance instead because a composed member would force `ZeroDimSolver` to *separately* be
`detail::Configured` and `Observable` and forward the entire solve/accessor surface — and, worse,
would duplicate the configuration state (two `Configured` bases: which one does the solve read?).
With inheritance there is exactly one `Configured`/`Observable`/tracker/endgame, owned by the engine
base, and `ZeroDimSolver`'s setup writes into it. The reuse goal ("no duplicated logic; the engine
is independently usable for monodromy/regen") is fully met: `HomotopySolver` stands alone and does
not depend on `ZeroDimSolver`.

### Where the new behaviors live

- **Classification stays in `HomotopySolver`** (finite/real/multiplicity/singular).
- **`ZeroDimSolver` adds the feasibility behaviors** the algorithm owns:
  - `ConsistencyCheck` rejects path-variable/non-polynomial/under-determined targets (the last with a
    friendly "positive-dimensional, not zero-dimensional" message).
  - `SquareUp` randomizes an over-determined target down to square (`System::Randomize`), keeping the
    original; after the solve it overrides the now-virtual `PostEGAction` to re-evaluate the original
    system at each finite endpoint **in double precision** and flag the extraneous solutions the
    squaring introduces. These are marked with a dedicated **`is_nonsolution`** metadata flag
    (orthogonal to `is_finite`, so they stay geometrically finite rather than masquerading as
    at-infinity); they are excluded from the finite/real/singular accessors, surfaced by
    `Nonsolutions()`, and counted in `SolveReport.num_nonsolutions`. `is_nonsolution` is load-bearing
    for the regeneration cascade, which must identify and discard nonsolutions. Exposed as
    `was_randomized()` / `randomization_matrix()`.
  - `RankCheck` rejects a square-by-count-but-positive-dimensional target via a generic-point
    Jacobian rank test.

To support the filter, `PostEGAction` is protected + virtual and the per-endpoint metadata/endpoints
are protected.

The Python surface gains a `solutions(**flags)` getter on both classes (finite genuine by default,
with `singular`/`nonsingular`, `real`/`nonreal`, `infinite`, `nonsolution` toggles); `all_solutions()`
remains the raw per-path list, so this is additive.

## Consequences

- `common/policies.hpp` (`CloneGiven` / `RefToGiven` / `CloneTarget` / `SysMgmtPolicy`) is **deleted**.
  `StartSystemFactory` / `MakeStartFactory` were not ownership policy — they are start-system
  *selection* — so they were rehomed from `policy::` to `start_system::` in `system/start_base.hpp`.
  NID (placeholder scaffolding whose `Solve()` throws) dropped its policy base and owns a cloned
  target inline.
- The **`SystemView`** idea (a System + optional bound time, considered as a way to shed the
  System/Policy template params) is **obviated**: the concrete two-class split achieves the
  engine/algorithm separation directly.
- The bigger deferred engine — a `RunBatch<Tracker>`-only kernel with type-erased `StartPoints` /
  `PathHandler` interfaces collapsing the instantiation lattice — remains future work. This split is
  the first concrete extraction of the continuation primitive, not that erasure.
