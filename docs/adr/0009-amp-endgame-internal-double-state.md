# ADR-0009: AMP endgame needs internal double-precision state, mirroring the tracker

**Status:** Proposed (design conclusion; implementation deferred)
**Date:** 2026-06-09
**Related:** ADR-0002 (endgame `Run` single-arg API), ADR-0007 (rational-config precision trap)

## Context

A performance investigation (PR with `SetStartPrecision` + explicit endgame-boundary precision,
"Fix 1 + Fix 4") removed the dominant cost in adaptive (AMP) zero-dim solves: the pre-endgame
homotopy track now runs in `std::complex<double>` for well-conditioned paths, giving ~5–7× on the
`medium`/`large` benchmarks. That left a residual: the **endgame** is still far slower than a
fixed-double solve of the same problem (`medium`: ~8.3 s adaptive vs ~0.5 s `mptype:0`).

We chased two *lighter* alternatives to a full endgame change and both dead-ended at the same wall:

1. **A `variant<Vec<dbl>,Vec<mpfr_complex>>` as the AMP *public contract* type.** Feasible and
   moderate in the C++ core, but it leaks to the boundaries: Boost.Python/eigenpy have no
   `std::variant` support (exposing it needs ~100 LOC of custom converters and breaks
   `solutions()[i][j]` / `final_approximation()`), MPI needs variant serialization, and it does
   nothing for *accuracy*. It also runs against the project's deliberate design that an AMP
   algorithm carries **one authoritative multiprecision type** (Bertini 1 experience: threading two
   types through every signature is miserable to reason about).

2. **"Make precision-16 mpfr honestly mean double" / lower `LowestMultiplePrecision()` 20→16.** An
   audit + a direct measurement disproved the premise:
   - `LowestMultiplePrecision()=20` is the *escalation floor* (skip the useless 16–19-digit mpfr
     dead zone), not the widening precision; asserts (amp_tracker.hpp:1507,1551) forbid mp ≤16.
   - The tracker's double output is *already* honest mpfr@16 (`CopyFinalSolution`, amp_tracker.hpp:618,
     at thread precision 16).
   - **Measurement** (instrumented `medium_seeded`, probes in `base_endgame::Run`,
     `CauchyEndgame::CircleTrack`, `AMPEndgame::RefineSampleImpl`, since reverted): the endgame
     **already runs at precision 16** — 229/243 boundary points at 16; 1859/1975 `CircleTrack`s and
     7436/7900 `RefineSampleImpl`s at `tracker_current_prec=16` (the remainder at 40 are genuinely
     ill-conditioned paths that correctly escalated). So the endgame's *tracker calls already execute
     in `std::complex<double>` internally*. There is no precision-20 contamination to remove, and no
     speedup to be had from this idea.

**Root of the residual.** Even though the tracker work inside the endgame runs in double, the
endgame *carries every value as `mpfr_complex@16`*. Across ~10k circle-track/refine operations it
pays mpfr allocate/convert churn at every tracker boundary, plus its own extrapolation / norm /
dehomogenization arithmetic in mpfr@16. That cost is intrinsic to the endgame's base type being
`mpfr_complex`. The only way to remove it is for the endgame to *execute in
`std::complex<double>`* on the regimes where double suffices — which is precisely what the **tracker
already does internally** and the endgame does not.

The tracker solved this years ago: `current_space_`/`tentative_`/`temporary_` are a
`tuple<Vec<dbl>,Vec<mpfr_complex>>` (via `NeededTypes`/`detail::TypeList`, config.hpp:333,
typelist.hpp:45); `current_precision_` selects the live arm; `TrackerIteration()` dispatches **per
step** to `TrackerIteration<dbl>` vs `<mpfr_complex>`; and `MultipleToDouble`/`DoubleToMultiple`
convert the one current point on a precision transition. All of this is *internal* — the tracker's
public contract is `mpfr_complex`.

An earlier attempt to give the endgame a double path by templating `Run`/`RunImpl` on the complex
type and adding a `RunImpl<dbl>` cascaded and was reverted: the endgame drives the **tracker** for
its sample circles (`GetTracker().TrackPath(...)`, cauchy.hpp:1040), and the AMP tracker's public
`TrackPath` is mpfr-only. A `RunImpl<dbl>` would therefore require a *public double `TrackPath`* on
the tracker — a second public type, contradicting the single-authoritative-type design.

## Decision

The AMP endgame should gain **internal double-precision state**, mirroring the tracker's existing
internal variant/tuple pattern, **behind an unchanged `mpfr_complex` public contract**.

Two distinct ideas must not be fused (an earlier draft did fuse them):

- **Precision as honest provenance (prerequisite cleanup).** Today precision *shadows* the regime —
  a widened double can land at 20 (via the parser floor and the double-meaning of
  `LowestMultiplePrecision`), indistinguishable from a number that genuinely needed 20 digits. That
  is a *bug*, not an inherent limit. Once the shadowing is fixed — decouple the two meanings of
  `LowestMultiplePrecision` (dead-zone/escalation boundary vs. parser storage floor) and stop
  widening double-regime values past `DoublePrecision()` — **a number's precision faithfully records
  the precision of the process that computed it.** Precision then *is* provenance (of precision), and
  serves as the regime signal directly, exactly as the tracker's `current_precision_` tag already
  does. (It still does **not** record *accuracy* — how many of those digits are correct — which
  remains separate metadata.)
- **The double arm is for speed, not signaling.** Even with precision honest, doing arithmetic in
  `mpfr@16` is ~100× slower than `std::complex<double>`. The dbl arm exists to *execute* the double
  regime in native double; precision already tells you *which* regime you're in. The two are
  orthogonal and reinforce each other: a dbl-arm value, widened to the contract, is an `mpfr@16` that
  honestly reports "computed in double."

Concretely:

- Set the endgame's `NeededTypes` to `{dbl, mpfr_complex}` (as the AMP tracker does), instead of
  `TypeList<BCT>`. The sample/time/derivative containers (`cauchy_samples_`, `pseg_samples_`, …) are
  **already** `TupleOfSamps`/`TupleOfTimes` over `NeededTypes`, so they gain a `dbl` slot
  automatically. Promote the remaining single-type state — `final_approximation_`,
  `previous_approximation_` (currently `Vec<BCT>`, base_endgame.hpp:124), and `start_time_` /
  `target_time_` — to the same tuple shape.
- Track a **current precision regime** in the endgame and **dispatch the main loop per iteration**
  on it, the way `Tracker::TrackerIteration()` does — rather than the present monolithic
  `RunImpl<ComplexT>` (cauchy.hpp:1076) whose type is fixed for the whole run and so cannot change
  regime mid-loop. The loop body becomes an outer driver owning the regime + tuple state, calling
  templated per-iteration step functions.
- On a double→multiprecision escalation (precision is expected to *rise* as the path converges
  toward a singularity), **migrate the accumulated sample history** (the deques), not just one
  point. Extend the existing deque-walking precision machinery — `ChangePrecision`,
  `adaptive::EnsureAtUniformPrecision` (adaptive_precision_utilities.hpp), which today re-precision
  mpfr samples to a uniform level — to convert the `dbl` arm into the `mpfr` arm at the regime
  boundary (the endgame analog of `DoubleToMultiple`). Template `RefineSampleImpl` (currently
  hardcoded `Vec<mpfr_complex>`, amp_endgame.hpp:94).
- Prefer the **tuple + regime-tag** representation (the tracker's choice, and the containers are
  already tuples) over a `std::variant` of the whole state. A variant is more honest that only one
  arm is live, but forces `std::visit` at every state access for no functional gain here.

Because all of this is internal, the **public contract stays `mpfr_complex`**: the Python bindings
and the MPI serialization (`parallel/path_result.hpp`) are **unaffected**, and the design stays
inside the one-authoritative-type philosophy — it is the *same* compromise already accepted for the
tracker, applied one level up.

## Consequences

**Positive**
- Removes the endgame's intrinsic mpfr@16 cost on well-conditioned regimes by executing in
  `std::complex<double>` there; this is the remaining lever toward the ~88× that `mptype:0` shows is
  achievable.
- Philosophically consistent: internal double, mpfr contract — identical posture to the tracker.
- No Python/MPI blast radius (contract unchanged), unlike the variant-as-contract alternative.

**Negative / cost**
- The endgame's Cauchy/PowerSeries main loops must be **restructured** from a monolithic
  `RunImpl<ComplexT>` into a per-iteration regime-dispatched driver — the bulk of the work, and
  numerically delicate (this is the most sensitive code in the library).
- New complexity the tracker does *not* have: migrating the **sample-history deques** across a
  mid-run precision transition, not just a single current point.
- `RefineSampleImpl` and any other hardcoded-mpfr endgame internals must be templated; the fixed
  (double / multiple) precision endgames must continue to compile with single-type `NeededTypes`.

**Does not address**
- **Accuracy** (how many digits are *correct*) remains orthogonal: it is carried separately in
  `SolutionMetaData` (`accuracy_estimate`, `condition_number`, `newton_residual`; computed in
  cauchy.hpp). This internal-double change neither solves nor worsens it. A unified
  precision-and-accuracy-aware solution type is a separate, larger decision.

## Status / next step

Recorded as the agreed *direction*; implementation is **deferred** (the loop restructuring is
substantial and delicate). The shipped Fix 1 + Fix 4 stand independently. If/when pursued, the first
implementation milestone is hoisting the Cauchy main loop into a per-iteration regime-dispatched
driver over tuple state, with a double→mp deque-migration step, validated against the full C++ +
pytest suites and the seeded benchmarks.
