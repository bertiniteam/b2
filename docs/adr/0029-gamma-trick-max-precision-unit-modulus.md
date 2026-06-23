# ADR-0029: The gamma-trick constant is a unit-modulus complex generated at maximum precision

**Status:** Accepted
**Date:** 2026-06-23

## Context

The gamma trick forms the homotopy `H = (1-t)*target + gamma*t*start`. A random
`gamma` moves the path-crossing/singular bad set off the real `t` interval; the
solve tracks `t: 1 -> 0`. `gamma` only needs to be *generic* (the bad set is
measure zero), so historically it was `node::Rational::Rand()` — a complex
rational of **arbitrary** norm.

Two problems with that choice surfaced:

1. **Arbitrary norm ill-conditions the homotopy.** When `|gamma|` is far from 1,
   the two terms `(1-t)*target` and `gamma*t*start` differ wildly in magnitude
   along the path, degrading conditioning and path separation. The standard gamma
   trick uses a unit-modulus `gamma` precisely to keep the two terms comparably
   scaled.

2. **A `node::Float` constant caps at its creation precision.** A `Float` stores a
   single `highest_precision_value_` (see `function_tree/symbols/number.hpp`) and
   can only *serve* values up to that stored precision — it cannot synthesize more
   digits on demand. `node::Rational`, by contrast, holds exact rationals and
   converts to `mpfr` at the current thread precision on every evaluation, so a
   Rational `gamma` was effectively unlimited-precision. Once `gamma` is a `Float`,
   if it is created at the current default precision it becomes a hard ceiling: an
   adaptive-precision (AMP) tracker that climbs past that precision can get no more
   correct digits of `gamma`, so the homotopy constant silently bottlenecks the
   achievable accuracy.

Separately, the unit-modulus generator itself was wrong: `rand_unit()` /
`RandomUnitAssign()` normalized by `sqrt(abs(z))` instead of `abs(z)`, so
`bertini.random.complex_unit()` — documented "magnitude 1" — actually returned
magnitudes ~0.76–1.09.

## Decision

- **`gamma` is a unit-modulus complex (`|gamma| = 1`)**, drawn from
  `bertini::multiprecision::RandomUnit()`, which now normalizes by `abs(z)`.
- **It is generated at `MaxPrecisionAllowed()`** (currently 1000 digits;
  `num_traits.hpp`), i.e.
  `node::Float::Make(RandomUnit(MaxPrecisionAllowed()))`, at both homotopy sites
  (`MakeHomotopy`, `MakeMovingHomotopy` in `core/src/system/system.cpp`). This
  mirrors how patch coefficients are generated (`system/patch.hpp` also uses
  `MaxPrecisionAllowed()`): a structural constant of the homotopy must carry enough
  digits to never limit an adaptive tracker.

`|gamma| = 1` is a **conditioning** choice, not a genericity one — a uniformly
random point on the unit circle is just as generic as one off it (the bad set is
still measure zero). The unit circle is chosen only to keep `(1-t)*target` and
`gamma*t*start` comparably scaled.

## Consequences

- Do **not** revert `gamma` to `node::Rational::Rand()`, and do **not** generate
  the `Float` `gamma` at the ambient/default precision — either reintroduces a
  precision ceiling or a magnitude imbalance. If a future representation lets a
  constant re-precision on demand, the max-precision generation can be revisited.
- `gamma`'s precision tracks `MaxPrecisionAllowed()`; raising that cap raises
  `gamma`'s digits automatically.
- This **complements, does not replace, ADR-0017** ("flakiness is a tolerance
  problem, not a seeding problem"). Improving `gamma`'s conditioning is legitimate
  numerical hygiene; it is not a substitute for fixing a flaky solve with
  tolerances/precision, and correctness must still hold across random `gamma` draws
  (count *distinct* solutions, per ADR-0015/0017). Empirically the unit-modulus
  `gamma` reduced cyclic-5 crossings under crossing-provoking settings, but that is
  a conditioning improvement, not a correctness guarantee.
- `complex_unit()` now genuinely has modulus 1 (regression test
  `python/test/random/complex_unit_test.py`).
