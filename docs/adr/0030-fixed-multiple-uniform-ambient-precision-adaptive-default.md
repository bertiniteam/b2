# ADR-0030: Fixed-multiple solves use one uniform ambient precision; ZeroDim defaults to adaptive

**Status:** Accepted
**Date:** 2026-06-23

## Context

A trivial fixed-multiple zero-dim solve threw out of the box:

```python
ZeroDim(sys).solve()   # sys = x**2 - 1
# RuntimeError: start point for fixed multiple precision tracker has differing
# precision from default (20!=16), tracking cannot start
```

Three precision sources must agree for a `MultiplePrecisionTracker` solve, and they
did not:

- the `MultiplePrecisionTracker` fixes its working precision to `DefaultPrecision()`
  (= 20) at construction, and its loop initialization *requires*
  start-point precision == thread precision == tracker precision;
- but `ZeroDimConfig::initial_ambient_precision` — which drives the thread precision
  and the start-point precision — defaulted to `DoublePrecision()` (= 16).

So the start points and thread sat at 16 while the tracker sat at 20: `16 != 20`,
throw. The fixed-multiple tracker also *cannot* re-precision itself
(`PrecisionSetup(FixedPrecisionConfig)` is a no-op; the class is documented as
fixed), so the precision genuinely has to be made uniform up front.

Separately, `ZeroDim(...)` defaulted to **fixed multiple** precision — the mode most
likely to surprise a new user with exactly this kind of precision-matching
requirement.

## Decision

1. **In a fixed-multiple solve the ambient precision is uniform and config-sourced.**
   `ZeroDimConfig<ComplexT>::initial_ambient_precision` now defaults per complex type
   via `DefaultInitialAmbientPrecision<ComplexT>()`: `DoublePrecision()` for
   `dbl_complex`, `DefaultPrecision()` for the multiprecision type. For a default
   fixed-multiple solve the tracker, thread, start points, and config precision are
   then **one value** — `DefaultPrecision()` — and that value follows whatever
   `default_precision(n)` was set to *before the solver is constructed*. The rule:
   in fixed multiple, the precision is the same everywhere, and it is set by the
   config (ultimately by `default_precision`).

2. **The `ZeroDim` Python factory defaults to adaptive (AMP), not fixed multiple.**
   `ZeroDim(system)` is now Cauchy endgame + adaptive precision + inferred start
   system. Adaptive is the robust general-purpose default that "just works"; fixed
   precision is an opt-in for users who want to control the cost/precision tradeoff
   themselves.

**Out of scope (deliberately not built):** letting a user *override*
`initial_ambient_precision` to a value different from the tracker's
construction-time `DefaultPrecision()` and have the fixed tracker re-precision to
match. That needs a real precision setter on `MultiplePrecisionTracker`. The fix
above makes the precision uniform and `default_precision`-driven, which is what the
bug required.

## Consequences

- For a fixed-multiple solve at precision `n`, call `default_precision(n)` **before**
  constructing the solver; the whole solve then runs uniformly at `n`.
- `ZeroDim(sys).solve()` works out of the box (adaptive); the previous
  precision-mismatch throw is gone.
- The fixed-multiple `ZeroDimConfig` is shared with the adaptive tracker, so this
  also nudges AMP's *starting* precision from 16 to `DefaultPrecision()`; AMP adapts
  upward from there, so this is benign.
- Related to ADR-0007 (rational config constants and the default-construct precision
  trap): both are about keeping multiprecision state internally consistent rather
  than letting a stray default leak in.
