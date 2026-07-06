# ADR-0007: Rational config constants, and the DefaultConstruct static-precision trap

**Status:** Accepted
**Date:** 2026-06-08

## Context

Stepping and endgame config structs hold numeric knobs — `SteppingConfig`'s
`initial_step_size`, `max_step_size`, `step_size_success_factor`,
`step_size_fail_factor`, and `EndgameConfig::sample_factor`. The type chosen for
these fields flip-flopped three times on this branch before settling, and the
reason is subtle enough to record so it is not re-litigated.

### The decision history

1. Originally `mpq_rational` (exact rationals), as a deliberate caution over
   Bertini 1's doubles.
2. `a434145d` — changed them to `mpfr_float`. Rationale at the time: the
   classic-input parsers already laundered every value through `double`
   (`NumTraits<double>::FromString`) before storing, so the "exact" guarantee was
   already void; and `mpfr_float` (over `double`) preserves exponent range, which
   `min_step_size` needs below double's ~1e-308 floor in high-precision AMP runs.
3. `d06448e7` — changed `sample_factor` to `double` as an expedient crash fix
   (see below).
4. `6abf488b` / `040429b3` — reverted to `mpq_rational` for the stepping and
   sample-factor constants, with an exact decimal-string parser. **This is the
   final state.**

### The trap that forced the revert

`core/include/bertini2/detail/enable_permuted_arguments.hpp` defines a non-local
static used to supply default-constructed config objects:

```cpp
template< typename T > struct DefaultConstruct { static const T value; };
template< typename T > const T DefaultConstruct< T >::value {};
```

This static is initialized **once at program startup**, before any
`DefaultPrecision()` call — when BMP's global default precision is still its
startup value (20 on this machine). It backs default construction such as
`TestedEGType my_endgame(tracker);` with no explicit config.

If a config field is `mpfr_float`, the field inside that static object is born at
precision 20 and stays there forever. With BMP **expression templates on** and
the `preserve_related_precision` policy, that stale precision propagates through
every arithmetic touch:

```
times[1] = end_time * sample_factor   // sample_factor prec=20 → times[1] prec=20
... → endtime_ prec=20 → delta_t_ prec=20 → space vector prec=20
→ System::SetVariables: vector at prec 20, system at 16 → runtime_error
```

`TrackerLoopInitialization` resets member precision before assigning, which is
the right pattern — but the contamination enters earlier, in the endgame's
`ComputeInitialSamples`, not via the tracker loop, so it isn't saved there.

## Decision

Store these rational config constants as **`mpq_rational`**, and convert to
`mpfr_float` lazily at each use site at the current precision.

- `mpq_rational` is **precision-free at rest** — no MPFR `prec` field, so it is
  safe inside `DefaultConstruct<T>::value`. It is **exact at use** — BMP converts
  the rational to `mpfr_float` at the *current* `DefaultPrecision()` exactly.
- Use the two-argument constructor `mpq_rational{1, 2}`, **never**
  `mpq_rational(1)/2` — the latter is integer division and evaluates to `0` when
  `p < q` (this latent bug was masked because `test_settings` had been comparing
  against zero).
- AMP call sites use the precision-carrying ctor `mpfr_float(field, prec)`, not
  `NumTraits<mpfr_float>::FromRational(field, prec)`.
- Classic-input config files are parsed to exact rationals via
  `decimal_str_to_rational` in
  `core/include/bertini2/io/parsing/settings_parsers/endgames.hpp` (so
  `SampleFactor: 0.647;` is exactly `647/1000`, not the nearest double). Note the
  GMP base-auto-detect trap: strip leading zeros before constructing `mpz_int`
  from a digit string (`mpz_int("08")` throws — leading `0` means octal).
- The Python bindings expose these `mpq_rational` fields as `mpfr_float`
  properties via an exact round-trip (`python/bertini/config.py`,
  `SteppingVisitor` in `python_bindings/include/tracker_export.hpp`), so the
  user-facing string/Float setter API is unchanged.

## Consequences

- **Rule:** never put `mpfr_float` (or any precision-carrying type) in a config
  struct that appears as a `DefaultConstruct<T>::value` static. Use `mpq_rational`
  (exact for all rationals, lazily converted) or, where exactness across all
  decimals does not matter, `double`.
- `mpq_rational` remains correct for the other exact-constant use it already had:
  Butcher tableau coefficients in `explicit_predictors.hpp`.
- The contamination is only observable with **expression templates on +
  `preserve_related_precision`**, which is the configured policy
  (`DefaultPrecisionPolicy()` / `thread_default_variable_precision_options` in
  `mpfr_complex.hpp`). The ET build is the default (`BMP_EXPRESSION_TEMPLATES=ON`);
  the ET-off build does not exhibit it but is not the shipping configuration.
- A `double`-typed knob (`min_step_size`) is retained where exact-decimal
  representation is unnecessary and only a threshold magnitude matters.
