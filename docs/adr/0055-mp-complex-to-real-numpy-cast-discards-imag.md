# ADR-0055: mp complex → real numpy cast discards the imaginary part (never throws)

**Status:** Accepted
**Date:** 2026-07-12

## Context

eigenpy registers numpy casts between the multiprecision dtypes and the builtin scalar
types. The `complex_mp → double` cast (registered in `mpfr_export.cpp`, so that a
`complex_mp` array can be stored into or `.astype()`-ed to a real array) went through the
generic `eigenpy::cast<complex_mp, To>` specialization in
`python_bindings/include/eigenpy_interaction.hpp`, whose body was
`return static_cast<To>(from);`.

`static_cast<double>(complex_mp)` routes through boost.multiprecision's complex→scalar
conversion, which **throws** `std::runtime_error("Could not convert imaginary number to
scalar.")` whenever the imaginary part is nonzero. numpy invokes the registered cast from an
internal C loop (`internal::cast::run`) that is not exception-safe: a C++ exception thrown
*out of* that loop unwinds through numpy's C frames and reaches an implicitly-`noexcept`
boundary, so `std::terminate()` fires and the whole interpreter dies with SIGABRT — not a
catchable Python exception.

This is trivially and naturally hit. Any real value produced by the solver carries tiny
(~1e-13) imaginary noise, and eval always returns a `complex_mp`, so the very common

```python
M = np.zeros((r, c))          # float64
M[i, j] = solver.real_solutions()[k][0]   # complex_mp with nonzero imag  → SIGABRT
```

hard-crashes. It surfaced while building a Macaulay matrix for a local-dimension test.

numpy's own builtin `complex128 → float64` cast does not throw: it **discards the imaginary
part** (with a `ComplexWarning`). A well-behaved user dtype should match that.

## Decision

The `complex_mp → To` (real target) cast converts the **real component** —
`return static_cast<To>(from.real());` — mirroring numpy's builtin complex→real cast, which
discards the imaginary part. It never throws, so it can never cross numpy's C cast loop and
abort the interpreter.

The only real-target cast registered for `complex_mp` is `complex_mp → double`
(`complex_mp → complex128` has its own specialization that preserves both components;
`complex_mp → integer` is deliberately unregistered), so this change affects exactly that
one cast.

## Consequences

- Storing / `.astype()`-ing a `complex_mp` (any imaginary part) into a real numpy array now
  behaves like numpy's builtin complex→real cast (real part kept, imaginary part dropped)
  instead of crashing.
- **Do not restore `static_cast<To>(from)` in `cast<complex_mp, To>`.** It looks cleaner and
  is symmetric with the `real_mp` specialization, but for a complex source with nonzero
  imaginary part it reintroduces the boost throw and the SIGABRT. The `.real()` conversion is
  load-bearing.
- The imaginary part is dropped **silently** (no `ComplexWarning`): the cast is a per-element
  callback with no clean once-per-operation hook, and per-element warning spam is worse than
  none. Callers who must not lose the imaginary part should cast to `complex128` (`.astype(
  complex)`), which preserves it.
- Regression coverage lives in `python/test/classes/eigenpy_numpy_test.py` (scalar setitem,
  vectorized `.astype`, broadcast assignment, parity with the numpy builtin, and the
  `complex_mp → complex128` no-regression check). This is part of the eigenpy numpy-slot
  hardening family — see ADR-0006 (uninitialized slots) and ADR-0051 (owned getitem).
