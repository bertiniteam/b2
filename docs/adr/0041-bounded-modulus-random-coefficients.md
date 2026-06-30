# ADR-0041: Random linear-form / scaling coefficients use a bounded-modulus draw (matching Bertini 1)

**Status:** Accepted
**Date:** 2026-06-30

## Context

Numerical homotopy continuation is sprinkled with *random* coefficients: the linear forms of the
**start systems** (total-degree binomial + linear-product, and multihomogeneous), the normalization
equations of projective **patches**, and the linear forms of **slices**. The numbers themselves are
arbitrary — genericity is all that is required — but their *scale* is not free: a wildly large, wildly
small, or near-zero coefficient badly conditions the object it parameterizes (a start-system linear
form, the n×n start-point solve and start Jacobian, a patch equation, a slice's forms).

The draws in use were poorly controlled:

- **Start systems** drew each coefficient as `complex_mp(real_mp(RandomRat()), real_mp(RandomRat()))`.
  `RandomRat()` is a ratio of two independent ~`[-10^50, 10^50]` integers, i.e. a Cauchy-like,
  **heavy-tailed** distribution: median modulus ~1, but fat log-tails that occasionally emit a `~10^20`
  or `~10^-20` coefficient.
- **Patches and slices** drew with `RandomComplexAssign` / `RandomRealAssign` — box-uniform
  (complex: real and imaginary parts each in `[0,1]`, so first-quadrant-biased *and* free to sit near
  0; real: `[-1,1]`, same near-0 risk).

On cyclic-5 the heavy-tailed total-degree coefficients drove the Cauchy endgame near `t = 0` into a
double-precision corrector roundoff floor, and the real linear-product start **hung**. This was
initially misdiagnosed as an *endgame stall* — even as a conditioning mystery, since the linear-product
start "should" be better-conditioned than the binomial one. The real cause was the coefficient scale,
not the endgame and not runaway precision.

## Decision

**All random linear-form / scaling coefficients are drawn with a bounded-modulus draw.** Draw a
box-uniform complex `z` (real, imag each in `[-1,1]`) and divide by `sqrt(|z|)`, so the result has
modulus `sqrt(|z|)`: pulled toward 1, bounded away from both 0 and infinity, but **not** collapsed onto
the unit circle (that would be `z/|z|`). This is exactly how **Bertini 1** generates its coefficients
(confirmed from b1 source). It lives in `core/include/bertini2/random.hpp` as
`multiprecision::RandomComplexBoundedModulus` / `RandomComplexBoundedModulusAssign`.

For objects that must stay **real** — a real patch exists precisely to keep a real path real (a complex
patch would complexify it), and likewise a real slice — there is a real-line analog
`RandomRealBoundedModulus` / `…Assign`: the same recipe with imaginary part 0.

Applied to:

- **Start systems** — `TotalDegreeBinomial` (coefficients *and* the binomial constants),
  `TotalDegreeLinearProduct`, and `MHom`.
- **Patches** — the complex `Patch` constructor and `Patch::RandomReal`.
- **Slices** — `Slice::RandomComplex` and `Slice::RandomReal`.

## Consequences

- cyclic-5 real linear-product: **hang → completes**, 70 finite, finite-solution endgames fast
  (~3.5 ms median, on par with the binomial start). The total-degree "stall" disappears.
- The change is **seed-affecting** (coefficient sequences differ), but the test suite is
  count/behavior-based, so it stays green (C++ ctest 10/10, pytest 543 passed/2 skipped).
- The draw compresses the modulus toward 1 but does **not** strictly bound it away from 0
  (`sqrt(|z|) → 0` as `|z| → 0`); the redraw-on-exactly-0 guard only removes the measure-zero
  degenerate case. This matches Bertini 1 and is sufficient in practice; a strict annulus draw was not
  adopted, as it would diverge from b1.

## Still open (follow-ups, not done here)

- **`System::Randomize`** — the `RandomizationBlock` squaring-up matrix (ADR-0025) still draws its
  entries with `multiprecision::RandomComplex` (box-uniform) at `core/src/system/system.cpp:1187,1197`.
  It is the same class of draw and should adopt the bounded-modulus version.
- **Slice generation, deeper pass** — a slice's constant column is drawn directly (never orthogonalized)
  and the `orthogonal=false` path uses the raw coefficients; whether these want different treatment
  (and exactly how b1 builds slices) is still to be worked out.
