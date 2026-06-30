# ADR-0041: Random coefficients are well-conditioned — bounded-modulus scalars, conjugate-orthonormal matrices (matching Bertini 1)

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

Two mechanisms, by what is being drawn. Both come from the same principle Bertini 1 follows: a random
object should be generic *and* well-scaled, never relying on luck to avoid a wild or near-degenerate
value.

### Individual scalar coefficients — bounded modulus

A scalar coefficient (one entry of a linear form, a patch/slice constant) is drawn **bounded-modulus**:
draw a box-uniform complex `z` (real, imag each in `[-1,1]`) and divide by `sqrt(|z|)`, so the result
has modulus `sqrt(|z|)` — pulled toward 1, bounded away from both 0 and infinity, but **not** collapsed
onto the unit circle (that would be `z/|z|`). This is how Bertini 1 draws its scalars. It lives in
`core/include/bertini2/random.hpp` as `multiprecision::RandomComplexBoundedModulus` /
`…Assign`. For objects that must stay **real** — a real patch exists precisely to keep a real path real
(a complex patch would complexify it), likewise a real slice — there is a real-line analog
`RandomRealBoundedModulus` / `…Assign`: the same recipe with imaginary part 0.

### Whole random matrices — conjugate-orthonormal

Whenever an entire random *matrix* is needed (a slice's coefficient block, a randomization tail), it is
drawn **conjugate-orthonormal** (unitary rows/columns), not entry-by-entry. Bertini 1 builds every
random complex matrix this way, and when it needs a non-square shape it generates a *square* one and
truncates. We mirror that in `bertini::RandomConjugateOrthonormalMatrix(rows, cols)`
(`core/include/bertini2/eigen_extensions.hpp`): draw a square seed of the larger dimension, QR-factor it
to a unitary `Q`, return the leading `rows × cols` block. The QR launders the seed draw away, so the
result is conjugate-orthonormal regardless of how the seed was drawn — and perfectly conditioned
(condition number 1), which is the whole point.

Applied to:

- **Start systems** (scalars) — `TotalDegreeBinomial` (coefficients *and* the binomial constants),
  `TotalDegreeLinearProduct`, and `MHom`.
- **Patches** (scalars) — the complex `Patch` constructor and `Patch::RandomReal`.
- **Slices** — the constant column and the `orthogonal=false` coefficients are **bounded-modulus
  scalars**; the `orthogonal=true` coefficient matrix is **conjugate-orthonormal**
  (`RandomConjugateOrthonormalMatrix`, which also retired the old transpose-dance QR).
- **Randomization** (matrix) — the `RandomizationBlock` squaring-up matrix in `System::Randomize`
  (ADR-0025) is **conjugate-orthonormal**: the dense multi-group `R`, and the random tail `C` of the
  single-affine-group `R = [I | C]`. Only `C` is randomized; the identity block is left exact, which
  keeps the degree-optimal structure — after the descending-degree sort every tail function is lower
  degree than row `i`'s leading `f_i`, so `target_md[i]` stays `d_i` even though `C` is dense.

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

- **Slice generation, deeper pass** — a slice's constant column is drawn directly (never orthogonalized)
  and the `orthogonal=false` path uses the raw coefficients; whether these want different treatment
  (and exactly how b1 builds slices) is still to be worked out. Bertini 1's slice generation is hard to
  read; this is deferred until that's understood.
