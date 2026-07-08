# ADR-0051: Full numpy ufunc coverage for the mp dtypes; getitem returns owned copies

**Status:** Accepted
**Date:** 2026-07-08

## Context

The eigenpy-registered numpy dtypes for `real_mp` / `complex_mp` covered only the
arithmetic core (add/subtract/multiply/divide, equality + real orderings,
negative/square/sqrt, matmul, and the hardened `dotfunc`).  Everything else —
`np.abs`, `np.conj`, the transcendental family, `power`, `sign`, `minimum/maximum`,
rounding, the `isnan` predicates — raised `ufunc ... not supported`; `np.sort` /
`np.argmax` failed on empty dtype slots ("type does not have compare function" /
"data type not ordered"); and the docs declared the identity-seeded reductions
(`np.sum`/`np.prod`/`np.mean`) permanently unsupported after a build-dependent
`SystemError` (documented 2026-06-29, "not something Bertini can patch").

Investigating the reductions on current numpy (2.3.2 and 2.4.6) showed the
`SystemError` no longer reproduces — but exposed something worse hiding behind them:

**eigenpy's `getitem` returns `boost::ref(slot)` — a Python scalar aliasing numpy
array storage.**  Our heal-on-read `getitem` specializations (ADR-0006) had copied
that behavior faithfully.  Consequences:

- A scalar extracted from a *temporary* array dangles once the array is freed.
  `s = np.sum(v)` is exactly that: the reduce result is a temporary 0-d array,
  `s` aliased its buffer, and reading `s` later gave zeros, garbage precision, or an
  MPFR assertion SIGABRT inside `str()`.
- `np.mean` returned a **silently wrong `0`** through the same mechanism.
- This is the root of the ADR-0031 / #259 hazard class ("indexing an eigenpy Vec
  returns an aliasing view"), previously mitigated consumer-by-consumer with a
  copy-at-extraction rule in Python code.

## Decision

1. **Register guarded loops for the full ufunc set** on both dtypes
   (`python_bindings/include/eigenpy_interaction.hpp`, `registerGuardedUfunct`):

   - both dtypes: `absolute` (complex → `real_mp` output), `conjugate`, `sign`
     (numpy-2 semantics: complex sign is `z/|z|`), `positive`, `reciprocal`,
     `power`, `exp`, `log`, `log10`, full trig/hyperbolic + inverses,
     `isnan`/`isinf`/`isfinite` (→ bool);
   - real only (ordering- or domain-dependent): `greater`/`less`/... , `fabs`,
     `exp2`, `log2`, `expm1`, `log1p`, `cbrt`, `floor`, `ceil`, `trunc`, `rint`,
     `signbit`, `arctan2`, `hypot`, `copysign`, `fmod`, `remainder`,
     `floor_divide`, `minimum`, `maximum`, `fmin`, `fmax`.

   Every loop reads through `value_or_zero` (the ADR-0006 uninitialized-slot
   doctrine) and calls the same boost::multiprecision free function the
   `bertini.multiprec` scalar function binds, so `np.f(a)[i] == mp.f(a[i])`
   exactly.  Deliberate semantic choices, matching numpy's float64 behavior:
   `rint` rounds half-to-even (direct `mpfr_rint` in `MPFR_RNDN`; boost's `rint`
   rounds half away), `remainder`/`mod` takes the sign of the divisor (`fmod`
   keeps C semantics), `minimum`/`maximum` propagate nan while `fmin`/`fmax`
   ignore it.

2. **Fill the `compare`/`argmax`/`argmin` dtype slots for `real_mp`**
   (`HardenCompare`, `HardenArgMinMax`, same install pattern as `HardenDotfunc`).
   Enables `np.sort`/`argsort`/`searchsorted`/`unique`/`median`/`argmax`/`argmin`.
   Complex stays unordered on purpose — numpy's lexicographic complex ordering is
   historical baggage we do not reproduce.

3. **`getitem` returns an owned copy, never `boost::ref`.**  An indexed element is
   a durable value, as numpy users expect.  This kills the ADR-0031/#259 hazard
   class at the source (that ADR's copy-at-extraction rule in Python remains good
   hygiene but is no longer load-bearing), and it is what makes the reductions
   *actually* safe rather than accidentally readable.  Cost: one mp copy per
   element read.

4. **Reductions are supported and regression-tested**, on numpy ≥ 2.3 (verified
   2.3.2 and 2.4.6; `python/test/classes/numpy_ufuncs_test.py::TestReductions`
   pins them in CI on all three platforms).  The docs note the historical
   `SystemError` and keep the `initial=` idiom as the fallback for older numpy.
   `pyproject.toml` keeps `numpy` unpinned.

5. **The float64 boundary stays closed for VALUES, open for tolerance
   comparisons** (both decided 2026-07-08).  `double → mp` casts remain
   registered *unsafe*, so float64 scalars/arrays do not silently promote into
   mp arrays: the conversion itself is bit-exact, but a promoted float64 `0.1`
   is not the decimal `0.1` the user typed — the user has to think.  However,
   mixed mp-vs-float64 **ordering** loops (`<`, `<=`, `>`, `>=`, both operand
   orders, real only) ARE registered: `np.abs(a - b) < 1e-10` is safe — the
   result is a bool, no float flows into an mp value, the comparison is exact
   (boost compares the number against the double directly), and it matches the
   C++ solvers' double `ToleranceT` and the scalar `GreatLessVisitor<T,double>`
   precedent.  Mixed EQUALITY stays unregistered (exact equality against a
   float literal is the 0.1-intent trap; not bound at scalar level either), as
   does mixed arithmetic; `np.isclose`/`np.allclose` still raise (they *compute*
   with float64 tolerances internally).  `mp → complex128` casts are registered
   unsafe alongside the pre-existing `mp → double`, so `arr.astype(complex)` /
   `astype(float)` are the explicit, conscious truncations.  Registration-order
   note: casts must be registered BEFORE the ufunc loops — registering the
   mixed loops makes numpy query the mp↔double casts, and a cast first
   registered after being queried is permanently ignored (numpy
   RuntimeWarning).

6. **Component access on complex arrays** goes through new array overloads of
   `multiprec.real`/`imag`/`arg` (returning `real_mp` arrays).  numpy cannot know
   a legacy user dtype is complex-like, so the ndarray `.real`/`.imag` attributes
   (and `np.real`/`np.imag`/`np.angle`) return silently wrong values — `.real`
   gives the complex values, `.imag` gives zeros.  Not hookable from the bindings;
   documented, plus a pinning test.  (En route, fixed a copy-paste bug: the scalar
   `mp.imag` was bound to `boost::multiprecision::real` and returned the real
   part.)

## Consequences

- The "Known gotchas" docs page shrinks to the real, permanent edges: the float64
  boundary (by design), complex component access (numpy limitation), complex
  ordering (by design), no mp→int casts.  The reductions section becomes a
  historical note.
- `python/test/classes/numpy_ufuncs_test.py` pins: element-wise agreement with the
  scalar functions, the numpy-semantics corners above, precision preservation
  through every loop shape (including `sign`/`reciprocal`, which cross the mixed
  real/complex division path with the known boost precision-mis-tagging hazard —
  loops re-tag via `at_precision_of`), unwritten-slot safety per loop shape,
  sorting/arg-extrema (nan-wins semantics), reductions, and the dangling-scalar
  regressions.
- Not upstreamed to eigenpy (the guarded loops already diverge; see ADR-0006).
  Upstreaming the owned-copy `getitem` would fix the aliasing class for all
  eigenpy user types and may be worth an issue later.

## Relation to prior ADRs

- **ADR-0006** — the slot-guard doctrine these loops follow; its "known gaps" list
  shrinks (compare/argmax slots now filled and guarded).
- **ADR-0031 / #259** — root-caused and fixed at the binding level by the owned-copy
  `getitem`; the Python-side copy rule is now belt-and-suspenders.
- **ADR-0001/0008** — unrelated eigenpy hazards, unchanged.
