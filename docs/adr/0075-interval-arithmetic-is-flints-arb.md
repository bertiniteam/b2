# ADR-0075: Interval arithmetic is FLINT's Arb, through a thin wrapper, with a Boost.Interval rectangle for double polynomials

**Status:** Accepted
**Date:** 2026-11-03

## Context

A point's accuracy is a scalar carried beside it, produced by whichever step last touched it, and nothing propagates it; measured on the Whitney umbrella against exact ground truth, the claimed error bars were wrong by 19 to 30 orders of magnitude (#412).  Every downstream same-point decision consumes the bar, never the value.  A rigorous enclosure removes the whole category: two points are the same, or provably not, with no constant to tune.  The concrete need is Krawczyk certification of points in the cellular decomposition work, including singular points certified through their deflated systems, which needs multiprecision.

Interval evaluation of a system and its Jacobian needs an interval scalar type that the SLP can evaluate in.  b2's scalar universe was closed at `{complex_dbl, complex_mp}` (ADR-0014).  The candidates for the multiprecision type:

- **Our own interval type over MPFR's directed rounding.**  Rejected outright by silviana: she wrote her own complex type early, and ripped it out the moment Boost.Multiprecision gained complex numbers; the rounding and the transcendental enclosures would be ours to get right, forever.
- **MPFI through Boost.Multiprecision's `mpfi_float` and `complex_adaptor`.**  Boost 1.88's `complex_adaptor` chooses formulas by the sign of the real or imaginary part (`eval_sqrt`, `eval_tan`, and `eval_abs` computes `re*re + im*im`, which for an interval straddling zero is not an enclosure of a square), so for a box across a branch cut it returns a box that does not enclose the image: wrong, not loose.  Fixing that is writing complex interval arithmetic ourselves.  MPFI also has no Windows build on conda-forge.
- **C-XSC.**  A complex double interval type with rigorous functions, but among the slowest interval libraries in published comparisons, with no conda-forge package.
- **FLINT 3, which carries Arb** (`arb_t` real balls, `acb_t` complex rectangles of balls).  Rigorous enclosures for every operation the SLP emits, branch cuts included (a box across a cut comes back wide or non-finite, never wrong); on conda-forge for all five wheel platforms, in Homebrew and in apt; Johansson's measurements put a real ball multiply at about one MPFR multiply, and a complex ball multiply at a quarter to a half of MPC's point multiply, since MPC pays for correct rounding.  There is no maintained C++ wrapper (arbxx was archived in 2024 when Arb merged into FLINT) and no Boost.Multiprecision backend.

The spike (draft PR #516, results in the b2_fun notes of 2026-11-01) built FLINT into every CI route and measured: on five platforms, `acb` at equal precision costs 0.23 to 0.83 of b2's non-rigorous `complex_mp`, mostly 0.3 to 0.6, so interval evaluation in mp is faster than today's point evaluation in mp.  A complex rectangle of `Boost.Interval<double>` with the rounding mode set once per evaluation costs 7 to 10 times `complex_dbl` on Linux and macOS Intel, 19 on Windows and 30 on Apple Silicon, and 4 to 6 times less than `acb` at 53 bits; with Boost's default per-operation rounding switch it costs 55 to 591 times `complex_dbl`.

## Decision

Two interval types, mirroring b2's dbl/mp pair.

- `real_iv_mp` / `complex_iv_mp` wrap FLINT's `arb_t` / `acb_t` in a thin C++ class of ours that manages lifetime and forwards every operation to FLINT.  No arithmetic, rounding or transcendental enclosure is ours.  Any precision, every SLP operation.
- `real_iv_dbl` / `complex_iv_dbl`: `boost::numeric::interval<double>` for the real part, and a complex rectangle of two of them by the textbook formulas, every step a Boost real-interval operation.  Polynomial arithmetic only: the transcendental functions are not defined (Boost's rely on libm, which does not promise directed rounding), and an SLP containing one refuses double-interval evaluation and points to the mp type.  The rounding mode is set once around a whole evaluation, never per operation.

FLINT is therefore a dependency of the library: `FindFLINT.cmake`, linked PUBLIC into `bertini2`; the Linux deps image builds GMP (6.3.0, fat on x86_64; FLINT needs >= 6.2.1 and AlmaLinux 9 ships 6.2.0), MPFR, MPC and FLINT from source into `/usr/local`, stripped, with the versions pinned in `.github/ci-deps-versions.env` and `-flint<version>` in the image tag; macOS takes Homebrew's `flint`, Windows conda-forge's `libflint`.  `FlintVersion()` and its export `info.flint_version` are the reference that keeps FLINT linked into the Python extension (a static `libbertini2` drops what nothing references, and the wheel repair tools then bundle nothing).  `test_flint` asserts, on every platform, the facts b2 relies on: enclosure of products and quotients, a non-finite ball from division by a ball around zero, exact `mpfr` <-> `arf` conversion, one GMP shared with b2, a zero-filled `acb` being an exact zero, and the Boost rectangle enclosing Arb's evaluation.

The work lives on `integration/intervals` (ADR-0074) until it has earned its merge; before it reaches `develop`, FLINT becomes a build option, default off.

## Consequences

- Do not write interval rounding or transcendental enclosures of our own; do not adopt MPFI through `complex_adaptor`; do not use `std::complex<interval<double>>` (unspecified by the standard, and libc++'s Annex G code branches on sign).
- Do not switch the FPU rounding mode per operation; one `save_state` scope around an SLP evaluation.
- Do not give the double interval type transcendental functions; the mp type is the answer for those.
- Do not run Eigen's pivoting decompositions on interval matrices (pivot choice compares magnitudes, meaningless for intervals); Krawczyk needs only products, with the approximate inverse computed in point arithmetic.
- Do not drop `info.flint_version` or the `FlintVersion()` reference: the wheels would silently stop carrying FLINT on Windows and Linux.
- Do not ship an unstripped `libflint.so` in the image: it adds about 25 MB to every Linux wheel.
- FLINT declares its `mpz`/`mpfr` conversions only when `<gmp.h>`/`<mpfr.h>` are included before its headers; include them first.
- The 4.0.0 precedent for the numpy dtypes still holds: `slot_write` for every output of the mp types (mantissas above 128 bits live on the heap, so numpy's memcpy of element 0 aliases them); a zero-filled slot needs no guard for a FLINT-backed type.
