# ADR-0052: The ExplicitRKPredictor Butcher tables live inline in the header

**Status:** Accepted
**Date:** 2026-07-10

## Context

`ExplicitRKPredictor` (`core/include/bertini2/trackers/explicit_predictors.hpp`) carries the
Runge–Kutta **Butcher tableaux** (Euler, Heun-Euler, RK4, RKF45, Cash-Karp45, Dormand-Prince56,
Verner67) as exact `mpq_rational` constants — 52 members total: 26 `mpq_rational[]` coefficient
arrays and the 26 `Eigen::Matrix<mpq_rational,…>` tables built from them. `mpq_rational` is a
non-literal type, so these can never be `constexpr`; they require dynamic initialization.

Historically they were **out-of-line static-const members**, defined in a *data-only*
`explicit_predictors.cpp`. That triggered a Windows-only link failure (issue #287):
`CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS` exports **functions but not data symbols**, and, more
concretely, `lld-link` will not pull a *data-only* object out of a static archive unless
something in it is referenced. So `explicit_predictors.obj`'s Butcher data went missing at link
time, and `core/CMakeLists.txt` compensated by forcing
`$<LINK_LIBRARY:WHOLE_ARCHIVE,bertini2>` + `/FORCE:MULTIPLE` onto **all nine test executables** —
dragging the entire library into each one and inflating the Windows build/link (a large part of
that job's wall-clock, the CI long pole).

## Decision

Define the Butcher tables as **C++17 `inline` static-const members, in the header**. Each using
translation unit gets a COMDAT copy (folded at link), so there is no out-of-line data definition
to miss and no export/archive-pull dependency on any single object.

Consequently:

- **`explicit_predictors.cpp` is retained, but holds ONLY the explicit template instantiations**
  (`Predict`, `FullStep`, `EvalRHS`, `SetNormsCond`, `SetErrorEstimate`, `SetSizeProportion`,
  `FillButcherTable`) that pair with the header's `extern template` declarations — the
  closed-universe pattern of **ADR-0014**. It is *not* deleted, and it is *not* the data TU.
- The Windows `/WHOLEARCHIVE` + `/FORCE:MULTIPLE` block in `core/CMakeLists.txt` is **removed**.
  The `/FORCE:MULTIPLE` existed only to tolerate the duplicate Boost.Serialization
  (`void_cast_register` / `oserializer`) strong symbols that whole-archive exposed from multiple
  TUs; with normal archive pull + COMDAT folding restored, there are no duplicates to force past.

This is deliberate placement, chosen over the alternatives of `__declspec(dllexport)` on the data
(fragile with `WINDOWS_EXPORT_ALL_SYMBOLS`, and wrong for the static-lib case) or out-of-line
accessor *functions* (works — functions export — but leaves the data in a `.cpp`; the inline
form removes the platform split entirely). The header cost is acceptable: it is included only via
`ode_predictors.hpp` → `base_tracker.hpp`, the added text is 52 short initializers (parse-only
unless ODR-used; codegen is COMDAT-deduped and setup-time, never a hot loop), and **these numbers
never change**, so the usual "definitions in headers force rebuilds" objection does not apply.

## Consequences

- **Do not move the Butcher tables back into a `.cpp`.** Seeing large `mpq_rational` tables in a
  header and relocating them "for cleanliness" reintroduces the exact Windows link failure this
  ADR retired. If you must, you also owe back the `/WHOLEARCHIVE` workaround.
- The tables are now lazily-safe COMDAT statics instead of eager namespace-scope globals,
  removing a latent static-initialization-order hazard.
- In-class initialization of a class-type member must use `= T(ptr)` (not `T name(ptr);`, which
  parses as a member-function declaration, nor `{…}`, which risks an `initializer_list` match).
  Each `…Ptr_` array is declared before its matrix so the matrix initializer, which reads the
  array, is sequenced after it within any TU.
- Values are unchanged: the migration preserved all 363 `mpq_rational(n,d)` tokens byte-for-byte
  and in order; `test_tracking_basics` and `test_endgames` pass. The Windows-specific claim (test
  exes link with no whole-archive) is proven by CI, not locally (the flag was MSVC-only).
- Relates to **ADR-0014** (the `.cpp` remains the single instantiating TU) and **ADR-0007** (the
  Butcher coefficients stay exact `mpq_rational`).
