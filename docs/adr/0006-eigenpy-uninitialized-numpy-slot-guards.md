# ADR-0006: Guard numpy mpfr/mpc dtypes against uninitialized slots

**Status:** Accepted
**Date:** 2026-06-08

## Context

`mpfr_float` and `mpfr_complex` are registered as custom numpy dtypes through
eigenpy (`registerNewType<T>()` in `python_bindings/src/mpfr_export.cpp`).
eigenpy registers them with `NPY_NEEDS_INIT`, which asks numpy to zero-fill
freshly allocated array buffers for the dtype.

The trap: **an all-zero `mpfr_t`/`mpc_t` is not a valid value — it is
Boost.Multiprecision's "uninitialized" sentinel** (`_mpfr_d == nullptr`). BMP's
own assignment operators check that sentinel and initialize the limb storage
before assigning, so *writes* into a fresh, zeroed slot are safe. But any path
that *reads* a never-written slot and hands the raw `mpfr_t`/`mpc_t` to
libmpfr/libmpc operates on a null limb pointer and crashes inside the C library.

Two distinct failure modes were observed:

1. **Read of an uninitialized slot → SIGSEGV.** eigenpy's stock ufunc loops
   (`add`, `subtract`, `multiply`, `divide`, comparisons, `matmul`, `negative`,
   `square`, `sqrt`) and its registered cast loops read input slots directly.
   `np.zeros(n, dtype=mpfr_complex) + ...`, `A @ A`, and `arr.astype(...)` on
   never-written slots segfaulted. **This reproduces locally**, not only in CI —
   it is not specific to any container or MPFR version.

2. **Write onto malloc-dirty memory → SIGABRT.** When the `NPY_NEEDS_INIT`
   zero-fill guarantee is violated (observed in the `manylinux_2_28` build
   container, see ADR-0003), a slot can contain garbage with a non-null
   `_mpfr_d`. That defeats BMP's null check on the write path: `mpc_set` runs on
   garbage limbs and MPFR's `MPFR_ASSERTN` calls `abort()`.

eigenpy provides a customization point for `getitem` (a template specialization)
but **none for `setitem` or for its ufunc/cast loop bodies** — those are fixed
implementations baked into `SpecialMethods<T, NPY_USERDEF>` and the
`EIGENPY_REGISTER_*_UFUNC` macros.

## Decision

Install three defenses in `python_bindings/include/eigenpy_interaction.hpp`,
all keyed on a single `internal::mpfr_slot<T>::uninitialized()` trait that
detects the zeroed-`_mpfr_d` sentinel:

1. **getitem heals on read.** The existing `getitem<mpfr_float>` /
   `getitem<mpfr_complex>` specializations replace an uninitialized slot with
   `T(0)` in place before handing it to Python.

2. **setitem zero-inits before assign.** `internal::zeroinit_setitem<NumT>`
   `memset`s the destination slot to zero, then delegates to eigenpy's original
   `setitem`. Zeroing forces `_mpfr_d == nullptr`, so BMP's `operator=` takes its
   init-before-set path regardless of what the allocator delivered. Because
   eigenpy has no setitem customization point, the guard is installed by
   **patching the registered dtype's `PyArray_ArrFuncs::setitem` function pointer
   immediately after `registerNewType<T>()`** — `eigenpy::HardenSetitem<T>()`,
   called for both `mpfr_float` and `mpfr_complex` in `mpfr_export.cpp`.

3. **Guarded ufunc and cast loops.** `eigenpy::registerGuardedUfunct<Scalar,
   WithOrderingComparitors>()` registers replacement loops that read each input
   through `internal::value_or_zero` (substituting an exact `T(0)` for an
   uninitialized slot), and `cast<mpfr_float, To>` / `cast<mpfr_complex, To>`
   specializations do the same on the cast path. The ordering comparitors are a
   **compile-time** parameter (`if constexpr`), not a runtime bool, because
   instantiating `op_greater` etc. for `mpfr_complex` is a hard compile error
   (no ordering on complex).

Regression tests live in
`python/test/classes/numpy_uninitialized_slots_test.py`.

## Consequences

- **No customization-point reliance for setitem.** The function-pointer patch is
  applied once at registration; it is robust to numpy 1.x/2.x (the
  `PyDataType_GetArrFuncs` accessor is shimmed by eigenpy for both ABIs).

- **A bounded leak on overwrite.** `zeroinit_setitem` `memset`s without first
  freeing any mpfr allocation the slot already held, so overwriting an
  initialized element leaks its limbs. numpy never destructs user-dtype elements
  anyway (every discarded array of these dtypes already leaks), so this trades a
  crash-on-dirty-memory for a small, bounded leak. Acceptable.

- **Upgrading eigenpy does not fix this.** eigenpy has set `NPY_NEEDS_INIT`,
  relied on it in `setitem`, and shipped unguarded ufunc/cast loops continuously
  from at least v3.1.0 through v3.13.0 (changelogs checked). The guards must live
  on our side.

- **Resolves the ADR-0003 blocker.** This is the fix ADR-0003 was waiting for
  before full Linux pytest can be restored in CI. See ADR-0003.

- **Related but distinct from ADR-0001.** ADR-0001 is about eigenpy's *writable
  `Eigen::Ref` + adjacent `const&` scalar* converter corruption. This ADR is
  about *uninitialized dtype slots*. Both are eigenpy-interaction hazards with
  workarounds on our side, not upstream.
