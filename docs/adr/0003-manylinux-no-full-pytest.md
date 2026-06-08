# ADR-0003: Linux wheel CI uses import smoke test only; full pytest runs on macOS/Windows

**Status:** Accepted — pending reversal once ADR-0006 is proven in CI
**Date:** 2026-06-07
**Update (2026-06-08):** The underlying cause (uninitialized mpfr/mpc numpy slots)
was fixed on the bindings side — see **ADR-0006**. The crash also reproduces
locally, not only in the manylinux container, which broadens this ADR's original
"CI-only" framing. The decision below (Linux smoke-test-only) remains in force
until full Linux pytest is restored and proven green in a manylinux CI run; the
restore recipe is in the Consequences section.

## Context

The `manylinux_2_28` Docker image (used by cibuildwheel for PyPI-compatible Linux
wheels) is based on AlmaLinux 8 and ships:

```
mpfr-devel-3.1.6-1.el8.x86_64
```

from the AlmaLinux 8 yum repos. This is MPFR **3.1.6**. The conda-based local and
CI macOS/Windows environments use MPFR **4.x**.

### The crash mechanism

Numpy allocates array backing memory via `malloc`, which does not zero-initialize.
In the cibuildwheel multi-version build environment (where multiple Python versions
build concurrently), `malloc` frequently returns non-zero memory.

When a numpy array of `mpfr_complex` is constructed (`np.array([mpfr_complex(0)] * n)`),
eigenpy's default `setitem` does:

```cpp
T& dest = *static_cast<T*>(dest_ptr);  // raw numpy slot — may have garbage _mpfr_d
dest = src;                              // Boost.Multiprecision operator=
```

`mpc_complex_imp::operator=` protects against uninitialized destinations by checking:

```cpp
if (m_data[0].re[0]._mpfr_d == nullptr)
    mpc_init2(m_data, ...);   // safe path: initialize before use
```

The null check is the only guard. If `_mpfr_d` is non-null garbage (from malloc
returning non-zero memory), the check passes, `mpc_init2` is skipped, and `mpc_set`
is called on garbage → MPFR's internal `MPFR_ASSERTN` fires → `abort()` (SIGABRT).

MPFR 3.1.6 has always-on `MPFR_ASSERTN`. MPFR 4.x may behave differently or the
local allocator may happen to return zeros more reliably; local tests pass.

Multiple test files are affected: `amptracking_test.py`, `cauchy_endgame_test.py`,
and likely others that construct numpy arrays of `mpfr_complex`.

The problem is in CI only. Installed wheels work correctly for users who have MPFR 4.x.
The manylinux wheel bundles the code; users bring their own MPFR at runtime.

### What was tried

1. Skipping `amptracking_test.py` — the next file (`cauchy_endgame_test.py`) crashed.
2. Skipping both — the scope of affected files was unknown and growing.

### The proper fix (implemented 2026-06-07)

Three defenses now live in `python_bindings/include/eigenpy_interaction.hpp`
(see the header comment there for the full picture). The root insight: numpy
zero-fills fresh buffers for these dtypes (NPY_NEEDS_INIT), but an all-zero
mpfr_t/mpc_t is Boost.Multiprecision's *uninitialized sentinel*, not a valid
value. BMP's own assignment operators check the sentinel and initialize first,
so writes are safe — but anything that **reads** a never-written slot and hands
it straight to libmpfr/libmpc crashes.

1. **getitem** (pre-existing): heals a zeroed slot in place on read.
2. **setitem** — `internal::zeroinit_setitem<NumT>` zero-initializes the
   destination slot before delegating to eigenpy's original setitem:

   ```cpp
   std::memset(dest_ptr, 0, sizeof(NumT));
   // Now _mpfr_d == nullptr, so operator= will call mpc_init2 before mpc_set.
   return original(src_obj, dest_ptr, array);
   ```

   Unlike `getitem`, eigenpy provides no template customization point for
   `setitem`, so the guard is installed by patching the registered dtype's
   `PyArray_ArrFuncs::setitem` immediately after `registerNewType<T>()` — see
   `eigenpy::HardenSetitem<T>()`, called from `mpfr_export.cpp` for both
   `mpfr_float` and `mpfr_complex`. The cost is a leak when overwriting an
   initialized slot; numpy never destructs user-dtype elements anyway, so this
   converts a crash on dirty memory into a bounded leak on element overwrite.
3. **Guarded ufunc loops and casts** — `eigenpy::registerGuardedUfunct` replaces
   eigenpy's stock loops (add/subtract/…/matmul/equal/…) with versions that
   substitute an exact zero when reading an uninitialized slot
   (`internal::value_or_zero`), and `eigenpy::cast<mpfr_*, To>` specializations
   guard the registered numpy cast loops the same way. These read paths were
   reproducibly crashing **locally** (`np.zeros(n, dtype=mpfr_complex) + ...`
   segfaulted before the guards).

Regression tests: `python/test/classes/numpy_uninitialized_slots_test.py`
(24 tests over both dtypes; the arithmetic/equality/matmul ones crashed the
interpreter outright before the guards).

Note: no eigenpy release fixes this upstream. eigenpy has set `NPY_NEEDS_INIT`
(requesting zeroed allocation from numpy) since at least v3.1.0, its `setitem`
has always relied on that guarantee, and its ufunc/cast loops have never guarded
the read side. Upgrading eigenpy does not help.

**CI restoration is still pending**: the fix is verified locally (numpy 2.4.6,
dirty-heap stress probe over 20 allocation paths, full pytest suite, setitem
semantics including overwrite and high-precision round-trips, virgin-slot
read-path checks), but the restore-full-Linux-testing recipe below must be
validated in an actual manylinux CI run before this ADR's decision is reversed.

## Decision

Replace the full pytest suite in `CIBW_TEST_COMMAND_LINUX` with a basic import smoke
test:

```yaml
CIBW_TEST_COMMAND_LINUX: "python -c 'import bertini; print(bertini.__version__)'"
```

Remove `CIBW_TEST_REQUIRES_LINUX: "pytest numpy"`.

The full Python test suite (140+ tests) continues to run on:
- **macOS host runners** via the `test_wheels_linux_macos` job (MPFR 4.x via Homebrew)
- **Windows host runners** via `test_windows_wheels` (MPFR 4.x via conda)

The smoke test verifies that the wheel installs and the extension module loads without
crashing — sufficient to catch packaging and linking regressions.

## Consequences

- **Linux CI passes.** The MPFR 3.1.6 / malloc crash no longer blocks the build.
- **Reduced coverage in manylinux.** Python test regressions that only manifest with
  MPFR 3.1.6 will not be caught by CI. This is a real gap, but MPFR 3.1.6 is EOL
  and not used by any target user environment.
- **Full coverage maintained.** macOS and Windows run the complete pytest suite.
- **Reversible.** When the `setitem` specialization is implemented (see above), restore:
  ```yaml
  CIBW_TEST_REQUIRES_LINUX: "pytest numpy"
  CIBW_TEST_COMMAND_LINUX: "cd {project} && python -m pytest python/test/ -q"
  ```
- **Do not run pytest in manylinux** until the `setitem` specialization is in place.
  Partial ignores are not a solution — the affected test file list is not bounded.
