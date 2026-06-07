# ADR-0001: Pass scalar mpc_complex args by value when adjacent to writable Eigen::Ref

**Status:** Accepted  
**Date:** 2026-06-04

## Context

`track_path_wrap` in `python_bindings/include/tracker_export.hpp` is the Boost.Python
wrapper for `AMPTracker::TrackPath`. Its signature originally was:

```cpp
SuccessCode track_path_wrap(TrackerT const& self,
    Eigen::Ref<Vec<ComplexT>> result,
    ComplexT const& start_time,    // ← const ref
    ComplexT const& end_time,      // ← const ref
    Vec<ComplexT> const& start_point)
```

All four AMP tracker Python tests crashed with an MPFR assertion in `set_prec.c:32`
immediately on entry to the wrapper — before any tracker logic ran. GDB showed that
`start_time` already had a garbage `_mpfr_prec` and an invalid limb pointer at wrapper
entry.

The cause is an eigenpy implementation detail: the from-Python converter for a
writable `Eigen::Ref<Vec<mpc_complex>>` writes into a static slot in Boost.Python's
rvalue-converter infrastructure. When the writable Ref is converted alongside adjacent
`const&` scalar arguments, the static-slot write clobbers the rvalue-converter storage
holding those scalars. The result is a garbage `mpc_complex` — reading it triggers
`mpc_set_prec` with an out-of-range precision → `MPFR_ASSERTN` → `abort()`.

This only affects bindings that combine a **writable** `Eigen::Ref` with adjacent
**`const&`** scalar `mpc_complex` arguments. Read-only `Vec<T> const&` and single-
argument bindings are unaffected (confirmed by targeted experiments).

## Decision

Pass the scalar time arguments **by value**, not by `const&`:

```cpp
SuccessCode track_path_wrap(TrackerT const& self,
    Eigen::Ref<Vec<ComplexT>> result,
    ComplexT start_time,     // by value
    ComplexT end_time,       // by value
    Vec<ComplexT> const& start_point)
```

By-value arguments are copied before eigenpy's Ref converter runs. They are not stored
in the rvalue-converter slot, so the static-slot write cannot corrupt them.

The writable `Eigen::Ref<Vec<ComplexT>> result` is kept — the caller's numpy array is
written back to in-place, which is the intended API behavior.

The fix applies to `track_path_wrap` only. Other bindings that use `Eigen::Ref` are
single-argument (no adjacent scalar) and are unaffected.

## Consequences

- **Safer:** Eliminates the class of corruption crashes caused by eigenpy's static-slot
  converter interacting with adjacent `const&` scalars.
- **Minor copy overhead:** Each `track_path` call copies two `mpc_complex` values. For
  these multi-second path tracking calls, this is unmeasurable.
- **Rule to follow:** Any future binding that places a writable `Eigen::Ref<Vec<mpc_complex>>`
  alongside scalar `mpc_complex` arguments must pass the scalars by value, not by `const&`.
  The root cause is in eigenpy's static-slot converter; it is not fixed upstream.
- **Endgame bindings:** The same pattern was present in `endgame_export.hpp`. It was
  resolved differently — the endgame API was refactored so `run()` takes only the start
  point vector (no adjacent scalar). See ADR-0002.
