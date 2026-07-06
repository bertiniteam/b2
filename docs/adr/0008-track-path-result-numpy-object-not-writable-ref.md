# ADR-0008: Take writable-Ref binding outputs as numpy objects (supersedes ADR-0001)

**Status:** Accepted
**Date:** 2026-06-08
**Supersedes:** ADR-0001

## Context

ADR-0001 found that eigenpy's from-Python converter for a **writable
`Eigen::Ref<Vec<mpc_complex>>`** overruns its rvalue-converter storage and clobbers an
*adjacent* `mpc_complex` argument, leaving it as garbage (invalid limb pointer /
precision 0) → `mpc_set_prec` on garbage → `MPFR_ASSERTN` → `abort()`. ADR-0001's
mitigation was to pass the adjacent scalar time arguments **by value** instead of by
`const&`.

That mitigation proved **insufficient on the x86_64 manylinux build**. The full Linux
wheel pytest SIGABRTed deterministically in `amptracking_test.py::test_tracker_linear`
→ `track_path`. Instrumenting `track_path_wrap` (env-gated `BERTINI_DIAG` prints) on a
scoped CI run showed the smoking gun on x86_64:

```
[DIAG] start_time.precision=30 end_time.precision=0   ← end_time corrupted, by-value
```

`start_time` (first by-value scalar) and `start_point` were fine; `end_time` (second
by-value scalar) arrived with **precision 0**, driving `mpc_set_prec(0)` inside
`AMPTracker::TrackerLoopInitialization`. So the by-value copy is *still* taken from
storage the writable-Ref converter has corrupted. The corruption is ABI/layout-sensitive:
it fires on x86_64 manylinux (gcc-toolset-14) but **never reproduced on aarch64** under
any condition (native, `MALLOC_PERTURB_`, valgrind, Eigen 3.4.0 *and* 5.0.1, eigenpy
3.13.0 from source). It is still present in eigenpy 3.13.0 (the newest pin) — the
converter has no Eigen-version-specific handling; the layout assumption in its
`rvalue_from_python_data`/`referent_storage` specialization simply doesn't hold for this
scalar/ABI combination.

A second attempt — taking the scalar times as `boost::python::object` and extracting them
*inside* the wrapper — made it **worse**: the clobber relocated onto the tracker handle
`self` (SIGSEGV in the first `self.GetSystem()` call). Removing the scalar converters
just moves the victim. **The writable `Eigen::Ref` itself is the corruptor**, not the
scalars.

## Decision

**Remove the writable `Eigen::Ref` from the binding entirely.** `track_path_wrap` now
takes the output vector as a `boost::python::object` (the caller's numpy array), tracks
into a local `Vec<ComplexT>`, and writes the result back **element-wise** through the
registered to-Python + (hardened, ADR-0006) numpy `setitem` path:

```cpp
SuccessCode track_path_wrap(TrackerT const& self,
    boost::python::object result_obj,          // was Eigen::Ref<Vec<ComplexT>>
    ComplexT start_time, ComplexT end_time,    // by value — now safe, no Ref converter
    Vec<ComplexT> const& start_point)
{
    Vec<ComplexT> temp_result(self.GetSystem().NumVariables());
    auto code = self.TrackPath(temp_result, start_time, end_time, start_point);
    for (Eigen::Index i = 0; i < temp_result.size(); ++i)
        result_obj[i] = temp_result(i);        // in-place writeback
    return code;
}
```

With no writable-Ref converter present, the by-value scalar times are no longer
corrupted. This mirrors **ADR-0002**, which resolved the same hazard in the endgame
bindings by refactoring the API so no scalar is adjacent to a writable Ref. A backstop
guard (gated to multiprecision trackers) rejects a precision-0 `mpc` input with a labelled
`std::runtime_error` rather than letting it reach libmpfr.

The Python API is unchanged: `track_path(result, start_time, end_time, start_point)` still
takes the result array first and updates it in place.

**Superseded rule.** ADR-0001's "pass adjacent scalars by value" is *not* a sufficient
defense. The rule going forward: **a binding must not expose a writable
`Eigen::Ref<Vec<mpfr/mpc>>` at all when mpfr/mpc data is involved** — take the output as a
numpy `boost::python::object` and write back element-wise (or return it).

## Consequences

- **Robust across eigenpy version, compiler, and ABI** — it removes the trigger condition
  (no writable-Ref converter), rather than depending on storage layout.
- Verified on x86_64 (CI run 27140927026): `end_time.precision=30`, `amptracking_test.py`
  9 passed; full suite 219 passed locally (writeback correctness confirmed).
- **Small per-element writeback cost** — negligible for multi-second path tracking.
- **Other bindings still carry the ADR-0001 shape and must get the same treatment.** The
  `refine` bindings (`tracker_export.hpp` `Refine3`/`Refine4`) take a writable `Vec<T>&`
  *plus* an adjacent `mpc` scalar and are bound as raw member-fn pointers with no wrapper.
  They are **not exercised by any test**, so they did not surface in CI, but they are a
  latent crash for users on x86_64. Fix them with `refine_wrap` functions (output as numpy
  object + writeback, scalars by value) **and add a `refine` test**. The single-arg writable
  Refs (`eval_wrap_1`, `rescale_wrap_inplace_mpfr`, `get_precision_vector`) have no adjacent
  scalar and are *not* in this corruption class (ADR-0001's scope note).
- **Worth an upstream eigenpy bug report.** Root cause is eigenpy's writable-`Eigen::Ref`
  rvalue converter (`eigen-from-python.hpp`) overrunning adjacent argument storage; related
  to stack-of-tasks/eigenpy#365. A minimal reproducer (a `probe(Eigen::Ref<VectorXcd>,
  mpc_complex a, mpc_complex b)` that returns the scalars' precision; one comes back 0 on
  x86_64 manylinux) would let this be fixed upstream and the workaround eventually dropped.
