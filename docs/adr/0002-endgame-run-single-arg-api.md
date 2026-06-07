# ADR-0002: Endgame run() takes only the start point; boundary times set at construction

**Status:** Accepted  
**Date:** 2026-06-05

## Context

The original Python endgame binding exposed:

```python
eg = FixedDoubleCauchyEG(tracker)
code = eg.run(complex(0.1, 0), sample)   # (start_time, start_point)
```

The C++ wrapper was:

```cpp
.def("run", &EndgameBaseVisitor::WrapRunDefaultTime,
     (arg("self"), arg("start_time"), "start_point"), ...)
```

This placed an `mpc_complex` scalar (`start_time`) adjacent to a `Vec<BCT> const&`
(`start_point`). The eigenpy from-Python converter for the vector argument writes into
a static rvalue-converter slot, potentially corrupting the adjacent scalar — the same
mechanism described in ADR-0001. On macOS CI this produced a SIGABRT in
`FixedDoubleCauchyEG::run`.

The `track_path_wrap` fix (ADR-0001) passed scalars by value to sidestep the issue.
That option was evaluated here too, but the endgame has two time arguments
(`start_time` and `target_time`) and a cleaner API was possible.

## Decision

Refactor so that boundary times are supplied at construction, not at each `run()` call:

```python
eg = FixedDoubleCauchyEG(tracker, boundary_time=complex(0.1, 0))
code = eg.run(sample)   # only the start point
```

In C++:
- `base_endgame.hpp`: added `start_time_` / `target_time_` members, `SetBoundaryTime()`,
  and `Run(Vec<BCT> const&)` that reads the stored times internally.
- `endgame_export.hpp`: replaced `WrapRunDefaultTime` / `WrapRunCustomTime` with a
  single `WrapRun(self, Vec<BCT> const& s)` — no adjacent scalar argument.
- `python/bertini/endgame/__init__.py`: endgame classes now take `(tracker, boundary_time)`
  at construction; `run(point)` takes only the point.

`WrapRun` has `self` (C++ reference, no conversion) and `Vec<BCT>` (numpy array).
There is no adjacent scalar for eigenpy to corrupt.

## Consequences

- **Eliminates the corruption crash** on macOS (and any platform where malloc returns
  non-zero memory). No writable Ref + adjacent scalar combination exists in the endgame
  binding.
- **API change:** `eg.run(t, point)` → `eg.run(point)` with `boundary_time` at
  construction. All Python tests updated. Any user code using the old API must update.
- **Cleaner semantics:** An endgame naturally has a fixed boundary time for its lifetime.
  Storing it at construction better reflects the object's identity.
- **Rule reinforced:** Never combine a writable `Eigen::Ref<Vec<mpc_complex>>` with
  adjacent `const&` scalar mpc_complex args in a Boost.Python binding. If you need both,
  either pass scalars by value (ADR-0001 approach) or restructure the API to eliminate
  the adjacency (this ADR's approach).
