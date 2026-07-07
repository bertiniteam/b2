# ADR-0039: The Python `ZeroDim` binding must keep the given `System` alive (custodian-and-ward)

**Status:** Accepted

## Context

A segfault was found from idiomatic Python usage: building a `System` in a helper function and
passing it straight to the solver, then solving with an observer attached:

```python
def target_system():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System(); s.add_function(...); s.add_variable_group(pb.VariableGroup([x, y]))
    return s

zd = ZeroDim(target_system(), mptype='adaptive')   # System passed as a temporary
coll = nag_algorithm.observers.SolutionPathCollector(); zd.add_observer(coll)
zd.solve()                                          # <-- SIGSEGV
```

Holding a Python reference to the system (`s = target_system(); zd = ZeroDim(s, ...)`) made it
disappear. The faulthandler C stack pinned the cause:

```
ZeroDim::ExecuteDuringEG -> System::EvalInPlace -> PolynomialBlock::Differentiate
  -> ~SumOperator -> ~NaryOperator
  -> boost::python::converter::shared_ptr_deleter::operator()  -> _Py_Dealloc   (on a freed object)
```

Function-tree nodes built in Python (via operator overloads like `x*y - 3*x + 2`) are held by
`shared_ptr` whose **deleter is tied to the Python wrapper object** (boost.python's
`shared_ptr_deleter`). The `ZeroDim` solver's internal `System` is a **shallow** copy that shares
those node `shared_ptr`s. When the caller passes a temporary `System` and keeps no reference, Python
destroys the System and its node wrappers; the solver, however, still holds the shared nodes. During
the endgame the solver differentiates the system, a shared node's refcount reaches zero, and its
deleter dereferences an already-freed Python object → crash. (An observer is what reliably drives the
per-path endgame differentiation that trips it; a plain solve happened to dodge it, which is why this
went unnoticed.)

## Decision

Make the binding own the lifetime dependency instead of relying on the caller. In
`python_bindings/include/zero_dim_export.hpp`, add boost.python **`with_custodian_and_ward`** to the
`ZeroDim` constructors so the given system(s) are kept alive as long as the solver:

- `ExportZeroDimSpecific` (total-degree solvers, `init<SystemT>()`):
  `init<SystemT>()[with_custodian_and_ward<1, 2>()]` — ward = the System (arg 2), custodian = the
  solver (self, arg 1).
- `ExportZeroDimUserHomotopy` (`RefToGiven`, holds references to all three systems): chain
  `with_custodian_and_ward<1, 2, with_custodian_and_ward<1, 3, with_custodian_and_ward<1, 4>>>()`
  over target / start / homotopy.

This makes the documented `ZeroDim(sys)` usage safe regardless of whether the caller retains a
reference to `sys`.

## Consequences

- The crash is gone (the minimal repro and the classic-cartoon tutorial both run clean). No behavior
  change for callers who already kept a reference.
- It is purely additive (a lifetime link), so it cannot shorten any lifetime or break existing code.
- Deeper alternative (not taken): make the system-management policy deep-`Clone()` so the solver owns
  pure-C++ nodes with no Python-tied deleters. That is a larger change to the C++ ownership model;
  the binding keep-alive is the minimal, correct fix and matches what a careful caller already does.
- Other bindings that take a `System`/start system by reference and outlive the call (e.g. NID) should
  be audited for the same keep-alive need; out of scope here.

## Update (ADR-0040)

The `ZeroDim` class was split into `ZeroDimSolver` (the algorithm) and `HomotopySolver` (the engine).
This keep-alive is unchanged in spirit and still applied: `ZeroDimSolver(system)` keeps the system
alive (`with_custodian_and_ward<1,2>`), and `HomotopySolver(target, start, homotopy)` — which holds
all three by reference — keeps all three alive via chained custodian-and-ward.
