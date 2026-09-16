# ADR-0057: Precision follows the point, not the System

**Status:** Accepted
**Date:** 2026-08-27

## Context

Three objects carried a precision and had to be kept in agreement by hand: the point, the
System, and the System's SLP. Evaluation *refused* when they disagreed, so every precision
boundary became a crash site and every caller grew the same boilerplate:

```cpp
auto target = max(Precision(point), system.precision(), DefaultPrecision());
system.precision(target);
auto lifted = lift(point, target);
```

Worse than unergonomic, it was **inescapable in one direction**. An SLP's `Memory` takes its
precision from the ambient `DefaultPrecision()` when the program is *lazily compiled*, while
the System keeps whatever it was told. Anything that moved the ambient default — an AMP
tracker or endgame, routinely — split them, and then *both* directions threw:

```
eval at the system's OWN precision  ->  "variable_values and SLP must be of same precision: 30 16"
eval at the ambient precision       ->  "input point (16) must match the system (30)"
```

and `System::precision(n)` could not repair it, because it and
`StraightLineProgram::precision(n)` both short-circuit when handed the value they already
hold. The System was wedged. See #377.

Downstream, this produced a small ecosystem of coping helpers, a precision *ratchet* on
shared systems (repeated elevated-precision work compounding 60 → 90 → 120), and ultimately
a two-tier workaround where sharpened results were handed back **downcast** to ambient
because delivering full precision would have meant touching every consumer.

## Decision

**A System has no precision to manage. Evaluation complies with the precision of the point
it is handed.**

1. `SetVariableValues` (SLP) and `System::SetVariables` **align** to the incoming point
   instead of throwing. Alignment re-materializes: constants are rebuilt from their exact
   recipes and block/patch coefficients from their highest-precision masters, so accuracy is
   reconstructed rather than padded with zeros.
2. `SetPathVariable` may only ever *raise* the precision — the variables are already in
   memory by then, and re-tagging downward would truncate them. Memory therefore ends an
   evaluation at the **max of its arguments' precisions**.
3. **A System carries no precision at all** — no setter, no getter, no `precision_` member,
   in C++ or in Python. There is also **nothing to prepare and nothing to fan out**: a
   System has no way to be told a precision.
4. **Every evaluable type self-aligns**, under one name and one body:

       template <typename T>
       void SyncPrecision(Vec<T> const& vars) const   // no-op for double
       {
           if (vars.size() && Precision(vars(0)) != precision_)
               Precision(Precision(vars(0)));
       }

   called first thing in each `EvalInPlace` / `JacobianInPlace` / `TimeDerivInPlace`. Six
   types have it: `linear_forms_block`, `products_of_linears_block`, `randomization_block`,
   `blend_block`, `patch`, and the SLP. Each keeps its own "materialized at" tag and
   short-circuits, so the steady state is one integer compare.

   This is not a new invention. `blend_block` and `randomization_block` already did exactly
   this independently — the change makes the pattern *and the name* uniform. The patch tells
   the same story from the other side: it carried a **commented-out assert** demanding that
   callers match its precision, disabled because it could not be honoured. It now aligns
   itself instead.
5. Result precision is **derived**: `CoerceBlockOutputPrecision` takes its target from the
   staged point rather than from a stored field.
6. **Coefficient generation is a construct-time concern.** The start systems (`mhom`,
   `total_degree_linear_product`) build their blocks' exact masters under
   `DefaultPrecision(MaxPrecisionAllowed())` — which is all that matters, since the ambient
   default governs the precision new mp values are *born* at. They no longer touch any
   working precision: setting the precision of working coefficients is an evaluation
   concern only.

### Who holds precision truth

| holder | role |
|---|---|
| **the point** | **the authority.** Evaluation happens at its precision. |
| ambient `DefaultPrecision()` | the precision new objects are *born* at. Not consulted during evaluation. |
| System / block / patch / SLP `precision_` | a *cache tag* — "what my working copies are currently materialized at". Derived, never authoritative. |
| the exact masters (`constant_recipes_`, `coefficients_highest_precision_`) | the source of truth for the values themselves, at any precision. |

The masters are what make this sound: every precision-carrying object can be
**re-materialized, never re-drawn**. Re-drawing a random patch or slice at a new precision
would silently change the system.

## Consequences

- **Do not re-introduce a precision setter on System**, in C++ or in Python, and do not make
  evaluation refuse a precision mismatch again. Both were the bug, not the guard rail.
- Callers must not compute a target precision and push it into a System. Evaluate at the
  precision you want; the System follows.
- Do not re-introduce a "prepare"/"materialize at" entry point on System. An earlier draft of
  this change kept one, on the theory that evaluations consuming *staged* values need the
  system prepared beforehand. That theory is false: every such evaluation is preceded by an
  `Eval(point)` that stages the point and aligns as a side effect. Removing all **22** of its
  callers — 13 in the trackers/endgames/`zero_dim_solve`, 5 pushing into sub-Systems from
  blocks, 4 in start systems — changed no test result. It was pure ceremony, and a long
  self-documenting name made it look deliberate.
- **Archive format changed.** `precision_` is no longer serialized. It was transient
  evaluation state that should never have been archived — the same `serialize` already
  excludes `current_variable_values_` for exactly that reason. Boost archives here carry
  Systems between MPI ranks of one run and are not a durable format (the forever-contract
  lives in the records/JSON path, ADR-0042), but a System archived by an older build will
  not load into a newer one.
- Do not restore a System-level precision field "for convenience", and do not re-introduce
  a "prepare"/"materialize at" entry point on System. An earlier draft of this change kept
  one, on the theory that evaluations consuming *staged* values need the system prepared
  beforehand. That theory is false: every such evaluation is preceded by an `Eval(point)`
  that stages the point and aligns as a side effect. Removing all 22 of its callers — 13 in
  the trackers/endgames/`zero_dim_solve`, 5 pushing into sub-Systems from blocks, 4 in start
  systems — changed no test result. It was pure ceremony, and a long self-documenting name
  made it look deliberate.
- If some code needs to know the working precision, it is the precision of the point being
  evaluated; ask the point.
- Threading is unaffected: `ZeroDimSolver` already gives every worker a deep-cloned System,
  precisely because "residual evaluation mutates System precision state".
- Regression tests fix this in place: a system whose ambient default moved underneath it must
  still evaluate; one system must serve a whole ladder of precisions; a constant inexact in
  decimal (1/3) compiled at 16 digits must return residual **0** at 100 digits, which holds
  only if re-tagging rebuilds rather than pads.

## What this does NOT extend to

`current_variable_values_` is superficially the same kind of thing — mutable per-evaluation
state on a shared object, and the reason `ZeroDimSolver` deep-clones a System per worker
thread. It stays, deliberately.

It earns its place: it backs the no-argument `Eval()` / `Jacobian()` entry points, so the
SLP's tape runs once and a caller can read function values now and the Jacobian later
without re-running it. That is a real optimization, not a leftover.

Moving it into a separate evaluation handle (`sys.At(point).Values()`) was considered and
**rejected**: it would tax the common case — `sys.Eval(point)`, one call, no ceremony — to
buy thread-safety that per-thread cloning already provides. Direct evaluation is the
interface worth protecting. Do not re-propose it.

The distinction that matters: `precision_` was state a CALLER had to reason about and keep
in sync, and it had no functional role left. `current_variable_values_` is state the caller
never sees, and it does real work.

Closes #377.
