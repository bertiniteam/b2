# ADR-0011: Differentiation emits already-simplified trees; held functions are never mutated

**Status:** Accepted
**Date:** 2026-06-11 (PR #7)

## Context

Derivative trees came out full of junk (`(0*x^2*1+2*x*1*3*1+0*3*x^2)*y+...` for fxx of
`x^3*y`). Two root causes:

1. The ad-hoc pruning inside `Differentiate()` was type-narrow:
   `MultOperator::Differentiate` cast candidate derivatives to **`Float`** only, while
   the chain rules actually emit `Integer(0)`/`Integer(1)` — so the zero/one terms
   sailed through. (`SumOperator` cast to `Number` and worked.)
2. The cleanup path was worse than the disease: `Simplify()` / `EliminateZeros` /
   `EliminateOnes` / `ReduceDepth` all **mutate trees in place** (swapping
   `operands_`/`signs_` vectors), and `DefaultAutoSimplify()` was `true`, so plain
   `System::Differentiate()` ran that machinery over derivative trees — which **share
   subtrees with the user's original functions**. Holding `f` while differentiating a
   System could restructure `f` itself (and System copies share whole trees, making
   this thread-hostile as well).

A user holding `f` (or any subexpression of it) must never observe it change because
something was differentiated; and a derivative must never sprout variables that were
not in `f`.

## Decision

**Simplify at construction time inside `Differentiate()`; never rewrite after the fact.**

- `Node` gains exact `IsLiteralZero()` / `IsLiteralOne()` (overridden by
  Integer/Float/Rational on their stored literal values — no double-eval round trips,
  no underflow surprises; non-literal expressions that happen to evaluate to zero are
  deliberately NOT detected).
- Free factories `SimplifiedNegate` / `SimplifiedSum` / `SimplifiedMult`
  (arithmetic.hpp/.cpp) build **fresh** nodes from (operand, flag) pairs: literal
  zeros vanish; a multiplied literal zero collapses the product; literal ones drop;
  singletons unwrap; nested Mult factors flatten (reading, never writing, their
  operands — flag-inverting through divisions); exact Integer/Rational constants fold
  via mpq arithmetic (`6*x*y`, not `3*2*x*y`). **Floats are never folded** (precision
  semantics stay untouched) and **a literal-zero divisor is left visible** (no
  silent x/0 rewrites).
- Every operator `Differentiate` (arithmetic + trig) routes through the factories.
  `Differential` leaves (the no-arg Jacobian form) are not Numbers and are never
  pruned.
- **`DefaultAutoSimplify()` now returns `false`.** The post-hoc mutating pass is
  redundant for derivatives and was the channel by which held functions changed.
  `System::Simplify()` / `AutoSimplify(true)` remain available as *explicit* opt-ins —
  mutation is then the user's deliberate choice.

## Consequences

**Positive**
- Derivative degrees are exact (`fxx` of `x^3*y` has degree 2 — algebra, not an upper
  bound); printed derivative forms are canonical and pinned by tests
  (`d/dx sin(x)` is `cos(x)`).
- The user-visible contract is enforced by `differentiation_safety_test.py`: f and
  held subexpressions byte- and value-identical after every flavor of
  differentiation; `gather_variables(f')` never a superset of f's; shared Variable
  objects stay live (one `set_current_value` drives f and f').

**Negative / to watch**
- Behavior change: code relying on auto-simplify of derivatives must now call
  `Simplify()` explicitly (and accept its in-place mutation of shared subtrees).
- The in-place Eliminate*/ReduceDepth machinery still exists and still mutates; a
  copy-on-write replacement is a possible future ADR. Do not wire it back into any
  default path.
- New `Differentiate` implementations must route through the factories; building
  `Sum`/`Mult` nodes directly reintroduces junk.
