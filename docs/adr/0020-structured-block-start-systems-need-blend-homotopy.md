# ADR-0020: A structured-block start system must be coupled into a homotopy with a blend, not node arithmetic

**Status:** Accepted
**Date:** 2026-06-16

## Context

A user can now author their own start system from Python as a product of linear forms
(`linalg.add_products_of_linears` → `System.add_products_of_linears_block` →
`blocks::ProductsOfLinearsBlock`) and solve a target with it via the user-homotopy path
(`nag_algorithm.user_homotopy`, ADR-0013/0014 machinery). To do that the user needs a *homotopy*
`H = (1-t)·target + γ·t·start` that they hand to `user_homotopy`.

The obvious way to build that homotopy from Python is the existing
`nag_algorithm.coefficient_parameter_homotopy`, which is `(1-t)*target + t*generic` using the
bound `System` arithmetic operators. **That silently produces a wrong homotopy when the start
system carries a structured block.** `System::operator+=` / `operator*=`
(`core/src/system/system.cpp`) operate only on `PolyBlock().Functions()` — the function-tree
(polynomial) block. A products-of-linears start system has *no* polynomial functions; its rows
live in a `ProductsOfLinearsBlock`. So the node-arithmetic combination drops the start system's
actual equations: at `t=1` the "homotopy" does not vanish at the start points, and the solve is
meaningless (no error is raised — hence "silently").

This is not a bug to fix in the operators: blending a non-function-tree block into a homotopy is
exactly what `blocks::BlendBlock<System>` exists for (see the block-composed-System design), and
the zero-dim solver's own generated start systems (MHom) already ride it. The internal
`policy::CloneGiven::FormHomotopy` (`nag_algorithms/common/policies.hpp`) already branches on
`start.HasStructuredBlocks()`: structured ⇒ build a `BlendBlock` homotopy; otherwise ⇒ node
arithmetic. But that logic lived inside a policy method, unreachable from the user-homotopy
(`RefToGiven`) path, where the user supplies the homotopy.

## Decision

**Factor the homotopy-forming logic into one reusable free function and expose it; do not
re-derive it, and do not route a structured-block start system through node arithmetic.**

- `bertini::MakeHomotopy(target, start, path_variable="t", gamma=nullptr)`
  (`system.hpp`/`system.cpp`) holds the single copy of the `(1-t)·target + γ·t·start`
  construction, branching on `HasStructuredBlocks()` (BlendBlock vs node arithmetic) and
  generating a random rational `γ` when none is given.
- `policy::CloneGiven::FormHomotopy` now *calls* `MakeHomotopy` (one home for the logic, per the
  reuse-don't-duplicate principle).
- It is bound as `system.make_homotopy` and wrapped as
  `nag_algorithm.straight_line_homotopy(target, start, *, path_variable, gamma)`, the Python entry
  point for turning a user-authored start system (structured or not) into a trackable homotopy to
  pass to `user_homotopy`.

That one wrapper is therefore the correct tool whether or not the start system carries a
structured block, and whether or not the deformation wants a gamma.

> **Naming note (3.5, b2#371).** This wrapper was called `blend_homotopy`, and a second wrapper
> `coefficient_parameter_homotopy` did the same thing with `γ` fixed at 1. "Blend" named the
> implementation, not the mathematics, and the second name claimed a specialization that was only
> a choice of `γ`. Both were replaced by the single `straight_line_homotopy`, whose `gamma=1` is
> the coefficient-parameter case. The decision recorded here is unchanged: a structured-block
> start must be coupled by a blend, never by `System` node arithmetic.

## Consequences

- A products-of-linears (or any structured-block) start system solves end to end from Python:
  validated affine, single-variable-group, with no homogenization needed — target
  `{x²+y²−1, y−x²}`, start `{(x−1)(x+1), (y−1)(y−2)}`, four hand-written hyperplane-intersection
  start points tracked to the four target roots
  (`python/test/zero_dim/products_of_linears_homotopy_test.py`).
- The blend branch is exercised by both the generated path (MHom, via `FormHomotopy`) and the
  user path (via `straight_line_homotopy`), so they cannot drift apart.
- `γ` is optional: omit it for a random rational (genericity), or pass an exact node
  (`linalg.coefficient(...)`) off the real axis for a reproducible straight-line path.
- The coefficient-parameter case is `gamma=1`, and it goes through `MakeHomotopy` like every
  other, so a structured-block `generic` is blended rather than silently dropped by node
  arithmetic. It keeps its `(1-t)/t` semantics; it is no longer a footgun.
- Authoring scope (validated): a user-authored products-of-linears start solves end to end for a
  single affine group and for multiple affine groups; projective (homogeneous) factors construct
  and evaluate. The block is evaluated as authored (its `Homogenize` is a no-op), so for a
  projective group you write homogeneous factors (constant column 0).
