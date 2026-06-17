# ADR-0023: Moving homotopies blend only the moving rows; fixed equations stay sibling blocks

**Status:** Accepted
**Date:** 2026-06-17

## Context

Regeneration and moving-slice continuation need a homotopy in which most equations are **fixed** (the
polynomial system, plus "below" linear slices that cut dimension) and only a small part **moves** with
the path variable `t`. The fixed equations must be evaluated **once** per step — never duplicated,
never scaled by the path coefficient, never differentiated in `t`. Two operations are wanted:

1. move one or several linear slices (slide hyperplanes), and
2. deform a product-of-linears into a polynomial — the regeneration "add a degree" step
   `(1−t)·f + γ·t·∏Lᵢ` for one row.

The existing whole-system homotopy (`MakeHomotopy` → `BlendBlock`) blends *whole* systems
`(1−t)·target + γt·start`. If the polynomial part is shared between the operands it is evaluated in
**both** and comes out as `((1−t)+γt)·f ≠ f` — both wrong and wasteful for the moving-rows case.

## Decision

**No new block types.** `BlendBlock<System>` already linearly combines operand systems with shared
path-variable coefficient nodes — it *is* a "moving wrapper". The fix is the construction: blend only
the **moving** rows, and keep every fixed equation as its own sibling block in the homotopy.

`MakeMovingHomotopy(fixed, start_moving, end_moving, t, gamma)` (`system.cpp`; Python
`nag_algorithm.moving_homotopy`, plus the `make_moving_homotopy` binding) builds

    H = [ fixed's blocks (unchanged) ;  BlendBlock( end_moving, start_moving ; (1−t), γ·t ) ]

- `fixed`'s blocks are kept (not `ClearBlocks`); being autonomous, the block loop
  (`system.hpp` `EvalBlocksInPlace`/`TimeDerivBlocksInPlace`) evaluates them once and their
  `TimeDerivInPlace` writes zeros — the fixed rows are exactly zero in `dH/dt`.
- only the moving operands enter the blend, so the fixed equations are never duplicated or scaled.
- at `t=1` the moving rows are `γ·start_moving` (start points = roots of `fixed` ∧ `start_moving`);
  at `t=0` they are `end_moving`. Row order is fixed rows then moving rows; the matching target is
  `concatenate(fixed, end_moving)`.

Both wanted operations are **output-level** blends and need nothing more: a moving linear slice is a
blend of linear-forms operands (= interpolating their coefficients = sliding hyperplanes), and a
∏L→polynomial deformation is a blend of a products-of-linears operand and a polynomial operand. The
tracked homotopy is affine, matching the existing user-homotopy pipeline (whose `SystemSetup` is a
no-op — it tracks the supplied homotopy as-is).

Supporting changes: bound `System::TimeDerivative` to Python as `eval_time_derivative` (so the
`dH/dt`-is-zero-on-fixed-rows invariant is checkable, and shown in the tutorial); taught
`System::NaturalFunctionsAsNodes` to expand a `LinearFormsBlock` (affine and homogeneous), so the
`ExpandToFunctionTree` oracle covers slice-bearing homotopies.

## Consequences

- The fixed system is provably evaluated once: `dH/dt` is exactly zero on every fixed-block row,
  asserted in `moving_homotopy_test.cpp` (also: `t=0`/`t=1` endpoints, expansion-oracle agreement in
  `dbl`/`mpfr`, a static linear slice + a moving slice both fixed, and a products-of-linears→polynomial
  deform) and in `python/test/zero_dim/moving_homotopy_test.py` (three hand-checkable demos).
- Tutorial `tutorials/moving_slice.rst` drives it (doctested), framed as regeneration.
- **"Factor-sliding" a product-of-linears** (each factor's hyperplane moving independently) is a
  *different*, deliberately out-of-scope semantics; the regeneration need is the output-level blend.
- Deferred: the higher-level regeneration cascade driver that sequences these moving homotopies, and
  a projective/patched variant if affine tracking proves insufficient for some target.
