# ADR-0025: Randomization of overdetermined systems is a first-class block

**Status:** Accepted
**Date:** 2026-06-17

## Context

An overdetermined system (N functions, n variables, N > n) with isolated solutions cannot be fed
to a total-degree start system, which requires a square target (`system/start/total_degree.hpp`
throws otherwise). The standard fix is **randomization**: replace the N functions with n generic
combinations `G = R·F` whose isolated solutions still contain the target's; solve the square
system and discard the extraneous solutions. This is the squaring-up the regeneration roadmap
needs, and ADR-0021 already names a future `RandomizationBlock`.

We implement it as a first-class, multidegree-native evaluation block rather than by inflating the
combinations into function-tree nodes ("everything is a block").

## Decision

### The block wraps a System and applies a constant matrix with homogenizing-variable powers

`RandomizationBlock<System>` (templated like `BlendBlock<System>` so `System` stays a dependent
name and the `Block` variant closes) holds the overdetermined system as a `shared_ptr<System>`
operand plus a constant coefficient matrix, and emits

    g_i(x) = sum_j  c_ij · f_j(x) · prod_g  h_g^{(D_{i,g} − d_{j,g})}

where `f_j` is the operand's j-th natural function (homogenized to its own multidegree `d_j`), `D_i`
is row i's target multidegree, and `h_g` is group g's homogenizing variable. The h-power factors
are what let a **constant** combination of functions of *differing* degree homogenize correctly:
before homogenization every `h_g = 1` and it collapses to `g_i = Σ c_ij f_j`. The block is
multidegree-native, so single-affine-group and multihomogeneous randomization share one evaluator,
and (per ADR-0021) it fully defines its own rows.

### Construction (`System::Randomize`)

- **Single affine group → optimal.** Sort the functions by descending degree and use `R = [I | C]`
  with random `C`. Then `deg g_i = d_i` and the total-degree path count is the product of the n
  largest degrees — the minimum. The sort is *required*, not cosmetic: it makes every deficit
  `D_{i,g} − d_{j,g} ≥ 0` (a negative deficit has no h-power representation).
- **Several variable groups → correct.** Multidegrees are only partially ordered, so use a dense
  random `R` with a common (componentwise-max) target multidegree. Correct, and optimal when the
  functions share a multidegree (e.g. a bilinear system).
- **User-supplied matrix.** `Randomize(R)` uses `R` verbatim, functions in their given order.

`Randomize` returns a **new** system and never mutates the caller's (the descending sort happens
on an internal copy). The matrix is retrievable via `RandomizationMatrix()` (toward future
`R · [f]` block stringification).

### Same homogenizing variables across the operand boundary

The h-power factors must reference the *actual* homogenizing-variable nodes, not coincidentally
aligned ordering positions. `RandomizationBlock::Homogenize` threads the owning system's hom var
into the operand by homogenizing it through the new `System::Homogenize(VariableGroup const&)`
overload, which reuses supplied hom vars instead of minting fresh ones. Operand evaluation and the
block's deficit padding then live in one coordinate system.

### Load-bearing fix: `MakeHomotopy` must blend whenever EITHER side is structured

`MakeHomotopy` previously chose the blend-block path only when the **start** system carried a
structured block, falling back to node arithmetic otherwise. But `operator+` / `operator*` on
Systems only combine the `PolynomialBlock` functions and **silently ignore structured blocks** — so
a structured-block *target* (exactly what a randomized system is) was dropped from the homotopy,
which then tracked from the start system to nothing (empty / singular endpoints). `MakeHomotopy`
now blends whenever **either** target or start has a structured block, and the homotopy shell uses
`ClearBlocks()` (not `ClearFunctions()`, which only removes a `PolynomialBlock`) so a structured
target's rows are not double-counted alongside the blend. This generalizes the MHom case (ADR-0020)
to a structured target.

## Consequences

- Overdetermined systems solve end to end: randomize → solve the square system → filter extraneous
  solutions by re-evaluating the original. Covered by `randomization_block_test.cpp` (eval/Jacobian
  vs the function-tree expansion oracle, affine and homogenized, single- and multi-group, plus an
  ADR-0021 dirty-buffer check) and `python/test/zero_dim/randomize_test.py` (solve + filter).
- A merged tutorial (`tutorials/randomize.rst`) drives it, written as executable `.. testcode::`
  blocks; the docs CI now runs `sphinx-build -b doctest` so the tutorials are verified code.
- **Deferred:** provably-minimal target-multidegree *selection* for unequal multidegrees across
  several groups (the partial-order optimization, no clean `[I|C]` analog); the constant-`C` common
  target is correct meanwhile. Also deferred: deep-cloning a block operand for per-thread tracking
  (shared with `BlendBlock`), and `R · [f]` stringification.
