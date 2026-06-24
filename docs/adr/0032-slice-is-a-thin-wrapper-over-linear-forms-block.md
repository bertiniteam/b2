# ADR-0032: A witness set's Slice is a thin wrapper over LinearFormsBlock, not a parallel type

**Status:** Accepted
**Date:** 2026-06-24

## Context

A witness set is a triple (system, slice, points). Exposing `WitnessSet` to Python (to *enable*
a numerical irreducible decomposition (NID) being prototyped in Python, to be ported to C++ later)
forced a decision about what the **slice** should be.

The slice has three faces at once: it is a *variety* (it evaluates, it has a Jacobian, it can be
added to a system), it is a *matrix of coefficients* (you want its rows, the first/last m of them),
and it carries *affine/projective* concerns we mostly hide from the user. We want it easy to reason
about, write code with, and carry around for further work — and, for regeneration, to compose with
products of linear forms.

There were already **two** representations of "a stack of linear forms" in the tree:

- `bertini::LinearSlice` (`core/include/bertini2/system/slice.hpp`) — the original slice type. It
  was **orphaned**: not a `System` evaluation block, used by no algorithm (its only "consumer",
  `nag_algorithms/trace.hpp`, is a stub returning 0), and **not exposed to Python at all**. It
  carried its own mpfr-master / per-type-working-copy / `Precision()` machinery, its own
  QR-orthogonal random-slice factory, and a separate constants vector.

- `bertini::blocks::LinearFormsBlock` (`core/include/bertini2/system/blocks/linear_forms_block.hpp`)
  — the modern, integrated equivalent, written later. A member of the `Block` variant, evaluated
  in-System, homogenization-aware, unit-tested, already authored from Python via
  `bertini.linalg.add_linear` / `add_linear_forms`. Its augmented coefficient matrix `M` (rows =
  forms, columns = num_vars+1, the last column the constant term) is the **same layout** that
  `ProductsOfLinearsBlock` (the m-homogeneous start system and regeneration) uses for its factors.
  Its own header even notes it "mirrors the pattern of `bertini::LinearSlice`."

So `LinearFormsBlock` had already obviated `LinearSlice`; the block was simply written without going
back to delete its predecessor. Keeping both means maintaining two parallel "linear forms with an
mpfr master" types forever — exactly the kind of duplication that invites drift.

## Decision

**Retire `LinearSlice`. A witness set's slice is `bertini::Slice`, a thin wrapper that *holds* a
`LinearFormsBlock` plus the slice-level metadata the block does not carry (the sliced
`VariableGroup`, and whether the forms were authored homogeneous).**

- `Slice` delegates eval / Jacobian / precision / coefficients to its `LinearFormsBlock`, and
  salvages `LinearSlice`'s QR-orthogonal random factories (`RandomReal` / `RandomComplex` / `Make`).
- The coefficient matrix is exposed directly (`Coefficients()`), and row-subsetting (`Head` / `Tail`
  / `Rows`, plus Python `[]`) returns a new `Slice` — so the matrix face and the variety face are
  the same object.
- `Slice` is **not** added to the `Block` variant. Instead `AddTo(System&)` / `AsSystem()` hand its
  block to a System, keeping the closed-set variant clean and avoiding the recursive-type machinery.
- Because a slice's rows and a `ProductsOfLinearsBlock`'s factor rows share the augmented
  `num_vars+1` layout, composing slices into a product-of-linears (the regeneration bridge,
  `linalg.add_slices_as_products`) is matrix-row assembly, not new math.

`WitnessSet`'s slice member is now a `Slice`. The wrapper is also where serialization,
`Concatenate` (Python `+`), and a readable `repr` live.

## Consequences

- One representation of linear forms in the evaluation path; `LinearSlice` and its duplicate
  precision/eval machinery are gone.
- The slice is numpy-native (coefficients round-trip through eigenpy) and authored by the existing
  `bertini.linalg` layer, so the Python NID prototype has an ergonomic, composable slice today, and
  the C++ `Slice` is the forward-investment for the later C++ NID port.
- "Homogeneous" is currently **authoring metadata** (a zero constant column), not full patch
  awareness. Genuine affine/projective patch handling for slices remains future work.
- A one-row coefficient matrix returns from eigenpy as a 1-D array (a general eigenpy convention);
  helpers that treat it as a matrix use `np.atleast_2d`.

See `core/include/bertini2/system/slice.hpp`, `core/include/bertini2/nag_datatypes/witness_set.hpp`,
and `python_bindings/src/nid_datatypes_export.cpp`. Related: ADR-0025 (randomization as a first-class
block), ADR-0001 (eigenpy writable-Ref hazard).
