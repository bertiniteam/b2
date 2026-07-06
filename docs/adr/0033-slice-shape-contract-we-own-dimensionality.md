# ADR-0033: Slice shape is our accessors' contract, not eigenpy's (the 1-D collapse)

**Status:** Accepted
**Date:** 2026-06-24

## Context

eigenpy converts a dynamic `Eigen::Matrix` to numpy by **collapsing "vector-shaped" matrices to
1-D**. The decision is hardcoded in `eigenpy/eigen-to-python.hpp` (`eigen_to_py_impl_matrix::convert`):

```cpp
if ( ( (!(C==1) != !(R==1)) && !MatrixDerived::IsVectorAtCompileTime )
     || MatrixDerived::IsVectorAtCompileTime )
{ npy_intp shape[1] = { C==1 ? R : C }; ... }   // 1-D
else
{ npy_intp shape[2] = { R, C }; ... }            // 2-D
```

`(!(C==1) != !(R==1))` is **`XOR(rows==1, cols==1)`** — *exactly one* extent is 1. Measured:

| C++ matrix (runtime) | numpy shape |
| --- | --- |
| 2×3 | `(2, 3)` |
| 1×3 | `(3,)` |
| 2×1 | `(2,)` |
| 1×1 | `(1, 1)` (not collapsed — both extents 1, XOR false) |

So a `Mat<T>`'s Python shape depends on the **runtime data**: a one-row matrix is 1-D, a many-row
matrix is 2-D, and 1×1 is a 2-D corner that contradicts the rule. This is **silent and
data-dependent** — code written/tested against multi-row data (`M[i,j]`, `M.shape[1]`, `list(M)` for
rows, `M @ v`) breaks the moment the data is a single row/column, with no error at the boundary.

**For NID this is not a corner case.** A dimension-`d` component is cut by `d` hyperplanes, so a
**curve (d=1) has a single-form slice** that collapses every time; curves are the most common
positive-dimensional component. The d=2 surface you test on won't collapse and the d=1 curve in
production will.

There is **no eigenpy knob.** The condition above consults nothing mutable (no flag, template
option, env, or macro); the only global state in eigenpy is `NumpyType::sharedMemory()` (copy vs
view, not dimensionality); and the old np.matrix mode (which forced 2-D) has been removed in this
version. Investigated and confirmed dead — see ADR-0001/0006/0031 for the other eigenpy hazards;
this is the same *family* (eigenpy's Python view isn't the faithful object you'd assume) but a
*shape/semantics* kind rather than a memory-corruption kind.

A footgun within the footgun: the obvious `np.atleast_2d` fix is **orientation-specific** — it
re-expands a collapsed *row* `(N,)` to `(1, N)` correctly, but a collapsed *column* `(N,)` to
`(1, N)` *wrongly* (should be `(N, 1)`). Only the side that produced the matrix knows its orientation.

## Decision

**The shape contract belongs to *our* accessors, not to eigenpy (nor nanoeigenpy, nor any binding
library). A `Slice` is a Python sequence of linear forms, with list semantics:**

- `Slice.coefficients()` → **always 2-D**, `(num_forms, num_variables+1)`. We reshape using the
  dimensions known in C++ (`dimension()`, `num_variables()`), so it is orientation-correct by
  construction and never collapses.
- `slice[i]` (integer) → an **element**: the i-th form's coefficient **vector** (1-D). Iterating a
  slice yields these form vectors. "A single linear form *is* a vector" — but you get it by indexing
  an element, explicitly, not because the form count happened to be 1.
- `slice[i:j]` (and a list/tuple of indices) → a **sub-collection**: a new (sub-)`Slice`.

Mental model: a slice is a *stack of form-vectors*; the matrix is the stack, a form is a vector. The
vector-vs-matrix choice is **named** (element index vs collection slice vs matrix accessor), never an
emergent property of the data size. The collapse you might *want* (one line = one vector) is still
there — it is `slice[i]` — but it is a choice you ask for, not a surprise the data hands you.

Implemented in `python/bertini/nag_algorithm/__init__.py` (the `Slice` ergonomics): a `coefficients`
wrapper that reshapes to the known 2-D shape, and a `__getitem__` that returns a form vector for an
integer index and a sub-`Slice` for a Python slice / index list.

## Consequences

- **Migration-proof.** When boost.python + eigenpy is eventually replaced (nanobind + nanoeigenpy,
  per the team's plan), whatever collapse rule the new layer uses is absorbed by our wrapper. Slice
  consumers — including the NID code — see no change. Owning the contract is what makes that swap a
  non-event.
- The raw hazard still lives in **un-wrapped** matrix-returning accessors: `System.eval_jacobian`
  (a 1-function or 1-variable system gives a collapsed Jacobian), the randomization matrix, and patch
  coefficients. These are **deliberate future retrofits** — changing their shape is an interface
  change with its own test fallout. Until then, defend at the call site, remembering `atleast_2d` is
  only correct for row-oriented results (`reshape(rows, cols)` with the known dims otherwise).
- Supersedes the working note that prompted this (`eigenpy returns a 1-row matrix as 1-D`).
