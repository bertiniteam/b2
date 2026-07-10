# ADR-0052: Multiprecision linear algebra by instantiating eigenpy's own decomposition visitors

**Status:** Accepted
**Date:** 2026-07-10

## Context

Building the codimension-1 step of NID in Python needs to solve `A x = b` where `A`, `b` are
`complex_mp` (slice coefficients). Two obvious tools both fail:

- **numpy** — `np.linalg.solve` / `lu` / `qr` / `svd` route into LAPACK, which only accepts
  `float32/64` and `complex64/128`. On a user dtype (`complex_mp` / `real_mp`) they raise; numpy's
  linear algebra never dispatches to the mp ufunc loops. numpy is the right tool for *double*, and
  the mp path exists precisely because numpy has none.
- **A stock `import eigenpy`** — its `PartialPivLU` etc. are C++ template instantiations baked into
  the PyPI shared library at *its* build time, for the standard scalars only. Registering our mp
  types with eigenpy (`registerNewType`, the Eigen⇄numpy converters) teaches eigenpy how to *marshal*
  a `Mat<complex_mp>` ⇄ numpy array, but does **not** create any decomposition instantiation for the
  mp scalar. So `eigenpy.PartialPivLU(A_mp)` finds no wrapper and bails.

We want real multiprecision LU / QR / SVD in Python, without reimplementing linear algebra and
without the caller having to type-if on dtype at every call site.

## Decision

### Instantiate eigenpy's own visitor templates on the mp matrix types, inside `_pybertini`

eigenpy exposes each decomposition through a **header-only `boost::python::def_visitor` template
parameterized on the matrix (scalar) type** — `eigenpy::PartialPivLUSolverVisitor<MatrixType>`,
`ColPivHouseholderQRSolverVisitor<…>`, `JacobiSVDVisitor<Eigen::JacobiSVD<…>>`, … — each with a
static `expose(name)`. Those templates instantiate at the *including* library's compile time.
`_pybertini` is compiled with knowledge of **both** eigenpy's headers **and** `complex_mp` /
`real_mp`, so it can do what the stock eigenpy `.so` cannot: instantiate those same visitors on
`Mat<complex_mp>` / `Mat<real_mp>` (`python_bindings/src/linalg_export.cpp`, one `expose("Name")`
line each). This **reuses eigenpy's binding code and Eigen's algorithm — it is not a reimplementation
and not a re-export.** Eigen's dense decompositions already compile and run on the mp scalars via the
`Eigen::NumTraits<complex_mp / mpfr_real>` specialization (`core/.../eigen_extensions.hpp`), proven
in the start-system code (`start/mhom.cpp`, `total_degree_linear_product.cpp`).

Use the explicit `expose(const std::string&)` overload — there is no `eigenpy::scalar_name`
specialization for the mp types, so the no-arg `expose()` (which builds a name from `scalar_name`)
would not compile.

### Home: a `bertini.linalg` submodule, because that is the honest home

The classes register only when `_pybertini` is imported, so they belong under `bertini.*`, not
injected into the `eigenpy` module namespace (which is fragile: import-order coupled, collision-prone,
and still requires importing bertini). The name is free: the old `bertini.linalg` (system-building
sugar, unrelated to linear algebra) was deleted in #64. Adding a new decomposition later is one
`expose()` line plus the header include (eigensolvers, `FullPivLU`, `BDCSVD`, `CompleteOrth…`).

### One dtype-agnostic surface (no type-iffing), with numpy adapters for double

The pure-Python `bertini.linalg` module routes by dtype so a single call site handles every scalar
type: `complex_mp` / `real_mp` → the mp path (eigenpy visitors); `float64` / `complex128` →
`numpy.linalg` (the right, fast tool there); integers → promoted to double. `solve` / `lstsq` return
a vector either way. The `lu` / `qr` / `svd` factories return an object with a **common method
surface for every dtype**: for mp it is the full eigenpy decomposition (`solve`, `determinant`,
`inverse`, `rank`, `singularValues`, `matrixU/V`, …); for double it is a thin numpy-backed adapter
(`_DoubleLU` / `_DoubleQR` / `_DoubleSVD`). The SVD adapter follows eigenpy's `A = U S Vᴴ` convention
(`matrixV()` returns `V`, not numpy's `Vh`). The double path uses numpy adapters rather than stock
eigenpy's double decompositions on purpose — eigenpy's double coverage is uneven (notably it does
**not** expose `JacobiSVD` for double), so numpy gives a uniform, predictable double surface.

### Companion: `bertini.precision(A, n)` is a functional (copy) setter

The related bulk-precision surface re-casts a whole vector/matrix's precision in one call, over
`{vector, matrix} × {complex_mp, real_mp}`, delegating to the existing C++
`Precision(container[, digits])`. It **returns a new array** (`A = bertini.precision(A, 20)`) rather
than mutating in place, because eigenpy marshals mp arrays by copy (see ADR-0031, ADR-0051): a
writable-`Eigen::Ref` setter would not write back to the caller's numpy buffer — the same
copy-semantics reasoning behind ADR-0008 (take writable-Ref binding outputs as numpy objects).

## Consequences

- Full-precision dense linear algebra in Python: `solve`, `lstsq`, `lu`, `qr`, `svd` on `complex_mp`
  / `real_mp`, verified in `python/test/classes/mp_linalg_test.py` (a 100-digit dyadic-rational solve
  is exact; QR/SVD residuals ~1e-20; singular-value product equals `|det|`; `U S Vᴴ = A`; exact QR
  rank on rank-deficient input) and `precision_bulk_test.py`. Interface-only tests, per the standing
  "correctness is a C++ gate" rule — the arithmetic is Eigen's and the mp scalars are exercised in
  the C++ tracker/start-system tests.
- Callers never branch on dtype. `bertini.linalg.solve(A, b)` works for an mp or a double `A`.
- A binding-only change (new `linalg_export.cpp`; `mpfr_export.cpp` for precision), no core C++ change.
- The double `lu`/`qr`/`svd` adapters expose only the common method subset; the mp objects carry the
  full eigenpy API (`matrixLU`, `permutationP`, `rcond`, `matrixQR`, …). Acceptable asymmetry — the
  shared surface is uniform.
- A 1-row matrix returned from `bertini.precision(M, n)` comes back 1-D `(N,)` — the known eigenpy
  1-row-Matrix→1-D behavior, not specific to this change.
