# ADR-0031: Indexing an eigenpy Vec returns an aliasing view — copy before storing

**Status:** Accepted
**Date:** 2026-06-24

## Context

`ZeroDim.to_dataframe()` (the pandas "database of solutions") built each row by reading the
solution out of the eigenpy solution vector and stashing it in a Python container that pandas
reads later:

```python
points = self.all_solutions(user_coords)      # an eigenpy container of Vec<mpc_complex>
row = {'x{}'.format(k): points[i][k] for k in range(len(points[i]))}   # WRONG — stores views
```

The resulting DataFrame was **silently wrong**: every row showed the *same* coordinate, as if one
solution had been copied over all of them. The underlying `all_solutions()` was always correct
(reading it directly, or re-reading it afterward, gave the right points); only the values *stored*
into the DataFrame were corrupt.

The cause is an eigenpy/NumPy detail. Indexing an eigenpy `Vec<mpc_complex>` from Python
(`pt[k]`, and likewise `list(pt)`) does **not** return an independent value — it returns a scalar
that **aliases a reused internal buffer**. While you read each element immediately the value is
correct, but if you *store the reference and read it later*, by then the buffer has been reused by
subsequent indexing, so all the stored references resolve to whatever was written last. Two
solutions never tripped it (no buffer reuse between two reads); three or more reliably did.

The same aliasing has a rarer, fatal face. When the reused buffer is not merely stale but *freed*,
the stored `mpc` carries a garbage limb pointer, and the first real arithmetic on it makes GMP read
a garbage limb count:

```
GNU MP: Cannot allocate memory (size=575897802350002184)   →  SIGABRT
```

This is the crash filed as **bertiniteam/b2 #259**. It was first mis-attributed to an
uninitialized-numpy-slot residual (the [ADR-0006](0006-eigenpy-uninitialized-numpy-slot-guards.md)
family) and to import order (numpy/pandas/matplotlib before bertini). Neither was the cause: import
order only changed what garbage occupied the buffer, and the slot guards were already in place. The
real cause is this aliasing, and #259 is closed by the fix below.

This is the same *class* of eigenpy hazard as [ADR-0001](0001-eigenpy-writable-ref-scalar-by-value.md)
(writable `Eigen::Ref` clobbering an adjacent scalar) and ADR-0006 (uninitialized `mpfr`/`mpc`
numpy slots): a Python-visible `mpc` value that is not the stable, fully-owned object it appears to
be.

## Decision

**Copy a value out of an eigenpy vector at the moment you read it, before you store it or read the
next element.** Never store the result of indexing (or `list()`-ing) an eigenpy vector and read it
later.

Two equivalent copies, both proven correct and precision-preserving (the copy keeps the native
element type — a Python `complex` for a double solve, a `bertini.multiprec.Complex` for a
multiprecision one, so nothing is truncated):

```python
# per element — make an independent copy of the right type, immediately
value = type(c)(c)            # where c = pt[k]

# per vector — snapshot the whole vector's buffer in one C-level copy
solution = points[i].copy()
```

`to_dataframe` uses the whole-vector form: each row's `solution` cell is `points[i].copy()`.

## Consequences

- **Correct and import-order independent.** Verified: pre-fix 15/15 runs silently corrupt and one
  observed SIGABRT under the worst import order; post-fix 60/60 clean. Import order no longer
  matters, because the hazard was never really about import order.
- **No precision loss.** Copying via `type(c)(c)` / `Vec.copy()` keeps the multiprecision element
  type; a multiprecision solve's coordinates round-trip exactly.
- **Negligible cost.** One copy per coordinate (or per solution vector) at table-build time.
- **Rule to follow:** any new Python code that pulls scalars or vectors out of an eigenpy
  `Vec`/`Matrix` and **keeps** them — builds a DataFrame, caches a list, returns them past the
  call — must copy at extraction. Reading-and-using-immediately (e.g. `complex(pt[k])` in a loop,
  or summing on the spot) is fine; deferring the read is the trap.
- **Regression test:** `python/test/zero_dim/solution_access_test.py::test_to_dataframe_solutions_are_independent_per_row`
  uses four distinct roots `(±1, ±1)` — an aliased frame would collapse them to one — and asserts
  the DataFrame reproduces `all_solutions()` exactly.
- **Relation to prior ADRs:** ADR-0001/0006/0008 hardened the *C++ binding* side of eigenpy `mpc`
  handling. This ADR is the *Python consumer* counterpart: even with correct bindings, a Python
  caller must not treat an indexed eigenpy element as a durable owned value.
