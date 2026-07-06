# ADR-0021: A block (and the patch) must fully define its own output rows

**Status:** Accepted
**Date:** 2026-06-16

## Context

In the block-composed `System`, evaluation results are assembled into a caller-allocated buffer:
`System::Eval`/`Jacobian`/`TimeDerivative` allocate `Vec<T>(NumTotalFunctions())` /
`Mat<T>(NumTotalFunctions(), NumVariables())` and then let each block write its contiguous slice
(`EvalBlocksInPlace` / `JacobianBlocksInPlace` hand each block a `segment`/`block` of the buffer).

**Eigen does not zero-initialize** `Vec<T>(n)` / `Mat<T>(r,c)`. So every output entry is garbage
until something writes it. The contract is therefore: **each block must fully define every entry
of the rows it owns** — not just the "interesting" entries.

This was violated by the patch. `Patch::JacobianInPlace` wrote only its *sparse* per-variable-group
coefficient entries and left the rest of its rows untouched, relying on the caller having handed it
a pre-zeroed matrix. The block-composed Jacobian path does not pre-zero. The unwritten patch-row
entries were then read uninitialized:

- benign on Linux/macOS, where fresh allocations usually sit on zeroed pages;
- **garbage on Windows**, where a degree-2 MHom blend homotopy's Jacobian "evaluated" to ~`4.45e252`,
  failing `mhom_homotopy_block_matches_function_tree` and (the strong hypothesis) inflating the AMP
  condition-number estimate enough to drive the multi-hour Windows MHom precision "grind".

It presented as a gamma/precision problem and a Windows-only flake; it was neither. It was
uninitialized memory — undefined behaviour on every platform, exposed only where the heap was dirty.

## Decision

**Every block's `EvalInPlace` / `JacobianInPlace` / `TimeDerivInPlace` must fully define all entries
of its slice** — sparse contributors zero their slice first, then write their non-zero entries.

- `Patch::JacobianInPlace` now does `jacobian.block(offset, 0, NumVariableGroups(), cols).setZero()`
  before writing the per-group coefficients (commit 574ddb56).
- `BlendBlock` and `LinearFormsBlock` already `setZero()` / fully assign their targets; this ADR
  makes that an explicit, named contract rather than an accident of each block's implementation.
- Callers must **not** assume the buffer is pre-zeroed, and conversely must not skip writing "the
  obviously-zero" entries.

A regression test (`patch_test.cpp::patch_jacobian_fully_defines_its_rows_into_a_dirty_buffer`)
hands `JacobianInPlace` a buffer poisoned with `1e300` and requires every owned entry to be correct
— catching this class of bug deterministically on any platform, not just where the heap is dirty.

## Consequences

- Windows CI's `C++ tests on Windows` job is green again; the fixed-double MHom test was hardened to
  assert it solves (the prior `try/except RuntimeError` guard had been masking this exact bug).
- This is the load-bearing reason to fold the patch into a `PatchBlock`: the block contract's
  value-in/value-out "fully defines its own rows" rule makes the bug structurally impossible, rather
  than something each ad-hoc post-block writer must remember. See the planned patches-as-blocks work.
- The same rule applies to any future block (slice, randomization, …): zero your slice, then fill it.
- Do **not** "fix" a future incarnation of this by zeroing the whole result buffer in `System`
  before the block loop — that hides the contract violation and pays an O(rows·cols) memset on every
  evaluation. Each block owning its rows is both correct and cheaper.
