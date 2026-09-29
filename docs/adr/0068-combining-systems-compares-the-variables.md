# ADR-0068: Combining systems compares the variables, and a loaded system holds the live ones

**Status:** Accepted
**Date:** 2026-09-29

## Context

Variables are canonical by name: `Variable::Make("x")` interns, so the `x` made here and the
`x` made there are one object.  That is what makes the symbolic engine cheap, and it gives a
plain meaning to "these two systems are over the same variables": the same objects, in the
same groups, in the same order.

Two things contradicted it.

**An archive handed back variables that were not interned.**  Boost deserialization builds a
node through its default constructor, never through `Make`.  `copy.deepcopy`, `copy.copy`,
`pickle`, and the MPI broadcast all load from an archive, so each returned a system
content-equal to the original and sharing none of its nodes.  Such a system could not be
concatenated with the one it was copied from ("differing variable orderings").  The repair
existed -- `System::ReinternNodes`, used by the records loader -- and nothing else called it.

**The operations that combine systems asked two different questions.**  `Concatenate`
compared the variable orderings.  `MakeHomotopy`, `MakeMovingHomotopy` and `+=` compared
counts: of variables, of homogenizing variables, of groups.  The homotopy builders evaluate
their operands as whole systems, each fed the point by position, so they never needed the
variables to be shared -- and never noticed when they were not the same.  Measured: moving
rows over `(y, x)` blended with a fixed system over `(x, y)` were accepted, and `y - 3` at
`(x, y) = (7, 11)` evaluated to 4.

The two met in `straight_line_homotopy`, which calls the moving-row builder and then
concatenates to assemble `.target` and `.start`: it refused operands the builder underneath
it accepted.  The cellular decomposition port, which deep-copies the three systems of its
projective move, found it.

## Decision

1. **Loading re-interns.**  `System::serialize` and `Slice::serialize` re-intern what they
   hold as the last step of a load.  It is done in the load, not asked of each caller, so a
   system that reaches anyone is over the live variables.  Re-interning changes no content
   (the digest is invariant) and drops the derived caches, which recompute on demand.

2. **One check, on the variables.**  `CheckVariableStructuresMatch` compares variable by
   variable: the ordering, the affine and projective groups, the ungrouped variables, the
   homogenizing variables.  `Concatenate`, `MakeHomotopy`, `MakeMovingHomotopy` and `+=` all
   make it.  The refusal names both structures and says which kind of difference it is: a
   different order, a different grouping, different variables.

3. **A system that declares no variables takes the other's** (`AdoptVariableStructure`).
   Functions with nothing said about variables are read as being over the variables of what
   they are combined with.  A system that declares any variables is compared, never adopted.

4. **In Python a node compares by identity**: `__eq__`, `__ne__` and `__hash__` on the node
   wrappers compare the node, not the wrapper.  Each call across the binding makes a new
   wrapper, and the default `==` said two handles on one node were different.

## Consequences

- A deep copy, a shallow copy and an unpickled system combine with their original, and with
  each other.  `clone` remains the direct way; it needs no archive.
- Systems over the same variables in a different order, or grouped differently, are refused
  by the homotopy builders and `+=`, which used to accept them.  That is a breaking change for
  code that relied on matching by position, and such code was computing with its coordinates
  exchanged.
- **Do not go back to counting** in `CheckVariableStructuresMatch`, and do not give an
  operation its own comparison.  Counts agree for `(x, y)` and `(y, x)`.
- **Do not remove the load hooks** on the grounds that a caller can re-intern: the callers
  that did not are the reason for this record.  A new class that holds nodes and can be
  archived needs the same last step.
- **Do not build a node outside `Make`.**  The protected default constructors exist for the
  archive and for nothing else.
- Left alone on purpose: a function that mentions a variable no group holds.  It is accepted
  at construction, the degree questions read the stray variable as a constant, and
  evaluation refuses ("a function references the variable 'y', which is not in the system's
  variable ordering").  The rule for it -- every variable in a function has a declared role
  in its system -- belongs with parameters, which are a role that does not exist yet.  Users
  may add functions and groups in any order, so that check cannot run at `AddFunction`.
- Pinned in `core/test/classes/system_identity_test.cpp` (the "combining systems" cases) and
  `python/test/classes/clone_concatenate_test.py`.
