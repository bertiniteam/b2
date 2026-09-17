# ADR-0012: Precedence-aware printing; printed trees must re-parse to the same values

**Status:** Accepted; amended 2026-09-17 by ADR-0059 (printing is a family of dialects, the
complex form is `(re+im*I)`, and the reparse invariant is per dialect)
**Date:** 2026-06-11 (PR #7)

## Context

Printing wrapped every operator unconditionally — `(((x^2)+((2*x)*y))-1)` — which made
trees unreadable and derivative output absurd. Dropping parentheses is only safe if it
can never change meaning: the printed form is fed back through the classic-format
parser (round-trip tests, user-written input files derived from printed systems).

## Decision

- Nodes report a printing precedence (`PrecSum < PrecNegate < PrecMult < PrecPower <
  PrecAtom`, node.hpp). A printing parent wraps a child **only when its precedence is
  too low for the position it occupies**: subtracted/divided operands group their own
  kind (`x-(y+z)`, `x/(y*z)`); `^` wraps any non-atom on either side (keeping
  `x^y^z` unambiguous); negation wraps sums and other leading-`-` printers (never
  emits `--`). Function-call operators (`sin(...)` etc.) and complex pairs
  self-delimit as atoms.
- **Literal constants participate in precedence by their printed shape, not their
  type.** Negative literals report `PrecNegate` (they print a leading `-`).
  Real-valued Rationals print bare as `p/q` — textually a division — so they report
  `PrecMult`: `x/(1/3)`, never `x/1/3`, which would re-parse as `x/9`. Real-valued
  Floats print bare; genuinely complex constants print as the self-delimiting
  `(re+im*I)` / `(re-im*I)` -- the one spelling Bertini 1 reads (amended by ADR-0059; the
  original `(re,im)` pair form is gone from every dialect).

## Consequences

- The load-bearing invariant, enforced by `printing_test.py`'s round-trip test
  (print → `bertini.parse` → identical evaluation): **removing parentheses must never
  change the parsed value.** Any printer change must keep that test green.
- A node whose printed form is textually an expression (contains an operator
  character) must report the corresponding precedence, not `PrecAtom` — the rational
  literal is the canonical example; a future literal type with a composite print form
  needs the same care.
- The exact printed forms are now a documented surface (tests assert them); cosmetic
  printer changes are API-visible.
- (ADR-0059) The invariant holds *per dialect*: Classic text (`to_classic()`,
  `to_classic_input()`) reparses through the classic parser; the exact Python text (`repr`)
  rebuilds through `eval`.  Printed text is never part of the content digest.
