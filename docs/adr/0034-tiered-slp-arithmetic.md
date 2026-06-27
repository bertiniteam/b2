# ADR-0034: Tiered numeric arithmetic in the SLP (integer / real / complex)

**Status:** Accepted (design; implementation staged and measurement-gated)
**Date:** 2026-06-27

> **Implementation status (2026-06-27):** **Stage 1 (real tier) implemented and landed.** The
> `Float`→`Complex` node rename (PR #37) and the scalar-name standardization (PR #38) are the
> foundation. Stage 2 (integer tier) is **deferred** — see the measured outcome below. JIT (the lever
> for competing with HomotopyContinuation.jl) is captured separately and deferred unless it proves
> straightforward.

## Measured outcome (Stage 1)

The real tier is correct (tiered eval is bit-identical to all-complex; verified by the full C++ suite
plus `SLP_tiered_numtype`) and a **modest, no-regression win**: ~10% at mpfr (256 digits) and ~8% at
double, with no fully-complex regression even at 4096 digits.

A raw-op probe (`mpfr_raw_op_microbench`) settles *why* the win is only modest — and rules out two
suspected culprits:
- Boost does **not** promote: `real_mp * complex_mp` costs 0.33–0.45× a complex×complex multiply, and
  `real_mp * real_mp` 0.14–0.21×. So the tier genuinely makes the ops it touches 2–3× cheaper.
- The hot loop is **not** allocation-bound: the register banks are pre-allocated and expression
  templates evaluate into the destination slot.

The ceiling is **structural**: in a polynomial over *complex* variables, the monomial arithmetic
(`x⁵`, `x·y·z`) is irreducibly complex and dominates per-eval cost, while the real-valued coefficients
are frozen constants computed once in the prologue. Only the coefficient×monomial multiplies (~1 in
`degree`) go cheap, so the aggregate per-eval win is bounded. The genuine "much faster" levers are
therefore **reducing the count of complex multiplies** (Horner / stronger CSE) and **JIT** (removing
interpreter overhead, biggest at double) — not more tiers. Hence Stage 2 (integer tier, expected to
hit the same ceiling for even less of the work) is deferred.

## Context

The SLP is the sole evaluator of polynomial systems (ADR-0027). It evaluates **everything in
complex**: the register file is `tuple<vector<dbl_complex>, vector<mpfr_complex>>` — one bank per
precision, both complex — and even an integer or rational coefficient is widened to a complex slot
(`ConstantRecipe::Produce`) before any arithmetic touches it.

`mpfr_complex` arithmetic is the most expensive thing the library does: a complex multiply is ~4 real
multiplies + 2 adds, each an arbitrary-precision mpfr operation. Yet polynomial systems are full of
cheaper structure — integer/rational/real coefficients, and coefficient×variable products where one
operand is narrow. Evaluating those in the **narrowest sufficient numeric type** and widening to
complex only when an operation demands it should be a large speedup at high precision, where the
`mpfr_complex` cost dominates everything else.

There is precedent: the existing **`IntPower`** opcode reads its exponent from a separate `integers_`
bank and dispatches to integer-power `pow` — a cross-tier operation, which this design generalizes.

The counter-risk, raised explicitly: a mixed-type tape adds register banks and per-operation tier
dispatch, which can **thrash cache and lose locality**, potentially eating the arithmetic savings.
So this is **measured, not assumed** — every tier must demonstrate a net speedup before it ships.

## Decision

Introduce a **tier lattice** `Integer ⊂ Real ⊂ Complex`, **orthogonal to precision**.

- **Precision stays the `Eval<NumT>` template parameter** (double vs mpfr), preserving the AMP
  template invariant: the numeric type is not smeared into runtime state. **Tier is a per-slot
  runtime property *within* a precision.** So `Eval<mpfr_complex>` works over mpfr banks
  `{mpz_int, mpfr_float, mpfr_complex}` and `Eval<dbl_complex>` over `{int, double, dbl_complex}`.

1. **Per-slot tier + per-tier banks.** The single complex bank becomes three banks per precision
   (integer/real/complex). Every slot carries a tier tag; slots are allocated per-tier in the
   compiler.

2. **Compile-time tier inference** in `SLPCompiler`. A subexpression's tier is the *join* of its
   operands' tiers, with op-specific escapes:
   - `+ - *` → join(operands).
   - `/` → escapes to at least real (integer/integer is not integer).
   - `sqrt, log, ^non-integer` → complex unless the operand is provably non-negative real (default
     complex, to stay correct on negative reals).
   - `IntPower` → base tier. `exp, trig` → operand tier.
   Leaves seed the lattice: `Kind::Integer→Integer`; `Rational→Real` (Integer if unit denominator);
   `Complex→Real` when `imag==0` else `Complex`; `Pi/E→Real`; variables and the path variable →
   `Complex` (homotopy unknowns are complex).

3. **Explicit promotion (widening) instructions.** Where an operation needs an operand at a higher
   tier than its slot, the compiler emits a `Promote` instruction (int→real, real→complex) into a
   target-tier slot. Widening is explicit and cheap, so eval stays branch-simple.

4. **Tier-aware eval dispatch.** `Add/Subtract/Multiply/Divide` carry their operand+result tiers,
   compactly encoded so the combinatorial blow-up stays bounded; the hot mixed case (narrow×complex)
   gets a dedicated fast path (`real*complex` = 2 real muls, not 4). The pure-complex path is kept
   byte-identical to today, so a fully-complex system never regresses.

5. **`IntPower` is folded in** as the Integer-tier-operand special case it already is.

## Staging (each stage gated by tests + benchmark)

- **Stage 0 — this ADR + scalar-name standardization.** One real + one complex type per precision,
  with an explicit `using` for the double real type (changeable in one place); the vocabulary the
  banks and inference are written against.
- **Stage 1 — Real tier.** Real bank, real↔complex promotion, tier-aware `+ - * /` for `{real,
  complex}`. The common, high-value, lowest-combinatorial case. **Gate:** equivalence holds and the
  benchmark shows a net mpfr speedup with no fully-complex regression.
- **Stage 2 — Integer tier.** Integer bank, int→real/complex promotion, integer fast paths, `IntPower`
  folded in. **Same gate.** If a tier does not pay off in measurement, it does not ship.

## Consequences

- **Speedup** on integer/rational/real-coefficient systems at high precision (the goal).
- **Correctness invariant:** tiered eval must equal all-complex eval to rounding — verified by a
  battery of systems at double and several mpfr precisions, for functions and Jacobian (the analogue
  of "threaded == serial").
- **Serialization format bump:** the tape and memory layout change; the SLP archive version is bumped
  and old dev archives will not load (acceptable pre-release; consistent with the node rename).
- **Complexity cost:** a larger eval switch and per-slot tier bookkeeping. Mitigated by keeping the
  complex path unchanged, by explicit (not implicit) promotion, and by the measurement gate — a tier
  that loses to cache effects is dropped rather than shipped.
- **Builds on** ADR-0027 (immutable Program + per-thread Memory; the freeze-set prologue) and the
  `Float`→`Complex` rename; precedent from the existing `IntPower` cross-tier op.
