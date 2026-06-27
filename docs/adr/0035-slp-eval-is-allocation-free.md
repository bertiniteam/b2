# ADR-0035: The SLP eval loop is allocation-free (lower integer powers to multiplies; no by-value temporaries)

**Status:** Accepted
**Date:** 2026-06-27

## Context

The straight-line program (SLP) is the sole evaluator of polynomial systems (ADR-0027), and its
inner `Eval<NumT>` loop is the hottest code in the library — it runs once per Newton step, per path,
per homotopy. At multiprecision the cost is dominated by `mpfr`/`mpc` arithmetic, and for years the
working assumption was that the arithmetic itself was the floor.

It was not. A counting GMP allocator (`mp_set_memory_functions`) over the eval loop showed the real
cost was **heap allocation churn**: a *single* evaluation of function + Jacobian for a tiny dense
4-variable, degree-≤5 system did **~1132 `malloc`/`free` pairs** at 256 digits (and ~2874 at 1024),
with `malloc` exactly balanced by `free` — pure churn. Per-operation probes isolated three sources,
none of them the arithmetic:

1. **`pow` for integer powers.** `pow(mpc, int)` heap-allocates ~14 temporaries per call; worse, the
   parser builds `x^k` as a general `PowerOperator`, so it was emitting `pow(mpc, mpc)` — the
   *transcendental* `exp(k·log x)` — at ~40 µs/call at 256 digits (and ~720 µs at 1600). A complex
   multiply into a preallocated slot, by contrast, is **0 allocations**.
2. **A differentiation bug.** `PowerOperator::Differentiate` built the new exponent as an
   *unevaluated* `exponent - 1` node, so `d(x^4)/dx` was `4·x^(4-1)` — a power with a *computed*
   exponent that the SLP could not recognize as integer, falling back to the allocation-heavy general
   `pow`. Jacobians were the worst offenders.
3. **By-value temporaries in the hot loop.** The eval's `Assign` and `Negate` were implemented with a
   by-value `return x` / `return -x` lambda, which materializes a fresh `mpc` temporary (an
   allocation) instead of writing in place. (Arithmetic `+ - * /` were already allocation-free via
   Boost expression templates evaluating into the destination slot.)

A further trap was measured and must be remembered: **aliased multiplication allocates.** `c = a*a`
(same operand twice) costs ~6 heap ops because `mpc` defensively allocates a temporary for the
squaring, whereas `c = a*b` (distinct operands) is 0.

## Decision

**The SLP eval loop must perform zero heap allocations.** Concretely:

1. **Lower integer powers to multiplications at compile time.** `Visit(IntegerPowerOperator)` and
   `Visit(PowerOperator)`-with-integer-exponent emit a chain of `Multiply` instructions using
   **exponentiation by squaring** (O(log n) multiplies). The squarings are kept allocation-free by
   copying the running square into a distinct slot before each square, so every multiply has distinct
   operands. The SLP's hash-consing shares repeated powers across terms. No `IntPower`/`Power`
   instruction is emitted for an integer exponent; the general `Power` opcode remains only for
   genuinely symbolic/non-integer exponents.
2. **Differentiation folds constant integer exponents to literals.** `PowerOperator::Differentiate`
   emits `n·x^(n-1)` with `n-1` a literal `Integer` when the exponent is an integer literal, so
   derivatives lower like any other integer power.
3. **No by-value temporaries in eval.** Results are written in place into the preallocated register
   slots: arithmetic uses expression templates; `Assign`/`Negate` are direct slot operations
   (`cplx[o] = cplx[i]` / `-cplx[i]`), not a value-returning lambda. Promotions and the rare
   transcendentals are the only remaining temporary constructions, and they must construct at the
   working precision (see ADR-0034's thread-precision pin).

No transcendental fallback for large exponents is warranted: the crossover sweep (30–1600 digits)
shows `pow(complex,complex)` always loses to the lowered multiplies, and by a wider margin at higher
precision (transcendentals scale worse than multiplies).

## Consequences

- **~10× faster multiprecision eval**, and **0 allocations per eval at ≤256 digits** (1132 → 0).
  This is what makes b2's evaluation competitive at high precision.
- **A standing invariant to defend.** The three traps — calling `pow`, multiplying aliased operands
  (`a*a`), and returning mpc/mpfr *by value* — each silently reintroduce per-op allocation. New eval
  code and new function-tree → SLP lowerings must avoid them. The probes that caught this live in
  `core/test/classes/slp_test.cpp` (`mpfr_alloc_churn_per_eval`, `mpfr_alloc_per_raw_op`,
  `power_method_crossover`, opt-in via `BERTINI_SLP_BENCH`).
- **Residual at extreme precision.** At ~1024+ digits some `realloc`/`malloc` remains (≈198/eval at
  1024) from `mpfr` limb growth; a pooling GMP allocator would mop that up but is process-global and
  thread-/wheel-sensitive, so it is deferred.
- **Follow-ups:** `Assign` *elision* (copy propagation / slot coalescing) would cut instruction count
  further (the lowering and output wiring emit copies). JIT (ADR-0034 / task) removes the interpreter
  entirely.
- **Builds on** ADR-0027 (immutable Program + preallocated per-thread Memory) and ADR-0034 (tiered
  arithmetic; the thread-precision pin and the per-slot bank dispatch).
