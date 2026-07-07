# ADR-0038: AMP Criterion B (cost model) must use the latest Newton residual, not `size_proportion`

**Status:** Accepted

## Context

Benchmarking the b2 CLI against Bertini 1.7 on adaptive-multiprecision (AMP) zero-dim solves showed
b2 ~55× slower on cyclic5 (11.0s vs 0.20s), *after* the SecurityLevel fix (ADR-0037). Per-path
inspection found **39 of 70 finite paths needlessly escalating from double to 30–40 digit
multiprecision in the endgame**, despite benign condition numbers (7–30). The same 70 solutions are
found in pure double (`mptype 0`) in 0.36s, so double is numerically sufficient — the escalation was
spurious.

Reconstructing the *byte-identical* homotopy (same gamma, same start points) in Bertini 1 via a
`UserHomotopy:1` input confirmed b1 tracks **all 70 paths entirely in double** (`Maximum precision
utilized: 52` bits, "first precision increase: 0.0" for every path). So b1 — the authoritative port
reference — never escalates here.

Instrumenting every precision change and attributing it to its trigger showed the dominant cause was
the per-step AMP cost model (`AMPTracker::AdjustAMPStepSuccess` → `MinimizeTrackingCost`) raising
precision on *successful* steps, because its Criterion B requirement `digits_B` was ~30 even for
well-conditioned matrices. Decomposing `B_RHS`:

```
B_RHS = safety_digits_1 + D(normJ, normJinv) + (-log10(tol) + log10(residual)) / max_newton_its
```

with measured `D ≈ 4` (innocent), the ~30 came entirely from the residual term, because the value
passed as the "latest Newton residual" was **`last_step_.size_proportion ≈ 1e51`**.

`size_proportion = err_est / |delta_t|^(p+1)` (the predictor's error *proportionality constant*).
In the endgame, `delta_t` is tiny and `err_est` is floored at double roundoff, so `size_proportion`
blows up to 1e51–1e53. Passed raw into Criterion B, `log10(size_proportion)/N ≈ 25` dominated the
RHS → `digits_B ≈ 30` → `min_precision = 30` → escalate.

Per AMP3 (Bates–Hauenstein–Sommese, *Adaptive multiprecision path tracking*, ODE refinement) and the
existing in-corrector usage, Criterion B's residual argument is the **latest Newton residual**: during
correction that is ‖Δz‖; during prediction it is the local truncation error estimate (= `size_proportion
* stepsize^(p+1)`, NOT `size_proportion` raw). The corrector's own `CriterionB` call
(`newton_corrector.hpp`) already passes `meta.norm_delta_z` correctly. Only the tracker cost-model
`B_RHS` substituted `size_proportion`.

## Decision

In `AMPTracker::B_RHS` (`core/include/bertini2/trackers/amp_tracker.hpp`), pass
`last_step_.norm_delta_z` (the latest Newton residual ‖Δz‖) to `amp::CriterionBRHS`, instead of
`last_step_.size_proportion`. This makes the cost-model Criterion B consistent with the corrector's
Criterion B and with AMP3.

## Consequences

- cyclic5 (serial, `mptype 2`, `randomseed 1`): **11.0s → ~3.5s**, paths escalating **39/70 → 2/70**,
  all 70 finite solutions unchanged (bijection vs the Bertini 1 oracle, worst coordinate difference
  2.6e-12). C++ tests (`test_endgames`, `test_nag_algorithms`, `test_tracking_basics`) pass.
- Pure double (`mptype 0`) is unaffected: `B_RHS` is only used by the AMP tracker.
- The remaining ~3.5s-vs-0.20s gap is separate follow-ups, not this fix: (1) AMP-at-double per-step
  overhead (condition probe solve + criteria every step even when staying double); (2) the endgame
  sample-refinement transient churn (refining each sample to `final_tolerance × 1e-2 = 1e-13`); and
  (3) two genuinely-deep-endgame paths (track to t≈1e-6/1e-7 where conditioning transiently degrades).
- Note (silvi): the b2 arithmetic cost model is calibrated higher for multiprecision than b1's
  (`ArithmeticCost`: ~`101.47 + 1.59 P` vs AMP2's `10.35 + 0.13 P`), so b2's preference for double
  should be at least as strong as b1's; this fix removes a spurious floor that defeated that preference.
