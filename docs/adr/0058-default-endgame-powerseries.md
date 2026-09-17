# 0058 — The default endgame is power series, matching Bertini 1

**Status:** Accepted
**Date:** 2026-09-16 (decided 2026-07-15; landed with 3.5.0)

## Context

b2 defaulted to the Cauchy endgame everywhere (the `EndgameChoiceConfig` default, the
classic-input parser's no-`endgamenum` fallback, and the Python `ZeroDimSolver` /
`HomotopySolver` / `solve` / `parameter_sweep` factories).  Bertini 1's **documented
default is the power-series endgame** (`EndgameNum: 1`); the Cauchy default here was an
acknowledged mistake, and the in-code parser comment even mis-stated B1's default as
Cauchy.

The cost of the divergence was measured, not hypothetical.  On a regeneration-cascade
NID workload (the manydims system, library defaults, adaptive tracking):

- Cauchy endgame: **190 s**.  Power series endgame: **2.5 s** — equal to Bertini 1's
  wall time on the same system, seed, and (adaptive) tracking mode.
- Root cause of the gap: slowly-diverging paths (fractional-order blowup, e.g.
  ‖x‖ ~ t^(−1/3)) are a Cauchy-specific pathology: the at-infinity truncation is an
  absolute norm threshold reached only after many decades of t, the stepsize hits its
  floor and stops scaling with t (one instrumented path: 153,389 endgame steps, 70
  digits), and the Cauchy security valve is gated to the pole-mass operating zone,
  which a slow diverger never enters — so it never arms.  The power-series endgame's
  time-sequence sampling truncates the same paths early.
- Bertini 1's own Cauchy endgame is ~11× slower than its power-series endgame on
  cyclic-6 (68 s vs 6.2 s), while both find all 156 solutions.  Even for B1, Cauchy is
  the specialist tool, not the default.

## Decision

Power series is the default endgame at every choice point: `EndgameChoiceConfig`
(hence the classic parser's no-setting fallback and the blackbox CLI), and the Python
factory defaults (`ZeroDimSolver`, `HomotopySolver`, `user_homotopy`, `parameter_sweep`,
`bertini.solve`).  The Cauchy endgame remains fully available by explicit request
(`endgamenum: 2`, `endgame='cauchy'`).

## Consequences

- Do **not** flip the default back to Cauchy for a workload where Cauchy performs
  better — pass `endgame='cauchy'` there instead; the default follows the
  B1-defaults doctrine.
- Runs relying on library defaults get different ask digests than before (the
  encoded `EndgameChoice` value changed, and the endgame kind is part of the ask).
  This is correct: it is a different computation.  No encoding-version bump is
  involved — the canonical encoding functions are unchanged, and the config keyspace
  test pins its `EndgameChoiceConfig` recipe to Cauchy so it keeps measuring the encoder
  rather than the default.
- The Cauchy slow-diverger pathology is still worth fixing on its own (planned:
  B1-style truncation on two consecutive loop-ROOT norms above `max_norm`,
  zone-free; the operating-zone gate is under re-evaluation with cyclic-6-via-B1
  as evidence).  This ADR only stops it from being the *default* cost.
