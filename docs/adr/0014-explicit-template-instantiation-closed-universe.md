# ADR-0014: Explicit template instantiation of the endgame/ZeroDim universe

**Status:** Accepted
**Date:** 2026-06-12

## Context

The tracker/endgame/algorithm stack is header-only, so every consumer TU
re-instantiated the entire dependency cone: each of the six zero_dim/endgame
python-binding TUs (14-line stubs!) cost up to 100 s / 5.9 GB peak RSS to
compile, and blackbox's `algorithm_builder.cpp` (which instantiates the
tracker×endgame switch ladder from `blackbox/switches_zerodim.hpp`) cost
141 s / 7.1 GB. This is the reason ADR-0004 pins CI wheel builds to
`CMAKE_BUILD_PARALLEL_LEVEL=2` and local builds run reduced parallelism.

The instantiation universe is **closed and small**: number types
`{dbl, mpfr_complex}`; trackers `{DoublePrecisionTracker,
MultiplePrecisionTracker, AMPTracker}`; endgames `{PowerSeries, Cauchy}` ×
those trackers' precision policies; ZeroDim = those six with `System`,
`TotalDegree`, and the default `CloneGiven` policy — exactly what the python
bindings and blackbox construct.

Precedent: `ExplicitRKPredictor` already used this pattern
(extern block at the bottom of explicit_predictors.hpp + one instantiating TU).

## Decision

Instantiate the universe **once, in libbertini2**:

- `core/src/eti/endgames_eti.cpp` — the 6 endgame types AND their
  `EndgameBase<Flavor, PrecT>` bases (explicitly instantiating a derived class
  does NOT instantiate base members).
- `core/src/eti/zero_dim_eti.cpp` — the 6 ZeroDim combos.
- Matching `extern template` declarations live at the bottom of the umbrella
  headers consumers actually include: `bertini2/endgames.hpp` and
  `bertini2/nag_algorithms/zero_dim_solve.hpp`.

**The rule when extending:** a new tracker, endgame flavor, or
production ZeroDim combo must be added to BOTH lists (the eti cpp and the
extern block). Forgetting the cpp → loud linker error. Forgetting the extern →
silently slow consumer TUs again (no breakage). Combos outside the universe
(other start systems, `RefToGiven`) instantiate implicitly as before —
graceful degradation.

ETI instantiates *every* member of a class, including never-called ones —
which immediately exposed three latent defects in dead endgame API
(`EndgameBase::ChangePrecision` calling flavor methods that never existed —
deleted; `SetFinalTolerance` assigning through a const accessor with a
mismatched parameter type — fixed; `prec_base` calling a member on a
`reference_wrapper` without `.get()` — fixed). This is a feature: the ETI TUs
now compile-check the whole public surface of the universe on every build.

## Measurements (aarch64 dev machine, GCC, Release)

Consumer-TU compile cost, wall / peak RSS, before → after:

| TU | before | after |
|---|---|---|
| zero_dim_amp_export | 101 s / 5.9 GB | 39 s / 3.9 GB |
| zero_dim_double_export | 79 s / 5.4 GB | 30 s / 3.4 GB |
| zero_dim_mp_export | 89 s / 5.7 GB | 32 s / 3.6 GB |
| endgame_{amp,double,mp}_export | 17–34 s / 2.4 GB | 16–24 s / 2.4 GB |
| blackbox algorithm_builder | 141 s / 7.1 GB | 100 s / 6.8 GB |

One-time ETI TU costs (parallelize within the libbertini2 build; consumers
compile concurrently and only wait at link): endgames_eti 29 s / 2.2 GB,
zero_dim_eti 109 s / 5.1 GB, zero_dim_blackbox_eti 116 s / 5.6 GB.

Reading: the binding zero_dim TUs drop ~2.6x in time and ~2 GB each in peak
RSS — the memory ceiling is what gated build parallelism.  The endgame export
TUs' residual cost is Boost.Python machinery, not endgames.
algorithm_builder's residual is its own switch-ladder breadth and io/parsing
includes.  Total CPU for a cold full build is roughly neutral (the work moved,
once, into the library); incremental rebuilds touching consumers are the big
winners.

## Consequences

- Bindings/blackbox/test TUs become "declare + link"; total build time and the
  peak-memory ceiling drop. Revisit ADR-0004's `CMAKE_BUILD_PARALLEL_LEVEL=2`
  once CI numbers confirm (follow-up; compounds with the ccache work, PR #15).
- Candidates deliberately not ETI'd yet: tracker member-template family
  (consolidated into the eti TUs as a side effect; method-level ETI only if
  future measurement warrants), NID (framework-only, kernels throw),
  blackbox user-homotopy start variants.
- `extern template` does not suppress what the optimizer chooses to inline;
  gains are concentrated in the big out-of-line loop bodies, which is where
  the cost was.
