# ADR-0015: Solution metadata classification (finite/real/singular) matches Bertini 1

**Status:** Accepted
**Date:** 2026-06-14

## Context

`SolutionMetaData` declared `is_finite`, `is_real`, and `is_singular`, but
`ComputePostTrackMetadata()` only computed multiplicities — the three flags were
never assigned, so they were always `false`. `PostProcessingConfig::endpoint_finite_threshold`
was parsed from classic input but **never applied** anywhere. So a caller asking
"is this solution finite?" got `false` for genuinely finite points, and "use the
configured cutoff" was impossible — there was nothing computed to use.

Two further problems surfaced:

- The default `endpoint_finite_threshold` was **inverted**: `1e-5`. An endpoint is
  at infinity when its norm is *large*, so the cutoff must be large (Bertini 1's
  default is `1e5`).
- `ComputeMultiplicities` compared the **internal** (homogenized, on-patch)
  coordinates with the **2-norm**. Internal coordinates carry the homogenizing
  variable and the patch scaling, so two projectively-identical endpoints can
  differ there; the comparison must be on **dehomogenized** coordinates (see
  ADR-0013), and the endgames consistently use the **infinity norm**.

Bertini2 is a re-implementation of Bertini 1, so classification must match Bertini 1's
behaviour except where a deviation is necessary and documented.

## Decision

In `ComputePostTrackMetadata()`, classify each **endgame-successful** endpoint,
dehomogenizing it once:

- **`is_finite`**: `infNorm(DehomogenizePoint(pt)) <= endpoint_finite_threshold`.
  A path the endgame already flagged divergent (`SuccessCode::GoingToInfinity` or
  `SecurityMaxNormReached`) is taken as infinite with **no recomputation**, so the
  metadata can never contradict the endgame's own verdict.
- **`is_real`**: `infNorm(imag(DehomogenizePoint(pt))) < real_threshold`.
- **`is_singular`**: `multiplicity > 1` **or** the spectral-norm condition-number
  estimate `> condition_number_threshold`.
- **`multiplicity`**: cluster endpoints whose dehomogenized coordinates agree to
  `final_tolerance * same_point_tolerance_multiplier` in the **infinity norm**.

The dehomogenize-then-infinity-norm measurement is factored into
`System::InfinityNormOfDehomogenized` and used by **both** the classifier and the
endgames' `Security::max_norm` divergence check (`cauchy.hpp`, `powerseries.hpp`),
so the two can never drift on what "going to infinity" means.

`PostProcessingConfig` is corrected to Bertini 1:

| Bertini 1 name | bertini2 field | default | meaning |
|---|---|---|---|
| `ImagThreshold` | `real_threshold` | 1e-8 | real if `infNorm(imag(dehom)) <` this |
| `EndpointFiniteThreshold` | `endpoint_finite_threshold` | **1e5** | at infinity if `infNorm(dehom) >` this |
| `EndpointSameThreshold` | `same_point_tolerance_multiplier` | **10** | a **multiplier** (≥1) on `final_tolerance` |
| `CondNumThreshold` | `condition_number_threshold` | **1e8** | singular if cond# `>` this |

`same_point_tolerance` was renamed to `same_point_tolerance_multiplier` because it
is a multiplier on `final_tolerance`, *not* an absolute tolerance: the same-point
test must stay a fixed factor looser than the accuracy you tracked to, even when
`final_tolerance` changes.

## Consequences

- **The thresholds are settings — never hardcode a finite/real/singular cutoff** in
  client code (Python examples, tutorials, downstream tools). Read `is_finite` /
  `is_real` / `is_singular` / `multiplicity` from the metadata, or the threshold
  from `PostProcessingConfig`.
- **Compare endpoints in dehomogenized (user) coordinates with the infinity norm.**
  Never compare raw internal coordinates (the homogenizing-variable + patch-scaling
  trap).
- When adding a new `SuccessCode` that means "diverged," add it to the
  `is_finite` short-circuit list so the metadata keeps agreeing with the endgame.
- Tests: `core/test/nag_algorithms/zero_dim.cpp` (`zero_dim_solution_metadata`
  suite) and `python/test/zero_dim/metadata_test.py` cover the classifications, that
  each threshold is actually applied, the config round-trip, and the B1 defaults
  (which guards against the inverted-default regression).
- Follow-up: the classic-input parser does not yet read `CondNumThreshold` (the
  field default is used).
