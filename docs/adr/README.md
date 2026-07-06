# Architecture Decision Records

This directory contains ADRs for Bertini2 — short committed documents capturing
significant design decisions: what was decided, why, and what to watch for.

Each ADR follows the template:

```
# ADR-NNNN: <title>
**Status:** Accepted | Superseded by ADR-XXXX | Deprecated
**Date:** YYYY-MM-DD
## Context
## Decision
## Consequences
```

## Index

| ADR | Title | Area |
|-----|-------|------|
| [0001](0001-eigenpy-writable-ref-scalar-by-value.md) | Pass scalar mpc_complex args by value when adjacent to writable Eigen::Ref _(superseded by 0008)_ | Python bindings |
| [0002](0002-endgame-run-single-arg-api.md) | Endgame run() takes only the start point; boundary times set at construction | Python bindings / API |
| [0003](0003-manylinux-no-full-pytest.md) | Linux wheel CI uses import smoke test only; full pytest runs on macOS/Windows | CI |
| [0004](0004-ci-ubuntu-parallel-level-2.md) | Ubuntu wheel builds use CMAKE_BUILD_PARALLEL_LEVEL=2 | CI |
| [0005](0005-github-actions-matrix-exclude-not-if.md) | Use matrix exclude to filter CI jobs, not matrix.os in job-level if | CI |
| [0006](0006-eigenpy-uninitialized-numpy-slot-guards.md) | Guard numpy mpfr/mpc dtypes against uninitialized slots | Python bindings / numerics |
| [0007](0007-rational-config-constants-default-construct-precision-trap.md) | Rational config constants, and the DefaultConstruct static-precision trap | Numerics / config |
| [0008](0008-track-path-result-numpy-object-not-writable-ref.md) | Take writable-Ref binding outputs as numpy objects (supersedes 0001) | Python bindings |
| [0009](0009-amp-endgame-internal-double-state.md) | AMP endgame needs internal double-precision state, mirroring the tracker _(proposed; implementation deferred)_ | Numerics / performance |
| [0010](0010-function-tree-single-inheritance-capability-classes.md) | function_tree uses single inheritance; cross-cutting traits are capability classes | Core / function_tree |
| [0011](0011-differentiation-simplified-construction-no-mutation.md) | Differentiation emits already-simplified trees; held functions are never mutated | Core / function_tree |
| [0012](0012-precedence-aware-printing-reparse-invariant.md) | Precedence-aware printing; printed trees must re-parse to the same values | Core / function_tree |
| [0013](0013-solutions-in-user-coordinates.md) | Solutions reported in user coordinates by default; internal coords explicit; HomogenizePoint lift | Core / API |
| [0014](0014-explicit-template-instantiation-closed-universe.md) | Explicit template instantiation of the endgame/ZeroDim universe | Build / core |
| [0015](0015-solution-metadata-classification-bertini1.md) | Solution metadata classification (finite/real/singular) matches Bertini 1; post-processing config corrected | Core / nag_algorithms |
| [0016](0016-distributed-zerodim-broadcast-homotopy-seed.md) | Distributed ZeroDim broadcasts the homotopy seed so all ranks form the identical homotopy | MPI / nag_algorithms |
| [0017](0017-flakiness-is-tolerances-not-seeding.md) | Solve flakiness is a tolerance problem, not a seeding problem | Numerics / tracking |
| [0018](0018-bound-every-test-with-a-timeout.md) | Bound every test with a timeout; pin seeds and skip known platform grinds | CI / testing |
| [0019](0019-omit-debug-info-from-release-builds.md) | Omit debug info (-g) from Release builds; per-TU peak RSS caps build parallelism | Build / core |
| [0020](0020-structured-block-start-systems-need-blend-homotopy.md) | A structured-block start system must be coupled into a homotopy via a blend (MakeHomotopy / blend_homotopy), not System node arithmetic | Core / nag_algorithms |
| [0021](0021-block-evaluations-fully-define-their-own-rows.md) | A block (and the patch) must fully define its own output rows; System allocates result buffers uninitialized | Core / system |
| [0022](0022-path-crossings-are-detected-and-retracked.md) | Path crossings at the endgame boundary are detected and re-tracked, not tolerated | Core / nag_algorithms |
| [0023](0023-unified-speculative-full-path-parallelism.md) | ZeroDim uses one per-path primitive and a speculative-full-path parallel model | Core / nag_algorithms |
| [0024](0024-condition-probe-tracker-owned-per-path.md) | The condition-number random direction is refreshed once per path, by the tracker | Numerics / tracking |
| [0025](0025-randomization-block.md) | Randomization of overdetermined systems is a first-class block | Core / system |
| [0026](0026-moving-homotopy-blend-moving-rows-only.md) | Moving homotopies blend only the moving rows; fixed equations stay sibling blocks | Core / nag_algorithms |
| [0027](0027-slp-management-output-sets-program-memory-freeze-partition.md) | SLPs belong to output-sets, not nodes; Program/Memory split with a freeze-set tape partition | Core / function_tree |
| [0028](0028-named-node-taxonomy.md) | Named-node taxonomy — NameHolder, three named kinds, NamedExpression replaces Handle, Jacobian deleted | Core / function_tree |
| [0029](0029-gamma-trick-max-precision-unit-modulus.md) | The gamma-trick constant is a unit-modulus complex generated at maximum precision | Numerics / system |
| [0030](0030-fixed-multiple-uniform-ambient-precision-adaptive-default.md) | Fixed-multiple solves use one uniform ambient precision; ZeroDim defaults to adaptive | Numerics / nag_algorithms |
| [0031](0031-eigenpy-vec-index-returns-aliasing-view-copy-before-storing.md) | Indexing an eigenpy Vec returns an aliasing view — copy values before storing them (fixes #259) | Python bindings |
| [0032](0032-slice-is-a-thin-wrapper-over-linear-forms-block.md) | A witness set's Slice is a thin wrapper over LinearFormsBlock; LinearSlice retired | Core / nag_datatypes |
| [0033](0033-slice-shape-contract-we-own-dimensionality.md) | Slice shape is OUR accessors' contract, not eigenpy's: coefficients() always 2-D, slice[i] a form vector, slice[i:j] a sub-Slice | Python bindings |
| [0034](0034-tiered-slp-arithmetic.md) | SLP evaluates in the narrowest sufficient numeric tier (integer<real<complex), orthogonal to precision; widen only when an op demands it; measurement-gated | Core / SLP |
| [0035](0035-slp-eval-is-allocation-free.md) | The SLP eval loop is allocation-free — lower integer powers to multiplies (exp-by-squaring), fold differentiated exponents to literals, no by-value temporaries; ~10x faster mp eval, 0 allocs/eval ≤256 digits | Core / SLP / performance |
| [0036](0036-cli-emits-bertini1-compatible-solution-files.md) | The zero-dim CLI emits Bertini 1.7-compatible machine-readable solution files (finite/real_finite/nonsingular/singular/raw_solutions, count-led); machine-readable parity, not byte-equality | Core / blackbox / IO |
| [0037](0037-cauchy-endgame-securitylevel-guard-inverted.md) | The Cauchy endgame's SecurityLevel guard was inverted (diverging paths ground at escalating precision) | Core / endgames |
| [0038](0038-amp-criterion-b-uses-newton-residual-not-size-proportion.md) | AMP Criterion B (cost model) must use the latest Newton residual, not size_proportion | Numerics / trackers |
| [0039](0039-zerodim-binding-keeps-system-alive.md) | The Python solver binding must keep the given System alive (custodian-and-ward) | Python bindings |
| [0040](0040-split-zerodimsolver-homotopysolver.md) | Split ZeroDim into HomotopySolver (continuation engine) + ZeroDimSolver (algorithm); retire the system-management policy | Core / nag_algorithms |
| [0041](0041-bounded-modulus-random-coefficients.md) | Well-conditioned random coefficients: bounded-modulus scalars, conjugate-orthonormal matrices (start systems, patches, slices, randomization), matching Bertini 1 | Numerics / core |
| [0042](0042-system-content-identity-seal-and-interning.md) | System content identity: versioned canonical encoding + SHA-256 digest, Seal(), System/Program interning, re-intern-on-load (ADR-0027 E4); append-only version registry enforces the bump | Core / system |
| [0043](0043-config-canonical-encoding-and-digests.md) | Config canonical encodings + digests (b2cfgenc/&lt;n&gt;): bit-exact doubles, fixed enum names, SettingsDigest composition; seed and num_threads excluded from identity; append-only version registry enforces the bump | Core / records |
| [0044](0044-seed-rooted-randomness.md) | Seed-rooted randomness (b2rand/1): every identity-relevant draw derives from a pinned SHA-256 counter stream; same seed = digest-identical homotopies, cross-platform, forever | Core / records |
| [0045](0045-structured-output-directory-b2rec1.md) | The structured output directory: record schema b2rec/1 + C++ OutputDirectory (JSONL history, content-addressed definitions, derived views, cross-impl tested) | Core / records |
| [0046](0046-solver-records-seam-ensure-answered.md) | The solver records seam: solve() is ensure-answered — track records ARE serialized FullPathResults, recall replays them through StoreFullPathResult, manager is sole writer, resume = memoization | Core / records |
| [0047](0047-casual-records-surface.md) | The casual records surface: bertini.solve/save/load, Solution = points that remember, CLI records-on-by-default beside b1 files (no flag) | Python + CLI / records |
| [0048](0048-cauchy-endgame-security-and-operating-zone.md) | Cauchy divergence handling: security check watches the ENDPOINT, truncates only in the operating zone, no pole-growth truncation (the acceptance gate alone cures junk-success); restores cyclic-6's 156 solutions (refines #70) | Core / endgames |
