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
