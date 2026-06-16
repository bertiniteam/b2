# ADR-0004: Ubuntu wheel builds use CMAKE_BUILD_PARALLEL_LEVEL=2

**Status:** Accepted  
**Date:** 2026-06-07

## Context

Linux wheels are built inside a manylinux Docker container on GitHub's `ubuntu-latest`
runners (7 GB RAM each). The build includes: Boost from source, eigenpy from source,
and the bertini2 Python bindings. All six Python versions (3.9–3.14) build in parallel
as separate matrix jobs.

`CMAKE_BUILD_PARALLEL_LEVEL=1` was originally set (commit `f4c58c1b` context) to
prevent OOM during the wheel build. After TU splitting (`f4c58c1b`) reduced per-TU
peak RAM from ~154 MB to ~119 MB, the one-job-at-a-time restriction was revisited.

### What happened with CMAKE_BUILD_PARALLEL_LEVEL=4

Three of six ubuntu wheel builds were killed with exit code 143 (SIGTERM), not 137
(SIGKILL / OOM killer). The GitHub Actions message was:

```
Killed    cibuildwheel . --output-dir wheelhouse
The runner has received a shutdown signal. This can happen when the runner
service is stopped, or a manually started runner is canceled.
```

The kills happened at nearly the same wall-clock time across different matrix jobs
(all while compiling heavy TUs near step 58/72). This pattern — multiple independent
jobs dying at the same moment, with SIGTERM not SIGKILL — indicates the runner HOST
processes (or Docker daemon) ran out of memory on shared GitHub infrastructure, causing
the Docker containers to be terminated.

GitHub's `ubuntu-latest` runners share physical hosts. With six matrix jobs each
running four compilation threads: 6 × 4 = 24 simultaneous compilations, each using
up to ~119 MB, for a theoretical peak of ~2.8 GB of compilation RAM distributed
across however many physical hosts hold these six containers. On a shared host, the
total resident memory can exceed what is available.

### Why 1-parallel was too conservative

With 1-parallel: build time was ~40 min per ubuntu wheel job. Total CI time was
gated by this. TU splitting had already paid its RAM cost; the 1-parallel restriction
provided no additional safety benefit after that optimization landed.

### Empirical results

| `CMAKE_BUILD_PARALLEL_LEVEL` | Build time | Outcome |
|-----|---|---|
| 1 | ~40 min | always stable |
| 2 | ~31 min | always stable (tested on 2026-06-07) |
| 4 | ~15 min | 3/6 SIGTERM killed (tested on 2026-06-07) |

## Decision

Set `CMAKE_BUILD_PARALLEL_LEVEL=2` in `CIBW_ENVIRONMENT_LINUX`.

```yaml
CIBW_ENVIRONMENT_LINUX: "... CMAKE_BUILD_PARALLEL_LEVEL=2"
```

This cuts ubuntu wheel build time from ~40 min to ~31 min without triggering the
shared-host memory pressure that 4-parallel caused.

## Consequences

- **~9 min saved** per ubuntu build cycle compared to 1-parallel. With 6 matrix jobs,
  this is the dominant contributor to total CI wall-clock time.
- **Stable.** 2-parallel has been tested and all 6 ubuntu builds pass consistently.
- **If builds start failing with SIGTERM again**, revert to 1. Do not go above 2
  without evidence that the shared-host memory situation has changed (e.g., GitHub
  moves to dedicated runners or increases runner RAM).
- **Do not confuse this with local builds.** On this machine (16 GB RAM, dedicated),
  4 parallel is the right setting via the `CMAKE_BUILD_PARALLEL_LEVEL` env var.
  The CI value is separate and exists only in `CIBW_ENVIRONMENT_LINUX`.
- **The 4-parallel failures were SIGTERM (exit 143), not SIGKILL (exit 137).** If
  future failures show exit 137, that is the OOM killer — a different problem with
  a different fix (further TU splitting or 1-parallel).

## Update (2026-06-16): raised to 4 on Linux

The two conditions this ADR named for going above 2 have both been met:

1. **Per-TU peak memory dropped.** ADR-0019 removed `-g` from Release builds, cutting the
   heaviest TU from ~6.5 GB to ~4.0 GB.
2. **Runner RAM increased.** GitHub's standard Linux runners are now 4-core / 16 GB (they were
   ~7 GB when this ADR was written).

So `CMAKE_BUILD_PARALLEL_LEVEL` is raised **2 → 4** for the Linux wheel build
(`CIBW_ENVIRONMENT_LINUX`) and the Linux C++ test build (`cmake --build --parallel`, made
OS-conditional via `runner.os == 'Linux'`).  4 = the full core count; worst-case concurrent
memory is bounded by the few heavy TUs (~4 GB each, and they rarely align), and ccache keeps most
rebuilds compile-free.  **macOS stays at 2** (3-core / 7 GB runner).  If a *cold* Linux build shows
exit 137 (OOM killer), drop Linux back to 3 (one heavy-TU slot of headroom).
