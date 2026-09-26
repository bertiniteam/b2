# ADR-0004: Linux wheel build parallelism is bounded by runner memory

**Status:** Accepted (amended 2026-06-16 and 2026-09-26 -- see the updates at the end; the
current values are 4 on x86_64 Linux and 3 on aarch64 Linux)  
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

## Update (2026-09-26): aarch64 Linux wheels build at 3

Linux wheels are now also built for aarch64 (#474), on GitHub's native `ubuntu-24.04-arm`
runners, which have the same 4 cores / 16 GB as `ubuntu-latest`.  At 4, three of the five
cold aarch64 wheel builds were killed while compiling the Python bindings: two with exit 143
(the pattern described above) and one with "the hosted runner lost communication with the
server", which GitHub attributes to a runner starved of CPU or memory.  The aarch64 C++ test
job, which builds no bindings, passed at 4, and the x86_64 wheels passed at 4 in the same run.

So `CMAKE_BUILD_PARALLEL_LEVEL` in `CIBW_ENVIRONMENT_LINUX` is now per architecture:
**4 on x86_64, 3 on aarch64** (`runner.arch == 'ARM64'`).  With 3, all five cold aarch64
wheel builds passed.  Why aarch64 needs the margin was not measured; the likely reason is a
larger per-TU peak from the aarch64 compiler, against x86_64's already zero-margin
4 x ~4 GB.

Consequences, in addition to the ones above:

- **Do not raise aarch64 wheels back to 4** to match x86_64 without measuring the heaviest
  binding TU's peak memory on aarch64 first.  A warm build will pass at 4 regardless (it is
  mostly ccache hits), so a green warm run is not evidence; only a cold build is.
- **The Linux C++ test job stays at 4 on both architectures.**  It builds no bindings, and
  passed cold at 4 on aarch64.
- A cold aarch64 wheel build takes somewhat longer than it would at 4.  Only cold builds
  (dependency bumps, cache eviction) pay this.
