# ADR-0019: Omit debug info (-g) from Release builds

**Status:** Accepted
**Date:** 2026-06-15

## Context

Cold rebuilds of the template-heavy core are constrained by **per-TU peak resident memory**,
not wall-clock alone (ccache keeps warm rebuilds cheap). On a ~60 GB host the heaviest TUs set
how high `-j` can go before the box OOMs, so lowering the per-TU peak is a compounding win:
it raises usable parallelism.

`core/CMakeLists.txt` added `-g` to *Release* builds (`$<$<CONFIG:Release>:-g>`), i.e. the
default `-O3` build also emitted full debug info. That had two costs:

- **Disk.** Debug info bloats objects ~3–4x. Example (this tree): `algorithm_builder.cpp.o`
  177 MB with `-g`. A full build (including `python_bindings`) exhausted a 45 GB root
  partition mid-compile (`as: No space left on device`).
- **Peak RSS.** `-g` raised per-TU peak memory enough to lower the achievable `-j`.

Measured per-TU peak RSS (standalone `-O3 -DNDEBUG` compile, GCC, this host):

| TU | with `-g` | without `-g` |
|---|---|---|
| `endgames_eti.cpp` | 2.07 GB | 1.78 GB |
| `zero_dim_eti.cpp` | 5.48 GB | 3.32 GB |
| `zero_dim_blackbox_eti.cpp` | 5.88 GB | 3.53 GB |
| `algorithm_builder.cpp` (cap) | **6.51 GB** | **4.07 GB** |

CI does not ship these symbols: wheels are stripped during the build, so Release `-g` bought
nothing for released artifacts — only local crash backtraces, which a developer can re-enable
on demand.

## Decision

**Do not emit `-g` in Release.** Remove `$<$<CONFIG:Release>:-g>` from `core/CMakeLists.txt`.
Debug builds (`-DCMAKE_BUILD_TYPE=Debug`) still get `-g`; a developer who wants symbols on an
optimized build can use `RelWithDebInfo` or add `-g` locally.

This single change drops the per-TU memory cap **6.51 → 4.07 GB** (≈ `-j8` → `-j14` headroom
on a 60 GB host), shrinks objects ~3–4x (resolving the disk exhaustion), and speeds compiles
(no debug-info emission). The full C++ test suite stays green (`ctest` 10/10).

**Splitting the heavy TUs was investigated and deliberately *not* done.** The idea (split the
ETI files / `algorithm_builder.cpp` by tracker × start × endgame to lower peak further) runs
into a floor: instantiating the ZeroDim template cone (trackers + endgames + system) costs
~2.9 GB *fixed*, with only ~0.08 GB marginal per added combo (measured: 2 combos → 3.00 GB,
6 combos → 3.32 GB). Consequences:

- Splitting an ETI file from 6 combos into 3 TUs barely lowers its peak (3.32 → 3.00 GB) while
  ~tripling that file's compile CPU (each TU re-pays the fixed cone cost).
- After the `-g` removal the ETI files (≤ 3.53 GB) are already *below* the cap
  (`algorithm_builder.cpp`, 4.07 GB), so splitting them does not lower the global cap at all.
- Restructuring `algorithm_builder.cpp` into per-tracker factory TUs would cut only ~0.5 GB
  (cap → ~3.5 GB, limited then by the ETI files) at the cost of new headers/TUs and edits
  rippling into `switches_zerodim.hpp`, `user_homotopy.hpp`, and two tests.

The cheap one-liner captures the bulk of the available win; the splitting is a poor trade.

## Consequences

- Release builds no longer carry debug symbols. Local optimized debugging uses
  `RelWithDebInfo` or a manual `-g`. Released wheels are unaffected (CI already strips).
- The per-TU memory cap is now `algorithm_builder.cpp` at ~4 GB. If a future change needs to
  push `-j` higher, the per-tracker factory split (above) is the next lever — but weigh it
  against the ~3 GB cone floor before spending the effort.
- This does not revisit ADR-0014: explicit instantiation of the closed ZeroDim/endgame
  universe stays as-is (one home in `libbertini2`, `extern template` in the umbrella headers).
