# ADR-0076: With ccache, clang-cl writes dependency files through /showIncludes

**Status:** Accepted
**Date:** 2026-11-04

## Context

Both Windows CI jobs compile with clang-cl through ccache, and for as long as the stats step has existed it has reported every call uncacheable: `Cacheable calls: 0 / 149`, `Unsupported compiler option: 149 / 149`.  The saved cache was 13 KB, and the C++ test job compiled cold on every run, 32 minutes.  The history was three commits that never read that number: ccache was switched off for the C++ job "as a precaution" (87bb4c612), the wheel job's cache directory was persisted "to measure the hit rate" (b4884c790), and ccache was switched back on because the wheel job had "run clang-cl through ccache for months without issue" (0d4d2b0e8, #322).  Nothing broke; nothing was cached either.

The ccache log names the refused option: `-clang:-MD`, with `-clang:-MT<obj>` and `-clang:-MF<depfile>`.  CMake's `Platform/Windows-Clang.cmake` sets those for clang-cl under Ninja on purpose ("Prefer clang-cl's gcc-style depfile over cl-style /showIncludes"), overriding the `/showIncludes` form its MSVC base module sets.  ccache 4.14.1 (2026-09-27) made ccache reject every `/clang:` pass-through option by design, to avoid caching compilations whose pass-through options it cannot see, and 4.14.2 keeps that.  The two meet on every translation unit, so nothing is cacheable.  Precompiled headers and debug information, the usual suspects, were not involved: PCH is used by only two targets, and the refusal was universal.

## Decision

When ccache drives clang-cl, the root `CMakeLists.txt` sets `CMAKE_DEPFILE_FLAGS_<lang>` to `/showIncludes` and `CMAKE_<lang>_DEPFILE_FORMAT` to `msvc` for C and C++, the values CMake's own Windows-MSVC module uses for cl.  ccache supports `/showIncludes` (it stores the include notes with the result and replays them on a hit, so Ninja's `deps = msvc` parsing still works).  The Windows jobs also set `CCACHE_SLOPPINESS=pch_defines,time_macros` so the two targets that use a precompiled header are cacheable, `CCACHE_COMPILERTYPE=clang-cl`, and print `ccache --show-stats -v`, which lists uncacheable calls by reason.

## Consequences

- Do not remove the `/showIncludes` override because "CMake prefers the gcc-style depfile": with ccache, that preference costs every cache hit on Windows.  If ccache ever accepts `-clang:` options again, the override can go, and the verbose stats will show it.
- Do not switch ccache off on Windows as a precaution, and do not conclude it works from the absence of failures; read `Cacheable calls` in the stats step.  A cache that saves kilobytes is a cache that saves nothing.
- Keep the verbose stats step: it is the measurement the three earlier commits lacked.
