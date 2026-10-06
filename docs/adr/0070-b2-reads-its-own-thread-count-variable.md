# ADR-0070: b2 reads its own thread-count variable, and defaults to every available CPU

**Status:** Accepted
**Date:** 2026-10-05

## Context

b2 tracks paths on a pool of `std::thread` workers.  It uses no OpenMP, but until 4.0 it took
its thread count from `OMP_NUM_THREADS`, because that variable is familiar and because SLURM
sets it from `--cpus-per-task`.  The default differed by entry point: a standalone solve used
every CPU, and an MPI rank used one thread.

`OMP_NUM_THREADS` is not b2's alone.  numpy's OpenBLAS reads it when `OPENBLAS_NUM_THREADS` is
unset, and so does any OpenMP library in the same process.  Two things followed.

- Pinning b2 to one thread for a reproducible run pinned numpy's linear algebra as well, and
  pinning numpy pinned b2.  An experiment that wanted one and not the other could not say so.
- A user who set the variable for one library changed the other without knowing it.

## Decision

b2 reads **`BERTINI_NUM_THREADS`** and ignores `OMP_NUM_THREADS`.  One function,
`parallel::EffectiveThreadCount`, decides the count for a standalone solve and for each MPI
rank, in this order:

1. `BERTINI_NUM_THREADS`, when set to a positive integer;
2. the solver's `ZeroDimConfig::num_threads`, when positive;
3. `parallel::AvailableCpuCount()`: the CPUs the process may run on.  On Linux that is the
   affinity mask (`sched_getaffinity`), so `taskset`, cpusets and a launcher's binding are
   respected; elsewhere it is `std::thread::hardware_concurrency()`.

The default is the one GCC's OpenMP runtime uses, so a user gets threaded runs without
configuring anything.

## Consequences

- **Do not read `OMP_NUM_THREADS` again**, not even as a fallback.  A fallback brings back the
  coupling with numpy for everyone who sets it for numpy's sake.
- **Do not give MPI ranks a different default.**  One rule for both entry points is the point;
  a rank that should be serial is told so by the variable or the config.
- Under MPI, several ranks on one machine can each see most of its cores, and together start
  more threads than there are cores.  The CLI help and the "Solving at scale" tutorial tell
  users to set `BERTINI_NUM_THREADS` when they run more than one rank per machine.  That is the
  same burden OpenMP places on hybrid MPI programs.
- A job script that set `OMP_NUM_THREADS` for b2 must change.  The 4.0.0 CHANGELOG says so.
- `AvailableCpuCount` counts CPUs, not CPU quota.  A container limited by a CFS quota
  (`docker --cpus`) without a cpuset still reports every CPU it can be scheduled on.
- Tools that need a serial, reproducible run (`tools/refresh_doc_artifacts.py`) set
  `BERTINI_NUM_THREADS=1` for b2 and `OMP_NUM_THREADS=1` for numpy, each for its own library.
- Pinned by `effective_thread_count_is_sane`, `bertini_num_threads_overrides_the_configured_count`
  and `omp_num_threads_does_not_reach_b2` in `core/test/nag_algorithms/threaded_solve.cpp`.
