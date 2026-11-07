# ADR-0079: A thread frees its own MPFR caches

**Status:** Accepted
**Date:** 2026-11-07

## Context

`test_nag_algorithms` crashed on Windows with heap corruption (0xc0000374) about two seconds in, deterministically, when PR #527 moved the C++ job to clang 22 and `link.exe`.  Rerun under cdb, a worker thread faulted reading `0xfeeefeee`, the debug heap's fill for freed memory, inside `mpfr_bernoulli_freecache`; under full page heap the same stack faulted one frame deeper, in libgcc's `_emutls_get_address`:

```
libgcc_s_seh_1!_emutls_get_address
libmpfr_6!mpfr_bernoulli_freecache
libmpfr_6!mpfr_free_cache2
test_nag_algorithms!__dyn_tls_dtor
ntdll!LdrpCallTlsInitializers
ntdll!LdrShutdownThread
```

Boost.Multiprecision frees MPFR's per-thread caches from a `thread_local` object (`mpfr_cleanup<true>::thread_initializer`), which every destroyed `real_mp` registers.  conda-forge's MPFR for Windows is built by gcc, so its thread-local variables live in libgcc's emulated TLS, whose per-thread memory goes when the DLLs detach from the thread.  In an executable, C++ `thread_local` destructors run from the executable's own TLS callback, which the stack above shows running inside thread shutdown, and by then the emulated TLS block has been freed.  MPFR reads its cache pointers out of freed memory and frees what they appear to point at.

The wheels did not crash: there the destructor belongs to the extension DLL, which sits on top of MPFR in the load order, and every Windows wheel's pytest run passed with the same compiler and linker.  That reasoning about the order is inferred from the stack and the wheels' results, not measured.

The read of freed memory does not depend on the compiler, so it was present before #527.  What made it fatal there is not established: whether freed memory still looks harmless depends on heap layout, and the C++ job passed with clang 19 and lld-link the day before.

## Decision

- On Windows, `bertini2/mpfr_extensions.hpp` explicitly specializes `mpfr_cleanup<true>::force_instantiate` to skip Boost's `thread_local` object and keep its process-wide cleanup.
- A thread b2 starts frees its own caches with `ReleaseThreadMpfrCaches()` as the last thing it does, while its TLS is alive.  `WorkerThreadPool` is the only place b2 starts threads, and it does this on every platform, so the behavior has one description.
- A thread b2 did not start, on Windows, keeps its caches until the process ends.  They are small and bounded by the thread count; an uncatchable crash is not.

## Consequences

- Do not remove the specialization because `ReleaseThreadMpfrCaches` looks sufficient, or the reverse.  Without the specialization Boost's destructor still runs after the thread's own cleanup, and without the release the caches leak per solve.
- Do not add a `thread_local` destructor that calls an MPFR function reading MPFR's thread-local state, in b2 or in a test.  On Windows that state is already gone when the destructor runs.
- A new place that starts a thread that uses `real_mp` calls `ReleaseThreadMpfrCaches()` last.
- `bertini2/mpfr_extensions.hpp` stays the first b2 header to reach Boost's MPFR backend in any translation unit that instantiates it, or the specialization comes too late.  If a Boost release renames `mpfr_cleanup` or `force_instantiate`, the specialization stops compiling, which is the signal to look again.
- This can go if conda-forge's Windows MPFR stops using emulated TLS, or Boost stops freeing the caches from a `thread_local`.
- `thread_pool_workers_that_filled_mpfr_caches_exit_cleanly` fills the Bernoulli and pi caches on pool workers and exits them.  On Windows, the C++ job's crash step (`cdb`, then full page heap) is how to read a failure like this one; it runs only when a suite failed.
