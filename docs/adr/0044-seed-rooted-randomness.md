# ADR-0044: Seed-rooted randomness — every identity-relevant draw derives from a pinned, versioned stream

**Status:** Accepted
**Date:** 2026-07-03

## Context

Seeds in numerical algebraic geometry are theoretically irrelevant (probability-1
genericity) and practically decisive: a draw determines conditioning, path crossings,
and which run you can reproduce (ADR-0041, ADR-0017, ADR-0018).  The
structured-output-directory arc (rung 2) makes the seed half of a computation's ask
identity — `solve(sys, seed=42)` must mean the same coefficients forever, across runs,
platforms, standard libraries, and Boost versions, so that recorded runs are exactly
replayable and "what seed 42 looks like" is a stable, shareable object.

The engine layer was already sound: a fully-specified `std::mt19937` with domain-
separated streams (setup / per-path / worker; `random.cpp`), giving same-machine
reproducibility.  The gap was the **distribution layer**: `boost::random` /
`std::uniform_*_distribution` algorithms carry no stability contract — a Boost upgrade
or a different standard library can silently change what a seed draws.  A recorded
run's replayability cannot depend on that.

## Decision

A pinned, versioned draw stream (`b2rand/1`) in `bertini::records`
(`records/derive.hpp`, ADR-0043's namespace), through which every identity-relevant
draw routes:

- **`DrawStream`**: SHA-256 in counter mode (the vendored, FIPS-pinned hash).  A stream
  is keyed by `SHA-256(b2rand/1 || master || domain || index)` (each 8 big-endian
  bytes), preserving the existing setup/path/worker domain separation;
  `SetGlobalSeed` / `ReseedThisThread` reseed it alongside the legacy engine.
- **Fully specified draw algorithms** (documented in the header, exact forever):
  `Uint64` (8 stream bytes, big-endian); `UnitDouble` = `(Uint64 >> 11) * 2^-53`;
  `SymmetricDouble` = `2*UnitDouble() - 1`; `Bits(n)`; `IntSymmetric(B)` by rejection;
  `UnitRealMp(digits)` = `ldexp(Bits(k), -k)` with `k = ceil(digits*log2(10)) + 1`,
  materialized at the target precision (mpz→mpfr conversion correctly rounded,
  power-of-two scaling exact — bit-identical on every platform).
- **Converted call sites**: `RandomInt`, `RandomRat` (which also gains a deterministic
  redraw of a zero denominator — a latent div-by-zero in the legacy draw),
  `RandomMp<digits>` (and thus the whole multiprecision cascade: bounded-modulus draws,
  patches, slices, TD/MHom coefficients, gamma), the double draws in `num_traits`
  (`rand_complex`, `RandomUnit`) and `double_extensions` (`RandReal`).  A tiny
  dependency-free header (`records/draw_functions.hpp`) breaks the include cycle
  `num_traits → random.hpp → derive.hpp → num_traits`.
- The `mt19937` engine and its seeding remain for any non-identity consumer, but no
  identity-relevant draw touches it.

## Consequences

- **Same seed ⇒ digest-identical homotopies** (tested: the full seeded construction
  cascade — TD start, gamma, patch — digests equally on rebuild, differently under a
  different seed).  This is the rung-1/rung-2 layering paying off: ADR-0042/0043
  digests are the test oracle for draw determinism.
- **Cross-platform enforcement is free via CI**: the pinned raw draws and the pinned
  seed-42 homotopy digest in `seeded_randomness_test.cpp` run on Linux/macOS/Windows —
  a platform-dependent draw or rounding difference fails the build.  The pinned values
  ARE the `b2rand/1` contract: a failing pin after a draw change means bump the
  version, never edit the expectation.
- **What a seed means changed once, now**: pre-`b2rand/1` seeds drew different values
  (they were never a contract; the suite's seed-pinned tests are count/behavior-based).
  From this ADR forward, seed meaning changes only by version bump.
- Per-path draws (the ADR-0024 condition probe, via `ReseedThisThread(path index)`)
  are deterministic per (seed, path) and domain-separated from setup draws — the
  order-independence parallel tracking requires.
- The v0 records prototype's pickled-homotopy stopgap is now removable: a homotopy
  instance is a pure function of (ask, seed).

## Addendum (2026-07-04): draw COMPOSITION order is part of the contract

The cross-platform seeded-homotopy fixture caught an architecture split (x86_64 vs
aarch64 computing different seed-42 homotopy digests) whose root cause was not the
draws — the raw `b2rand/1` pins passed everywhere — but the ORDER two draws were
composed into one complex number: expressions like `complex(RandomMp(...),
RandomMp(...))` leave the evaluation order of the two draws unspecified in C++, and
gcc really does order them differently on the two architectures (every (re, im) pair
of every seeded coefficient was transposed between them).

Contract, from this addendum forward: **component draws are sequenced explicitly,
real part first, then imaginary — never two draws inside one full-expression.**  All
eight such sites (`random.hpp`, `num_traits.hpp`) were rewritten with named
`re`/`im` locals; the pinned seed-42 digest (minted on aarch64, which already
evaluated left-to-right) is unchanged, and x86_64 now agrees with it.  Any future
draw-composition helper must follow the same rule; the cross-platform fixture is the
enforcement.
