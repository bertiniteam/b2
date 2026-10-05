# ADR-0072: Seed 0 is a seed; the classic input format keeps Bertini 1's meaning

**Status:** Accepted
**Date:** 2026-10-05

## Context

`SetGlobalSeed(0)`, and so `bertini.set_random_seed(0)`, meant "draw a seed from entropy".  A
user who set seed 0 to make a run reproducible got a different run every time (#475).  0 is the
first value people try, and in numpy (`default_rng(0)`, `np.random.seed(0)`) and in Python's
`random` (`random.seed(0)`) it is an ordinary seed.  Those libraries ask for entropy with `None`
or no argument.

Bertini 1's input file means something else.  There `randomseed: 0` is the default, and it means
"choose a seed".  b2's command-line program reads Bertini 1 input files, and a file that says 0,
or says nothing, expects a fresh seed.

## Decision

- **In the library API, every value is a seed, 0 included.**  `SetGlobalSeed(0)` sets seed 0.
  In Python, `set_random_seed(0)` sets seed 0.
- **Entropy is its own request.**  `SetGlobalSeedFromEntropy()` in C++, and `set_random_seed()`
  or `set_random_seed(None)` in Python, draw a seed and return it, so the run can be reproduced
  from that value.  `GetGlobalSeed()` draws one on first use if none was set, and keeps it.
- **The classic input format keeps Bertini 1's meaning.**  `ApplyClassicRandomSeed` reads
  `randomseed: 0` as "draw from entropy" and any other value as the seed.  The command-line
  program applies a classic input's seed through it.

## Consequences

- Seed 0 does not carry between the two interfaces.  `set_random_seed(0)` is reproducible, while
  `randomseed: 0` in a classic input file is not.  This is the one value on which they disagree,
  and the docstrings of both say so.
- **Do not make 0 mean entropy in the library API again**, and do not make the classic format's 0
  an ordinary seed: each interface follows the convention its users bring to it.
- "Not set yet" is a flag of its own in `random.cpp`, not the value 0.
- The derived seeds (`DeriveSolveSeed`, `DerivedWorkerSeed`) still never return 0.  That way a
  seed recorded from a library solve replays through a classic input file too.
- Pinned by `seed_zero_is_an_ordinary_seed`,
  `a_seed_from_entropy_is_reported_and_reproduces_the_run` and
  `classic_randomseed_zero_draws_from_entropy` (`core/test/classes/seeded_randomness_test.cpp`),
  and in Python by `python/test/random/seeding_test.py`, which includes the reproducer from #475.
