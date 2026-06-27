# ADR-0034: The CLI emits Bertini 1.7-compatible machine-readable solution files

**Status:** Accepted
**Date:** 2026-06-27

## Context

Bertini 2 already produces Bertini 1.7-compatible *input* (`System::to_classic_input`, the
classic writer) so the same problem can be solved in either tool. The *output* side had no such
parity: the CLI (`core/src/blackbox/main_mode_switch.cpp`) wrote only `main_data` and `raw_data`,
in a Bertini-2-specific layout. Neither carries a solution-count header, neither uses Bertini 1's
`-1` raw terminator, and **none** of Bertini 1's companion solution files were written at all.

That broke backwards compatibility for any tool (or person, or benchmark) that consumes Bertini 1
output: there was no `finite_solutions` to read, and `main_data`'s first line is the *variable*
count, not the solution count — a silent trap (the old `benchmark/run_benchmark.py` parsed it as a
solution count and got the wrong number).

Backwards compatibility was the whole point of the classic writer; the reading side should match.
Note: **byte-for-byte** equality with Bertini 1 is impossible and is explicitly *not* a goal — the
implementations differ, and even the random gamma/seed differ, so the actual solution *values* and
ordering will never coincide. The contract is **machine-readable parity**: same file names, same
count-led block format, so a parser written for Bertini 1.7 reads Bertini 2 unchanged.

## Decision

The zero-dim CLI additionally writes the Bertini 1.7 solution files, in the classic format
(first line = count, then one block of `NumVariables` `"re im"` coordinate lines per solution, in
user/dehomogenized coordinates, blocks blank-separated):

- `finite_solutions`, `real_finite_solutions`, `nonsingular_solutions`, `singular_solutions`
- `raw_solutions` (each successful endpoint preceded by its path number)

Implementation is additive and reuses what already exists:

- `output::Classic<ZeroDim>` (`core/include/bertini2/nag_algorithms/output.hpp`) gains the writers,
  built on the existing classified accessors (`FiniteSolutions()`, `RealSolutions()`,
  `NonsingularSolutions()`, `SingularSolutions()`, `FinalSolutionMetadata()`; see ADR-0015 for the
  classification thresholds) and the existing `generators::Classic` number formatter.
- `AnyZeroDim` (`zero_dim_solve.hpp`) gains five pure-virtual `Write*Solutions` methods; `ZeroDim`
  is the only implementer, so this does not affect any other type.
- The CLI opens and writes the five files after `main_data`/`raw_data`.

`main_data` and `raw_data` are left as-is: they are informational, and chasing byte-compatibility
there has no payoff given the no-byte-equality reality above.

## Consequences

- Tools and scripts that parse Bertini 1.7 solution files now read Bertini 2 CLI output unchanged.
- A single parser serves both solvers — used by `benchmark/comparison/` to time Bertini 2 against
  Bertini 1 on the same emitted input and cross-check that they agree on the solution count.
- Verified: on cyclic-5 the CLI writes 70 finite solutions whose *set* matches the trusted
  Bertini 1.7 oracle (`python/test/zero_dim/data/cyclic5_finite_solutions.txt`) to ~2.5e-12. A C++
  test (`core/test/nag_algorithms/zero_dim.cpp::classic_solution_file_output`) pins the file format
  and counts.
- The five extra files are written only by the manager rank, alongside the existing output writes.
