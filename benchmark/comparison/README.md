# Bertini 2 vs other solvers — serial & MPI comparison benchmark

This benchmark measures how fast the **bertini2 CLI** solves zero-dimensional polynomial systems
compared with other numerical-algebraic-geometry solvers, under **identical settings**, and at the
same time checks that the solvers **agree on the answer**.

Today the only other solver wired up is **Bertini 1**. The design is solver-agnostic: a future
phase will add [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/) simply by
adding one more adapter in `solvers.py` — nothing else changes.

These numbers are intended to be reproducible and citable in a paper, so the methodology below is
spelled out in full.

---

## Why this is a fair comparison

Each test system is built once in Python with pybertini and emitted to a **single classic Bertini
input file** via `System.to_classic_input(...)`. That one file — same equations, same precision
mode, same predictor, same tolerances, same step-size cadence — is handed to *every* solver. So
"same settings in both" is not a promise we have to police; it is true by construction. (The
emitter is exercised by `python/test/classes/classic_writer_test.py` and round-trip fidelity is
pinned in C++ by `classic_parsing_test::classic_writer_round_trips_a_system`.)

Two phases:

1. **Serial** — one process, one thread. The MPI-built binaries are simply run directly (a single
   MPI process), which is the serial baseline.
2. **MPI, single thread** — `mpirun -n N` with `OMP_NUM_THREADS=1`, sweeping `--ranks`.

## Correctness/agreement comes for free — and is *not* part of the timing

Because all solvers read the same input, every solver should return the **same number of
solutions**, equal to the system's known count. The driver reports this as a `matches_expected`
column and flags any disagreement.

Both solvers are parsed with **one** parser: the Bertini 2 CLI now writes Bertini 1.7-compatible
machine-readable solution files (`finite_solutions`, `real_finite_solutions`,
`nonsingular_solutions`, `singular_solutions`, `raw_solutions` — each count-led; see ADR-0036), so
the count comes from `finite_solutions` for both. Byte-for-byte equality with Bertini 1 is *not*
expected (different implementations and RNG); only the file format and the solution *count* match.

This check is performed **after** the timed region, by reading each solver's output files. The
stopwatch (`time.perf_counter()`) brackets **only** the solve subprocess; parsing solution counts
and comparing them adds nothing to any reported wall time.

---

## Prerequisites

- An **MPI-built `bertini2`** (`cmake -DENABLE_MPI=ON ...`; target `bertini2_exe`, binary
  `./build/core/bertini2`).
- For the Bertini 1 column, an **MPI-built Bertini 1** binary (you supply it; it is not part of
  this repo). Point `--bertini1` at it.
- `mpirun` on `PATH` for any `--ranks` greater than 1.
- `pybertini` importable (the `bertini` Python package, matching your built `_pybertini`).

> Both Bertini 1 and the bertini2 CLI write their output files (`main_data`, `raw_data`, `output`,
> `failed_paths`, …) into the current directory. The driver runs **every** invocation in its own
> fresh temp directory and deletes it afterward, so runs never collide and nothing is left behind.

---

## How to run

Serial only, bertini2 vs Bertini 1, on the default system subset:

```bash
python benchmark/comparison/run_comparison.py \
    --bertini2 ./build/core/bertini2 \
    --bertini1 /path/to/bertini \
    --output comparison_results.csv
```

Pick systems explicitly and add an MPI single-thread sweep:

```bash
python benchmark/comparison/run_comparison.py \
    --bertini2 ./build/core/bertini2 \
    --bertini1 /path/to/bertini \
    --systems cyclic6 cyclic7 katsura5 \
    --ranks 1 2 4 8
```

bertini2 only (no Bertini 1 installed) — still times and checks counts:

```bash
python benchmark/comparison/run_comparison.py --bertini2 ./build/core/bertini2 --systems all
```

Key options: `--systems NAME ... | all`, `--ranks N ...`, `--repeats N` (keep fastest),
`--timeout SECS`, `--mptype {0,1,2}`, `--predictor {0,2,5}`, `--mpirun`, `--mpirun-args`,
`--output FILE`. Run with `-h` for the full list.

For publishable numbers: use a quiet, non-throttling machine and `--repeats 3` or more (the fastest
time is kept).

---

## Historical data (tracking b2's performance over time)

Every run **appends** one row per measurement to a committed history file (default
`benchmark/comparison/history.csv`; `--history FILE` to redirect, `--no-history` to skip). Each row
records full provenance so results stay interpretable years later:

```
timestamp, host, cpu, os, b2_commit, mptype, predictor,
solver, solver_version, system, ranks, threads,
wall_time_s, solutions_found, expected_count, matches_expected, note
```

Commit this file. As you improve Bertini 2, new rows accumulate and you can plot wall time for a
given `(system, solver, ranks)` over `b2_commit`/`timestamp` to *see* the speedups land. Use
`--note` to tag a row (e.g. `--note "after SLP rework"`).

**Wall-times are only comparable within the same machine.** `host`+`cpu` are recorded precisely so
you filter to one machine before comparing across time; a `b2_commit` ending in `-dirty` means the
working tree had uncommitted changes (not a clean, citable point). Bertini 1 will likely never move,
but HomotopyContinuation.jl is actively developed, so its `solver_version` matters too.

---

## Test bed

All systems are genuinely zero-dimensional, so the solution count is a fixed known value used as
the correctness ground truth. Generators live in `systems.py`; add a family by registering a
builder in the `SYSTEMS` dict.

| Name(s)            | Family    | Vars | Solutions | Notes |
|--------------------|-----------|------|-----------|-------|
| `diag3/5/6`        | diagonal  | 3/5/6 | 27/243/729 | `xᵢ³ − cᵢ`; trivial warmup/sanity tier |
| `cyclic5`          | cyclic-n  | 5    | 70        | |
| `cyclic6`          | cyclic-n  | 6    | 156       | |
| `cyclic7`          | cyclic-n  | 7    | 924       | |
| `katsura3..6`      | Katsura-n | 4..7 | 2ⁿ (8..64)| |

**cyclic-n caveat:** cyclic-n is zero-dimensional only when *n* is squarefree (Backelin). The
squarefree small cases are n ∈ {5, 6, 7, 10, 11}; n ∈ {4, 8, 9, 12, …} have positive-dimensional
components and are **not** valid zero-dim test cases, so they are not registered. The `cyclic(n)`
generator accepts any n, but only register squarefree ones.

---

## Results

Fill these in from your runs (machine, date, and `bertini2`/Bertini 1 versions matter — record
them). `b1/b2` is the wall-time ratio (>1 means bertini2 was faster).

### Serial (ranks=1, threads=1)

| System  | bertini2 (s) | bertini1 (s) | b1/b2 | counts agree |
|---------|--------------|--------------|-------|--------------|
| cyclic6 |              |              |       |              |
| cyclic7 |              |              |       |              |
| katsura5|              |              |       |              |

### MPI, single thread (sweep ranks)

| System  | ranks | bertini2 (s) | bertini1 (s) | b1/b2 |
|---------|-------|--------------|--------------|-------|
| cyclic7 | 1     |              |              |       |
| cyclic7 | 2     |              |              |       |
| cyclic7 | 4     |              |              |       |

Environment to record alongside results: CPU model, core count, OS, `bertini2` version/commit, the
Bertini 1 version, MPI implementation, and the `--mptype`/`--predictor` used.

---

## Files

| File | Purpose |
|------|---------|
| `run_comparison.py` | Driver: emit shared input, run each solver, time, parse, tabulate, write CSV |
| `systems.py`        | Python generators for the test bed (cyclic-n, Katsura-n, diagonal) |
| `solvers.py`        | Solver adapters (`bertini2`, `bertini1`) behind one `run() -> RunResult` interface |

## Future work

Add a HomotopyContinuation.jl adapter to `solvers.py` (and an entry to `ADAPTERS`). It writes its
own problem from the same system rather than the classic input, so the "same settings" guarantee
will need an equivalent mapping of tracking knobs — note that explicitly when added.
