# Bertini2 Parallel Benchmark

Scripts for measuring MPI + thread parallel speedup of the `bertini2` solver.

## Prerequisites

- `bertini2` built with MPI support (`cmake -DENABLE_MPI=ON ...`)
- `mpirun` on `PATH`
- Python 3.7+, no external packages needed

## Quick Start

```bash
# From the repo root
python benchmark/run_benchmark.py \
    --bertini2 ./build/core/bertini2 \
    --input    benchmark/inputs/small.b2 \
    --ranks    1 2 \
    --threads  1 2
```

This runs the `small.b2` system (27 paths) with 4 combinations of ranks × threads,
prints a speedup table, and writes `benchmark_results.csv`.

## Sample Input Files

| File | Variables | Degree | Paths |
|------|-----------|--------|-------|
| `inputs/small.b2`   | 3 | 3 |    27 |
| `inputs/medium.b2`  | 5 | 3 |   243 |
| `inputs/large.b2`   | 6 | 3 |   729 |
| `inputs/xlarge.b2`  | 7 | 3 |  2187 |
| `inputs/huge.b2`    | 9 | 3 | 19683 |

Each system is diagonal (`xi^3 - ci = 0`) with distinct prime constants, so solutions
are known analytically and correctness is easy to verify (solution count = number of
paths for these fully real systems).

## Options

```
--bertini2 PATH    Path to bertini2 executable (default: ./build/core/bertini2)
--input    FILE    Bertini input file to solve (required)
--ranks    N ...   MPI rank counts to sweep (default: 1 2 4)
--threads  N ...   OMP thread counts per rank (default: 1)
--output   FILE    CSV output path (default: benchmark_results.csv)
--timeout  SECS    Per-run timeout (default: 600)
--repeats  N       Timed repeats per combo; records minimum (default: 1)
--mpirun   CMD     mpirun command (default: mpirun)
```

## Cluster Submission (SLURM)

Edit `submit.slurm.sh` (paths, module loads, rank/thread sweep) and submit:

```bash
sbatch benchmark/submit.slurm.sh
```

Results are written to `benchmark/results_<JOBID>.csv`.

## Notes

- `--bind-to none` is passed to `mpirun` automatically so that threaded workers
  can use all cores on a node without affinity conflicts.
- Each (ranks, threads) combo runs in its own temp directory to avoid file conflicts.
- The serial baseline (ranks=1, threads=1) is always run first; speedup is computed
  relative to that measurement.
- `OMP_NUM_THREADS` is set per-run by the script; do not set it externally when sweeping.
