# Bertini2 parallel benchmark

Measure how much faster `bertini2` solves a system when you give it more compute. There are
two independent kinds of parallelism, and this benchmark can sweep either or both:

- **Shared-memory threads** — one process, many worker threads, all on one machine. Needs **no
  MPI**: it works with a plain `pip install` build or any `bertini2` compiled without MPI. This
  is the everyday "use all the cores on my laptop/workstation" case.
- **MPI ranks** — multiple processes, typically spread across the nodes of a cluster, each rank
  itself optionally running a pool of threads. This is the "scale across a supercomputer" case
  and requires a `bertini2` built with MPI plus `mpirun`.

The driver, `run_benchmark.py`, times each configuration, prints a speedup table, writes a CSV,
and (optionally) fails if parallelism didn't actually help.

---

## 1. Run it on a single machine (no MPI, just threads)

This is the simplest case and the one most people want. You only need a `bertini2` executable
and Python 3 (no extra packages).

```bash
# from the repo root
python benchmark/run_benchmark.py \
    --bertini2 ./build/core/bertini2 \
    --input    benchmark/inputs/medium.b2 \
    --no-mpi \
    --threads  1 2 4 8
```

`--no-mpi` runs the solver directly (no `mpirun`) and sweeps the thread counts you list. The
run at `threads=1` is the serial baseline; every other run's speedup is reported relative to it.
The thread count is passed to the solver through `OMP_NUM_THREADS`, which the script sets for
you — don't set it yourself.

Example output:

```
 ranks  threads  workers    time(s)  solutions   speedup
------------------------------------------------------------
     1        1        1    12.6141          6    1.0000
     1        2        2    11.9876          6    1.0523
     1        4        4     3.4350          6    3.6722
     1        8        8     2.5991          6    4.8533

Best parallel speedup: 4.85x vs serial baseline.
```

The `solutions` column should be identical across every row: threading is a speed feature, never
a correctness change. If it ever differs, that's a bug — please report it.

### Make the run *fail* if threads don't help

Handy for CI or an acceptance check on a new machine:

```bash
python benchmark/run_benchmark.py --bertini2 ./build/core/bertini2 \
    --input benchmark/inputs/medium.b2 --no-mpi --threads 1 4 \
    --assert-speedup 1.5
```

This exits non-zero unless the best multi-thread run beats serial by at least 1.5×. Pick a
threshold that suits the machine: a 2-core CI runner will not hit 4×, while a 32-core
workstation should comfortably exceed it. Leave `--assert-speedup` off to just measure.

---

## 2. Run it with MPI (multiple ranks)

If your `bertini2` was built with MPI (`cmake -DENABLE_MPI=ON ...`) and `mpirun` is on your
`PATH`, you can sweep MPI ranks as well as threads. Drop `--no-mpi` and add `--ranks`:

```bash
python benchmark/run_benchmark.py \
    --bertini2 ./build/core/bertini2 \
    --input    benchmark/inputs/large.b2 \
    --ranks    1 2 4 \
    --threads  1 2 4
```

This launches each combination under `mpirun -n <ranks>` with `OMP_NUM_THREADS=<threads>` per
rank. `--bind-to none` is passed to `mpirun` automatically so threaded workers can spread across
a node's cores without affinity conflicts. The total worker count in the table is
`ranks × threads`.

You do **not** need a cluster for this — a single multi-core machine with MPI installed runs the
rank sweep locally just fine.

---

## 3. Run it on an HPC cluster

### With SLURM (`sbatch`)

`submit.slurm.sh` is a ready-to-edit batch script. Open it and adjust:

- the `#SBATCH` directives — `--nodes`, `--ntasks-per-node`, `--cpus-per-task` (this is the
  cores-per-rank, i.e. the max threads), `--time`;
- the **module loads / environment activation** near the top (load your MPI, your Python, and/or
  activate the conda env that has `bertini2`);
- the `INPUT`, `RANKS`, and `THREADS` sweep near the bottom.

Keep two consistency rules in mind:

- `RANKS` should not exceed `--nodes × --ntasks-per-node`.
- the largest value in `THREADS` should not exceed `--cpus-per-task`.

Then submit:

```bash
sbatch benchmark/submit.slurm.sh
```

Results land in `benchmark/results_<JOBID>.csv`, and stdout/stderr in `benchmark_<JOBID>.log` /
`.err`. The script discovers the repo root relative to its own location, so you can submit it
from anywhere.

### On a cluster without SLURM (or interactively)

There's nothing SLURM-specific about the measurement — `submit.slurm.sh` just sets up an
environment and calls `run_benchmark.py`. On a cluster with a different scheduler (PBS, LSF,
...), or in an interactive allocation, load your modules / activate your environment by hand and
run the driver directly, exactly as in section 2:

```bash
python benchmark/run_benchmark.py --bertini2 /path/to/bertini2 \
    --input benchmark/inputs/large.b2 --ranks 1 2 4 8 --threads 1 4 8 \
    --mpirun srun        # if your site launches MPI jobs with srun instead of mpirun
```

Use `--mpirun` to point at whatever launcher your site uses (e.g. `srun`, or a full path to a
specific `mpirun`), and `--mpirun-args` to pass site-specific flags. For example, on a host
whose CPU topology `hwloc` cannot read (some VMs/containers), OpenMPI's default mapper fails
with "failed to map" / "all nodes already filled"; pass

```bash
--mpirun-args "--map-by slot:OVERSUBSCRIBE --bind-to none"
```

to map by slot instead of by core. On a normal cluster the default (`--bind-to none`) is fine.

---

## Sample input files

Each system is diagonal (`xi^3 - ci = 0`) with distinct prime constants, so the answers are
known analytically and easy to check. "Paths" is the total-degree path count the solver tracks
(the work), which is what scales with parallelism.

| File | Variables | Degree | Paths |
|------|-----------|--------|-------|
| `inputs/small.b2`   | 3 | 3 |    27 |
| `inputs/medium.b2`  | 5 | 3 |   243 |
| `inputs/large.b2`   | 6 | 3 |   729 |
| `inputs/xlarge.b2`  | 7 | 3 |  2187 |
| `inputs/huge.b2`    | 9 | 3 | 19683 |

Bigger systems show parallel scaling more clearly (more paths to spread over workers) but take
longer. Start with `small`/`medium` to sanity-check your setup, then move up.

---

## All options

```
--bertini2 PATH      Path to the bertini2 executable (default: ./build/core/bertini2)
--input    FILE      Bertini input file to solve (required)
--threads  N ...     Thread counts to sweep (default: 1)
--ranks    N ...     MPI rank counts to sweep (default: 1 2 4); ignored with --no-mpi
--no-mpi             Run the solver directly, no mpirun: a pure thread sweep for a build
                     without MPI. Forces ranks=1 and sweeps --threads only.
--assert-speedup F   Exit non-zero unless the best parallel run beats serial by >= F (e.g. 1.5).
                     Off by default.
--repeats  N         Timed repeats per combination; the fastest is recorded (default: 1).
--timeout  SECS      Per-run timeout (default: 600).
--mpirun   CMD       Launcher for MPI runs (default: mpirun; e.g. srun on some clusters).
--mpirun-args STR    Extra flags for mpirun, one quoted string (default: '--bind-to none').
--output   FILE      CSV output path (default: benchmark_results.csv).
```

## Notes

- The serial baseline (`ranks=1, threads=1`) always runs first; all speedups are relative to it.
- Each run executes in its own temp directory, so concurrent files never collide.
- Timings on a busy or thermally-throttling machine are noisy. Use `--repeats 3+` (the minimum
  time is kept) and a quiet machine for numbers you intend to publish.
- The `*_seeded` and `baseline_*` CSVs checked in here are reference results from past runs; your
  numbers will differ by machine.
```
