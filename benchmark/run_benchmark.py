#!/usr/bin/env python3
"""
Benchmark parallel speedup of the bertini2 CLI solver.

Sweeps MPI rank counts and OMP thread counts, times each run, and emits
a CSV showing wall time and speedup relative to the serial (1 rank, 1 thread) baseline.

Usage:
    python run_benchmark.py \
        --bertini2 ./build/core/bertini2 \
        --input benchmark/inputs/medium.b2 \
        --ranks 1 2 4 8 \
        --threads 1 2 4 \
        --output results.csv

Requirements:
    - mpirun must be on PATH
    - bertini2 must be compiled with MPI support (BERTINI2_HAVE_MPI)
    - OMP_NUM_THREADS controls threads per MPI rank
"""

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time


def parse_args():
    p = argparse.ArgumentParser(description="Benchmark bertini2 parallel speedup")
    p.add_argument("--bertini2", default="./build/core/bertini2",
                   help="Path to bertini2 executable (default: ./build/core/bertini2)")
    p.add_argument("--input", required=True,
                   help="Bertini input file to solve")
    p.add_argument("--ranks", nargs="+", type=int, default=[1, 2, 4],
                   help="MPI rank counts to sweep (default: 1 2 4)")
    p.add_argument("--threads", nargs="+", type=int, default=[1],
                   help="OMP thread counts to sweep (default: 1)")
    p.add_argument("--output", default="benchmark_results.csv",
                   help="CSV output file (default: benchmark_results.csv)")
    p.add_argument("--timeout", type=float, default=600.0,
                   help="Per-run timeout in seconds (default: 600)")
    p.add_argument("--repeats", type=int, default=1,
                   help="Number of timed repeats per (ranks, threads) combo (default: 1)")
    p.add_argument("--mpirun", default="mpirun",
                   help="mpirun command (default: mpirun)")
    return p.parse_args()


def run_once(bertini2_path, input_file, ranks, threads, mpirun_cmd, timeout):
    """
    Run bertini2 with the given parallelism settings in a fresh temp directory.

    Returns (wall_time_s, solutions_found) or (float('nan'), -1) on failure.
    """
    tmpdir = tempfile.mkdtemp(prefix="b2_bench_")
    try:
        dest_input = os.path.join(tmpdir, "input")
        shutil.copy2(input_file, dest_input)

        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(threads)

        cmd = [mpirun_cmd, "-n", str(ranks), "--bind-to", "none",
               os.path.abspath(bertini2_path)]

        t0 = time.perf_counter()
        try:
            result = subprocess.run(
                cmd,
                cwd=tmpdir,
                env=env,
                capture_output=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT after {timeout}s", flush=True)
            return float("nan"), -1
        t1 = time.perf_counter()

        if result.returncode != 0:
            stderr_tail = result.stderr.decode(errors="replace")[-500:]
            print(f"  FAILED (exit {result.returncode}): {stderr_tail}", flush=True)
            return float("nan"), -1

        wall_time = t1 - t0
        solutions = _parse_solution_count(os.path.join(tmpdir, "main_data"))
        return wall_time, solutions

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _parse_solution_count(main_data_path):
    """Read solution count from the first line of main_data."""
    try:
        with open(main_data_path) as f:
            first_line = f.readline().strip()
        return int(first_line)
    except (OSError, ValueError):
        return -1


def validate_inputs(args):
    if not os.path.isfile(args.input):
        sys.exit(f"Error: input file not found: {args.input}")
    if not os.path.isfile(args.bertini2):
        sys.exit(f"Error: bertini2 executable not found: {args.bertini2}")
    if shutil.which(args.mpirun) is None:
        sys.exit(f"Error: {args.mpirun!r} not found on PATH")


def main():
    args = parse_args()
    validate_inputs(args)

    ranks_list = sorted(set(args.ranks))
    threads_list = sorted(set(args.threads))

    # Ensure serial baseline (1, 1) is always first
    if 1 not in ranks_list:
        ranks_list = [1] + ranks_list
    if 1 not in threads_list:
        threads_list = [1] + threads_list

    print(f"bertini2:  {args.bertini2}")
    print(f"input:     {args.input}")
    print(f"ranks:     {ranks_list}")
    print(f"threads:   {threads_list}")
    print(f"repeats:   {args.repeats}")
    print(f"output:    {args.output}")
    print()

    rows = []
    serial_time = None

    combos = [(r, t) for r in ranks_list for t in threads_list]
    # Run (1,1) first regardless of order so speedup can be computed incrementally
    combos = sorted(combos, key=lambda rt: (rt[0] != 1 or rt[1] != 1, rt[0], rt[1]))

    for ranks, threads in combos:
        label = f"ranks={ranks}, threads={threads}"
        times = []
        solutions = -1
        for rep in range(args.repeats):
            rep_label = f"  rep {rep+1}/{args.repeats}" if args.repeats > 1 else ""
            print(f"Running {label}{rep_label} ... ", end="", flush=True)
            t, sol = run_once(args.bertini2, args.input, ranks, threads,
                              args.mpirun, args.timeout)
            if math.isnan(t):
                print("failed")
                times.append(float("nan"))
            else:
                print(f"{t:.2f}s  ({sol} solutions)")
                times.append(t)
                solutions = sol

        valid = [t for t in times if not math.isnan(t)]
        wall_time = min(valid) if valid else float("nan")

        if ranks == 1 and threads == 1 and not math.isnan(wall_time):
            serial_time = wall_time

        if serial_time is not None and not math.isnan(wall_time):
            speedup = serial_time / wall_time
        else:
            speedup = float("nan")

        rows.append({
            "ranks": ranks,
            "threads": threads,
            "total_workers": ranks * threads,
            "wall_time_s": f"{wall_time:.4f}" if not math.isnan(wall_time) else "nan",
            "solutions_found": solutions,
            "speedup_vs_serial": f"{speedup:.4f}" if not math.isnan(speedup) else "nan",
        })

    # Write CSV
    fieldnames = ["ranks", "threads", "total_workers", "wall_time_s",
                  "solutions_found", "speedup_vs_serial"]
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Print summary table
    print()
    print(f"{'ranks':>6}  {'threads':>7}  {'workers':>7}  {'time(s)':>9}  {'solutions':>9}  {'speedup':>8}")
    print("-" * 60)
    for row in rows:
        print(f"{row['ranks']:>6}  {row['threads']:>7}  {row['total_workers']:>7}  "
              f"{row['wall_time_s']:>9}  {row['solutions_found']:>9}  {row['speedup_vs_serial']:>8}")

    print(f"\nResults written to: {args.output}")


if __name__ == "__main__":
    main()
