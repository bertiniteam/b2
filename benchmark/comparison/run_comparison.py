#!/usr/bin/env python3
"""Benchmark the bertini2 CLI against other zero-dim solvers (Bertini 1 today) on the same systems.

Each test system is built in Python with pybertini and emitted ONCE to a classic Bertini-1 input
file (System.to_classic_input). That single file is fed to every solver, so the problem and all
tracking settings (precision mode, predictor, tolerances, step cadence) are identical across
solvers by construction. Two phases:

  * serial            -- one process, one thread, every system through every solver;
  * MPI single-thread -- sweep --ranks with OMP_NUM_THREADS=1, each solver under `mpirun -n N`.

We get a correctness/agreement check for free: every solver should report the same solution count,
equal to the system's known expected count. That check is done AFTER timing, from output files, and
is NEVER part of any reported wall time.

Usage:
    python benchmark/comparison/run_comparison.py \
        --bertini2 ./build/core/bertini2 \
        --bertini1 /path/to/bertini \
        --systems cyclic6 katsura5 \
        --ranks 1 2 4 \
        --output comparison_results.csv

Requirements:
    - an MPI-built bertini2 and (for the bertini1 column) an MPI-built Bertini 1
    - mpirun on PATH for the rank sweep (--ranks beyond 1)
    - pybertini importable (the `bertini` package)
"""

import argparse
import csv
import datetime
import math
import os
import platform
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))           # benchmark/comparison -> repo root
DEFAULT_HISTORY = os.path.join(HERE, "history.csv")          # committed, append-only

sys.path.insert(0, HERE)
import solvers          # noqa: E402
import systems          # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description="Compare bertini2 vs other solvers on the same systems")
    p.add_argument("--bertini2", default="./build/core/bertini2",
                   help="Path to the bertini2 executable (default: ./build/core/bertini2)")
    p.add_argument("--bertini1", default=None,
                   help="Path to a Bertini 1 executable. If omitted, only bertini2 is run.")
    p.add_argument("--systems", nargs="+", default=systems.DEFAULT_SYSTEMS,
                   help=f"System names, or 'all' (default: {' '.join(systems.DEFAULT_SYSTEMS)}). "
                        f"Known: {', '.join(sorted(systems.SYSTEMS))}")
    p.add_argument("--ranks", nargs="+", type=int, default=[1],
                   help="MPI rank counts (single-thread) to sweep (default: 1 = serial only)")
    p.add_argument("--repeats", type=int, default=1,
                   help="Timed repeats per (system, solver, ranks); the fastest is kept (default: 1)")
    p.add_argument("--timeout", type=float, default=600.0, help="Per-run timeout, seconds (default: 600)")
    p.add_argument("--mpirun", default="mpirun", help="MPI launcher for rank sweeps (default: mpirun)")
    p.add_argument("--mpirun-args", default="--bind-to none",
                   help="Extra flags before -n, one quoted string (default: '--bind-to none')")
    p.add_argument("--output", default="comparison_results.csv",
                   help="Per-run snapshot CSV (overwritten each run)")
    p.add_argument("--history", default=DEFAULT_HISTORY,
                   help=f"Append-only history CSV, committed to the repo so performance can be "
                        f"tracked over time (default: {os.path.relpath(DEFAULT_HISTORY, REPO_ROOT)})")
    p.add_argument("--no-history", action="store_true",
                   help="Do not append to the history file (just write the per-run --output CSV)")
    p.add_argument("--note", default="",
                   help="Free-text note recorded with each history row (e.g. 'after SLP rework')")
    # Tracking knobs, emitted into the shared classic input (identical for every solver).
    p.add_argument("--mptype", type=int, default=2, help="0 double, 1 fixed-multiple, 2 adaptive (default 2)")
    p.add_argument("--predictor", type=int, default=5, help="odepredictor: 0 Euler, 2 RK4, 5 RKF45 (default 5)")
    return p.parse_args()


def _git_commit():
    """Short b2 commit hash, with a -dirty suffix if the working tree has uncommitted changes."""
    try:
        commit = subprocess.run(["git", "-C", REPO_ROOT, "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", REPO_ROOT, "status", "--porcelain"],
                               capture_output=True, text=True).stdout.strip()
        return (commit or "unknown") + ("-dirty" if dirty else "")
    except OSError:
        return "unknown"


def _cpu_brand():
    """Human-readable CPU model -- timings are only comparable across rows with the same CPU."""
    try:
        if platform.system() == "Darwin":
            return subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                                  capture_output=True, text=True).stdout.strip() or "unknown"
        if platform.system() == "Linux":
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("model name"):
                        return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def _solver_version(exe):
    """First line of `<solver> --version`; 'unknown' if it doesn't support the flag."""
    try:
        out = subprocess.run([os.path.abspath(exe), "--version"],
                             capture_output=True, text=True, timeout=15)
        text = (out.stdout + out.stderr).strip()
        return text.splitlines()[0][:120] if text else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def build_metadata(solver_cols, mptype, predictor, note):
    """Per-invocation provenance recorded on every history row (machine, date, versions, settings)."""
    return {
        "timestamp": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "host": platform.node(),
        "cpu": _cpu_brand(),
        "os": platform.platform(),
        "b2_commit": _git_commit(),
        "mptype": mptype,
        "predictor": predictor,
        "note": note,
        "_versions": {name: _solver_version(exe) for name, exe in solver_cols},
    }


HISTORY_FIELDS = ["timestamp", "host", "cpu", "os", "b2_commit", "mptype", "predictor",
                  "solver", "solver_version", "system", "ranks", "threads",
                  "wall_time_s", "solutions_found", "expected_count", "matches_expected", "note"]


def append_history(path, rows, meta):
    """Append one history row per measurement; write the header if the file is new."""
    is_new = not os.path.exists(path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HISTORY_FIELDS)
        if is_new:
            w.writeheader()
        for r in rows:
            w.writerow({
                "timestamp": meta["timestamp"], "host": meta["host"], "cpu": meta["cpu"],
                "os": meta["os"], "b2_commit": meta["b2_commit"],
                "mptype": meta["mptype"], "predictor": meta["predictor"], "note": meta["note"],
                "solver": r["solver"], "solver_version": meta["_versions"].get(r["solver"], ""),
                "system": r["system"], "ranks": r["ranks"], "threads": r["threads"],
                "wall_time_s": r["wall_time_s"], "solutions_found": r["solutions_found"],
                "expected_count": r["expected_count"], "matches_expected": r["matches_expected"],
            })


def emit_input(system, mptype, predictor):
    """Build the single classic-Bertini input string shared by all solvers for this system."""
    return system.to_classic_input(mptype=mptype, odepredictor=predictor)


def timed(adapter, exe, input_text, ranks, mpirun, mpirun_args, timeout, repeats):
    """Run one (solver, system, ranks) cell `repeats` times; keep the fastest valid time.

    Solution count comes from the runs (it is constant across repeats); timing is the minimum.
    """
    best_time = float("nan")
    solutions = -1
    detail = ""
    for _ in range(repeats):
        r = adapter(exe, input_text, ranks=ranks, threads=1,
                    mpirun=mpirun, mpirun_args=mpirun_args, timeout=timeout)
        detail = r.detail
        if r.ok:
            solutions = r.solutions_found
            if math.isnan(best_time) or r.wall_time_s < best_time:
                best_time = r.wall_time_s
    return best_time, solutions, detail


def main():
    args = parse_args()

    names = sorted(systems.SYSTEMS) if args.systems == ["all"] else args.systems
    for n in names:
        if n not in systems.SYSTEMS:
            sys.exit(f"Error: unknown system {n!r}. Known: {', '.join(sorted(systems.SYSTEMS))}")
    def resolve(exe):
        """Accept either a path to an executable or a bare command resolved on PATH."""
        return exe if os.path.isfile(exe) else shutil.which(exe)

    bertini2 = resolve(args.bertini2)
    if not bertini2:
        sys.exit(f"Error: bertini2 not found: {args.bertini2}")

    if args.bertini1:
        bertini1 = resolve(args.bertini1)
        if not bertini1:
            sys.exit(f"Error: bertini1 not found: {args.bertini1}")
    else:
        # Auto-detect Bertini 1, conventionally installed as `bertini` on PATH.
        bertini1 = shutil.which("bertini")
        if bertini1:
            print(f"Auto-detected Bertini 1 on PATH: {bertini1}  (use --bertini1 to override)")

    ranks_list = sorted(set(args.ranks) | {1})              # always include the serial baseline
    needs_mpi = any(r > 1 for r in ranks_list)
    if needs_mpi and shutil.which(args.mpirun) is None:
        sys.exit(f"Error: {args.mpirun!r} not on PATH but --ranks includes >1")
    if bertini1 is None:
        print("NOTE: no Bertini 1 found; running bertini2 only (no cross-solver comparison).\n")

    solver_cols = [("bertini2", bertini2)]
    if bertini1:
        solver_cols.append(("bertini1", bertini1))

    print(f"systems:   {names}")
    print(f"solvers:   {[s for s, _ in solver_cols]}")
    print(f"ranks:     {ranks_list}   (threads fixed at 1)")
    print(f"settings:  mptype={args.mptype}, odepredictor={args.predictor}  (identical for all solvers)")
    print(f"repeats:   {args.repeats}\n")

    rows = []
    for name in names:
        system, expected = systems.build(name)
        input_text = emit_input(system, args.mptype, args.predictor)

        for ranks in ranks_list:
            for solver_name, exe in solver_cols:
                adapter = solvers.ADAPTERS[solver_name]
                label = f"{name:10s} {solver_name:9s} ranks={ranks}"
                print(f"Running {label} ... ", end="", flush=True)
                t, found, detail = timed(adapter, exe, input_text, ranks,
                                         args.mpirun, args.mpirun_args, args.timeout, args.repeats)
                match = "yes" if (found == expected) else "NO"
                if math.isnan(t):
                    print(f"failed ({detail})")
                else:
                    flag = "" if match == "yes" else f"  <-- count {found} != expected {expected}"
                    print(f"{t:8.3f}s  ({found} solutions){flag}")
                rows.append({
                    "system": name, "solver": solver_name, "ranks": ranks, "threads": 1,
                    "wall_time_s": f"{t:.4f}" if not math.isnan(t) else "nan",
                    "solutions_found": found, "expected_count": expected,
                    "matches_expected": match, "detail": detail,
                })

    fieldnames = ["system", "solver", "ranks", "threads", "wall_time_s",
                  "solutions_found", "expected_count", "matches_expected", "detail"]
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    if not args.no_history:
        meta = build_metadata(solver_cols, args.mptype, args.predictor, args.note)
        append_history(args.history, rows, meta)

    _print_summary(rows, names, ranks_list, solver_cols)
    print(f"\nPer-run results written to: {args.output}")
    if not args.no_history:
        print(f"Appended {len(rows)} row(s) to history:  {args.history}  "
              f"(commit this file to track performance over time)")

    # Correctness gate is reported, not timed: did everything that ran match its expected count?
    bad = [r for r in rows if r["wall_time_s"] != "nan" and r["matches_expected"] != "yes"]
    if bad:
        print(f"\nWARNING: {len(bad)} run(s) did not match the expected solution count "
              f"(see matches_expected column).")


def _print_summary(rows, names, ranks_list, solver_cols):
    """Side-by-side timing table per system & rank, plus the bertini2/bertini1 speed ratio."""
    solver_names = [s for s, _ in solver_cols]
    print()
    header = f"{'system':10s} {'ranks':>5s} " + "".join(f"{s+'(s)':>13s}" for s in solver_names)
    if {"bertini2", "bertini1"} <= set(solver_names):
        header += f"{'b1/b2':>9s}"
    print(header)
    print("-" * len(header))
    by_key = {(r["system"], r["solver"], r["ranks"]): r for r in rows}
    for name in names:
        for ranks in ranks_list:
            line = f"{name:10s} {ranks:5d} "
            times = {}
            for s in solver_names:
                r = by_key.get((name, s, ranks))
                t = r["wall_time_s"] if r else "nan"
                times[s] = t
                line += f"{t:>13s}"
            if {"bertini2", "bertini1"} <= set(solver_names):
                try:
                    ratio = float(times["bertini1"]) / float(times["bertini2"])
                    line += f"{ratio:>9.2f}"
                except (ValueError, ZeroDivisionError):
                    line += f"{'-':>9s}"
            print(line)


if __name__ == "__main__":
    main()
