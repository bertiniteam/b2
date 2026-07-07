#!/usr/bin/env python
"""Re-measure the distributed-solve timings behind the ``solving_at_scale`` tutorial.

The tutorial ``python/docs/source/tutorials/solving_at_scale/index.rst`` shows wall-clock +
speedup tables for the ``solve_cyclic.py`` and ``solve_eigenvalues.py`` example scripts at several
rank/thread layouts.  Those numbers are hardware- and version-specific and go stale -- especially
after performance work -- so this tool re-runs the scripts and regenerates the data behind the
tables.

The numbers do **not** live in the prose.  This tool writes them as *data*:

* one CSV per table (``cyclic_timings.csv``, ``eigen_timings.csv``, ``hybrid_timings.csv``) that
  the tutorial pulls in with ``.. csv-table:: :file:``; and
* ``_timing_data.txt``, a set of ``.. |tw-...| replace::`` substitution definitions the tutorial
  ``.. include::``s, carrying every scalar the prose quotes (paths, finite count, top speedups,
  the eigenvalue plateau, the host description, core counts, the date).

so re-running this tool updates every table cell *and* every number in the surrounding sentences at
once -- there is no number to hand-edit in the ``.rst``.

Run it **from the repository root, inside the bertini environment**::

    python tools/update_scaling_timings.py             # full sweep, rewrite the data files
    python tools/update_scaling_timings.py --dry-run   # measure + print, do NOT write
    python tools/update_scaling_timings.py --only eigen # refresh just one table's data
    python tools/update_scaling_timings.py --quick      # tiny problems, to smoke-test this tool

The host description (``|tw-host|``, e.g. "a 16-core Apple M3 Max (12 performance + 4 efficiency
cores)") and the performance-core count (``|tw-perf-cores|``) are **preserved** across runs unless
you override them with ``--host`` / ``--perf-cores`` -- they describe the reference machine, which
this tool cannot fully introspect.  The core count and date come from the machine and the clock.

Each configuration is run **alone** (sequentially), so wall-clock numbers are not polluted by
contention; the full default sweep is a few minutes on a fast 12+-core machine.  The tool **aborts
without writing** if any run fails its built-in correctness check, so it can never publish numbers
from a broken solve.  Prefer driving it through ``tools/refresh_doc_artifacts.py`` (the single
artifact entry point).
"""

import argparse
import csv
import datetime
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "python" / "docs" / "source" / "tutorials" / "parallelism" / "solving_at_scale"
SUBS_FILE = OUT / "_timing_data.txt"
CYCLIC = "python/examples/solve_cyclic.py"
EIGEN = "python/examples/solve_eigenvalues.py"

WALL_RE = re.compile(r"wall=([0-9.]+)s")
CYCLIC_PATHS_RE = re.compile(r"paths tracked=(\d+)")
CYCLIC_FINITE_RE = re.compile(r"finite solutions=(\d+)")

# Canonical order + full set of substitutions the tutorial references.  A partial run (--only ...)
# updates just the keys it measures and PRESERVES the rest from the existing file, so the page is
# never left half-defined.
SUBS_ORDER = [
    "tw-host", "tw-cores", "tw-perf-cores", "tw-date",
    "tw-cyclic-paths", "tw-cyclic-finite", "tw-cyclic-top-speedup",
    "tw-eigen-paths", "tw-eigen-floor-wall", "tw-eigen-floor-speedup", "tw-eigen-12w-speedup",
]

_cache = {}                       # (script, args, nprocs, omp, bind) -> (wall_seconds, output_text)


def launch(script, args, nprocs, omp, bind_none, timeout):
    """Run one configuration alone; return (wall_seconds, combined_output).

    Aborts the whole tool if the run fails or does not print its ``OK:`` correctness line --
    we never want to publish a number from a solve that did not verify.
    """
    env = dict(os.environ, OMP_NUM_THREADS=str(omp))
    if nprocs <= 1:
        cmd = [sys.executable, script, *args]
    else:
        cmd = ["mpirun", "-n", str(nprocs)]
        # rank 0 is a near-idle manager, so a worker-per-core pool wants nprocs = cores + 1;
        # let mpirun place that extra rank when we are short a slot.
        if nprocs > (os.cpu_count() or nprocs):
            cmd += ["--map-by", ":OVERSUBSCRIBE"]
        if bind_none:
            cmd += ["--bind-to", "none"]
        cmd += [sys.executable, script, *args]

    shown = (f"OMP_NUM_THREADS={omp} " if omp != 1 else "") + " ".join(cmd)
    print(f"    $ {shown}", flush=True)
    proc = subprocess.run(cmd, env=env, cwd=REPO, capture_output=True, text=True, timeout=timeout)
    out = proc.stdout + proc.stderr
    m = WALL_RE.search(out)
    if proc.returncode != 0 or "OK:" not in out or not m:
        return None, out      # caller decides whether to retry or abort
    wall = float(m.group(1))
    print(f"      -> {wall:.1f}s", flush=True)
    return wall, out


# A correct double-precision solve is occasionally subject to a rare path failure under an unlucky
# random homotopy (a fresh gamma each run -- see ADR-0017); each retry draws a new one.  A couple of
# retries absorb that tail without masking a systematic problem (3 strikes -> abort, never publish a
# number from a run that did not pass its own correctness check).
MAX_ATTEMPTS = 3


def measure(script, args, nprocs, omp=1, bind_none=False, timeout=7200):
    key = (script, tuple(args), nprocs, omp, bind_none)
    if key not in _cache:
        last_out = ""
        for attempt in range(MAX_ATTEMPTS):
            wall, last_out = launch(script, list(args), nprocs, omp, bind_none, timeout)
            if wall is not None:
                _cache[key] = (wall, last_out)
                break
            print(f"      (attempt {attempt + 1}/{MAX_ATTEMPTS} failed its correctness check; retrying)",
                  flush=True)
        else:
            sys.exit(f"\nABORT: a config failed {MAX_ATTEMPTS} times; refusing to write stale/garbage "
                     f"numbers.\n--- last output ---\n{last_out[-2500:]}")
    return _cache[key]


def speedup(serial, t):
    return f"{serial / t:.1f}x"


def ladder_rows(script, args):
    """The serial -> 2/4/8/12-worker ladder shared by the cyclic and eigenvalue tables.

    Returns ``(serial_wall, rows, by_workers)`` where ``rows`` are the CSV rows
    ``[launch, workers, wall, speedup]`` and ``by_workers`` maps worker count -> (wall, speedup).

    The top rung (12 workers, ``-n 13``) targets a 12-performance-core machine (e.g. Apple M3 Max);
    on a smaller box the wide rungs oversubscribe -- measure() adds ``:OVERSUBSCRIBE`` -- and the
    speedup simply plateaus, which is the honest result.
    """
    serial, _ = measure(script, args, 1)
    rows = [["serial", "--", f"{serial:.1f}", "1.0x"]]
    by_workers = {}
    for nprocs, workers in ((3, 2), (5, 4), (9, 8), (13, 12)):
        t, _ = measure(script, args, nprocs)
        sp = speedup(serial, t)
        rows.append([f"``-n {nprocs}``", str(workers), f"{t:.1f}", sp])
        by_workers[workers] = (t, sp)
    return serial, rows, by_workers


def write_csv(path, header, rows, dry_run):
    print(f"  -> {os.path.relpath(path, REPO)}")
    for row in rows:
        print("     " + " | ".join(row))
    if dry_run:
        return
    with open(path, "w", newline="", encoding="utf-8") as fh:
        # csv.writer quotes only the fields that need it (the hybrid 'layout' cells contain commas),
        # matching what csv-table's RFC4180 reader expects.
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def load_subs(path):
    """Parse an existing ``_timing_data.txt`` into {name: value} so a partial run preserves the rest."""
    subs = {}
    if not path.exists():
        return subs
    for line in path.read_text(encoding="utf-8").splitlines():
        m = re.match(r"\.\. \|([\w-]+)\| replace:: (.*)", line)
        if m:
            subs[m.group(1)] = m.group(2)
    return subs


def write_subs(path, subs, dry_run):
    lines = [
        ".. Machine-generated by tools/update_scaling_timings.py -- do not hand-edit.",
        "   These substitutions carry every measured number in solving_at_scale out of the prose.",
        "",
    ]
    for name in SUBS_ORDER:
        if name in subs:
            lines.append(f".. |{name}| replace:: {subs[name]}")
    body = "\n".join(lines) + "\n"
    print(f"  -> {os.path.relpath(path, REPO)}")
    for name in SUBS_ORDER:
        if name in subs:
            print(f"     |{name}| = {subs[name]}")
    if not dry_run:
        path.write_text(body, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cyclic-n", type=int, default=6, help="which cyclic-n (default 6)")
    parser.add_argument("--eigen-size", type=int, default=24, help="matrix size n (default 24)")
    parser.add_argument("--only", choices=["all", "cyclic", "eigen", "hybrid"], default="all",
                        help="refresh just one table's data (default all)")
    parser.add_argument("--quick", action="store_true",
                        help="tiny problems (cyclic-5, eigen-6) to smoke-test this tool")
    parser.add_argument("--seed", type=int, default=1,
                        help="RNG seed passed to both scripts (default 1).  A fixed seed makes "
                             "every rank/thread layout solve the IDENTICAL homotopy, so the "
                             "speedups compare like with like.")
    parser.add_argument("--host", default=None,
                        help="reference-machine description for |tw-host| (preserved if omitted)")
    parser.add_argument("--perf-cores", type=int, default=None,
                        help="performance-core count for |tw-perf-cores| (preserved if omitted)")
    parser.add_argument("--dry-run", action="store_true",
                        help="measure and print, but do not write the data files")
    args = parser.parse_args()

    if args.quick:
        args.cyclic_n, args.eigen_size = 5, 6

    cyclic_args = ["--n", str(args.cyclic_n), "--seed", str(args.seed)]
    eigen_args = ["--size", str(args.eigen_size), "--seed", str(args.seed)]
    cores = os.cpu_count()
    want = {args.only} if args.only != "all" else {"cyclic", "eigen", "hybrid"}

    # start from what is on disk so a partial run preserves the other keys / the descriptive host
    subs = load_subs(SUBS_FILE)
    subs["tw-date"] = datetime.date.today().isoformat()
    subs["tw-cores"] = str(cores)
    if args.host is not None:
        subs["tw-host"] = args.host
    subs.setdefault("tw-host", f"a {cores}-core machine")
    if args.perf_cores is not None:
        subs["tw-perf-cores"] = str(args.perf_cores)
    subs.setdefault("tw-perf-cores", str(cores))

    header4 = ["launch", "workers", "wall-clock (s)", "speedup"]

    if "cyclic" in want or "hybrid" in want:
        # cyclic serial + ladder; needed by the cyclic table and (for the speedup base) the hybrid one
        _, serial_out = measure(CYCLIC, cyclic_args, 1)
        subs["tw-cyclic-paths"] = CYCLIC_PATHS_RE.search(serial_out).group(1)
        subs["tw-cyclic-finite"] = CYCLIC_FINITE_RE.search(serial_out).group(1)

    if "cyclic" in want:
        serial, rows, by_workers = ladder_rows(CYCLIC, cyclic_args)
        subs["tw-cyclic-top-speedup"] = by_workers[12][1]
        write_csv(OUT / "cyclic_timings.csv", header4, rows, args.dry_run)

    if "eigen" in want:
        serial, rows, by_workers = ladder_rows(EIGEN, eigen_args)
        subs["tw-eigen-paths"] = str(args.eigen_size)
        # the floor is the fastest rung (highest speedup); the plateau the prose quotes
        floor_workers = max(by_workers, key=lambda w: float(by_workers[w][1].rstrip("x")))
        floor_wall, floor_sp = by_workers[floor_workers]
        subs["tw-eigen-floor-wall"] = f"{floor_wall:.1f}"
        subs["tw-eigen-floor-speedup"] = floor_sp
        subs["tw-eigen-12w-speedup"] = by_workers[12][1]
        write_csv(OUT / "eigen_timings.csv", header4, rows, args.dry_run)

    if "hybrid" in want:
        serial, _ = measure(CYCLIC, cyclic_args, 1)
        # 12 worker-cores split every way the factors of 12 allow: ranks x threads = 12.
        t12_1, _ = measure(CYCLIC, cyclic_args, 13, omp=1)                # == cyclic -n 13 row
        t6_2, _ = measure(CYCLIC, cyclic_args, 7, omp=2, bind_none=True)
        t4_3, _ = measure(CYCLIC, cyclic_args, 5, omp=3, bind_none=True)
        t2_6, _ = measure(CYCLIC, cyclic_args, 3, omp=6, bind_none=True)
        rows = [
            ["``-n 13``, 12 workers x 1 thread", f"{t12_1:.1f}", speedup(serial, t12_1)],
            ["``-n 7``, 6 workers x 2 threads", f"{t6_2:.1f}", speedup(serial, t6_2)],
            ["``-n 5``, 4 workers x 3 threads", f"{t4_3:.1f}", speedup(serial, t4_3)],
            ["``-n 3``, 2 workers x 6 threads", f"{t2_6:.1f}", speedup(serial, t2_6)],
        ]
        write_csv(OUT / "hybrid_timings.csv", ["layout", "wall-clock (s)", "speedup"], rows, args.dry_run)

    write_subs(SUBS_FILE, subs, args.dry_run)
    print("\n(dry run: nothing written)" if args.dry_run else f"\nWrote data files under {os.path.relpath(OUT, REPO)}")


if __name__ == "__main__":
    main()
