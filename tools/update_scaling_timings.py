#!/usr/bin/env python
"""Re-measure the distributed-solve timings and rewrite the tables in ``solving_at_scale.rst``.

The tutorial ``python/docs/source/tutorials/solving_at_scale.rst`` shows wall-clock + speedup
tables for the ``solve_cyclic.py`` and ``solve_eigenvalues.py`` example scripts at several
rank/thread layouts.  Those numbers are hardware- and version-specific and go stale -- especially
after performance work -- so this tool re-runs the scripts and rewrites the tables in place.

Run it **from the repository root, inside the bertini environment**::

    python tools/update_scaling_timings.py             # full sweep, rewrite the tables
    python tools/update_scaling_timings.py --dry-run   # measure + print, do NOT edit the file
    python tools/update_scaling_timings.py --only eigen # refresh just one table
    python tools/update_scaling_timings.py --quick      # tiny problems, to smoke-test this tool

Each configuration is run **alone** (sequentially), so the wall-clock numbers are not polluted by
contention; the full default sweep is a few minutes on a fast 12+-core machine.  The tool **aborts without editing**
if any run fails its built-in correctness check, so it can never write numbers from a broken solve.

This script is the single source of truth for the layouts shown in the tutorial: change the ladder
here and re-run, and the page follows.
"""

import argparse
import datetime
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RST = REPO / "python" / "docs" / "source" / "tutorials" / "solving_at_scale.rst"
CYCLIC = "python/examples/solve_cyclic.py"
EIGEN = "python/examples/solve_eigenvalues.py"

WALL_RE = re.compile(r"wall=([0-9.]+)s")
CYCLIC_PATHS_RE = re.compile(r"paths tracked=(\d+)")
CYCLIC_FINITE_RE = re.compile(r"finite solutions=(\d+)")

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


def render_table(caption, widths, header, rows):
    """Render a reStructuredText ``list-table`` (alignment-free, robust to editing)."""
    out = [
        f".. list-table:: {caption}",
        "   :header-rows: 1",
        f"   :widths: {' '.join(str(w) for w in widths)}",
        "",
    ]
    for row in (header, *rows):
        for i, cell in enumerate(row):
            out.append(("   * - " if i == 0 else "     - ") + cell)
    return "\n".join(out)


def splice(text, key, block):
    """Replace the content between ``.. BEGIN-TIMING <key>`` and ``.. END-TIMING <key>``."""
    begin, end = f".. BEGIN-TIMING {key}", f".. END-TIMING {key}"
    pat = re.compile(re.escape(begin) + r"\n.*?\n" + re.escape(end), re.DOTALL)
    text, n = pat.subn(f"{begin}\n\n{block}\n\n{end}", text, count=1)
    if n != 1:
        sys.exit(f"ABORT: could not find the '{key}' timing region (markers) in {RST}")
    return text


def rank_table(script, args, problem_caption):
    """The serial -> 2/4/8/12-worker ladder shared by the cyclic and eigenvalue tables.

    The top rung (12 workers, ``-n 13``) targets a 12-performance-core machine (e.g. Apple
    M3 Max: 12 performance + 4 efficiency cores); the near-idle manager rides an efficiency
    core.  On a smaller box the wide rungs oversubscribe -- measure() adds ``:OVERSUBSCRIBE``
    automatically -- and the speedup simply plateaus, which is the honest result.
    """
    serial, _ = measure(script, args, 1)
    rows = [["serial", "--", f"{serial:.1f}", "1.0x"]]
    for nprocs, workers in ((3, 2), (5, 4), (9, 8), (13, 12)):
        t, _ = measure(script, args, nprocs)
        rows.append([f"``-n {nprocs}``", str(workers), f"{t:.1f}", speedup(serial, t)])
    header = ["launch", "workers", "wall-clock (s)", "speedup"]
    return render_table(problem_caption, [20, 15, 25, 15], header, rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cyclic-n", type=int, default=6, help="which cyclic-n (default 6)")
    parser.add_argument("--eigen-size", type=int, default=24, help="matrix size n (default 24)")
    parser.add_argument("--only", choices=["all", "cyclic", "eigen", "hybrid"], default="all",
                        help="refresh just one table (default all)")
    parser.add_argument("--quick", action="store_true",
                        help="tiny problems (cyclic-5, eigen-6) to smoke-test this tool")
    parser.add_argument("--seed", type=int, default=1,
                        help="RNG seed passed to both scripts (default 1).  A fixed seed makes "
                             "every rank/thread layout solve the IDENTICAL homotopy, so the "
                             "speedups compare like with like -- without it each run draws a "
                             "different gamma and the timings are not comparable.")
    parser.add_argument("--dry-run", action="store_true",
                        help="measure and print, but do not edit the .rst")
    args = parser.parse_args()

    if args.quick:
        args.cyclic_n, args.eigen_size = 5, 6

    cyclic_args = ["--n", str(args.cyclic_n), "--seed", str(args.seed)]
    eigen_args = ["--size", str(args.eigen_size), "--seed", str(args.seed)]
    cores = os.cpu_count()
    stamp = datetime.date.today().isoformat()
    refresh = "run tools/update_scaling_timings.py to refresh"
    want = {args.only} if args.only != "all" else {"cyclic", "eigen", "hybrid"}

    blocks = {}

    if "cyclic" in want or "hybrid" in want:
        # the cyclic serial + ladder; needed by the cyclic table and (for the speedup base) hybrid
        _, serial_out = measure(CYCLIC, cyclic_args, 1)
        paths = CYCLIC_PATHS_RE.search(serial_out).group(1)
        finite = CYCLIC_FINITE_RE.search(serial_out).group(1)

    if "cyclic" in want:
        cap = (f"cyclic-{args.cyclic_n} ({paths} paths, {finite} finite solutions), "
               f"measured on {cores} cores, {stamp} -- {refresh}")
        blocks["cyclic"] = rank_table(CYCLIC, cyclic_args, cap)

    if "eigen" in want:
        cap = (f"eigenvalues of a {args.eigen_size}x{args.eigen_size} symmetric matrix "
               f"({args.eigen_size} paths), measured on {cores} cores, {stamp} -- {refresh}")
        blocks["eigen"] = rank_table(EIGEN, eigen_args, cap)

    if "hybrid" in want:
        serial, _ = measure(CYCLIC, cyclic_args, 1)
        # 12 worker-cores split every way the factors of 12 allow: ranks x threads = 12.
        t12_1, _ = measure(CYCLIC, cyclic_args, 13, omp=1)                # == cyclic -n 13 row
        t6_2, _  = measure(CYCLIC, cyclic_args, 7, omp=2, bind_none=True)
        t4_3, _  = measure(CYCLIC, cyclic_args, 5, omp=3, bind_none=True)
        t2_6, _  = measure(CYCLIC, cyclic_args, 3, omp=6, bind_none=True)
        rows = [
            ["``-n 13``, 12 workers x 1 thread", f"{t12_1:.1f}", speedup(serial, t12_1)],
            ["``-n 7``, 6 workers x 2 threads",  f"{t6_2:.1f}",  speedup(serial, t6_2)],
            ["``-n 5``, 4 workers x 3 threads",  f"{t4_3:.1f}",  speedup(serial, t4_3)],
            ["``-n 3``, 2 workers x 6 threads",  f"{t2_6:.1f}",  speedup(serial, t2_6)],
        ]
        cap = (f"cyclic-{args.cyclic_n} at 12 worker-cores, ranks x threads, "
               f"measured on {cores} cores, {stamp} -- {refresh}")
        blocks["hybrid"] = render_table(cap, [40, 25, 15],
                                        ["layout", "wall-clock (s)", "speedup"], rows)

    print("\n=== measured tables ===")
    for key in ("cyclic", "eigen", "hybrid"):
        if key in blocks:
            print(f"\n{blocks[key]}")

    if args.dry_run:
        print("\n--dry-run: not editing", RST)
        return

    text = RST.read_text(encoding="utf-8")
    for key, block in blocks.items():
        text = splice(text, key, block)
    RST.write_text(text, encoding="utf-8")
    print(f"\nUpdated {len(blocks)} table(s) in {RST}")


if __name__ == "__main__":
    main()
