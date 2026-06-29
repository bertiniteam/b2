"""Benchmark Bertini 2 against Bertini 1 on the same systems, and plot the wall-time comparison.

The fair way to compare two solvers is to hand them the *same problem*: we build each system once in
pybertini and emit it as a Bertini-1 classic input with ``System.to_classic_input(...)``, so both
solvers get identical equations and identical tracking settings (precision mode, predictor,
tolerances, step cadence).  Each solver then builds its own random start system and solves; we time
the solver subprocess only (never parsing/IO) and take the fastest of a few repeats.

This is the same idea as ``benchmark/comparison/run_comparison.py`` (the maintained harness, which also
sweeps MPI ranks and appends a committed history.csv); here it is a small self-contained script whose
crescendo is a plot.

Run (needs matplotlib; Bertini 1 optional)::

    python python/examples/b1_vs_b2_timing.py \
        --bertini2 ./build/core/bertini2 --bertini1 /usr/local/bin/bertini --out .
"""
import argparse
import datetime
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
import bertini as pb


def _first_nonempty_version_line(exe):
    if not exe:
        return "n/a"
    try:
        out = subprocess.run([os.path.abspath(exe), "--version"], capture_output=True, text=True, timeout=15)
        for line in (out.stdout + out.stderr).splitlines():
            if line.strip():
                return line.strip()[:60]
    except Exception:
        pass
    return "unknown"


def provenance(b2_exe, b1_exe):
    """Record WHAT was timed and WHERE, so the numbers can be refreshed and stay interpretable."""
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        commit = subprocess.run(["git", "-C", here, "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True).stdout.strip()
    except Exception:
        commit = ""
    cpu = platform.processor() or "unknown"
    try:
        if platform.system() == "Darwin":
            cpu = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                                 capture_output=True, text=True).stdout.strip() or cpu
    except Exception:
        pass
    b2 = _first_nonempty_version_line(b2_exe) + (" @ " + commit if commit else "")
    return {"date": datetime.date.today().isoformat(),
            "b2": b2,
            "b1": _first_nonempty_version_line(b1_exe),
            "machine": "{} / {} {}".format(cpu, platform.system(), platform.release())}


def diagonal(n):
    """x_i^3 - c_i = 0 : 3^n well-separated, well-conditioned roots (an easy warm-up family)."""
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23]
    x = pb.variables('x', n)
    s = pb.System()
    for i in range(n):
        s.add_function(x[i] ** 3 - primes[i])
    s.add_variable_group(pb.VariableGroup(x))
    return s


def cyclic(n):
    """The cyclic-n roots system (n=5: 70 finite solutions; the endgame-stressing case)."""
    x = pb.variables('x', n)
    w = list(x) + list(x)
    s = pb.System()
    for length in range(1, n):
        s.add_function(np.sum([np.prod(w[start:start + length]) for start in range(n)]))
    s.add_function(np.prod(x) - 1)
    s.add_variable_group(pb.VariableGroup(x))
    return s


SYSTEMS = [("diag-3", diagonal(3)), ("diag-5", diagonal(5)),
           ("diag-6", diagonal(6)), ("cyclic-5", cyclic(5))]


def solution_count(d):
    for name in ("finite_solutions", "nonsingular_solutions", "raw_solutions"):
        p = os.path.join(d, name)
        if os.path.exists(p):
            try:
                return int(open(p).readline().strip())
            except ValueError:
                pass
    return None


def time_solver(exe, input_text, repeats=3, timeout=600):
    """Fastest serial wall time (seconds) over `repeats` runs of `exe` on `input_text`; (time, count)."""
    best, count = float("inf"), None
    for _ in range(repeats):
        d = tempfile.mkdtemp()
        try:
            open(os.path.join(d, "input"), "w").write(input_text)
            env = dict(os.environ, OMP_NUM_THREADS="1")
            t0 = time.perf_counter()
            r = subprocess.run([os.path.abspath(exe), "input"], cwd=d, env=env,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=timeout)
            dt = time.perf_counter() - t0
            if r.returncode == 0:
                best = min(best, dt)
                count = solution_count(d)
        finally:
            shutil.rmtree(d, ignore_errors=True)
    return (best if best < float("inf") else None), count


def make_plot(rows, prov, out_stem):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = [r["name"] for r in rows]
    t2 = [r["t2"] for r in rows]
    t1 = [r["t1"] for r in rows]
    xs = np.arange(len(names))
    w = 0.38

    fig, ax = plt.subplots(figsize=(10, 6))
    b2bars = ax.bar(xs - w / 2, t2, w, label="Bertini 2", color="#1565c0")
    have_b1 = any(v is not None for v in t1)
    if have_b1:
        b1vals = [v if v is not None else 0.0 for v in t1]
        ax.bar(xs + w / 2, b1vals, w, label="Bertini 1", color="#ef6c00")

    ax.set_yscale("log")
    ax.set_ylabel("serial wall time (s), fastest of 3   —   lower is better")
    ax.set_xticks(xs); ax.set_xticklabels(names)
    ax.set_title("Bertini 1 vs Bertini 2 — adaptive-precision zero-dim solve\n"
                 "same system & settings to each solver (its own random start system)")
    ax.legend()
    ax.grid(True, axis="y", which="both", alpha=0.2)

    # annotate the b2/b1 slowdown factor over each pair
    if have_b1:
        for xi, r in zip(xs, rows):
            if r["t1"]:
                ax.annotate("{:.0f}x".format(r["t2"] / r["t1"]),
                            xy=(xi, max(r["t2"], r["t1"])), ha="center", va="bottom", fontsize=9)
    # solution counts under the x labels
    ax.set_xlabel(" ".join("{}={}".format(r["name"], r["n2"]) for r in rows if r["n2"] is not None))

    # provenance — so the plot stays interpretable and can be refreshed as the library improves
    caption = ("{date}   |   {b2}   |   {b1}   |   {machine}").format(**prov)
    fig.text(0.5, 0.005, caption, ha="center", va="bottom", fontsize=7.5, color="0.35")

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    svg, png = out_stem + ".svg", out_stem + ".png"
    fig.savefig(svg); fig.savefig(png, dpi=150)
    plt.close(fig)
    return svg, png


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bertini2", default="./build/core/bertini2")
    ap.add_argument("--bertini1", default=shutil.which("bertini"))
    ap.add_argument("--out", default=".")
    ap.add_argument("--mptype", type=int, default=2)
    args = ap.parse_args()

    prov = provenance(args.bertini2, args.bertini1)
    print("  ".join("{}={}".format(k, v) for k, v in prov.items()))

    rows = []
    for name, system in SYSTEMS:
        txt = system.to_classic_input(mptype=args.mptype)
        t2, n2 = time_solver(args.bertini2, txt)
        t1, n1 = time_solver(args.bertini1, txt) if args.bertini1 else (None, None)
        rows.append(dict(name=name, t2=t2, n2=n2, t1=t1, n1=n1))
        print("{:9s} b2={:7.3f}s ({}) b1={} ({})".format(
            name, t2 or float("nan"), n2,
            "{:.3f}s".format(t1) if t1 else "n/a", n1))

    os.makedirs(args.out, exist_ok=True)
    svg, png = make_plot(rows, prov, os.path.join(args.out, "b1_vs_b2_timing"))
    print("wrote", svg, "\n      ", png)


if __name__ == "__main__":
    main()
