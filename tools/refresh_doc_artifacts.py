#!/usr/bin/env python
"""Regenerate the machine-generated artifacts embedded in the tutorials.

Two very different kinds of artifact live in ``python/docs/source/tutorials/`` and go stale:

* **timing tables** -- pure numbers (wall-clock, speedups) rendered from data files, refreshed
  by :mod:`tools.update_scaling_timings`.  They are version-independent, so they are safe to
  refresh on every release and cheap to gate in CI.
* **plots** -- matplotlib figures (``.svg`` + ``.png``).  These are *doubly* fragile: a different
  matplotlib **version** reshapes the whole SVG, and even same-version runs churn on hashed
  element ids, an embedded ``<dc:date>``, and scatter point-order (data-order nondeterminism).
  So plots are regenerated **deliberately**, on a pinned matplotlib, never on every release.

Because of that split this tool defaults to **timings only**.  You must opt in to touch plots::

    python tools/refresh_doc_artifacts.py                 # timings only (the default)
    python tools/refresh_doc_artifacts.py --timings       # same, explicit
    python tools/refresh_doc_artifacts.py --plots         # ONLY regenerate the plot images
    python tools/refresh_doc_artifacts.py --all           # timings + plots
    python tools/refresh_doc_artifacts.py --plots --only real_points   # one plot
    python tools/refresh_doc_artifacts.py --list          # show what would run, do nothing

Run it **from the repository root, inside the bertini environment** (the b2-B venv on the /B
worktree).  Numbers/plots should come from the designated **reference machine** -- CI never
generates them (see the release runbook / CI staleness gate).

Plot determinism: this tool removes the churn it *can* -- it runs each plot script with a
temporary ``MATPLOTLIBRC`` pinning ``svg.hashsalt`` and then strips the ``<dc:date>`` line from
every SVG it produced.  The residual point-order churn is inherent to the plotting data and is
the reason plots stay opt-in.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TUT = REPO / "python" / "docs" / "source" / "tutorials"
SHOWPIECES = REPO / "python" / "docs" / "source" / "showpieces"
EXAMPLES = REPO / "python" / "examples"

# Where each plot's solves record (BERTINI_RECORDS_DIR), one folder per plot under the docs'
# gitignored build output, so records never land in the docs source.  Each folder is emptied
# before its plot runs: persistent records would RECALL identical solves instead of computing
# them, which draws a figure from an older library and leaves path observers with nothing to
# see.  Kept after the run, for looking at what a figure's solves did.
RECORDS_SCRATCH = REPO / "python" / "docs" / "build" / "refresh_records"

# A fixed salt makes matplotlib's SVG element ids deterministic run-to-run (same matplotlib
# version).  Any stable string works; keep it constant so ids don't move.
SVG_HASHSALT = "bertini2-docs"


# --- plot manifest -----------------------------------------------------------------------------
# Each entry: the regenerator script, the argv it needs, the image files it is expected to write
# (relative to `outdir`, which is the script's own directory unless noted), and an optional
# `needs` predicate for artifacts that require something extra (e.g. the built CLI binary).
# The script runs with `outdir` as its working directory, so a plain `savefig("name.svg")`
# lands where the docs expect it; the script needs no path handling of its own.
class Plot:
    def __init__(self, key, script, outputs, argv=None, outdir=None, needs=None, note=None):
        self.key = key
        self.script = script                      # Path
        self.outputs = outputs                    # list[str] basenames
        self.outdir = outdir or script.parent     # where the images land
        self.argv = argv or []                    # extra argv
        self.needs = needs                        # optional (label, Path) that must exist
        self.note = note


def _plots():
    B2 = REPO / "build" / "core" / "bertini2"   # optional CLI binary for the b1-vs-b2 benchmark
    return [
        Plot("real_points",
             TUT / "formulating_and_solving" / "real_points" / "real_points.py",
             ["real_points.png"]),
        Plot("parameter_homotopy",
             TUT / "formulating_and_solving" / "parameter_homotopy" / "parameter_homotopy.py",
             ["parameter_homotopy_circle.svg", "parameter_homotopy_circle.png"]),
        Plot("solution_dataframe",
             TUT / "formulating_and_solving" / "solution_dataframe" / "solution_dataframe.py",
             ["solution_dataframe_real_plane.svg", "solution_dataframe_real_plane.png",
              "solution_dataframe_complex_planes.svg", "solution_dataframe_complex_planes.png"]),
        Plot("plotting",
             TUT / "manipulating_solutions" / "plotting" / "plotting.py",
             ["plotting_solutions.svg", "plotting_solutions.png",
              "plotting_paths.svg", "plotting_paths.png"]),
        Plot("observers_and_path_data",
             TUT / "observing_metadata_more" / "observers_and_path_data" / "observers_and_path_data.py",
             ["observers_and_path_data.svg", "observers_and_path_data.png",
              "cyclic3_paths.svg", "cyclic3_paths.png",
              "griewank_osborn_endgame.svg", "griewank_osborn_endgame.png"]),
        Plot("classic_continuation_cartoon",
             TUT / "observing_metadata_more" / "classic_continuation_cartoon" / "classic_continuation_cartoon.py",
             ["classic_continuation_cartoon.svg", "classic_continuation_cartoon.png"]),
        Plot("seeing_precision_change",
             TUT / "observing_metadata_more" / "seeing_precision_change" / "amp_precision_plot.py",
             ["amp_precision_cyclic5.svg", "amp_precision_cyclic5.png"]),
        Plot("parallel_parameter_homotopy",
             EXAMPLES / "parallel_parameter_homotopy.py",
             ["parallel_parameter_homotopy.svg"],
             outdir=TUT / "parallelism" / "parallel_parameter_homotopy",
             argv=["--save", "parallel_parameter_homotopy.svg"]),
        Plot("bertini1_vs_bertini2_timing",
             TUT / "performance_benchmarking" / "bertini1_vs_bertini2_timing" / "b1_vs_b2_timing.py",
             ["b1_vs_b2_timing.svg", "b1_vs_b2_timing.png"],
             argv=["--bertini2", str(B2)],
             needs=("bertini2 CLI binary", B2),
             note="benchmark vs Bertini 1; needs the built CLI and (optionally) a `bertini` on PATH"),
        Plot("tracking_analytic",
             TUT / "doing_things_manually" / "tracking_analytic" / "tracking_analytic.py",
             ["tracking_analytic.svg", "tracking_analytic.png"]),
        Plot("critical_points",
             TUT / "formulating_and_solving" / "critical_points" / "critical_points_plot.py",
             ["critical_points.svg", "critical_points.png"]),
        Plot("chained_homotopies",
             EXAMPLES / "chained_homotopies.py",
             ["chain_progression.svg", "chain_progression.png"],
             outdir=TUT / "record_keeping" / "chained_homotopies",
             argv=["--plot", "chain_progression"]),
        Plot("monodromy_loom",
             SHOWPIECES / "monodromy_loom" / "monodromy_loom.py",
             ["monodromy_loom.png", "monodromy_loom_teaching.png"],
             note="showpiece: raster PNG only (a 3-D render has no meaningful SVG)"),
        Plot("flight_recorder",
             SHOWPIECES / "flight_recorder" / "flight_recorder.py",
             ["flight_recorder.png", "flight_recorder_setup.png"],
             note="showpiece: raster PNG only; tight-tolerance mult-35 solve -> slow (~1 min)"),
        Plot("homotopy_basins",
             SHOWPIECES / "homotopy_basins" / "homotopy_basins.py",
             ["homotopy_basins.png", "homotopy_basins_annotated.png",
              "homotopy_basins_teaching.png", "homotopy_basins_teaching_annotated.png"],
             note="showpiece: raster PNG only; ~17M tracked paths -> SLOW (~10-15 min on 12 cores)"),
    ]


# critical_points_fibers.svg/.png is a DELIBERATE hand-authored illustration: it has no in-repo
# source (no testcode, no script draws it) and depicts the projection-fiber cartoon by hand.  It is
# listed here so `--plots` reports it as intentionally uncovered rather than a silent gap.  (The
# other critical_points figure, critical_points.svg, now has a regenerator -- see _plots() above.)
_NO_REGENERATOR = {
    "critical_points": ["critical_points_fibers.svg"],
}


def _env_with_determinism():
    """A child-process env whose matplotlib pins svg.hashsalt (deterministic element ids)."""
    rc = tempfile.NamedTemporaryFile("w", suffix="matplotlibrc", delete=False)
    rc.write(f"svg.hashsalt: {SVG_HASHSALT}\n")
    rc.close()
    # serial b2 (BERTINI_NUM_THREADS) and serial numpy (OMP_NUM_THREADS, which OpenBLAS reads):
    # a figure must not depend on thread scheduling in either library
    env = dict(os.environ, MATPLOTLIBRC=rc.name, MPLBACKEND="Agg", BERTINI_NUM_THREADS="1",
               OMP_NUM_THREADS="1")
    return env, rc.name


_DCDATE = re.compile(rb"^\s*<dc:date>.*</dc:date>\n", re.MULTILINE)


def _strip_svg_date(path: Path):
    """Remove the per-run <dc:date> line so an unchanged plot stays byte-identical."""
    if path.suffix != ".svg" or not path.exists():
        return
    data = path.read_bytes()
    stripped = _DCDATE.sub(b"", data)
    if stripped != data:
        path.write_bytes(stripped)


def run_timings(extra_argv, dry_run):
    cmd = [sys.executable, str(REPO / "tools" / "update_scaling_timings.py"), *extra_argv]
    if dry_run:
        cmd.append("--dry-run")
    print(f"[timings] $ {' '.join(cmd)}", flush=True)
    return subprocess.run(cmd, cwd=REPO).returncode


def run_plot(plot: Plot, dry_run):
    if plot.needs:
        label, path = plot.needs
        if not Path(path).exists():
            print(f"[plots] SKIP {plot.key}: missing {label} ({path})", flush=True)
            return 0  # not a failure -- this artifact just can't be built here
    if not plot.outdir.is_dir():
        print(f"[plots] FAIL {plot.key}: its image directory does not exist ({plot.outdir})", flush=True)
        return 1
    cmd = [sys.executable, str(plot.script), *plot.argv]
    print(f"[plots] $ {' '.join(cmd)}", flush=True)
    if dry_run:
        return 0
    started = time.time()
    records = RECORDS_SCRATCH / plot.key
    shutil.rmtree(records, ignore_errors=True)
    env, rc_path = _env_with_determinism()
    env["BERTINI_RECORDS_DIR"] = str(records)
    try:
        proc = subprocess.run(cmd, cwd=plot.outdir, env=env)
    finally:
        os.unlink(rc_path)
    if proc.returncode != 0:
        return proc.returncode
    # a clean exit is not proof: the script may have saved somewhere else, or not at all
    stale = [name for name in plot.outputs
             if not (plot.outdir / name).exists() or (plot.outdir / name).stat().st_mtime < started]
    if stale:
        print(f"[plots] FAIL {plot.key}: not written to {plot.outdir}: {', '.join(stale)}", flush=True)
        return 1
    for name in plot.outputs:
        _strip_svg_date(plot.outdir / name)
    print(f"[plots] wrote: {', '.join(plot.outputs)}", flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--timings", action="store_true", help="refresh the timing tables (default)")
    ap.add_argument("--plots", action="store_true", help="regenerate the plot images (opt-in)")
    ap.add_argument("--all", action="store_true", help="timings + plots")
    ap.add_argument("--only", metavar="KEY", help="restrict --plots to one plot by key")
    ap.add_argument("--list", action="store_true", help="show what would run, then exit")
    ap.add_argument("--dry-run", action="store_true", help="print actions but do not write")
    ap.add_argument("rest", nargs=argparse.REMAINDER,
                    help="args after -- are forwarded to update_scaling_timings.py")
    args = ap.parse_args()

    do_timings = args.timings or args.all or not (args.plots or args.all)  # default: timings
    do_plots = args.plots or args.all

    plots = _plots()
    if args.only:
        plots = [p for p in plots if p.key == args.only]
        if not plots:
            sys.exit(f"no plot with key {args.only!r}; known: {', '.join(p.key for p in _plots())}")

    if args.list:
        print("timings:", "yes" if do_timings else "no")
        print("plots:", "yes" if do_plots else "no")
        for p in plots if do_plots else []:
            gated = f"  (needs {p.needs[0]})" if p.needs else ""
            print(f"  - {p.key}: {p.script.relative_to(REPO)}{gated}")
        if do_plots:
            for key, imgs in _NO_REGENERATOR.items():
                print(f"  ! {key}: NO regenerator for {', '.join(imgs)} (inline-testcode plot)")
        return

    forwarded = [a for a in args.rest if a != "--"]
    rc = 0
    if do_timings:
        rc |= run_timings(forwarded, args.dry_run)
    if do_plots:
        for p in plots:
            rc |= run_plot(p, args.dry_run)
        for key, imgs in _NO_REGENERATOR.items():
            print(f"[plots] NOTE: {key} has no regenerator for {', '.join(imgs)} "
                  f"(deliberate hand-authored illustration -- edit the image by hand if it changes)",
                  flush=True)
    sys.exit(rc)


if __name__ == "__main__":
    main()
