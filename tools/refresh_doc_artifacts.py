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
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TUT = REPO / "python" / "docs" / "source" / "tutorials"
EXAMPLES = REPO / "python" / "examples"

# A fixed salt makes matplotlib's SVG element ids deterministic run-to-run (same matplotlib
# version).  Any stable string works; keep it constant so ids don't move.
SVG_HASHSALT = "bertini2-docs"


# --- plot manifest -----------------------------------------------------------------------------
# Each entry: the regenerator script, the argv it needs, the image files it is expected to write
# (relative to `outdir`, which is the script's own directory unless noted), and an optional
# `needs` predicate for artifacts that require something extra (e.g. the built CLI binary).
class Plot:
    def __init__(self, key, script, outputs, argv=None, outdir=None, needs=None, note=None):
        self.key = key
        self.script = script                      # Path
        self.outputs = outputs                    # list[str] basenames
        self.outdir = outdir or script.parent     # where the images land
        self.argv = argv or []                    # extra argv (may reference {outdir})
        self.needs = needs                        # optional (label, Path) that must exist
        self.note = note

    def resolved_argv(self):
        return [a.format(outdir=str(self.outdir)) for a in self.argv]


def _plots():
    B1 = REPO / "build" / "core" / "bertini2"     # optional CLI binary for the b1-vs-b2 benchmark
    return [
        Plot("real_points",
             TUT / "real_points" / "real_points.py",
             ["real_points.png"]),
        Plot("solution_database",
             TUT / "solution_database" / "solution_database.py",
             ["solution_database_real_plane.svg", "solution_database_real_plane.png",
              "solution_database_complex_planes.svg", "solution_database_complex_planes.png"]),
        Plot("observers_and_path_data",
             TUT / "observers_and_path_data" / "observers_and_path_data.py",
             ["observers_and_path_data.svg", "observers_and_path_data.png",
              "cyclic3_paths.svg", "cyclic3_paths.png",
              "griewank_osborn_endgame.svg", "griewank_osborn_endgame.png"]),
        Plot("classic_continuation_cartoon",
             TUT / "classic_continuation_cartoon" / "classic_continuation_cartoon.py",
             ["classic_continuation_cartoon.svg", "classic_continuation_cartoon.png"],
             argv=["{outdir}"]),
        Plot("homotopy_cartoon_from_real_data",
             TUT / "homotopy_cartoon_from_real_data" / "amp_precision_cartoon.py",
             ["amp_precision_cartoon_cyclic5.svg", "amp_precision_cartoon_cyclic5.png"],
             argv=["{outdir}"]),
        Plot("parallel_parameter_homotopy",
             EXAMPLES / "parallel_parameter_homotopy.py",
             ["parallel_parameter_homotopy.svg"],
             outdir=TUT / "parallel_parameter_homotopy",
             argv=["--save", "{outdir}/parallel_parameter_homotopy.svg"]),
        Plot("bertini1_vs_bertini2_timing",
             TUT / "bertini1_vs_bertini2_timing" / "b1_vs_b2_timing.py",
             ["b1_vs_b2_timing.svg", "b1_vs_b2_timing.png"],
             argv=["--out", "{outdir}", "--bertini2", str(B1)],
             needs=("bertini2 CLI binary", B1),
             note="benchmark vs Bertini 1; needs the built CLI and (optionally) a `bertini` on PATH"),
        Plot("critical_points",
             TUT / "critical_points" / "critical_points_plot.py",
             ["critical_points.svg", "critical_points.png"]),
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
    env = dict(os.environ, MATPLOTLIBRC=rc.name, MPLBACKEND="Agg", OMP_NUM_THREADS="1")
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
    cmd = [sys.executable, str(plot.script), *plot.resolved_argv()]
    print(f"[plots] $ {' '.join(cmd)}", flush=True)
    if dry_run:
        return 0
    env, rc_path = _env_with_determinism()
    try:
        proc = subprocess.run(cmd, cwd=REPO, env=env)
    finally:
        os.unlink(rc_path)
    if proc.returncode == 0:
        for name in plot.outputs:
            _strip_svg_date(plot.outdir / name)
        print(f"[plots] wrote: {', '.join(plot.outputs)}", flush=True)
    return proc.returncode


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
