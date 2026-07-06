"""Solver adapters for the external-solver comparison benchmark.

Each adapter runs one CLI solver on a given classic-Bertini input *string*, in its own fresh
scratch directory (both Bertini 1 and the bertini2 CLI pollute their CWD with main_data,
raw_data, output, failed_paths, ... -- isolation keeps runs from colliding and keeps the repo
clean), times only the subprocess, and -- separately, after timing -- parses the solution count
from the files the run wrote.

Adding HomotopyContinuation.jl later is just a new adapter with the same ``run(...) -> RunResult``
shape; nothing else in the suite needs to change.

IMPORTANT: timing measures ONLY the solve subprocess (perf_counter brackets subprocess.run).
Parsing the solution count and checking agreement happen afterwards, from already-written output
files, and are never included in any reported wall time.
"""

import os
import re
import shlex
import shutil
import subprocess
import tempfile
import time
from collections import namedtuple

# wall_time_s: float (nan on failure); solutions_found: int (-1 if unparseable);
# ok: did the process exit 0 and produce a parseable count; detail: short human note.
RunResult = namedtuple("RunResult", ["wall_time_s", "solutions_found", "ok", "detail"])


def _launch(exe, ranks, mpirun, mpirun_args):
    """Build the argv to launch a solver: direct for serial, under mpirun for an MPI rank sweep."""
    exe = os.path.abspath(exe)
    if ranks <= 1:
        # Serial: invoke the binary directly. An MPI-built binary run without mpirun is a single
        # MPI process, which is exactly the serial baseline we want.
        return [exe]
    return [mpirun, *shlex.split(mpirun_args), "-n", str(ranks), exe]


def _time_run(argv, input_text, threads, timeout):
    """Run argv in a fresh scratch dir with the given input; time ONLY the subprocess.

    Returns (wall_time_s, returncode, scratch_dir, stderr_tail).  Caller parses output files in
    scratch_dir, then removes it.  wall_time_s is nan on timeout.
    """
    scratch = tempfile.mkdtemp(prefix="b2_cmp_")
    with open(os.path.join(scratch, "input"), "w") as f:   # both solvers read a file named "input"
        f.write(input_text)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)                  # single-thread for this benchmark phase

    t0 = time.perf_counter()
    try:
        # stdin=DEVNULL is essential: the solvers read no input, but if stdin is an inherited
        # pipe/tty they can block forever waiting on it (observed as a hang under capture).
        proc = subprocess.run(argv, cwd=scratch, env=env, stdin=subprocess.DEVNULL,
                              capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return float("nan"), None, scratch, "timeout"
    t1 = time.perf_counter()
    stderr_tail = proc.stderr.decode(errors="replace")[-500:]
    return t1 - t0, proc.returncode, scratch, stderr_tail


def _count_from_first_int_line(path):
    """A classic Bertini solution file begins with an integer count on its first line."""
    try:
        with open(path) as f:
            return int(f.readline().strip())
    except (OSError, ValueError):
        return -1


# --------------------------------------------------------------------------------------------
# Adapters
# --------------------------------------------------------------------------------------------

def _run_classic(exe, input_text, ranks, threads, mpirun, mpirun_args, timeout):
    """Run a solver that writes Bertini 1.7 classic output files; return a RunResult.

    Both Bertini 1 and (as of the classic-output work) the Bertini 2 CLI write the count-led
    `finite_solutions` file, so a single parser serves both.
    """
    argv = _launch(exe, ranks, mpirun, mpirun_args)
    wall, rc, scratch, err = _time_run(argv, input_text, threads, timeout)
    try:
        if wall != wall:                                   # nan -> timeout
            return RunResult(float("nan"), -1, False, "timeout")
        if rc != 0:
            return RunResult(float("nan"), -1, False, f"exit {rc}: {err}")
        count = _finite_solution_count(scratch)
        return RunResult(wall, count, count >= 0, "ok" if count >= 0 else "unparseable output")
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def _finite_solution_count(scratch):
    """Solution count from the classic count-led solution files (first line is the count).

    `finite_solutions` is the apples-to-apples cross-solver count; fall back to the other
    count-led files, then a regex scrape of main_data as a last resort.
    """
    for name in ("finite_solutions", "nonsingular_solutions", "real_finite_solutions",
                 "raw_solutions"):
        c = _count_from_first_int_line(os.path.join(scratch, name))
        if c >= 0:
            return c
    try:
        with open(os.path.join(scratch, "main_data")) as f:
            text = f.read()
        m = re.search(r"(\d+)\s+(?:finite\s+)?solutions", text, re.IGNORECASE)
        if m:
            return int(m.group(1))
    except OSError:
        pass
    return -1


def run_bertini2(exe, input_text, *, ranks=1, threads=1,
                 mpirun="mpirun", mpirun_args="--bind-to none", timeout=600.0):
    """Run the bertini2 CLI (writes Bertini 1.7-compatible solution files)."""
    return _run_classic(exe, input_text, ranks, threads, mpirun, mpirun_args, timeout)


def run_bertini1(exe, input_text, *, ranks=1, threads=1,
                 mpirun="mpirun", mpirun_args="--bind-to none", timeout=600.0):
    """Run Bertini 1.7 (assumed MPI-built)."""
    return _run_classic(exe, input_text, ranks, threads, mpirun, mpirun_args, timeout)


# Registry of available solver adapters, keyed by the column name used in output.
ADAPTERS = {
    "bertini2": run_bertini2,
    "bertini1": run_bertini1,
}
