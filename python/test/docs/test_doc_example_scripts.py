"""Run the example scripts that the tutorials pull in via ``.. literalinclude::``.

Two tutorials -- :doc:`crossed_paths` and :doc:`solving_at_scale` -- do not have inline,
doctestable code: their substantive code lives in standalone scripts under ``python/examples/``
(shown in the rendered docs with ``literalinclude``), and the inline blocks are MPI/shell
fragments that cannot run under ``sphinx -b doctest``.  So we verify those scripts here instead,
by running them **serially** (a lone MPI rank takes the no-communicator path) at small sizes.

This keeps the literalinclude'd code honest -- it is exactly how the bit-rot in the other
tutorials was caught, just at the script level rather than the doctest level.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

# Every one of these example scripts does `from mpi4py import MPI` at import time (they are
# MPI-aware, run serially when launched without mpirun).  The CI test environments do not install
# mpi4py, so skip the whole module there rather than fail.
pytest.importorskip("mpi4py")

# python/test/docs/this_file.py -> parents[2] == the `python/` dir, which holds `examples/`.
PYTHON_DIR = Path(__file__).resolve().parents[2]
EXAMPLES = PYTHON_DIR / "examples"


def _run(script, *args, timeout=170):
    """Run an example script serially, importing the same bertini this test process does."""
    env = dict(os.environ)
    # Resolve `bertini` in the subprocess exactly as this test process resolves it -- whether that
    # is an installed wheel (the CI wheel-test envs, where the compiled bertini._pybertini lives
    # only in site-packages, NOT in the source python/ tree) or a source checkout on PYTHONPATH
    # (local dev).  Propagate our own import path; do NOT hard-code the source python/ dir, which
    # would shadow an installed wheel with a _pybertini-less source package.
    env["PYTHONPATH"] = os.pathsep.join(p for p in sys.path if p)
    proc = subprocess.run(
        [sys.executable, str(EXAMPLES / script), *args],
        capture_output=True, text=True, timeout=timeout, env=env,
    )
    assert proc.returncode == 0, (
        f"{script} {' '.join(args)} exited {proc.returncode}\n"
        f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
    )
    return proc


# cyclic-5 in ADAPTIVE precision (the point of the example) is much slower than double and varies a
# lot by runner -- on a slow Windows CI box it exceeds the global 180s pytest-timeout, so override it
# here.  The inner subprocess budget (540s) sits just under the test budget (600s) so a genuine hang
# surfaces as a clean subprocess.TimeoutExpired rather than pytest-timeout killing the test.
@pytest.mark.timeout(600)
def test_solve_cyclic_runs():
    # cyclic-5 (120 paths, 70 finite) is the smallest *zero-dimensional* cyclic case;
    # cyclic-4 is positive-dimensional, so do not use it here.
    _run("solve_cyclic.py", "--n", "5", timeout=540)


def test_solve_eigenvalues_runs():
    _run("solve_eigenvalues.py", "--size", "5")


# Two full cyclic-5 solves with a re-track repair.  This used to be gated behind
# BERTINI_RUN_SLOW_DOC_EXAMPLES (~30s locally, minutes on slow CI), but the adaptive-numeric-type
# endgame + well-conditioned start coefficients made it run in a few seconds, so it now runs by
# default.  The timeout is a slow-CI ceiling, matching test_solve_cyclic_runs.
@pytest.mark.timeout(600)
def test_crossed_paths_runs():
    # No size knob: it always provokes a crossing on cyclic-5 and shows the re-track repair.
    _run("crossed_paths.py", timeout=540)
