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

# python/test/docs/this_file.py -> parents[2] == the `python/` dir, which holds `examples/`.
PYTHON_DIR = Path(__file__).resolve().parents[2]
EXAMPLES = PYTHON_DIR / "examples"


def _run(script, *args, timeout=170):
    """Run an example script serially, importing *this* worktree's bertini, and require success."""
    env = dict(os.environ)
    # Prepend this worktree's python dir so the subprocess imports the bertini we are testing,
    # not whatever a shared environment's site-packages might point at.
    env["PYTHONPATH"] = os.pathsep.join([str(PYTHON_DIR), env.get("PYTHONPATH", "")])
    proc = subprocess.run(
        [sys.executable, str(EXAMPLES / script), *args],
        capture_output=True, text=True, timeout=timeout, env=env,
    )
    assert proc.returncode == 0, (
        f"{script} {' '.join(args)} exited {proc.returncode}\n"
        f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
    )
    return proc


def test_solve_cyclic_runs():
    # cyclic-5 (120 paths, 70 finite) is the smallest *zero-dimensional* cyclic case;
    # cyclic-4 is positive-dimensional, so do not use it here.
    _run("solve_cyclic.py", "--n", "5")


def test_solve_eigenvalues_runs():
    _run("solve_eigenvalues.py", "--size", "5")


def test_crossed_paths_runs():
    # No size knob: it always provokes a crossing on cyclic-5 and shows the re-track repair.
    _run("crossed_paths.py")
