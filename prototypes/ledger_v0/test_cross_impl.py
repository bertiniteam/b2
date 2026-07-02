# This file is part of Bertini 2 (prototypes/ledger_v0 -- experimental, unshipped).
# GPL v3+; see the repository's licenses/ directory.

"""The cross-implementation bridge (arc rung 3): the Python prototype and the C++
OutputDirectory must read each other's directories.

- reads `build/core/cross_impl_cpp_output` (written by the C++ suite's
  output_directory/write_the_cross_impl_example case; skip if ctest hasn't run);
- writes `build/core/cross_impl_python_output` for the C++ suite's
  read_the_python_written_example_if_present case (rerun ctest after this to close
  the loop).

Run: PYTHONPATH=python pytest prototypes/ledger_v0/test_cross_impl.py
"""

import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

from ledger import Ledger

REPO = Path(__file__).resolve().parents[2]
CPP_DIR = REPO / "build" / "core" / "cross_impl_cpp_output"
PY_DIR = REPO / "build" / "core" / "cross_impl_python_output"


@pytest.mark.skipif(not CPP_DIR.exists(), reason="C++ example absent; run ctest first")
def test_prototype_reads_cpp_written_directory():
    lg = Ledger(CPP_DIR)
    records = lg.scan()
    kinds = [r["kind"] for r in records]
    assert "run" in kinds and "track" in kinds

    run = next(r for r in records if r["kind"] == "run")
    assert run["schema"] == "ledgerrec/1"
    done = lg.completed_paths(run["run"])
    assert len(done) == run["num_paths"] == 2
    # the endpoint decodes: sqrt(2) as recorded by C++
    assert abs(float(done[0]["endpoint"][0][0]) - 2 ** 0.5) < 1e-15

    # the definition is content-addressed identically in both implementations
    target = lg.get_object(run["target_object"])
    assert hashlib.sha256(target).hexdigest() == run["target_object"]
    assert b"x^2-2" in target

    # the C++-rendered results.json is one json.load away, as promised
    machine = json.loads((CPP_DIR / "results.json").read_text())
    assert isinstance(machine, dict)


def test_prototype_writes_directory_for_cpp_to_read():
    """Deterministic-enough example for the C++ read leg (rerun ctest to close)."""
    import shutil
    if PY_DIR.exists():
        shutil.rmtree(PY_DIR)
    PY_DIR.parent.mkdir(parents=True, exist_ok=True)
    lg = Ledger(PY_DIR)
    target_id = lg.put_object("function f;\nvariable_group x;\nf = x^3-2;\n")
    lg.append({"kind": "run", "schema": "ledgerrec/1", "when": "2026-07-03 10:00",
               "run": "pydemo00000000001", "op": "solve", "target_object": target_id,
               "num_paths": 1,
               "start_points": [[["1", "0"]]]})
    lg.append({"kind": "track", "run": "pydemo00000000001", "index": 0,
               "status": "success", "endpoint": [["1.2599210498948732", "0"]],
               "start": {"kind": "start_label", "index": 0}})
    assert len(lg.scan()) == 2
