# This file is part of Bertini 2 (prototypes/ledger_v0 -- experimental, unshipped).
# GPL v3+; see the repository's licenses/ directory.

"""Pilot tests for ledger v0: object store semantics, journal crash tolerance, and the
kill-and-rerun memoized solve -- the arc's acceptance test in miniature.

Run: PYTHONPATH=python pytest prototypes/ledger_v0/test_ledger.py
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))

import bertini as pb
from bertini import System, VariableGroup
from bertini.function_tree.symbol import Variable

from ledger import Ledger
from memo_solve import (ensure_solved, ensure_continued, provenance_chain,
                        annotate, annotations_for, SimulatedCrash)


def circle_line():
    x, y = Variable("x"), Variable("y")
    s = System()
    s.add_variable_group(VariableGroup([x, y]))
    s.add_function(x**2 + y**2 - 1)
    s.add_function(x - y)
    return s


# ---- object store ---------------------------------------------------------------------

def test_object_store_is_idempotent_and_self_verifying(tmp_path):
    lg = Ledger(tmp_path)
    oid1 = lg.put_object("hello ledger")
    oid2 = lg.put_object("hello ledger")
    assert oid1 == oid2                       # content-addressed: no conflict, no dup
    assert lg.get_object(oid1) == b"hello ledger"
    import hashlib
    assert oid1 == hashlib.sha256(b"hello ledger").hexdigest()  # self-verifying


def test_object_store_honors_external_ids(tmp_path):
    lg = Ledger(tmp_path)
    sys_ = circle_line()
    oid = lg.put_object(sys_.to_classic_input(), object_id=sys_.content_digest())
    assert oid == sys_.content_digest()
    assert b"function" in lg.get_object(oid)  # the human-readable rendering


# ---- journals -------------------------------------------------------------------------

def test_journal_appends_and_scan_round_trips(tmp_path):
    lg = Ledger(tmp_path)
    lg.append({"kind": "run", "run": "runA", "n": 1})
    lg.append({"kind": "track", "run": "runA", "index": 0})
    recs = lg.scan()
    assert [r["kind"] for r in recs] == ["run", "track"]
    # one session = one journal file, with a human-scannable date-stamped name
    journals = list((tmp_path / "journals").glob("*.jsonl"))
    assert len(journals) == 1
    assert journals[0].name[:8].isdigit()     # YYYYMMDD prefix: `ls` reads as history


def test_torn_final_line_is_tolerated(tmp_path):
    lg = Ledger(tmp_path)
    lg.append({"kind": "run", "run": "runB"})
    # simulate a kill mid-append: a truncated JSON line at EOF
    journal_file = next((tmp_path / "journals").glob("*.jsonl"))
    with open(journal_file, "a") as f:
        f.write('{"kind": "track", "ind')
    recs = lg.scan()
    assert len(recs) == 1                     # the torn tail is invisible to replay


def test_torn_interior_line_is_an_error(tmp_path):
    lg = Ledger(tmp_path)
    journal_file = tmp_path / "journals" / "bad-0-x.jsonl"
    journal_file.write_text('{"kind": "run"\n{"kind": "track", "index": 0}\n')
    with pytest.raises(ValueError, match="corrupt"):
        lg.scan()


# ---- the memoized solve + kill-and-rerun pilot ------------------------------------------

def test_fresh_solve_records_everything(tmp_path):
    lg = Ledger(tmp_path)
    result = ensure_solved(circle_line(), lg)
    assert result.num_computed == 2           # total degree: circle(2) x line(1) = 2 paths
    assert result.num_reused == 0
    assert len(result.statuses) == result.num_computed


def test_rerun_is_a_noop(tmp_path):
    lg = Ledger(tmp_path)
    first = ensure_solved(circle_line(), lg)
    again = ensure_solved(circle_line(), lg)  # independently built, equal content
    assert again.num_computed == 0            # memo hit end to end
    assert again.num_reused == len(first.statuses)
    # and the answers are the recorded ones
    for i, sol in again.solutions.items():
        assert np.allclose(
            np.array([complex(z) for z in sol]),
            np.array([complex(z) for z in first.solutions[i]]), atol=1e-25)


def test_kill_and_rerun_resumes(tmp_path):
    lg = Ledger(tmp_path)
    with pytest.raises(SimulatedCrash):
        ensure_solved(circle_line(), lg, crash_after=1)

    resumed = ensure_solved(circle_line(), lg)
    assert resumed.num_reused == 1            # the pre-crash path survived on disk
    assert resumed.num_computed == len(resumed.statuses) - 1

    # the resumed run agrees with a from-scratch solve on the correct solutions x=y=+-1/sqrt(2)
    values = sorted(complex(sol[0]).real for sol in resumed.solutions.values())
    assert np.allclose(values, [-(0.5 ** 0.5), 0.5 ** 0.5], atol=1e-10)


def test_different_systems_get_different_runs(tmp_path):
    lg = Ledger(tmp_path)
    a = ensure_solved(circle_line(), lg)

    x = Variable("x")
    other = System()
    other.add_variable_group(VariableGroup([x]))
    other.add_function(x**3 - 2)
    b = ensure_solved(other, lg)

    assert a.run_id != b.run_id
    assert b.num_computed == 3                # cube roots of 2


# ---- chains: the parameter-homotopy two-link lineage --------------------------------

def circle_family(a):
    """x^2 + y^2 = a intersect x = y: the family the parameter homotopy sweeps."""
    x, y = Variable("x"), Variable("y")
    s = System()
    s.add_variable_group(VariableGroup([x, y]))
    s.add_function(x**2 + y**2 - a)
    s.add_function(x - y)
    return s


def test_continuation_chains_and_answers_correctly(tmp_path):
    lg = Ledger(tmp_path)
    generic = circle_family(13)
    ensure_solved(generic, lg)

    target = circle_family(2)                     # solutions x = y = +-1
    cont = ensure_continued(target, generic, lg)
    assert cont.num_computed == 2

    values = sorted(complex(sol[0]).real for sol in cont.solutions.values())
    assert np.allclose(values, [-1.0, 1.0], atol=1e-10)

    # rerun of the whole chain is a no-op at every link
    again = ensure_continued(circle_family(2), circle_family(13), lg)
    assert again.num_computed == 0 and again.num_reused == 2


def test_continuation_requires_recorded_generic(tmp_path):
    lg = Ledger(tmp_path)
    with pytest.raises(ValueError, match="no recorded solve"):
        ensure_continued(circle_family(2), circle_family(13), lg)


def test_crash_mid_continuation_resumes(tmp_path):
    lg = Ledger(tmp_path)
    ensure_solved(circle_family(13), lg)
    with pytest.raises(SimulatedCrash):
        ensure_continued(circle_family(5), circle_family(13), lg, crash_after=1)
    resumed = ensure_continued(circle_family(5), circle_family(13), lg)
    assert resumed.num_reused == 1 and resumed.num_computed == 1


def test_provenance_walks_back_to_the_beginning(tmp_path):
    lg = Ledger(tmp_path)
    generic = circle_family(13)
    base = ensure_solved(generic, lg)
    cont = ensure_continued(circle_family(7), generic, lg)

    chain = provenance_chain(lg, cont.run_id, 0)
    # link 1: the continuation's own track record; link 2: the generic run's track
    # record it started from; terminator: the total-degree start label.
    assert chain[0]["run"] == cont.run_id
    assert chain[1]["run"] == base.run_id
    assert "start_label" in chain[-1]
    assert len(chain) == 3


def test_chains_nest_and_provenance_reaches_depth(tmp_path):
    """sample <- midpoint <- witness <- start label: continuation off a continuation."""
    lg = Ledger(tmp_path)
    ensure_solved(circle_family(13), lg)                                  # witness
    ensure_continued(circle_family(5), circle_family(13), lg)             # midpoint
    sample = ensure_continued(circle_family(3), circle_family(5), lg)     # sample off midpoint

    chain = provenance_chain(lg, sample.run_id, 0)
    assert len(chain) == 4                        # three track links + the start label
    assert "start_label" in chain[-1]


def test_annotations_round_trip(tmp_path):
    lg = Ledger(tmp_path)
    result = ensure_solved(circle_family(13), lg)
    annotate(lg, result.run_id, 0, "projection", 1.5)
    annotate(lg, result.run_id, 0, "edge", "top")
    assert annotations_for(lg, result.run_id, 0) == {"projection": 1.5, "edge": "top"}
    assert annotations_for(lg, result.run_id, 1) == {}


def test_ledger_is_plain_text(tmp_path):
    """The no-special-software property: grep-able journals, readable objects."""
    lg = Ledger(tmp_path)
    ensure_solved(circle_line(), lg)
    journal_text = "".join(p.read_text() for p in (tmp_path / "journals").glob("*.jsonl"))
    assert '"kind":"run"' in journal_text
    assert '"kind":"track"' in journal_text
    # the target system object is its classic input -- readable by a human or bertini 1
    import json
    run = next(r for r in lg.scan() if r["kind"] == "run")
    assert b"function" in lg.get_object(run["target_object"])
