# This file is part of Bertini 2.
#
# python/test/records/records_test.py is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Bertini 2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
# See the GNU General Public License for more details.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license, as well as COPYING.

"""The casual records surface (arc rung 5): solve / save / load, the Solution type,
and ensure-answered recall.  Correctness of the seam lives in C++
(test_nag_algorithms/zero_dim_records); these pin the Python feel."""

import json

import numpy as np
import pytest

import bertini as pb
from bertini import Variable, VariableGroup, System


def circle_line():
    x, y = Variable('x'), Variable('y')
    s = System()
    s.add_variable_group(VariableGroup([x, y]))
    s.add_function(x**2 + y**2 - 1)
    s.add_function(x - y)
    return s


def test_solve_records_and_rerun_recalls(tmp_path):
    d = str(tmp_path / 'records')
    first = pb.solve(circle_line(), seed=42, directory=d)
    assert len(first) == 2
    assert first.num_recalled == 0
    assert first.run_id

    again = pb.solve(circle_line(), seed=42, directory=d)
    assert again.num_recalled == 2          # ensure-answered: nothing recomputed
    assert again.run_id == first.run_id
    for a, b in zip(first, again):
        assert all(abs(complex(u) - complex(v)) < 1e-12 for u, v in zip(a, b))


def test_solutions_are_points_that_remember(tmp_path):
    r = pb.solve(circle_line(), seed=42, directory=str(tmp_path / 'records'))
    s = r[0]
    assert isinstance(s, pb.Solution)
    assert s.provenance == {'run': r.run_id, 'index': s.provenance['index']}
    # behaves like the array you expect
    assert abs(complex(s[0])) == pytest.approx(2 ** -0.5, abs=1e-8)
    # a derived point is a NEW point: provenance honestly absent
    doubled = 2 * s
    assert doubled.provenance is None


def test_save_and_load_round_trip(tmp_path):
    d = str(tmp_path / 'records')
    r = pb.solve(circle_line(), seed=42, directory=d)

    pb.save('my favorites', r, directory=d)
    pb.save('notes', {'count': 2, 'nice': True}, directory=d)

    everything = pb.load(directory=d)
    assert 'my favorites' in everything and 'notes' in everything
    assert 'solutions [run %s]' % r.run_id in everything   # solve auto-declares

    favorites = pb.load('my favorites', directory=d)
    assert len(favorites['points']) == 2
    pt = favorites['points'][0]
    assert pt['status'] == 'success'
    # endpoints are recorded in INTERNAL coordinates (homogenized: x, y + hom var),
    # the representation restart needs; user-coordinate rendering is a noted follow-up
    assert {'x', 'y'} <= set(pt['coordinates'])
    assert pb.load('notes', directory=d)['value'] == {'count': 2, 'nice': True}


def test_records_are_plain_json(tmp_path):
    """The no-special-software property, from Python's side: raw json suffices."""
    d = tmp_path / 'records'
    pb.solve(circle_line(), seed=7, directory=str(d))

    assert (d / 'README.txt').exists()
    assert (d / 'INDEX.txt').exists()
    machine = json.loads((d / 'results.json').read_text())
    assert isinstance(machine, dict) and machine   # one json.load away

    kinds = []
    for journal in (d / 'history').glob('*.jsonl'):
        for line in journal.read_text().splitlines():
            if line.strip():
                kinds.append(json.loads(line)['kind'])
    assert 'run' in kinds and 'track' in kinds and 'result' in kinds


def test_nameless_save(tmp_path):
    d = str(tmp_path / 'records')
    r = pb.solve(circle_line(), seed=42, directory=d)
    name = pb.save(r, directory=d)               # save(stuff): magic happens
    assert name.startswith('saved ')
    assert len(pb.load(name, directory=d)['points']) == 2


def test_ambient_directory_resolution(tmp_path, monkeypatch):
    import bertini.records as records_module
    monkeypatch.setattr(records_module, '_ambient', None)
    monkeypatch.setenv('BERTINI_RECORDS_DIR', str(tmp_path / 'ambient'))
    assert pb.records_dir() == str(tmp_path / 'ambient')
    monkeypatch.setattr(records_module, '_ambient', None)   # don't leak


def test_records_dir_is_the_ambient_switch_for_bare_solvers(tmp_path, monkeypatch):
    """bertini.records_dir(path): the one line that turns ambient recording on for the
    bare solver classes too -- the path is exported to BERTINI_RECORDS_DIR, which the
    solvers' ambient attach reads; recording(False) still wins."""
    import os
    import bertini.records as records_module
    d = tmp_path / 'ambient_records'
    monkeypatch.setattr(records_module, '_ambient', None)
    monkeypatch.setattr(records_module, '_recording_enabled', True)
    monkeypatch.delenv('BERTINI_RECORDS_DIR', raising=False)

    pb.records_dir(str(d))                                  # the one line
    assert os.environ['BERTINI_RECORDS_DIR'] == str(d)

    solver = pb.ZeroDimSolver(circle_line(), mptype='adaptive')
    solver.solve()                                          # bare solver, no directory named
    assert (d / 'history').exists()
    assert (d / 'results.json').exists()

    # the off switch still wins, and coming back restores the ambient export
    pb.recording(False)
    assert os.environ['BERTINI_RECORDS_DIR'] == ''
    pb.records_dir(str(tmp_path / 'elsewhere'))             # setting while off: no export
    assert os.environ['BERTINI_RECORDS_DIR'] == ''
    pb.recording(True)
    assert os.environ['BERTINI_RECORDS_DIR'] == str(tmp_path / 'elsewhere')


def test_annotate_renders_beside_the_point(tmp_path):
    """annotate(sol, key, value): margin notes land in results.json by the point."""
    d = str(tmp_path / 'records')
    r = pb.solve(circle_line(), seed=42, directory=d)
    sol = r.solutions[0]
    pb.annotate(sol, 'projection', 1.5, directory=d)
    pb.annotate(sol, 'the one i meant', True, directory=d)
    pb.save('picked', r, directory=d)

    entry = pb.load('picked', directory=d)
    by_index = {p['provenance']['index']: p for p in entry['points']}
    notes = by_index[sol.provenance['index']]['annotations']
    assert notes['projection'] == 1.5
    assert notes['the one i meant'] is True
    assert sol.annotations['projection'] == 1.5      # the live object learns it too

    # newest wins per (point, key)
    pb.annotate(sol, 'projection', 2.5, directory=d)
    entry = pb.load('picked', directory=d)
    by_index = {p['provenance']['index']: p for p in entry['points']}
    assert by_index[sol.provenance['index']]['annotations']['projection'] == 2.5


def test_annotate_without_provenance_is_an_error(tmp_path):
    with pytest.raises(ValueError):
        pb.annotate(np.array([1.0, 2.0]), 'key', 'value', directory=str(tmp_path / 'r'))


# --- chains and givens (rung 6) --------------------------------------------------------

def circle_line_r(r2):
    x, y = Variable('x'), Variable('y')
    s = System()
    s.add_variable_group(VariableGroup([x, y]))
    s.add_function(x**2 + y**2 - r2)
    s.add_function(x - y)
    return s


def test_chained_solve_records_point_refs(tmp_path):
    """solve(B, homotopy=H, start=r1): every new track's start is a point_ref into r1 --
    the provenance chain is walkable back to the beginning."""
    from bertini.nag_algorithm import blend_homotopy
    d = str(tmp_path / 'records')
    A = circle_line_r(1)
    r1 = pb.solve(A, seed=42, directory=d)

    B = circle_line_r(4)
    r2 = pb.solve(B, homotopy=blend_homotopy(B, A), start=r1, seed=42, directory=d)
    assert len(r2) == 2
    assert {abs(round(complex(s[0]).real, 10)) for s in r2} == {round(2 ** 0.5, 10)}

    starts = []
    for journal in (tmp_path / 'records' / 'history').glob('*.jsonl'):
        for line in journal.read_text().splitlines():
            rec = json.loads(line)
            if rec.get('kind') == 'track' and rec.get('run') == r2.run_id:
                starts.append(rec['start'])
    assert starts and all(
        st['kind'] == 'point_ref' and st['run'] == r1.run_id for st in starts)


def test_chained_solve_recalls_on_rerun(tmp_path):
    """A chained ask is an ask like any other: the identical rerun computes nothing."""
    from bertini.nag_algorithm import blend_homotopy
    d = str(tmp_path / 'records')
    A = circle_line_r(1)
    B = circle_line_r(4)

    pb.set_random_seed = None   # noqa -- explicit seeds below control everything
    r1 = pb.solve(A, seed=42, directory=d)
    first = pb.solve(B, homotopy=blend_homotopy(B, A), start=r1, seed=7, directory=d)
    assert first.num_recalled == 0

    r1b = pb.solve(A, seed=42, directory=d)      # recalls r1
    assert r1b.num_recalled == 2
    again = pb.solve(B, homotopy=blend_homotopy(B, A), start=r1b, seed=7, directory=d)
    assert again.num_recalled == 2


def test_raw_start_points_become_a_given(tmp_path):
    """Raw start points (no provenance) are archived as a given: provenance bottoms out
    honestly at data the user supplied."""
    from bertini.nag_algorithm import blend_homotopy
    d = str(tmp_path / 'records')
    A = circle_line_r(1)
    B = circle_line_r(4)
    raw = [np.array([2 ** -0.5 + 0j, 2 ** -0.5 + 0j]),
           np.array([-(2 ** -0.5) + 0j, -(2 ** -0.5) + 0j])]
    r = pb.solve(B, homotopy=blend_homotopy(B, A), start=raw, seed=7, directory=d)
    assert len(r) == 2

    givens, starts = [], []
    for journal in (tmp_path / 'records' / 'history').glob('*.jsonl'):
        for line in journal.read_text().splitlines():
            rec = json.loads(line)
            if rec.get('kind') == 'given':
                givens.append(rec)
            if rec.get('kind') == 'track' and rec.get('run') == r.run_id:
                starts.append(rec['start'])
    assert len(givens) == 1 and givens[0]['role'] == 'start_points'
    assert all(st['kind'] == 'given_ref' and st['given'] == givens[0]['source']
               for st in starts)
    # the given's definition is plain json, coordinates readable without bertini,
    # filed by kind with the FULL digest in the filename and an honest extension
    gid = givens[0]['source']
    body = json.loads((tmp_path / 'records' / 'definitions' / 'givens'
                       / gid[:2] / ('given-%s.json' % gid)).read_text())
    assert len(body['points']) == 2


def test_depth_4_provenance_walk(tmp_path):
    """Rung 6 acceptance: four chained solves; the final point's provenance walks back
    through three point_refs to a canonical start label -- all the way to the beginning."""
    from bertini.nag_algorithm import blend_homotopy
    d = str(tmp_path / 'records')

    radii = [1, 4, 9, 16]
    systems = [circle_line_r(r) for r in radii]
    results = [pb.solve(systems[0], seed=42, directory=d)]
    for prior, target in zip(systems, systems[1:]):
        results.append(pb.solve(target, homotopy=blend_homotopy(target, prior),
                                start=results[-1], seed=42, directory=d))
    assert all(len(r) == 2 for r in results)

    trail = pb.provenance(results[-1].solutions[0], directory=d)
    runs_walked = [hop['run'] for hop in trail if 'run' in hop]
    assert runs_walked == [r.run_id for r in reversed(results)]   # depth 4
    assert trail[-1]['kind'] == 'start_label'                     # the beginning

    # solutions_of reads the records cold (no solver object): chains from ANY run
    cold = pb.solutions_of(results[1].run_id, directory=d)
    assert len(cold) == 2
    assert {c.provenance['index'] for c in cold} == {0, 1}
    resumed = pb.solve(systems[2], homotopy=blend_homotopy(systems[2], systems[1]),
                       start=cold, seed=42, directory=d)
    assert resumed.num_recalled == 2      # identical ask as results[2]: pure memo


def test_recording_off_is_one_line(tmp_path, monkeypatch):
    """bertini.recording(False): solve runs bare -- nothing written, nothing consulted,
    provenance honestly absent; recording(True) restores everything."""
    import bertini.records as records_module
    d = tmp_path / 'records'
    monkeypatch.setattr(records_module, '_recording_enabled', True)   # isolate

    pb.recording(False)
    try:
        r = pb.solve(circle_line(), seed=42, directory=str(d))
        assert len(r) == 2
        assert r.run_id == '' and r.solutions[0].provenance is None
        assert not d.exists()
    finally:
        pb.recording(True)

    assert pb.recording() is True
    r2 = pb.solve(circle_line(), seed=42, directory=str(d))
    assert r2.run_id and d.exists()


# --- navigation / visualization tools ---------------------------------------------------

def _small_chain(tmp_path):
    from bertini.nag_algorithm import blend_homotopy
    d = str(tmp_path / 'records')
    members = [circle_line_r(r2) for r2 in (1, 4, 9)]
    results = [pb.solve(members[0], seed=42, directory=d)]
    for prev, tgt in zip(members, members[1:]):
        results.append(pb.solve(tgt, homotopy=blend_homotopy(tgt, prev),
                                start=results[-1], seed=42, directory=d))
    return d, results


def test_runs_and_tracks_dataframes(tmp_path):
    pytest.importorskip('pandas')
    d, results = _small_chain(tmp_path)

    r = pb.runs(directory=d)
    assert len(r) == 3
    assert set(r['run']) == {res.run_id for res in results}
    assert (r['num_paths'] == 2).all() and (r['seed'] == 42).all()
    assert r['producer_version'].notna().all()

    t = pb.tracks(directory=d)
    assert len(t) == 6 and (t['status'] == 'success').all()
    kinds = t.groupby('run')['start_kind'].agg(set)
    assert kinds[results[0].run_id] == {'start_label'}
    assert kinds[results[1].run_id] == {'point_ref'}

    # coordinates stay OUT of the table unless asked (scale guard)
    assert 'endpoint_user' not in t.columns
    tc = pb.tracks(run=results[0].run_id, directory=d, coordinates=True)
    assert len(tc) == 2 and len(tc['endpoint_user'].iloc[0]) == 2


def test_provenance_graph_and_chain_plot(tmp_path):
    pytest.importorskip('pandas')
    nx = pytest.importorskip('networkx')
    matplotlib = pytest.importorskip('matplotlib')
    matplotlib.use('Agg')
    d, results = _small_chain(tmp_path)

    g = pb.provenance_graph(directory=d)
    # 2 origins + 2 points x 3 runs; each run contributes 2 edges
    assert g.number_of_nodes() == 8 and g.number_of_edges() == 6
    # walk one lineage through the graph: final endpoint reaches an origin
    final = (results[-1].run_id, 0)
    ancestors = nx.ancestors(g, final)
    assert any(len(a) == 3 and a[0] == 'start_label' for a in ancestors)

    ax = pb.plot_chain(directory=d)                        # per-path drawing
    assert ax.figure is not None
    ax2 = pb.plot_chain(directory=d, max_paths_drawn=1)    # aggregation fallback
    assert 'aggregated' in ax2.get_title()
    import matplotlib.pyplot as plt
    plt.close('all')


@pytest.mark.xfail(reason="KNOWN BUG, fix owned by a DEDICATED PR (user decision "
                          "2026-07-03; see z_notes/junk-success-cauchy-pole-blindness.md): "
                          "the Cauchy endgame's circle-mean is blind to poles, so on raw "
                          "affine user homotopies toward deficient targets it can report "
                          "Success at a NON-ROOT (function_residual ~ 1; inherited from "
                          "Bertini 1, reproduction committed).  In the affine case infinity "
                          "is really infinity -- the endgame can never converge there.  The "
                          "records stay honest; the plain (patched) zero-dim pipeline is "
                          "correct.  Un-xfail when the endgame PR merges.",
                   strict=False)
def test_chained_deficient_target_paths_never_junk_success(tmp_path):
    from bertini.nag_algorithm import blend_homotopy
    from bertini import Variable
    x, y = Variable('x'), Variable('y')
    A = System(); A.add_variable_group(VariableGroup([x, y]))
    A.add_function(x**2 - 1); A.add_function(x * y - 1)
    B = System(); B.add_variable_group(VariableGroup([x, y]))
    B.add_function(x**2 - 1); B.add_function(y * (x + 1) - 1)   # deficient: 1 finite root

    d = str(tmp_path / 'records')
    r1 = pb.solve(A, seed=42, directory=d)
    r2 = pb.solve(B, homotopy=blend_homotopy(B, A), start=r1, seed=42, directory=d)
    for m in r2.solver.solution_metadata():
        # a path recorded successful must actually sit on a root
        assert float(m.function_residual) < 1e-6
