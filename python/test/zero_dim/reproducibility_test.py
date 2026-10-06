"""A path's result does not depend on what ran before it, nor on how many threads ran (#378).

Trackers and endgames are reused path after path: by the solver on each thread, and by anyone
driving one by hand.  These tests ask for bit-identical results -- every coordinate compared
through its full-precision repr, and every metadata field that is not a wall-clock time -- where
the older tests asked only for the same solution set within a tolerance, which this bug passed.

The C++ suite is the gate (core/test/nag_algorithms/threaded_solve.cpp,
core/test/tracking_basics/amp_tracker_test.cpp, the generic endgame tests); these confirm the
same holds through the Python interface.
"""

import numpy as np
import pytest

import bertini as pb
import bertini.multiprec as mp
from bertini import ZeroDimSolver
from bertini.tracking import AMPTracker, amp_config_from, Predictor, SteppingConfig, NewtonConfig


def _cyclic(n):
    xs = pb.variables('x', n)
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup(list(xs)))
    for k in range(1, n):
        total = None
        for i in range(n):
            term = None
            for j in range(k):
                v = xs[(i + j) % n]
                term = v if term is None else term * v
            total = term if total is None else total + term
        sys.add_function(total)
    prod = None
    for v in xs:
        prod = v if prod is None else prod * v
    sys.add_function(prod - 1)
    return sys


# every metadata field a path's computation determines; path_time_seconds is a wall-clock time
_FIELDS = ('endgame_success_code', 'pre_endgame_success_code', 'num_successful_steps',
           'num_failed_steps', 'max_precision_used', 'precision_changed',
           'time_of_first_prec_increase', 'final_time_used', 'condition_number',
           'accuracy_estimate_internal_coords', 'accuracy_estimate_user_coords', 'cycle_num',
           'function_residual', 'newton_residual', 'multiplicity', 'is_finite', 'is_real',
           'is_singular')


def _path_by_path(endgame, num_threads):
    pb.random.set_random_seed(20261005)
    solver = ZeroDimSolver(_cyclic(5), endgame=endgame)
    cfg = solver.get_config(pb.nag_algorithm.ZeroDimConfig)
    cfg.num_threads = num_threads
    solver.set_config(cfg)
    solver.solve()
    out = []
    for point, md in zip(solver.all_solutions(), solver.solution_metadata()):
        out.append(([repr(c) for c in point],
                    {f: repr(getattr(md, f)) for f in _FIELDS}))
    return out


@pytest.mark.parametrize("endgame", ["powerseries", "cauchy"])
def test_threaded_solve_is_bit_identical_to_serial(endgame, monkeypatch):
    # the variable would override num_threads and make every solve here serial
    monkeypatch.delenv('BERTINI_NUM_THREADS', raising=False)
    serial = _path_by_path(endgame, 1)
    assert len(serial) == 120
    for num_threads in (2, 4):
        threaded = _path_by_path(endgame, num_threads)
        assert len(threaded) == len(serial)
        for index, (a, b) in enumerate(zip(serial, threaded)):
            assert a[0] == b[0], 'path {} endpoint differs at {} threads'.format(index, num_threads)
            for field in _FIELDS:
                assert a[1][field] == b[1][field], \
                    'path {} {} differs at {} threads: {} vs {}'.format(
                        index, field, num_threads, a[1][field], b[1][field])


@pytest.mark.parametrize("endgame", ["powerseries", "cauchy"])
def test_max_precision_used_covers_every_track_of_a_path(endgame):
    """A path's max_precision_used is at least every precision its tracker stepped at, across the
    track to the endgame boundary and every endgame sub-track.  (On this solve it once reported 30
    for a path whose endgame had run at 40: the record started afresh with every sub-track.)"""
    pb.random.set_random_seed(20261005)
    solver = ZeroDimSolver(_cyclic(5), endgame=endgame)
    collector = pb.nag_algorithm.SolutionPathCollector()
    solver.add_observer(collector)
    solver.solve()

    metadata = solver.solution_metadata()
    assert len(collector.series) == len(metadata) == 120
    for path in collector.series:
        stepped = path.diagnostics()[:, 2]           # the precision column, at every accepted step
        if len(stepped) == 0:
            continue
        reported = int(metadata[path.path_index].max_precision_used)
        assert reported >= int(stepped.max()), \
            'path {}: max_precision_used {} below the {} it stepped at'.format(
                path.path_index, reported, int(stepped.max()))


def _double_root_tracker():
    """x^2 - t, y - x: toward a double root at t = 0, tracked at tolerance 1e-12 so that a track
    to t = 1e-20 ends in multiple precision.  The seed is set before building, so every tracker
    made here draws the same condition-number probe."""
    x, y, t = pb.Variable('x'), pb.Variable('y'), pb.Variable('t')
    s = pb.System()
    s.add_function(x**2 - t)
    s.add_function(y - x)
    s.add_variable_group(pb.VariableGroup([x, y]))
    s.add_path_variable(t)
    mp.default_precision(30)
    pb.random.set_random_seed(20261005)
    tracker = AMPTracker(s)
    tracker.setup(Predictor.RKF45, 1e-12, 1e5, SteppingConfig(), NewtonConfig())
    tracker.precision_setup(amp_config_from(s))
    return s, tracker


def _track(tracker, digits, t_end):
    mp.default_precision(digits)
    start = np.array([mp.complex_mp(1), mp.complex_mp(1)])
    end = np.array([mp.complex_mp(0)] * 2)
    code = tracker.track_path(end, mp.complex_mp(1), mp.complex_mp(t_end), start)
    return (code, tracker.num_total_steps_taken(), tracker.current_precision(),
            repr(tracker.current_stepsize()), repr(tracker.latest_condition_number()),
            [repr(c) for c in end])


@pytest.mark.parametrize("before_digits, digits", [(60, 16), (16, 30), (20, 30), (30, 60)])
def test_a_track_does_not_depend_on_the_track_before_it(before_digits, digits):
    """A hand-driven tracker, reused: the same track after another must match a fresh tracker's."""
    _, used = _double_root_tracker()
    used.precision_preservation(True)      # the previous track ends where it started
    before = _track(used, before_digits, '0.5')
    used.precision_preservation(False)
    assert before[0] == pb.SuccessCode.Success
    assert before[2] != digits

    after_another = _track(used, digits, '1e-20')
    _, fresh = _double_root_tracker()
    first_ever = _track(fresh, digits, '1e-20')

    assert first_ever[0] == pb.SuccessCode.Success
    assert first_ever[2] > 16              # ends in multiple precision, where the mp probe is used
    assert after_another == first_ever
