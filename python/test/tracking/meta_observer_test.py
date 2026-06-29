"""Tests for the path-collecting (meta-)observers in bertini.tracking.

PathCollectionObserver (A) attaches a fresh PathDataCollector (B) to the tracker
on every TrackingStarted and harvests it on TrackingEnded, exercising the
attach/detach-during-dispatch machinery from the observer refactor.
"""

import numpy as np
import pytest

import bertini as pb
import bertini.tracking as tk
from bertini.multiprec import Complex as mpfr_complex
from bertini import System, VariableGroup, Variable
from bertini.tracking import amp_config_from
from bertini.nag_algorithm import ZeroDimSolver


def _force_serial(solver):
    """Pin a ZeroDimSolver solver to a single thread (no pool).

    The default solve is multi-threaded, where each path runs on a thread-local tracker clone --
    so an observer attached to the solver's MEMBER tracker sees nothing.  Tests that exercise the
    member-tracker attachment pattern force serial; threaded collection is covered separately via
    SolutionPathCollector / event.tracker() in threaded_solve_test.py.
    """
    cfg = solver.get_config(pb.nag_algorithm.ZeroDimConfig)
    cfg.num_threads = 1
    solver.set_config(cfg)
    return solver


# ---------------------------------------------------------------------------
# Bare-tracker: deterministic, no endgame noise -> exactly one series per track
# ---------------------------------------------------------------------------

@pytest.fixture
def amp_tracker():
    """y - t = 0 on one variable y, with a configured AMP tracker."""
    y = Variable("y")
    tt = Variable("t")
    s = System()
    vg = VariableGroup()
    vg.append(y)
    s.add_function(y - tt)
    s.add_path_variable(tt)
    s.add_variable_group(vg)

    tracker = tk.AMPTracker(s)
    tracker.setup(tk.Predictor.Euler, 1e-5, 1e5, tk.SteppingConfig(), tk.NewtonConfig())
    tracker.precision_setup(amp_config_from(s))
    return tracker, s


def _track_once(tracker, s):
    end = np.zeros(s.num_variables(), dtype=mpfr_complex)
    tracker.track_path(end, mpfr_complex(1), mpfr_complex(0), np.array([mpfr_complex(1)]))
    return end


def test_meta_observer_one_series_per_track(amp_tracker):
    tracker, s = amp_tracker
    a = tk.observers.amp.PathCollectionObserver()
    tracker.add_observer(a)

    _track_once(tracker, s)   # path 1
    _track_once(tracker, s)   # path 2

    tracker.remove_observer(a)

    # one finished collector per track_path call, in order
    assert len(a.series) == 2
    for b in a.series:
        assert len(b) > 0
        # start time was captured from TrackingStarted
        assert abs(b.start_time) == pytest.approx(1.0, abs=1e-9)


def test_collector_array_shapes_and_dtypes(amp_tracker):
    tracker, s = amp_tracker
    a = tk.observers.amp.PathCollectionObserver()
    tracker.add_observer(a)
    _track_once(tracker, s)
    tracker.remove_observer(a)

    b = a.series[0]
    n = len(b)
    times = b.times()
    points = b.points()
    diag = b.diagnostics()

    assert times.shape == (n,)
    assert times.dtype == np.complex128
    assert points.shape == (n, s.num_variables())
    assert points.dtype == np.complex128
    assert diag.shape == (n, len(b.DIAGNOSTIC_COLUMNS))
    assert diag.dtype == np.float64
    # condition numbers and stepsizes are positive reals
    cond = diag[:, b.DIAGNOSTIC_COLUMNS.index("condition_number")]
    step = diag[:, b.DIAGNOSTIC_COLUMNS.index("stepsize")]
    assert np.all(cond > 0)
    assert np.all(step > 0)


def test_collector_as_dataframe(amp_tracker):
    pd = pytest.importorskip("pandas")
    tracker, s = amp_tracker
    a = tk.observers.amp.PathCollectionObserver()
    tracker.add_observer(a)
    _track_once(tracker, s)
    tracker.remove_observer(a)

    df = a.series[0].as_dataframe()
    n = len(a.series[0])
    assert len(df) == n
    assert "t" in df.columns
    assert "z0" in df.columns
    for name in a.series[0].DIAGNOSTIC_COLUMNS:
        assert name in df.columns


def test_no_observers_leak_after_run(amp_tracker):
    """After a path, the per-path collector is detached: a later observer-free
    run must not append to the old collector."""
    tracker, s = amp_tracker
    a = tk.observers.amp.PathCollectionObserver()
    tracker.add_observer(a)
    _track_once(tracker, s)
    first_len = len(a.series[0])
    tracker.remove_observer(a)

    # run again with no observers; the harvested collector must be frozen
    _track_once(tracker, s)
    assert len(a.series[0]) == first_len


# ---------------------------------------------------------------------------
# Whole ZeroDimSolver solve: the dream -- collect every path of a multi-path run.
# ---------------------------------------------------------------------------

def _circle_meets_line():
    # x^2 + y^2 - 1 and x - y -> total degree 2, two nonsingular solutions
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_function(x * x + y * y - 1)
    sys.add_function(x - y)
    return sys


def test_zerodim_lifecycle_events():
    """A nag observer attached to the ZeroDimSolver itself sees AlgorithmStarted once,
    AlgorithmComplete once, and a PathStarted/PathComplete per path. event.solver()
    resolves (via RTTI) to the concrete solver with its full API."""
    from bertini._pybertini import nag_algorithms as nag

    sys = _circle_meets_line()
    solver = ZeroDimSolver(sys, mptype='amp')

    started = []
    completed = []
    begins = []
    ends = []
    solver_nsol = []

    class Lifecycle(nag.observers.CustomObserver):
        def Observe(self, e):
            if isinstance(e, nag.observers.AlgorithmStarted):
                started.append(1)
            elif isinstance(e, nag.observers.AlgorithmComplete):
                completed.append(1)
                solver_nsol.append(len(e.solver().all_solutions()))   # concrete solver API
            elif isinstance(e, nag.observers.PathStarted):
                begins.append(e.path_index())
            elif isinstance(e, nag.observers.PathComplete):
                ends.append(e.path_index())

    solver.add_observer(Lifecycle())   # attached to the ZeroDimSolver, co-owned (no ref kept)
    solver.solve()

    assert len(started) == 1
    assert len(completed) == 1
    assert len(begins) == 2 and len(ends) == 2     # two total-degree paths
    assert sorted(begins) == [0, 1]
    assert solver_nsol == [2]                       # e.solver() gave the real solver


def test_nag_observer_rejected_on_tracker_and_vice_versa():
    """The attach guard distinguishes observables: a nag observer can't attach to a
    tracker, and a tracker observer can't attach to the solver."""
    from bertini._pybertini import nag_algorithms as nag

    sys = _circle_meets_line()
    solver = ZeroDimSolver(sys, mptype='amp')

    with pytest.raises(TypeError):
        solver.get_tracker().add_observer(nag.observers.CustomObserver())
    with pytest.raises(TypeError):
        solver.add_observer(tk.observers.amp.CustomObserver())


def test_solution_path_collector_one_series_per_solution_path():
    """The two-level meta-observer: attached to the SOLVER, it yields exactly one
    collector per solution path, each capturing the whole journey (main homotopy
    track AND endgame sub-tracks), with no start-time filtering."""
    from bertini.nag_algorithm import SolutionPathCollector

    sys = _circle_meets_line()
    solver = ZeroDimSolver(sys, mptype='amp')

    a = SolutionPathCollector()
    solver.add_observer(a)
    solver.solve()

    assert len(solver.all_solutions()) == 2
    assert len(a.series) == 2                       # exactly one per solution path
    for path in a.series:
        assert len(path) > 0
        assert hasattr(path, "path_index")
        # the tracked point lives in the homogenized space the start system uses
        assert path.points().shape[1] >= sys.num_variables()

    # the collected indices are the two total-degree paths
    assert sorted(p.path_index for p in a.series) == [0, 1]


def test_solution_path_collector_captures_more_than_the_main_track():
    """A SolutionPathCollector series (whole path incl. endgame) has at least as many
    steps as the bare main homotopy track alone -- it picks up the endgame sub-tracks
    that the tracker-level start-time filter discards."""
    from bertini.nag_algorithm import SolutionPathCollector

    sys = _circle_meets_line()

    # Pin serial: this test attaches a collector directly to the solver's MEMBER tracker
    # (solver2 below), which only runs the paths in serial mode.  The default solve is threaded,
    # where paths run on thread-local clones -- see threaded_solve_test.py for the threaded path,
    # which collects via event.tracker().
    solver = ZeroDimSolver(sys, mptype='amp')
    _force_serial(solver)
    a = SolutionPathCollector()
    solver.add_observer(a)
    solver.solve()
    whole_path_steps = sum(len(p) for p in a.series)

    # tracker-level collector keeps only the main tracks (|t| start > 0.5)
    solver2 = ZeroDimSolver(sys, mptype='amp')
    _force_serial(solver2)
    b = tk.observers.amp.PathCollectionObserver()
    solver2.get_tracker().add_observer(b)
    solver2.solve()
    main_only_steps = sum(len(s) for s in b.series if abs(s.start_time) > 0.5)

    assert whole_path_steps >= main_only_steps > 0


def test_zerodim_solve_collects_all_paths():
    sys = _circle_meets_line()
    solver = ZeroDimSolver(sys, mptype='amp')
    # Attaching to the solver's member tracker collects only in serial mode; the default solve is
    # threaded (paths run on clones).  For threaded collection use SolutionPathCollector /
    # event.tracker() -- see threaded_solve_test.py.
    _force_serial(solver)

    a = tk.observers.amp.PathCollectionObserver()
    solver.get_tracker().add_observer(a)
    solver.solve()

    assert len(solver.all_solutions()) == 2

    # The solver reuses one tracker for both the homotopy paths AND the endgame
    # sub-tracks, so we collect more than two series ...
    assert len(a.series) >= 2

    # ... but the *main* homotopy paths are exactly those that start at the
    # global start time (|t| = 1); the endgame sub-tracks start near the
    # endgame boundary.  There must be one main path per start-system path.
    main = [b for b in a.series if b.start_time is not None and abs(b.start_time) > 0.5]
    assert len(main) == 2
    # the tracked points live in the homogenized space the start system works in,
    # so there is (at least) one coordinate per original variable, consistent
    # across the paths.
    widths = {b.points().shape[1] for b in main}
    assert len(widths) == 1
    assert widths.pop() >= sys.num_variables()
    for b in main:
        assert len(b) > 0
