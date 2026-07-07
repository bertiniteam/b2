"""MPI-less shared-memory threading for the ZeroDimSolver solve, exercised from Python.

The contract: a threaded solve returns the SAME solutions as a serial one, Python observers
fire safely from worker threads (the GIL is released for the C++ solve and re-acquired in the
observer trampoline), and the meta-observer per-path collection keeps working because it attaches
to event.tracker() -- the thread-local clone -- not the member tracker.

These need no MPI and no free-threaded Python, so they run on every platform's wheel.
"""

import pytest

import bertini as pb
import bertini.tracking as tk
from bertini import ZeroDimSolver, SolutionPathCollector


def _two_cubics():
    """{x^3 - x, y^3 - y}: 9 total-degree paths, 9 finite real solutions, well separated."""
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_function(x**3 - x)
    sys.add_function(y**3 - y)
    return sys


def _solve_with(num_threads):
    solver = ZeroDimSolver(_two_cubics(), mptype='amp')
    cfg = solver.get_config(pb.nag_algorithm.ZeroDimConfig)
    cfg.num_threads = num_threads
    solver.set_config(cfg)
    solver.solve()
    return solver


def _solution_key_set(solver, tol=6):
    """A hashable, order-independent representation of the solution set."""
    keys = set()
    for v in solver.all_solutions():
        if len(v) == 0:
            continue
        keys.add(tuple(sorted((round(complex(c).real, tol), round(complex(c).imag, tol))
                              for c in v)))
    return keys


def test_num_threads_config_roundtrips():
    solver = ZeroDimSolver(_two_cubics(), mptype='amp')
    cfg = solver.get_config(pb.nag_algorithm.ZeroDimConfig)
    assert cfg.num_threads == 0          # default: auto
    cfg.num_threads = 3
    solver.set_config(cfg)
    assert solver.get_config(pb.nag_algorithm.ZeroDimConfig).num_threads == 3


@pytest.mark.parametrize("num_threads", [2, 4, 8])
def test_threaded_matches_serial(num_threads):
    """The whole point: more threads, identical answers."""
    serial = _solve_with(1)
    threaded = _solve_with(num_threads)

    assert len(serial.all_solutions()) == 9
    assert _solution_key_set(threaded) == _solution_key_set(serial)


def test_lifecycle_observer_fires_under_threads():
    """A Python observer attached to the solver sees the right event counts even though the
    per-path events fire from C++ worker threads -- proves the GIL is re-acquired (no crash) and
    notifications are serialized (no lost/garbled counts)."""
    from bertini._pybertini import nag_algorithms as nag

    solver = ZeroDimSolver(_two_cubics(), mptype='amp')
    cfg = solver.get_config(pb.nag_algorithm.ZeroDimConfig)
    cfg.num_threads = 4
    solver.set_config(cfg)

    starts, completes, begins = [], [], []

    class Lifecycle(nag.observers.CustomObserver):
        def Observe(self, e):
            if isinstance(e, nag.observers.AlgorithmStarted):
                starts.append(1)
            elif isinstance(e, nag.observers.AlgorithmComplete):
                completes.append(1)
            elif isinstance(e, nag.observers.PathStarted):
                begins.append(e.path_index())

    solver.add_observer(Lifecycle())
    solver.solve()

    assert len(starts) == 1
    assert len(completes) == 1
    assert len(begins) == 9                     # one per path, none lost to a race
    assert sorted(begins) == list(range(9))     # every index exactly once


def test_event_tracker_is_attachable_and_distinct_under_threads():
    """event.tracker() returns a concrete, attachable tracker; under threads it is NOT the
    solver's member tracker (it's the thread-local clone)."""
    from bertini._pybertini import nag_algorithms as nag

    solver = ZeroDimSolver(_two_cubics(), mptype='amp')
    cfg = solver.get_config(pb.nag_algorithm.ZeroDimConfig)
    cfg.num_threads = 4
    solver.set_config(cfg)

    member = solver.get_tracker()
    saw_member = []
    have_observers_attr = []

    class Probe(nag.observers.CustomObserver):
        def Observe(self, e):
            if isinstance(e, nag.observers.PathStarted):
                t = e.tracker()
                have_observers_attr.append(hasattr(t, "add_observer") and hasattr(t, "observers"))
                # identity by the underlying C++ object: boost.python compares by pointer
                saw_member.append(t == member)

    solver.add_observer(Probe())
    solver.solve()

    assert all(have_observers_attr)        # the concrete tracker API is present
    assert not any(saw_member)             # never the member tracker -> a clone


def test_solution_path_collector_under_threads():
    """The two-level meta-observer collects exactly one series per path under threading, because
    it now attaches to event.tracker() (the clone that actually runs the path)."""
    solver = ZeroDimSolver(_two_cubics(), mptype='amp')
    cfg = solver.get_config(pb.nag_algorithm.ZeroDimConfig)
    cfg.num_threads = 4
    solver.set_config(cfg)

    collector = SolutionPathCollector()
    solver.add_observer(collector)
    solver.solve()

    assert len(solver.all_solutions()) == 9
    assert len(collector.series) == 9                          # one per solution path
    assert sorted(p.path_index for p in collector.series) == list(range(9))
    for path in collector.series:
        assert len(path) > 0                                   # actually captured steps
