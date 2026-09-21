"""Printing a bertini object says something about it, never an address (#99).

The first test is the ratchet: it walks every module under ``bertini``, finds every public
class, and fails if any of them still inherits the default ``__repr__`` -- the one that
prints a type name and a memory address.  A binding added without a repr fails here rather
than reaching a user.

The rest check that the descriptions are true: that a tracker reports where the path
actually is, that a collector's count matches what it collected, and that the types whose
repr used to print pointers now print their contents.
"""

import inspect
import types

import numpy as np
import pytest

import bertini as pb
import bertini.tracking as tracking
import bertini.endgame as eg
import bertini.nag_algorithm as nag
from bertini.multiprec import complex_mp


DEFAULT_REPR_OWNERS = ('object', 'instance')


def _public_classes():
    """Every class reachable by a public name under ``bertini``, by class object."""
    seen_modules = set()
    found = {}

    def walk(module, prefix, depth=0):
        if depth > 4 or id(module) in seen_modules:
            return
        seen_modules.add(id(module))
        for name in dir(module):
            if name.startswith('_'):
                continue
            try:
                value = getattr(module, name)
            except Exception:
                continue
            if isinstance(value, types.ModuleType):
                if getattr(value, '__name__', '').startswith('bertini'):
                    walk(value, prefix + '.' + name, depth + 1)
            elif inspect.isclass(value) and getattr(value, '__module__', '').startswith('bertini'):
                found.setdefault(value, []).append('%s.%s' % (prefix, name))

    walk(pb, 'bertini')
    return found


def _defining_class(cls, method):
    for ancestor in cls.__mro__:
        if method in ancestor.__dict__:
            return ancestor.__name__
    return '-'


def test_no_public_class_prints_an_address():
    classes = _public_classes()
    assert len(classes) > 100, 'the walk found suspiciously few classes'

    stragglers = sorted(sorted(names)[0] for cls, names in classes.items()
                        if _defining_class(cls, '__repr__') in DEFAULT_REPR_OWNERS
                        and _defining_class(cls, '__str__') in DEFAULT_REPR_OWNERS)

    assert not stragglers, (
        '%d public classes still print a type name and a memory address; give each one a '
        '__repr__ in bertini/_repr.py:\n  %s' % (len(stragglers), '\n  '.join(stragglers)))


@pytest.fixture
def homotopy():
    """y^2 - 1 at t=1, deformed to y^2 - 4 at t=0.  Start it from y=1."""
    y, t = pb.Variable('y'), pb.Variable('t')
    system = pb.System()
    system.add_function((y * y - 1) * t + (y * y - 4) * (1 - t))
    system.add_path_variable(t)
    system.add_variable_group(pb.VariableGroup([y]))
    return system


@pytest.fixture
def tracker(homotopy):
    tracked = pb.AMPTracker(homotopy)
    tracked.setup(pb.Predictor.RKF45, 1e-5, 1e5,
                  tracking.SteppingConfig(), tracking.NewtonConfig())
    tracked.precision_setup(tracking.amp_config_from(homotopy))
    return tracked


def _track(tracker, homotopy):
    endpoint = np.zeros(homotopy.num_variables(), dtype=complex_mp)
    code = tracker.track_path(endpoint, complex_mp(1), complex_mp(0),
                              np.array([complex_mp(1)]))
    assert code == pb.SuccessCode.Success
    return endpoint


def test_a_tracker_describes_its_system_and_settings(tracker):
    described = repr(tracker)

    assert described.startswith('<AMPTracker:')
    assert 'homotopy in 1 variable' in described
    assert 'RKF45' in described
    assert 'adaptive precision' in described
    assert 'not started' in described


def test_a_tracker_reports_where_the_path_got_to(tracker, homotopy):
    _track(tracker, homotopy)
    described = repr(tracker)

    assert 'not started' not in described
    assert 'at t=0' in described                      # the path ran to t=0
    assert 'after %d steps' % tracker.num_total_steps_taken() in described


def test_a_fresh_tracker_has_taken_no_steps(homotopy):
    """The per-path counters are zeroed at construction, not only inside TrackPath."""
    assert pb.AMPTracker(homotopy).num_total_steps_taken() == 0


def test_an_event_says_what_changed_and_where(tracker, homotopy):
    described = []

    class Watch(tracking.observers.amp.CustomObserver):
        def Observe(self, event):
            described.append(repr(event))

    watch = Watch()
    tracker.add_observer(watch)
    _track(tracker, homotopy)
    tracker.remove_observer(watch)

    assert described, 'no events were observed'
    assert not any('0x' in text for text in described)

    precision_changes = [text for text in described if text.startswith('<PrecisionChanged:')]
    if precision_changes:
        assert 'digits' in precision_changes[0]
        assert '->' in precision_changes[0]


def test_a_collector_counts_what_it_collected(tracker, homotopy):
    collector = tracking.observers.amp.PathDataCollector()
    assert 'nothing recorded' in repr(collector)

    tracker.add_observer(collector)
    _track(tracker, homotopy)
    tracker.remove_observer(collector)

    assert repr(collector) == '<PathDataCollector: %d steps>' % len(collector.times())


def test_a_path_collector_counts_paths_and_steps():
    x, y = pb.Variable('x'), pb.Variable('y')
    system = pb.System()
    system.add_variable_group(pb.VariableGroup([x, y]))
    system.add_function(x * x + y * y - 1)
    system.add_function(y - x * x)

    paths = pb.SolutionPathCollector()
    solver = pb.ZeroDimSolver(system, mptype='adaptive')
    solver.add_observer(paths)
    solver.solve()

    described = repr(paths)
    assert '%d paths' % len(paths.series) in described
    assert 'steps' in described


def test_an_endgame_describes_its_interval(tracker):
    described = repr(eg.AMPPowerSeriesEndgame(tracker, complex_mp('0.1')))

    assert described.startswith('<AMPPowerSeriesEndgame:')
    assert 'from t=0.1' in described
    assert 'to t=0' in described
    assert 'not yet run' in described


def test_a_decomposition_gives_its_shape():
    matrix = np.array([[complex_mp(1), complex_mp(0)], [complex_mp(0), complex_mp(1)]])
    assert repr(pb.linalg.PartialPivLU(matrix)) == '<PartialPivLU: 2x2>'


def test_a_variable_group_prints_its_variables():
    """It used to print the addresses of the variables it holds."""
    x, y = pb.Variable('x'), pb.Variable('y')
    assert repr(pb.VariableGroup([x, y])) == 'VariableGroup([x, y])'


def test_a_list_of_variable_groups_prints_its_groups():
    x, y = pb.Variable('x'), pb.Variable('y')
    system = pb.System()
    system.add_variable_group(pb.VariableGroup([x, y]))
    system.add_function(x * y - 1)

    assert repr(system.variable_groups()) == '[VariableGroup([x, y])]'


def test_a_solve_result_names_only_what_it_has():
    x, y = pb.Variable('x'), pb.Variable('y')
    system = pb.System()
    system.add_variable_group(pb.VariableGroup([x, y]))
    system.add_function(x * x + y * y - 1)
    system.add_function(y - x * x)

    described = repr(pb.solve(system))

    assert described.startswith('SolveResult(')
    assert 'solutions' in described
    assert 'None' not in described          # it used to say "records at None"
    assert 'run ,' not in described         # ... and "run " with nothing after it


def test_the_types_with_nothing_to_report_still_name_themselves():
    assert repr(tracking.observers.amp.GoryDetailLogger()) == '<GoryDetailLogger>'
    assert repr(tracking.observers.amp.FirstPrecisionRecorder()) == '<FirstPrecisionRecorder>'


def test_no_repr_raises_on_a_freshly_built_object(homotopy):
    """A repr that throws turns a print in a debugging session into a traceback."""
    fresh = [
        pb.AMPTracker(homotopy),
        pb.DoublePrecisionTracker(homotopy),
        pb.MultiplePrecisionTracker(homotopy),
        tracking.observers.amp.PathDataCollector(),
        tracking.observers.amp.GoryDetailLogger(),
        tracking.observers.amp.CallbackObserver(),
        pb.SolutionPathCollector(),
        eg.observers.amp_pseg.SampleSequenceCollector(),
    ]

    for thing in fresh:
        text = repr(thing)
        assert text.startswith('<') and '0x' not in text
