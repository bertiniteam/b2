import numpy as np
import pytest

import bertini as b
import bertini.multiprec as mp
from bertini.multiprec import Complex as mpfr_complex
from bertini import System, VariableGroup, Variable
from bertini.tracking import amp_config_from
import bertini.tracking as t


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def linear_system():
    """y - t = 0, one variable y, path variable t."""
    y = Variable("y")
    tt = Variable("t")
    s = System()
    vg = VariableGroup()
    vg.append(y)
    s.add_function(y - tt)
    s.add_path_variable(tt)
    s.add_variable_group(vg)
    return s


@pytest.fixture
def amp_tracker(linear_system):
    s = linear_system
    tracker = t.AMPTracker(s)
    tracker.setup(t.Predictor.Euler, 1e-5, 1e5, t.SteppingConfig(), t.NewtonConfig())
    tracker.precision_setup(amp_config_from(s))
    return tracker, s


@pytest.fixture
def double_tracker(linear_system):
    s = linear_system
    tracker = t.DoublePrecisionTracker(s)
    tracker.setup(t.Predictor.Euler, 1e-5, 1e5, t.SteppingConfig(), t.NewtonConfig())
    return tracker, s


def _run_amp(tracker, s):
    t_start = mpfr_complex(1)
    t_end   = mpfr_complex(0)
    y_start = np.array([mpfr_complex(1)])
    y_end   = np.zeros(s.num_variables(), dtype=mpfr_complex)
    tracker.track_path(y_end, t_start, t_end, y_start)
    return y_end


def _run_double(tracker, s):
    t_start = complex(1)
    t_end   = complex(0)
    y_start = np.array([complex(1)])
    y_end   = np.zeros(s.num_variables(), dtype=complex)
    tracker.track_path(y_end, t_start, t_end, y_start)
    return y_end


# ---------------------------------------------------------------------------
# Tests — AMP tracker
# ---------------------------------------------------------------------------

def test_callback_observer_receives_started_and_ended(amp_tracker):
    tracker, s = amp_tracker
    obs = t.observers.amp.CallbackObserver()

    fired = []
    obs.on(t.observers.amp.TrackingStarted, lambda e: fired.append("started"))
    obs.on(t.observers.amp.TrackingEnded,   lambda e: fired.append("ended"))

    tracker.add_observer(obs)
    _run_amp(tracker, s)
    tracker.remove_observer(obs)

    assert "started" in fired
    assert "ended"   in fired


def test_callback_observer_isinstance(amp_tracker):
    """Verify isinstance discrimination works for events arriving in Observe."""
    tracker, s = amp_tracker

    seen_types = []

    class TypeCheckObserver(t.observers.amp.Abstract):
        def Observe(self, event):
            if isinstance(event, t.observers.amp.TrackingStarted):
                seen_types.append("started")
            elif isinstance(event, t.observers.amp.TrackingEnded):
                seen_types.append("ended")
            elif isinstance(event, t.observers.amp.PrecisionChanged):
                seen_types.append("precision_changed")
            elif isinstance(event, t.observers.amp.SuccessfulStep):
                seen_types.append("step")

    obs = TypeCheckObserver()
    tracker.add_observer(obs)
    _run_amp(tracker, s)
    tracker.remove_observer(obs)

    assert "started" in seen_types
    assert "ended"   in seen_types
    assert "step"    in seen_types


def test_callback_observer_precision_changed_data(amp_tracker):
    """PrecisionChanged events expose .previous() and .next() integers."""
    tracker, s = amp_tracker

    precision_changes = []
    obs = t.observers.amp.CallbackObserver()
    obs.on(t.observers.amp.PrecisionChanged,
           lambda e: precision_changes.append((e.previous(), e.next())))

    tracker.add_observer(obs)
    _run_amp(tracker, s)
    tracker.remove_observer(obs)

    for prev, nxt in precision_changes:
        assert isinstance(prev, int)
        assert isinstance(nxt,  int)
        assert prev != nxt


def test_callback_observer_tracker_accessor(amp_tracker):
    """event.tracker() returns the live AMPTracker."""
    tracker, s = amp_tracker

    tracker_refs = []
    obs = t.observers.amp.CallbackObserver()
    obs.on(t.observers.amp.TrackingStarted,
           lambda e: tracker_refs.append(e.tracker()))

    tracker.add_observer(obs)
    _run_amp(tracker, s)
    tracker.remove_observer(obs)

    assert len(tracker_refs) >= 1
    assert isinstance(tracker_refs[0], t.AMPTracker)


def test_tracker_step_diagnostics_accessors(amp_tracker):
    """During a SuccessfulStep, the tracker exposes the step diagnostics the
    path-visualization observers need: condition number, stepsize, delta_t,
    norm of step, error estimate (in addition to point/time/precision)."""
    tracker, s = amp_tracker

    rows = []
    obs = t.observers.amp.CallbackObserver()

    def grab(e):
        trk = e.tracker()
        rows.append((
            trk.current_time(),
            trk.current_point(),
            trk.current_precision(),
            trk.current_stepsize(),
            trk.delta_t(),
            trk.latest_condition_number(),
            trk.latest_norm_of_step(),
            trk.latest_error_estimate(),
        ))

    obs.on(t.observers.amp.SuccessfulStep, grab)

    tracker.add_observer(obs)
    _run_amp(tracker, s)
    tracker.remove_observer(obs)

    assert len(rows) >= 1
    time, point, prec, stepsize, dt, cond, norm_step, err = rows[-1]
    # everything must cast cleanly to plain python numbers for plotting
    assert complex(time) == complex(time)            # not NaN
    assert int(prec) > 0
    assert float(stepsize) > 0.0
    assert complex(dt) == complex(dt)
    assert float(cond) > 0.0
    assert float(norm_step) >= 0.0
    assert float(err) >= 0.0
    assert len(point) == s.num_variables()


def test_observer_self_unsubscribe_via_return(amp_tracker):
    """A python observer can drop itself by returning ObserveResult.Unsubscribe;
    it then receives no further events for the rest of the path."""
    tracker, s = amp_tracker

    counts = {"n": 0}

    class OneShot(t.observers.amp.Abstract):
        def Observe(self, event):
            counts["n"] += 1
            return t.ObserveResult.Unsubscribe

    obs = OneShot()
    tracker.add_observer(obs)
    _run_amp(tracker, s)
    # dropped after the first event, so it saw exactly one
    assert counts["n"] == 1


def test_remove_observer_stops_callbacks(amp_tracker):
    """Removing an observer before track_path means no callbacks fire."""
    tracker, s = amp_tracker

    fired = []
    obs = t.observers.amp.CallbackObserver()
    obs.on(t.observers.amp.TrackingStarted, lambda e: fired.append("started"))

    tracker.add_observer(obs)
    tracker.remove_observer(obs)
    _run_amp(tracker, s)

    assert fired == []


def test_multiple_callbacks_same_event(amp_tracker):
    """Multiple callbacks on the same event type all fire."""
    tracker, s = amp_tracker

    counts = [0, 0]
    obs = t.observers.amp.CallbackObserver()
    obs.on(t.observers.amp.TrackingEnded, lambda e: counts.__setitem__(0, counts[0] + 1))
    obs.on(t.observers.amp.TrackingEnded, lambda e: counts.__setitem__(1, counts[1] + 1))

    tracker.add_observer(obs)
    _run_amp(tracker, s)
    tracker.remove_observer(obs)

    assert counts[0] == 1
    assert counts[1] == 1


# ---------------------------------------------------------------------------
# Tests — DoublePrecisionTracker
# ---------------------------------------------------------------------------

def test_double_tracker_callback_observer(double_tracker):
    tracker, s = double_tracker

    fired = []
    obs = t.observers.double.CallbackObserver()
    obs.on(t.observers.double.TrackingStarted, lambda e: fired.append("started"))
    obs.on(t.observers.double.TrackingEnded,   lambda e: fired.append("ended"))

    tracker.add_observer(obs)
    _run_double(tracker, s)
    tracker.remove_observer(obs)

    assert "started" in fired
    assert "ended"   in fired


def test_double_tracker_isinstance(double_tracker):
    tracker, s = double_tracker

    seen = []

    class Obs(t.observers.double.Abstract):
        def Observe(self, event):
            if isinstance(event, t.observers.double.TrackingStarted):
                seen.append("started")
            elif isinstance(event, t.observers.double.TrackingEnded):
                seen.append("ended")

    obs = Obs()
    tracker.add_observer(obs)
    _run_double(tracker, s)
    tracker.remove_observer(obs)

    assert "started" in seen
    assert "ended"   in seen
