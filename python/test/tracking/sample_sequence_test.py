"""The SampleSequenceCollector from Python: attach, run, read the sequence.

Interface checks only -- the events, the per-run boundaries and the precision handling are
covered in C++ (generic_pseg_test.hpp, generic_cauchy_test.hpp).  Here: the collector is
reachable under bertini.endgame.observers.<flavor>, attaches to an endgame, and hands back
its buckets as Python lists that mean what they say.
"""

import numpy as np
import pytest

import bertini
from bertini import System, VariableGroup
from bertini.symbolics import Variable
from bertini import AMPTracker, DoublePrecisionTracker, Predictor, SuccessCode
from bertini.tracking import SteppingConfig, NewtonConfig, amp_config_from
from bertini.endgame import AMPPowerSeriesEndgame, FixedDoubleCauchyEndgame
from bertini.endgame import observers
from bertini.multiprec import complex_mp


@pytest.fixture
def cubic_homotopy():
    """f(x,t) = (x-1)^3*(1-t) + (x^3+1)*t: a triple root at t=0, cycle number 3."""
    x = Variable("x")
    t = Variable("t")
    s = System()
    vg = VariableGroup(); vg.append(x)
    s.add_variable_group(vg); s.add_path_variable(t)
    s.add_function((x - 1)**3 * (1 - t) + (x**3 + 1) * t)
    return s


def _times_shrink(times):
    mags = [abs(complex(t)) for t in times]
    return all(b < a for a, b in zip(mags, mags[1:]))


def test_sequence_collector_on_the_adaptive_power_series_endgame(cubic_homotopy):
    s = cubic_homotopy
    tracker = AMPTracker(s)
    tracker.setup(Predictor.HeunEuler, 1e-6, 1e5, SteppingConfig(), NewtonConfig())
    tracker.precision_setup(amp_config_from(s))
    eg = AMPPowerSeriesEndgame(tracker, complex_mp("0.1"))

    seq = observers.amp_pseg.SampleSequenceCollector()
    eg.add_observer(seq)
    code = eg.run(np.array([complex_mp("5.000000000000001e-01", "9.084258952712920e-17")]))
    assert code == SuccessCode.Success

    assert seq.num_runs() == 1
    assert seq.run_path_starts() == [0]
    assert seq.num_samples() == len(seq.path_samples()) == len(seq.path_times()) > 0
    assert _times_shrink(seq.path_times())                      # toward the target time
    assert abs(complex(seq.path_times()[0]) - 0.1) < 1e-12      # starting at the boundary
    assert len(seq.approximations()) == len(seq.approximation_errors()) == len(seq.cycle_numbers()) > 0
    assert seq.circle_samples() == [] and seq.advance_times() == []   # power series has no circle
    assert isinstance(seq.num_precision_increases(), int)

    seq.clear()
    assert seq.num_samples() == 0 and seq.num_runs() == 0 and seq.num_precision_increases() == 0


def test_sequence_collector_on_the_cauchy_endgame_keeps_circle_points_apart(cubic_homotopy):
    s = cubic_homotopy
    tracker = DoublePrecisionTracker(s)
    tracker.setup(Predictor.HeunEuler, 1e-5, 1e5, SteppingConfig(), NewtonConfig())
    eg = FixedDoubleCauchyEndgame(tracker, complex(0.1, 0))

    seq = observers.double_cauchy.SampleSequenceCollector()
    eg.add_observer(seq)
    code = eg.run(np.array([complex(5.000000000000001e-01, 9.084258952712920e-17)]))
    assert code == SuccessCode.Success

    assert seq.num_runs() == 1
    assert seq.num_samples() > 0 and _times_shrink(seq.path_times())
    assert len(seq.circle_samples()) == len(seq.circle_times()) > 0   # the loop points, kept apart
    assert len(seq.advance_times()) > 0
    root = seq.approximations()[-1]
    first = abs(complex(np.asarray(seq.path_samples()[0]).ravel()[0]) - complex(np.asarray(root).ravel()[0]))
    last = abs(complex(np.asarray(seq.path_samples()[-1]).ravel()[0]) - complex(np.asarray(root).ravel()[0]))
    assert last < first                                            # the approach is real
