# This file is part of Bertini 2.
#
# python/test/tracking/powerseries_endgame_test.py is free software: you can
# redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This file is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this file.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

# individual authors of this file include:
#   silviana amethyst

"""Power-series endgame tests across tracker types.

Mirrors C++ tests in:
  core/test/endgames/generic_pseg_test.hpp  (included by)
  core/test/endgames/fixed_double_powerseries_test.cpp
  core/test/endgames/fixed_multiple_powerseries_test.cpp
  core/test/endgames/amp_powerseries_test.cpp

The homotopy used throughout is  f(x,t) = (x-1)^3*(1-t) + (x^3+1)*t.
At t=1: x^3+1=0 (three simple roots).
At t=0: (x-1)^3=0 (triple root x=1, cycle number = 3).
"""

import numpy as np
import pytest

import bertini
from bertini import System, VariableGroup
from bertini.function_tree.symbol import Variable
from bertini.tracking import (
    AMPTracker, DoublePrecisionTracker, MultiplePrecisionTracker,
    Predictor, SuccessCode,
)
from bertini.tracking.config import SteppingConfig, NewtonConfig, amp_config_from
from bertini.endgame import AMPPSEG, FixedDoublePSEG, FixedMultiplePSEG

import bertini.multiprec as mp
from bertini.multiprec import Complex as mpfr_complex


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def cubic_homotopy():
    """f(x,t) = (x-1)^3*(1-t) + (x^3+1)*t  with path variable t."""
    x = Variable("x")
    t = Variable("t")
    s = System()
    vg = VariableGroup(); vg.append(x)
    s.add_variable_group(vg); s.add_path_variable(t)
    s.add_function((x - 1)**3 * (1 - t) + (x**3 + 1) * t)
    return s, x, t


@pytest.fixture
def quadratic_homotopy():
    """f(x,t) = (x-1)^2*(1-t) + (x^2-1)*t  (cycle number 1 at x=1)."""
    x = Variable("x")
    t = Variable("t")
    s = System()
    vg = VariableGroup(); vg.append(x)
    s.add_variable_group(vg); s.add_path_variable(t)
    s.add_function((x - 1)**2 * (1 - t) + (x**2 - 1) * t)
    return s, x, t


# ---------------------------------------------------------------------------
# FixedDoublePSEG
# ---------------------------------------------------------------------------

def test_fixed_double_pseg_full_run(cubic_homotopy):
    """Track to endgame boundary then run PSEG; should converge to x=1.

    Mirrors generic_pseg_test.hpp/pseg_full_run with DoublePrecisionTracker.
    Start point (0.5, ~0i) at t=0.1 is pre-tracked from the C++ reference.
    """
    s, x, t = cubic_homotopy

    tracker = DoublePrecisionTracker(s)
    tracker.setup(Predictor.HeunEuler, 1e-6, 1e5, SteppingConfig(), NewtonConfig())

    current_time = complex(0.1, 0)
    eg = FixedDoublePSEG(tracker, current_time)

    # Pre-computed boundary point from the C++ test
    current_space = np.array([complex(5.000000000000001e-01, 9.084258952712920e-17)])

    code = eg.run(current_space)

    fa = eg.final_approximation()
    assert code == SuccessCode.Success
    assert abs(fa[0] - complex(1, 0)) < 1e-11


def test_fixed_double_pseg_full_run_track_to_boundary(cubic_homotopy):
    """Track from t=1 to t=0.1 with DoublePrecisionTracker, then run PSEG.

    Mirrors generic_pseg_test.hpp/pseg_full_run but derives boundary point
    from an actual tracking step rather than using a hard-coded value.
    """
    s, x, t = cubic_homotopy

    tracker = DoublePrecisionTracker(s)
    tracker.setup(Predictor.HeunEuler, 1e-6, 1e5, SteppingConfig(), NewtonConfig())

    # At t=1: x^3+1=0 → x=-1 is a real root
    start = np.array([complex(-1, 0)])
    bdry = np.zeros(1, dtype=complex)
    code = tracker.track_path(bdry, complex(1, 0), complex(0.1, 0), start)
    assert code == SuccessCode.Success

    eg = FixedDoublePSEG(tracker, complex(0.1, 0))
    code = eg.run(bdry)

    fa = eg.final_approximation()
    assert code == SuccessCode.Success
    assert abs(fa[0] - complex(1, 0)) < 1e-11


def test_fixed_double_pseg_cycle_num_1(quadratic_homotopy):
    """PSEG with cycle number 1 (simple root).

    Mirrors generic_pseg_test.hpp/full_run_cycle_num_2 but the quadratic
    system has cycle number 1 for the x=1 component.
    """
    s, x, t = quadratic_homotopy

    tracker = DoublePrecisionTracker(s)
    tracker.setup(Predictor.HeunEuler, 1e-6, 1e5, SteppingConfig(), NewtonConfig())

    start = np.array([complex(1, 0)])
    bdry = np.zeros(1, dtype=complex)
    code = tracker.track_path(bdry, complex(1, 0), complex(0.1, 0), start)
    assert code == SuccessCode.Success

    eg = FixedDoublePSEG(tracker, complex(0.1, 0))
    code = eg.run(bdry)

    fa = eg.final_approximation()
    assert code == SuccessCode.Success
    assert eg.cycle_number() == 1
    assert abs(fa[0] - complex(1, 0)) < 1e-11


def test_fixed_double_pseg_multiple_variables():
    """PSEG on a 2-variable decoupled system.

    Mirrors generic_pseg_test.hpp/full_run_multiple_variables.
    f1 = (x-1)^3*(1-t) + (x^3+1)*t,  f2 = (y-1)^2*(1-t) + (y^2+1)*t
    Solution: (x,y) = (1, 1).
    """
    x = Variable("x"); y = Variable("y"); t = Variable("t")
    s = System()
    vg = VariableGroup(); vg.append(x); vg.append(y)
    s.add_variable_group(vg); s.add_path_variable(t)
    s.add_function((x - 1)**3 * (1 - t) + (x**3 + 1) * t)
    s.add_function((y - 1)**2 * (1 - t) + (y**2 + 1) * t)

    tracker = DoublePrecisionTracker(s)
    tracker.setup(Predictor.HeunEuler, 1e-6, 1e5, SteppingConfig(), NewtonConfig())

    start = np.array([complex(-1, 0), complex(0, 1)])
    bdry = np.zeros(2, dtype=complex)
    code = tracker.track_path(bdry, complex(1, 0), complex(0.1, 0), start)
    assert code == SuccessCode.Success

    eg = FixedDoublePSEG(tracker, complex(0.1, 0))
    code = eg.run(bdry)

    fa = eg.final_approximation()
    assert code == SuccessCode.Success
    assert abs(fa[0] - complex(1, 0)) < 1e-10
    assert abs(fa[1] - complex(1, 0)) < 1e-10


# ---------------------------------------------------------------------------
# FixedMultiplePSEG
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("precision", [30, 50], indirect=True)
def test_fixed_multiple_pseg_full_run(cubic_homotopy, precision):
    """FixedMultiplePSEG converges to x=1 at several precisions.

    Mirrors fixed_multiple_powerseries_test.cpp/pseg_full_run.
    """
    s, x, t = cubic_homotopy
    s.precision(precision)

    tracker = MultiplePrecisionTracker(s)
    tracker.setup(Predictor.HeunEuler, 1e-6, 1e5, SteppingConfig(), NewtonConfig())

    start = np.array([mpfr_complex(-1)])
    bdry = np.array(np.zeros(1, dtype=np.int64), dtype=mpfr_complex)
    code = tracker.track_path(bdry, mpfr_complex(1), mpfr_complex("0.1"), start)
    assert code == SuccessCode.Success

    eg = FixedMultiplePSEG(tracker, mpfr_complex("0.1"))
    code = eg.run(bdry)

    fa = eg.final_approximation()
    assert code == SuccessCode.Success
    # EndgameConfig.final_tolerance defaults to 1e-11; actual accuracy is
    # bounded by that, not by the floating-point precision of the numbers.
    assert mp.abs(fa[0] - mpfr_complex(1)) < 1e-11


# ---------------------------------------------------------------------------
# AMPPSEG
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("precision", [16, 30, 50], indirect=True)
def test_amp_pseg_full_run(cubic_homotopy, precision):
    """AMPPSEG converges to x=1 starting from three ambient precisions.

    Mirrors amp_powerseries_test.cpp/generic_tests_ambient_precision_*.
    """
    s, x, t = cubic_homotopy

    ampconfig = amp_config_from(s)
    tracker = AMPTracker(s)
    tracker.setup(Predictor.HeunEuler, 1e-6, 1e5, SteppingConfig(), NewtonConfig())
    tracker.precision_setup(ampconfig)

    start = np.array([mpfr_complex(-1)])
    bdry = np.array(np.zeros(1, dtype=np.int64), dtype=mpfr_complex)
    code = tracker.track_path(bdry, mpfr_complex(1), mpfr_complex("0.1"), start)
    assert code == SuccessCode.Success

    eg = AMPPSEG(tracker, mpfr_complex("0.1"))
    code = eg.run(bdry)

    fa = eg.final_approximation()
    assert code == SuccessCode.Success
    assert mp.abs(fa[0] - mpfr_complex(1)) < 1e-11
