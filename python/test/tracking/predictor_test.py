# This file is part of Bertini 2.
#
# python/test/tracking/predictor_test.py is free software: you can redistribute
# it and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
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

"""Tests that every exported Predictor enum value produces successful tracking.

Mirrors the intent of C++ tests in:
  core/test/tracking_basics/euler_test.cpp
  core/test/tracking_basics/heun_test.cpp
  core/test/tracking_basics/higher_predictor_test.cpp

The C++ tests exercise internal ExplicitRKPredictor::Predict directly (not
exported to Python).  Here we verify end-to-end tracking with each predictor
on the same circle-line homotopy used in the C++ suite.
"""

import numpy as np
import pytest

import bertini
from bertini import System, VariableGroup
from bertini.symbolics import Variable
from bertini.symbolics import *
from bertini import AMPTracker, Predictor, SuccessCode
from bertini.tracking import SteppingConfig, NewtonConfig, amp_config_from

import bertini.multiprec as mp
from bertini.multiprec import complex_mp as mpfr_complex


AMP_PREDICTORS = [
    Predictor.Euler,
    Predictor.HeunEuler,
    Predictor.RKCashKarp45,
    Predictor.RKDormandPrince56,
    Predictor.RKVerner67,
    Predictor.RKF45,
]

AMP_PREDICTOR_IDS = [p.name for p in AMP_PREDICTORS]

# RK4 works with the fixed-precision tracker but not AMPTracker on this homotopy.
DOUBLE_ONLY_PREDICTORS = [
    Predictor.RK4,
]


@pytest.fixture
def circle_line_homotopy():
    """Circle-line homotopy:
        f1 = t*(x^2 - 1) + (1-t)*(x^2 + y^2 - 4)
        f2 = t*(y - 1)   + (1-t)*(2x + 5y)
    Start: (x,y,t) = (2, 0, 1)  →  end at t=0.
    At t=1 the system is x^2=1, y=1, so (1,1) and (-1,1) are solutions.
    We track from (x,y)=(2,0) at t=1 to t=0.
    """
    x = Variable("x")
    y = Variable("y")
    t = Variable("t")

    s = System()
    vars = VariableGroup()
    vars.append(x)
    vars.append(y)
    s.add_variable_group(vars)
    s.add_path_variable(t)
    s.add_function(t * (x**2 - 1) + (1 - t) * (x**2 + y**2 - 4))
    s.add_function(t * (y - 1) + (1 - t) * (2 * x + 5 * y))

    return s, x, y, t


@pytest.mark.parametrize("predictor", AMP_PREDICTORS, ids=AMP_PREDICTOR_IDS)
def test_circle_line_amp_predictors(circle_line_homotopy, predictor):
    """Every embedded predictor should track the circle-line path to success with AMPTracker.

    Mirrors: euler_test, higher_predictor_test (circle_line_* cases).
    """
    s, x, y, t = circle_line_homotopy

    ampconfig = amp_config_from(s)
    tracker = AMPTracker(s)
    tracker.setup(predictor, 1e-5, 1e5, SteppingConfig(), NewtonConfig())
    tracker.precision_setup(ampconfig)

    t_start = mpfr_complex(1)
    t_end = mpfr_complex(0)
    # At t=1 the system is x^2-1=0, y-1=0, so (1,1) is a valid start point.
    start = np.array([mpfr_complex(1), mpfr_complex(1)])
    end = np.array(np.zeros(s.num_variables(), dtype=np.int64), dtype=mpfr_complex)

    code = tracker.track_path(end, t_start, t_end, start)

    assert code == SuccessCode.Success
    assert end.shape == (s.num_variables(),)


@pytest.mark.parametrize("predictor", DOUBLE_ONLY_PREDICTORS,
                         ids=[p.name for p in DOUBLE_ONLY_PREDICTORS])
def test_circle_line_double_only_predictors(circle_line_homotopy, predictor):
    """RK4 and RKF45 work with DoublePrecisionTracker but not AMPTracker.

    Mirrors: higher_predictor_test/circle_line_RK4_double (and RKF45 variant).
    """
    from bertini import DoublePrecisionTracker

    s, x, y, t = circle_line_homotopy

    tracker = DoublePrecisionTracker(s)
    tracker.setup(predictor, 1e-5, 1e5, SteppingConfig(), NewtonConfig())

    t_start = complex(1, 0)
    t_end = complex(0, 0)
    start = np.array([complex(1, 0), complex(1, 0)])
    end = np.zeros(s.num_variables(), dtype=complex)

    code = tracker.track_path(end, t_start, t_end, start)

    assert code == SuccessCode.Success
    assert end.shape == (s.num_variables(),)
