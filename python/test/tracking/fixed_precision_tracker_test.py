# This file is part of Bertini 2.
#
# python/test/tracking/fixed_precision_tracker_test.py is free software: you
# can redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.
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

"""Tests for DoublePrecisionTracker and MultiplePrecisionTracker.

Mirrors C++ tests in core/test/tracking_basics/fixed_precision_tracker_test.cpp:
  - fixed_precision_tracker_basics/double_tracker_track_linear
  - fixed_precision_tracker_basics/multiple_100_tracker_track_linear
"""

import numpy as np
import pytest

import bertini
from bertini import System, VariableGroup
from bertini.function_tree.symbol import Variable
from bertini.tracking import DoublePrecisionTracker, MultiplePrecisionTracker, Predictor, SuccessCode
from bertini.tracking import SteppingConfig, NewtonConfig, amp_config_from

import bertini.multiprec as mp
from bertini.multiprec import Complex as mpfr_complex


@pytest.fixture
def linear_system_y_minus_t():
    """y - t: trivial single-variable linear homotopy; solution at t=0 is y=0."""
    y = Variable("y")
    t = Variable("t")

    s = System()
    vars = VariableGroup()
    vars.append(y)
    s.add_function(y - t)
    s.add_path_variable(t)
    s.add_variable_group(vars)
    return s, y, t


def test_double_tracker_linear(linear_system_y_minus_t):
    s, y, t = linear_system_y_minus_t

    tracker = DoublePrecisionTracker(s)
    tracker.setup(Predictor.Euler, 1e-5, 1e5, SteppingConfig(), NewtonConfig())

    t_start = complex(1, 0)
    t_end = complex(0, 0)

    y_start = np.array([complex(1, 0)])
    y_end = np.zeros(s.num_variables(), dtype=complex)

    code = tracker.track_path(y_end, t_start, t_end, y_start)

    assert code == SuccessCode.Success
    assert y_end.shape == (s.num_variables(),)
    assert abs(y_end[0] - complex(0, 0)) < 1e-5


def test_multiple_precision_tracker_linear(linear_system_y_minus_t):
    s, y, t = linear_system_y_minus_t

    bertini.default_precision(100)
    s.precision(100)

    tracker = MultiplePrecisionTracker(s)
    tracker.setup(Predictor.Euler, 1e-5, 1e5, SteppingConfig(), NewtonConfig())

    t_start = mpfr_complex(1)
    t_end = mpfr_complex(0)

    y_start = np.array([mpfr_complex(1)])
    y_end = np.array(np.zeros(s.num_variables(), dtype=np.int64), dtype=mpfr_complex)

    code = tracker.track_path(y_end, t_start, t_end, y_start)

    assert code == SuccessCode.Success
    assert y_end.shape == (s.num_variables(),)
    assert mp.abs(y_end[0] - mpfr_complex(0)) < 1e-5
