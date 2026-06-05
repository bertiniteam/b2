# This file is part of Bertini 2.
#
# python/test/system_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/system_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/system_test.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
#
#   James Collins
#   West Texas A&M University
#   Spring 2016
#
#  silviana amethyst
#  UWEC
#  Spring, Summer 2018
#


__author__ = 'James Collins'

from bertini import *
from bertini.function_tree.symbol import *
from bertini.function_tree.root import *
from bertini.function_tree import *
from bertini.tracking import *
from bertini.tracking.config import *

import numpy as np
import pytest

import bertini.multiprec as mp
from bertini.multiprec import Float as mpfr_float
from bertini.multiprec import Complex as mpfr_complex


# global default precision is reset to 30 before every test by the autouse
# _reset_precision fixture in python/test/conftest.py.


@pytest.fixture
def xyzt():
    return Variable("x"), Variable("y"), Variable("z"), Variable("t")


def test_tracker_linear(xyzt):
    x, y, z, t = xyzt
    s = System()

    vars = VariableGroup()
    vars.append(y)
    s.add_function(y-t)
    s.add_path_variable(t)
    s.add_variable_group(vars)

    ampconfig = amp_config_from(s)

    tracker = AMPTracker(s)

    stepping_pref = SteppingConfig()
    newton_pref = NewtonConfig()

    tracker.setup(Predictor.Euler, 1e-5, 1e5, stepping_pref, newton_pref)
    tracker.precision_setup(ampconfig)

    t_start = mpfr_complex(1)
    t_end = mpfr_complex(0)

    y_start = np.array([mpfr_complex(1)])

    y_end = np.array(np.zeros(shape=(s.num_variables()), dtype=np.int64), dtype=mpfr_complex)

    tracker.track_path(y_end, t_start, t_end, y_start)

    assert y_end.shape == (s.num_variables(),)
    assert mp.abs(y_end[0]-mpfr_complex(0)) <= 1e-5


def test_tracker_quad(xyzt):
    x, y, z, t = xyzt
    s = System()

    vars = VariableGroup()
    vars.append(y)
    s.add_function(y-t**2)
    s.add_path_variable(t)
    s.add_variable_group(vars)

    s.precision(30)

    ampconfig = amp_config_from(s)

    tracker = AMPTracker(s)

    stepping_pref = SteppingConfig()
    newton_pref = NewtonConfig()

    tracker.setup(Predictor.Euler, 1e-5, 1e5, stepping_pref, newton_pref)
    tracker.precision_setup(ampconfig)

    t_start = mpfr_complex(1)
    t_end = mpfr_complex(-1)

    y_start = np.array([mpfr_complex(1)])

    y_end = np.array(np.zeros(shape=(s.num_variables()), dtype=np.int64), dtype=mpfr_complex)

    tracker.track_path(y_end, t_start, t_end, y_start)

    assert y_end.shape == (s.num_variables(),)
    assert mp.abs(y_end[0]-mpfr_complex(1)) <= 1e-5


def test_tracker_sqrt(xyzt):
    x, y, z, t = xyzt
    s = System()

    vars = VariableGroup()
    vars.append(y); vars.append(x)
    s.add_function(x-t)
    s.add_function(y**2 - x)
    s.add_path_variable(t)
    s.add_variable_group(vars)
    s.precision(30)
    ampconfig = amp_config_from(s)

    tracker = AMPTracker(s)

    stepping_pref = SteppingConfig()
    newton_pref = NewtonConfig()

    tracker.setup(Predictor.Euler, 1e-5, 1e5, stepping_pref, newton_pref)
    tracker.precision_setup(ampconfig)

    t_start = mpfr_complex(1)
    t_end = mpfr_complex(0)

    y_start = np.array([mpfr_complex(1), mpfr_complex(1)])

    y_end = np.array(np.zeros(shape=(s.num_variables()), dtype=np.int64), dtype=mpfr_complex)

    track_success = tracker.track_path(y_end, t_start, t_end, y_start)

    assert track_success == SuccessCode.Success
    assert y_end.shape == (s.num_variables(),)
    assert mp.abs(y_end[0]-mpfr_complex(0)) <= 1e-5
    assert mp.abs(y_end[1]-mpfr_complex(0)) <= 1e-5

    y_start = np.array([mpfr_complex(1), mpfr_complex(-1)])

    tracker.track_path(y_end, t_start, t_end, y_start)

    assert y_end.shape == (s.num_variables(),)
    assert mp.abs(y_end[0]-mpfr_complex(0)) <= 1e-5
    assert mp.abs(y_end[1]-mpfr_complex(0)) <= 1e-5

    y_start = np.array([mpfr_complex(-1), mpfr_complex(-1)])

    tracker.track_path(y_end, t_start, t_end, y_start)

    assert y_end.shape == (s.num_variables(),)
    assert mp.abs(y_end[0]-mpfr_complex(0)) <= 1e-5
    assert mp.abs(y_end[1]-mpfr_complex(0)) <= 1e-5

    y_start = np.array([mpfr_complex(-1), mpfr_complex(0, 1)])

    track_success = tracker.track_path(y_end, t_start, t_end, y_start)

    assert track_success == SuccessCode.Success
    assert y_end.shape == (s.num_variables(),)
    assert mp.abs(y_end[0]-mpfr_complex(0)) <= 1e-5
    assert mp.abs(y_end[1]-mpfr_complex(0)) <= 1e-5


def test_tracker_decic(xyzt):
    """Track y - t^10 from (t=1,y=1) to (t=-2,y=1024).

    Mirrors C++ AMP_tracker_basics/AMP_tracker_track_decic.
    """
    x, y, z, t = xyzt
    s = System()

    vars = VariableGroup()
    vars.append(y)
    s.add_function(y - t**10)
    s.add_path_variable(t)
    s.add_variable_group(vars)

    ampconfig = amp_config_from(s)
    tracker = AMPTracker(s)

    tracker.setup(Predictor.Euler, 1e-5, 1e5, SteppingConfig(), NewtonConfig())
    tracker.precision_setup(ampconfig)

    t_start = mpfr_complex(1)
    t_end = mpfr_complex(-2)

    y_start = np.array([mpfr_complex(1)])
    y_end = np.array(np.zeros(shape=(s.num_variables(),), dtype=np.int64), dtype=mpfr_complex)

    code = tracker.track_path(y_end, t_start, t_end, y_start)

    assert code == SuccessCode.Success
    assert y_end.shape == (s.num_variables(),)
    assert mp.abs(y_end[0] - mpfr_complex(1024)) <= 1e-5


def test_tracker_nonhomogeneous_prec16():
    """Two-variable nonhomogeneous system tracked from precision-16 start.

    Mirrors C++ AMP_simple_nonhomogeneous_system_trackable_initialprecision16.
    System: x^2 + (1-t)*x - 1 = 0, y^2 + (1-t)*x*y - 2 = 0
    Known solution: (x,y) ≈ (0.6180, 1.1386)
    """
    mp.default_precision(16)

    x = Variable("x")
    y = Variable("y")
    t = Variable("t")

    s = System()
    vars = VariableGroup()
    vars.append(x); vars.append(y)
    s.add_function(x**2 + (1 - t) * x - 1)
    s.add_function(y**2 + (1 - t) * x * y - 2)
    s.add_path_variable(t)
    s.add_variable_group(vars)

    ampconfig = amp_config_from(s)
    tracker = AMPTracker(s)
    tracker.setup(Predictor.Euler, 1e-5, 1e5, SteppingConfig(), NewtonConfig())
    tracker.precision_setup(ampconfig)

    t_start = mpfr_complex(1)
    t_end = mpfr_complex(0)

    start = np.array([mpfr_complex(1), mpfr_complex("1.41421356237309504880168872421")])
    end = np.array(np.zeros(s.num_variables(), dtype=np.int64), dtype=mpfr_complex)

    code = tracker.track_path(end, t_start, t_end, start)

    assert code == SuccessCode.Success
    assert end.shape == (s.num_variables(),)
    assert mp.abs(end[0] - mpfr_complex("6.180339887498949e-01")) < 1e-5
    assert mp.abs(end[1] - mpfr_complex("1.138564265110173e+00")) < 1e-5


def test_tracker_nonhomogeneous_prec30():
    """Same system as test_tracker_nonhomogeneous_prec16 but starting at precision 30.

    Mirrors C++ AMP_simple_nonhomogeneous_system_trackable_initialprecision30.
    Also verifies that precision_preservation(True) keeps default_precision at 30.
    """
    x = Variable("x")
    y = Variable("y")
    t = Variable("t")

    s = System()
    vars = VariableGroup()
    vars.append(x); vars.append(y)
    s.add_function(x**2 + (1 - t) * x - 1)
    s.add_function(y**2 + (1 - t) * x * y - 2)
    s.add_path_variable(t)
    s.add_variable_group(vars)

    ampconfig = amp_config_from(s)
    tracker = AMPTracker(s)
    tracker.setup(Predictor.Euler, 1e-5, 1e5, SteppingConfig(), NewtonConfig())
    tracker.precision_setup(ampconfig)
    tracker.precision_preservation(True)

    t_start = mpfr_complex(1)
    t_end = mpfr_complex(0)

    start = np.array([mpfr_complex(1), mpfr_complex("1.414")])
    end = np.array(np.zeros(s.num_variables(), dtype=np.int64), dtype=mpfr_complex)

    code = tracker.track_path(end, t_start, t_end, start)

    assert mp.default_precision() == 30
    assert code == SuccessCode.Success
    assert end.shape == (s.num_variables(),)
    assert mp.abs(end[0] - mpfr_complex("6.180339887498949e-01")) < 1e-5
    assert mp.abs(end[1] - mpfr_complex("1.138564265110173e+00")) < 1e-5


def test_tracker_nonhomogeneous_prec100():
    """Same system tracked from precision-100 start.

    Mirrors C++ AMP_simple_nonhomogeneous_system_trackable_initialprecision100.
    """
    mp.default_precision(100)

    x = Variable("x")
    y = Variable("y")
    t = Variable("t")

    s = System()
    vars = VariableGroup()
    vars.append(x); vars.append(y)
    s.add_function(x**2 + (1 - t) * x - 1)
    s.add_function(y**2 + (1 - t) * x * y - 2)
    s.add_path_variable(t)
    s.add_variable_group(vars)

    ampconfig = amp_config_from(s)
    tracker = AMPTracker(s)
    tracker.setup(Predictor.Euler, 1e-5, 1e5, SteppingConfig(), NewtonConfig())
    tracker.precision_setup(ampconfig)

    t_start = mpfr_complex(1)
    t_end = mpfr_complex(0)

    start = np.array([mpfr_complex(1), mpfr_complex("1.41421356237309504880168872421")])
    end = np.array(np.zeros(s.num_variables(), dtype=np.int64), dtype=mpfr_complex)

    code = tracker.track_path(end, t_start, t_end, start)

    assert code == SuccessCode.Success
    assert end.shape == (s.num_variables(),)
    assert mp.abs(end[0] - mpfr_complex("6.180339887498949e-01")) < 1e-5
    assert mp.abs(end[1] - mpfr_complex("1.138564265110173e+00")) < 1e-5


def test_tracker_fails_singularity_on_path(xyzt):
    """Path with a singularity at t=0.5 should not succeed.

    Mirrors C++ AMP_tracker_basics/AMP_tracker_fails_with_singularity_on_path.
    System: x^2 - (2t-1) = 0, y^2 - (2t-1) = 0
    The discriminant vanishes at t=0.5, making the path singular there.
    """
    x, y, z, t = xyzt

    s_val = -1 * (1 - t) + 1 * t  # = 2t - 1; zero at t=0.5

    s = System()
    vars = VariableGroup()
    vars.append(x); vars.append(y)
    s.add_function(x**2 - s_val)
    s.add_function(y**2 - s_val)
    s.add_path_variable(t)
    s.add_variable_group(vars)

    ampconfig = amp_config_from(s)
    tracker = AMPTracker(s)
    tracker.setup(Predictor.Euler, 1e-5, 1e5, SteppingConfig(), NewtonConfig())
    tracker.precision_preservation(True)
    tracker.precision_setup(ampconfig)

    t_start = mpfr_complex(1)
    t_end = mpfr_complex(0)

    start = np.array([mpfr_complex(1), mpfr_complex(1)])
    end = np.empty(s.num_variables(), dtype=mpfr_complex)

    code = tracker.track_path(end, t_start, t_end, start)

    assert code != SuccessCode.Success
    assert mp.default_precision() == 30


def test_tracker_singular_start(xyzt):
    x, y, z, t = xyzt
    s = System()

    vars = VariableGroup()
    vars.append(y); vars.append(x)
    s.add_function(x**2 + (1-t)*x)
    s.add_function(y**2 + (1-t)*y)
    s.add_path_variable(t)
    s.add_variable_group(vars)
    s.precision(30)
    ampconfig = amp_config_from(s)

    tracker = AMPTracker(s)

    stepping_pref = SteppingConfig()
    newton_pref = NewtonConfig()

    tracker.setup(Predictor.Euler, 1e-5, 1e5, stepping_pref, newton_pref)
    tracker.precision_setup(ampconfig)

    t_start = mpfr_complex(1)
    t_end = mpfr_complex(0)

    y_start = np.array([mpfr_complex(0), mpfr_complex(0)])

    y_end = np.empty(shape=(s.num_variables(),), dtype=mpfr_complex)

    track_success = tracker.track_path(y_end, t_start, t_end, y_start)

    assert track_success == SuccessCode.SingularStartPoint
    assert y_end.shape == (s.num_variables(),)
