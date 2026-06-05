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
