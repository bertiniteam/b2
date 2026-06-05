# This file is part of Bertini 2.
#
# python/test/tracking/cauchy_endgame_test.py is free software: you can
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

"""Cauchy endgame tests across tracker types.

Mirrors C++ tests in:
  core/test/endgames/generic_cauchy_test.hpp  (included by)
  core/test/endgames/fixed_double_cauchy_test.cpp
  core/test/endgames/fixed_multiple_cauchy_test.cpp
  core/test/endgames/amp_cauchy_test.cpp

Three homotopies are used:
  linear:    f(x,t) = (x-1)*(1-t) + (x+1)*t       cycle_num=1 at x=1
  quadratic: f(x,t) = (x-1)^2*(1-t) + (x^2+1)*t   cycle_num=2 at x=1
  cubic:     f(x,t) = (x-1)^3*(1-t) + (x^3+1)*t   cycle_num=3 at x=1
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
from bertini.endgame import AMPCauchyEG, FixedDoubleCauchyEG, FixedMultipleCauchyEG

import bertini.multiprec as mp
from bertini.multiprec import Complex as mpfr_complex


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def linear_homotopy():
    """(x-1)*(1-t) + (x+1)*t — simple root at x=1, cycle number 1."""
    x = Variable("x"); t = Variable("t")
    s = System()
    vg = VariableGroup(); vg.append(x)
    s.add_variable_group(vg); s.add_path_variable(t)
    s.add_function((x - 1) * (1 - t) + (x + 1) * t)
    return s, x, t


@pytest.fixture
def quadratic_homotopy():
    """(x-1)^2*(1-t) + (x^2+1)*t — double root at x=1, cycle number 2."""
    x = Variable("x"); t = Variable("t")
    s = System()
    vg = VariableGroup(); vg.append(x)
    s.add_variable_group(vg); s.add_path_variable(t)
    s.add_function((x - 1)**2 * (1 - t) + (x**2 + 1) * t)
    return s, x, t


@pytest.fixture
def cubic_homotopy():
    """(x-1)^3*(1-t) + (x^3+1)*t — triple root at x=1, cycle number 3."""
    x = Variable("x"); t = Variable("t")
    s = System()
    vg = VariableGroup(); vg.append(x)
    s.add_variable_group(vg); s.add_path_variable(t)
    s.add_function((x - 1)**3 * (1 - t) + (x**3 + 1) * t)
    return s, x, t


# ---------------------------------------------------------------------------
# FixedDoubleCauchyEG
# ---------------------------------------------------------------------------

def test_fixed_double_cauchy_cycle_num_1(linear_homotopy):
    """Cauchy EG on linear homotopy: cycle number 1, converges to x=1.

    Mirrors generic_cauchy_test.hpp/full_test_cycle_num_1 with DoublePrecisionTracker.
    Pre-computed sample at t=0.1 taken from the C++ reference test.
    """
    s, x, t = linear_homotopy

    tracker = DoublePrecisionTracker(s)
    tracker.setup(Predictor.HeunEuler, 1e-5, 1e5, SteppingConfig(), NewtonConfig())

    eg = FixedDoubleCauchyEG(tracker, complex(0.1, 0))
    # Pre-computed from C++ test at t=0.1
    sample = np.array([complex(7.999999999999999e-01, 2.168404344971009e-19)])

    code = eg.run(sample)

    fa = eg.final_approximation()
    assert code == SuccessCode.Success
    assert eg.cycle_number() == 1
    assert abs(fa[0] - complex(1, 0)) < 1e-5


def test_fixed_double_cauchy_cycle_num_2(quadratic_homotopy):
    """Cauchy EG on quadratic homotopy: cycle number 2, converges to x=1.

    Mirrors generic_cauchy_test.hpp/full_test_cycle_num_greater_than_1.
    Pre-computed sample at t=0.1 taken from the C++ reference test.
    """
    s, x, t = quadratic_homotopy

    tracker = DoublePrecisionTracker(s)
    nc = NewtonConfig()
    nc.max_num_newton_iterations = 2
    nc.min_num_newton_iterations = 1
    tracker.setup(Predictor.HeunEuler, 1e-5, 1e5, SteppingConfig(), nc)

    eg = FixedDoubleCauchyEG(tracker, complex(0.1, 0))
    # Pre-computed from C++ test at t=0.1
    sample = np.array([complex(9.000000000000001e-01, 4.358898943540673e-01)])

    code = eg.run(sample)

    fa = eg.final_approximation()
    assert code == SuccessCode.Success
    assert abs(fa[0] - complex(1, 0)) < 1e-5
    assert eg.cycle_number() == 2


def test_fixed_double_cauchy_track_to_boundary(cubic_homotopy):
    """Track then run Cauchy EG on cubic homotopy (cycle number 3).

    Tracks from t=1 to t=0.1 first, then runs the endgame.
    """
    s, x, t = cubic_homotopy

    tracker = DoublePrecisionTracker(s)
    tracker.setup(Predictor.HeunEuler, 1e-6, 1e5, SteppingConfig(), NewtonConfig())

    # At t=1: x^3+1=0 → x=-1 is a real root
    start = np.array([complex(-1, 0)])
    bdry = np.zeros(1, dtype=complex)
    code = tracker.track_path(bdry, complex(1, 0), complex(0.1, 0), start)
    assert code == SuccessCode.Success

    eg = FixedDoubleCauchyEG(tracker, complex(0.1, 0))
    code = eg.run(bdry)

    fa = eg.final_approximation()
    assert code == SuccessCode.Success
    assert abs(fa[0] - complex(1, 0)) < 1e-5


# ---------------------------------------------------------------------------
# FixedMultipleCauchyEG
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("precision", [16, 30, 50], indirect=True)
def test_fixed_multiple_cauchy_cycle_num_1(linear_homotopy, precision):
    """Cauchy EG on linear homotopy at several mp precisions, cycle number 1.

    Mirrors fixed_multiple_cauchy_test.cpp/full_test_cycle_num_1 at precision 16/30/50.
    """
    s, x, t = linear_homotopy
    s.precision(precision)

    tracker = MultiplePrecisionTracker(s)
    tracker.setup(Predictor.HeunEuler, 1e-5, 1e5, SteppingConfig(), NewtonConfig())

    eg = FixedMultipleCauchyEG(tracker, mpfr_complex("0.1"))
    # Pre-computed sample (converted to mpfr_complex at the current precision)
    sample = np.array([mpfr_complex("7.999999999999999e-01")])

    code = eg.run(sample)

    fa = eg.final_approximation()
    assert code == SuccessCode.Success
    assert eg.cycle_number() == 1
    assert mp.abs(fa[0] - mpfr_complex(1)) < 1e-5


@pytest.mark.parametrize("precision", [16, 30, 50], indirect=True)
def test_fixed_multiple_cauchy_cycle_num_2(quadratic_homotopy, precision):
    """Cauchy EG on quadratic homotopy at several mp precisions, cycle number 2.

    Mirrors fixed_multiple_cauchy_test.cpp/full_test_cycle_num_greater_than_1.
    """
    s, x, t = quadratic_homotopy
    s.precision(precision)

    tracker = MultiplePrecisionTracker(s)
    nc = NewtonConfig()
    nc.max_num_newton_iterations = 2
    nc.min_num_newton_iterations = 1
    tracker.setup(Predictor.HeunEuler, 1e-5, 1e5, SteppingConfig(), nc)

    eg = FixedMultipleCauchyEG(tracker, mpfr_complex("0.1"))
    sample = np.array([mpfr_complex("9.000000000000001e-01", "4.358898943540673e-01")])

    code = eg.run(sample)

    fa = eg.final_approximation()
    assert code == SuccessCode.Success
    assert mp.abs(fa[0] - mpfr_complex(1)) < 1e-5
    assert eg.cycle_number() == 2


# ---------------------------------------------------------------------------
# AMPCauchyEG
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("precision", [16, 30, 50], indirect=True)
def test_amp_cauchy_cycle_num_1(linear_homotopy, precision):
    """AMPCauchyEG on linear homotopy at three ambient precisions, cycle number 1.

    Mirrors amp_cauchy_test.cpp/generic_tests_ambient_precision_*.
    """
    s, x, t = linear_homotopy

    ampconfig = amp_config_from(s)
    tracker = AMPTracker(s)
    tracker.setup(Predictor.HeunEuler, 1e-5, 1e5, SteppingConfig(), NewtonConfig())
    tracker.precision_setup(ampconfig)

    eg = AMPCauchyEG(tracker, mpfr_complex("0.1"))
    sample = np.array([mpfr_complex("7.999999999999999e-01")])

    code = eg.run(sample)

    fa = eg.final_approximation()
    assert code == SuccessCode.Success
    assert eg.cycle_number() == 1
    assert mp.abs(fa[0] - mpfr_complex(1)) < 1e-5


@pytest.mark.parametrize("precision", [16, 30, 50], indirect=True)
def test_amp_cauchy_cycle_num_2(quadratic_homotopy, precision):
    """AMPCauchyEG on quadratic homotopy: cycle number 2.

    Mirrors amp_cauchy_test.cpp/generic_tests_ambient_precision_* (cycle>1 variant).
    """
    s, x, t = quadratic_homotopy

    ampconfig = amp_config_from(s)
    tracker = AMPTracker(s)
    nc = NewtonConfig()
    nc.max_num_newton_iterations = 2
    nc.min_num_newton_iterations = 1
    tracker.setup(Predictor.HeunEuler, 1e-5, 1e5, SteppingConfig(), nc)
    tracker.precision_setup(ampconfig)

    eg = AMPCauchyEG(tracker, mpfr_complex("0.1"))
    sample = np.array([mpfr_complex("9.000000000000001e-01", "4.358898943540673e-01")])

    code = eg.run(sample)

    fa = eg.final_approximation()
    assert code == SuccessCode.Success
    assert mp.abs(fa[0] - mpfr_complex(1)) < 1e-5
    assert eg.cycle_number() == 2
