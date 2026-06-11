# This file is part of Bertini 2.
#
# python/test/classes/hessian_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/classes/hessian_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/classes/hessian_test.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
#
#  silviana amethyst
#  University of Wisconsin - Eau Claire
#  Spring 2026
#

"""Tests for computing Hessians (second derivatives) of functions and systems.

There is no first-class Hessian API; the supported route is repeated
explicit-variable differentiation, f.differentiate(v).differentiate(w), which
yields a plain evaluable tree (no Differential leaves, unlike the no-argument
Jacobian form).  These tests pin down that route: exact values at integer
points, mixed-partial symmetry at complex points, transcendental entries
against independently computed values, and assembly of the full Hessian
tensor of a System.
"""

import numpy as np
import pytest

import bertini as pb
import bertini.multiprec as mp
from bertini.function_tree.symbol import Variable, Rational
from bertini.function_tree import sin

from bertini.multiprec import Float as mpfr_float
from bertini.multiprec import Complex as mpfr_complex


# global default precision is reset to 30 before every test by the autouse
# _reset_precision fixture in python/test/conftest.py.

TOL_D = 1e-14


@pytest.fixture
def tol_mp():
    return mpfr_float("1e-27")


def hessian(f, variables):
    """The matrix of second-derivative trees of f, H[i][j] = d2f/(dvi dvj)."""
    return [[f.differentiate(vi).differentiate(vj) for vj in variables]
            for vi in variables]


def set_value(v, re, im):
    # each variable carries both a double and a multiprecision current value
    # (eval_d uses the former, eval_mp the latter), so set both.
    v.set_current_value(complex(re, im))
    v.set_current_value(mpfr_complex(str(re), str(im)))


@pytest.fixture
def xy_integer_point():
    x, y = Variable('x'), Variable('y')
    set_value(x, 2, 0)
    set_value(y, 3, 0)
    return x, y


@pytest.fixture
def xy_complex_point():
    x, y = Variable('x'), Variable('y')
    set_value(x, -2.43, .21)
    set_value(y, 4.84, -1.94)
    return x, y


# --- exact entries at an integer point (values exactly representable) ---

def test_polynomial_hessian_entries(xy_integer_point, tol_mp):
    x, y = xy_integer_point
    f = x**3 * y + x * y**2
    H = hessian(f, [x, y])
    # at (2,3): fxx = 6xy = 36, fxy = fyx = 3x^2 + 2y = 18, fyy = 2x = 4
    expected = [[36, 18], [18, 4]]
    for i in range(2):
        for j in range(2):
            d = H[i][j].eval_d()
            assert abs(d.real - expected[i][j]) <= TOL_D
            assert abs(d.imag) <= TOL_D
            m = H[i][j].eval_mp()
            assert mp.abs(m.real - mpfr_float(expected[i][j])) <= tol_mp
            assert mp.abs(m.imag) <= tol_mp


def test_hessian_of_absent_variable_is_zero(xy_integer_point, tol_mp):
    x, y = xy_integer_point
    z = Variable('z')
    set_value(z, -6.48, -.731)
    f = x**2 * y
    assert abs(f.differentiate(z).differentiate(z).eval_d()) <= TOL_D
    assert abs(f.differentiate(x).differentiate(z).eval_d()) <= TOL_D
    assert mp.abs(f.differentiate(z).differentiate(z).eval_mp()) <= tol_mp


def test_hessian_degrees_exact(xy_integer_point):
    # differentiation emits already-simplified trees, so degree() on second
    # derivatives is exact -- this is algebra, not an upper bound.
    x, y = xy_integer_point
    f = x**3 * y  # total degree 4

    fxx = f.differentiate(x).differentiate(x)  # 6xy
    assert fxx.degree() == 2
    assert fxx.degree(x) == 1
    assert fxx.degree(y) == 1
    assert fxx.is_polynomial()

    fxy = f.differentiate(x).differentiate(y)  # 3x^2
    assert fxy.degree() == 2
    assert fxy.degree(x) == 2
    assert fxy.degree(y) == 0

    fyy = f.differentiate(y).differentiate(y)  # 0
    assert fyy.degree() == 0


# --- mixed partials commute (the trees differ; the values must not) ---

def test_mixed_partial_symmetry_complex_point(xy_complex_point, tol_mp):
    x, y = xy_complex_point
    f = x**3 * y + x * y**2 - Rational('1/3') * x * y + sin(x * y)
    fxy = f.differentiate(x).differentiate(y)
    fyx = f.differentiate(y).differentiate(x)
    assert str(fxy) != str(fyx)  # genuinely different trees...
    assert abs(fxy.eval_d() - fyx.eval_d()) <= 1e-12  # ...same value
    assert mp.abs(fxy.eval_mp() - fyx.eval_mp()) <= mpfr_float("1e-25")


# --- transcendental entries, against independently computed values ---

def test_transcendental_hessian(xy_complex_point, tol_mp):
    x, y = xy_complex_point
    f = sin(x * y)
    H = hessian(f, [x, y])

    # doubles: independent expected values via numpy
    x0, y0 = complex(-2.43, .21), complex(4.84, -1.94)
    expected_d = [
        [-y0**2 * np.sin(x0 * y0), np.cos(x0 * y0) - x0 * y0 * np.sin(x0 * y0)],
        [np.cos(x0 * y0) - x0 * y0 * np.sin(x0 * y0), -x0**2 * np.sin(x0 * y0)],
    ]
    for i in range(2):
        for j in range(2):
            got = H[i][j].eval_d()
            assert abs(got - expected_d[i][j]) / abs(expected_d[i][j]) <= 1e-12

    # multiprecision: independent expected values via the multiprec library
    x0m, y0m = mpfr_complex("-2.43", ".21"), mpfr_complex("4.84", "-1.94")
    s, c = mp.sin(x0m * y0m), mp.cos(x0m * y0m)
    expected_m = [
        [-y0m * y0m * s, c - x0m * y0m * s],
        [c - x0m * y0m * s, -x0m * x0m * s],
    ]
    for i in range(2):
        for j in range(2):
            got = H[i][j].eval_mp()
            assert mp.abs(got - expected_m[i][j]) / mp.abs(expected_m[i][j]) <= tol_mp


# --- the Hessian tensor of a System ---

def test_system_hessian_tensor(xy_integer_point, tol_mp):
    x, y = xy_integer_point
    sys = pb.System()
    sys.add_function(x**2 * y + y**3)
    sys.add_function(x * y)
    sys.add_variable_group(pb.VariableGroup([x, y]))

    variables = [x, y]
    tensor = [hessian(sys.function(i), variables) for i in range(2)]

    # at (2,3):
    # f1 = x^2 y + y^3: fxx = 2y = 6, fxy = 2x = 4, fyy = 6y = 18
    # f2 = x y:         fxx = 0,      fxy = 1,      fyy = 0
    expected = [
        [[6, 4], [4, 18]],
        [[0, 1], [1, 0]],
    ]
    for i in range(2):
        for j in range(2):
            for k in range(2):
                d = tensor[i][j][k].eval_d()
                assert abs(d.real - expected[i][j][k]) <= TOL_D
                assert abs(d.imag) <= TOL_D
                m = tensor[i][j][k].eval_mp()
                assert mp.abs(m.real - mpfr_float(expected[i][j][k])) <= tol_mp
