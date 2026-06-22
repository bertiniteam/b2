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
#  Spring 2018
#


__author__ = 'jcollins'

from bertini import *
from bertini.function_tree.symbol import *
from bertini.function_tree.root import *
from bertini.function_tree import *
import numpy as np
import pytest

import bertini as pb

import bertini.multiprec as mp
from bertini.multiprec import Float as mpfr_float
from bertini.multiprec import Complex as mpfr_complex


# global default precision is reset to a known baseline before every test by the
# autouse _reset_precision fixture in python/test/conftest.py; the `precision` fixture
# (used indirectly below) overrides it for the multiprecision-sensitive cases.

TOLDBL = 1e-15


@pytest.fixture
def variables():
    return Variable("x"), Variable("y"), Variable("z")


@pytest.fixture
def fg(variables):
    x, y, z = variables
    a = Float("4.897", "1.23")
    f = (x*y)
    g = (pow(x, 2)*y - a*z*x)
    return f, g


def test_system_create():
    x = Variable("x")
    y = Variable("y")
    f = (x*y)
    #
    s = System()

    vg = pb.VariableGroup()
    vg.append(x)
    vg.append(y)

    s.add_variable_group(vg)

    s.add_function(f)


def test_gather_variables_alphabetical(variables):
    x, y, z = variables
    # declared out of order, x used twice; expect distinct, sorted by name
    f1 = (pow(z, 2) + y*x)
    f2 = (x - y)
    found = gather_variables([f1, f2])
    assert [str(v) for v in found] == ['x', 'y', 'z']


def test_system_from_functions(variables):
    x, y, z = variables
    f1 = (x*y*z)
    f2 = (x + y + z)
    s = System([f1, f2])
    assert s.num_functions() == 2
    assert s.num_variable_groups() == 1
    assert s.num_variables() == 3
    assert str(s.variable_groups()[0]) == '[x,y,z]'


def test_set_variable_groups(variables):
    x, y, z = variables
    s = System([(x*y*z)])
    assert s.num_variable_groups() == 1
    s.set_variable_groups([pb.VariableGroup([x]), pb.VariableGroup([y, z])])
    assert s.num_variable_groups() == 2
    assert s.num_variables() == 3


def test_fix_variable(variables):
    x, y, z = variables
    s = System([(x + y)])
    assert s.num_variables() == 2
    assert s.fix_variable(y, complex(3.0))
    assert s.num_variables() == 1
    # x + y, with y pinned to 3, at x = 2  ->  5
    e = s.eval(np.array([complex(2.0)]))
    assert e[0] == complex(5.0)
    # fixing a variable not in the system returns False
    assert not s.fix_variable(Variable("w"), complex(1.0))


def test_system_eval_double(variables, fg):
    x, y, z = variables
    f, g = fg
    exact_real = (-32.841085, -150.5480559)
    exact_imag = (-26.66705, -258.97936865)
    #
    s = System()

    vg = pb.VariableGroup()
    vg.append(x)
    vg.append(y)
    vg.append(z)
    s.add_variable_group(vg)

    s.add_function(f)
    s.add_function(g)
    #
    v = np.array([complex(0), complex(0), complex(0)])
    v[0] = complex(3.5, 2.89); v[1] = complex(-9.32, .0765); v[2] = complex(5.4, -2.13)
    #
    e = s.eval(v)
    #
    assert np.abs(e[0].real - exact_real[0]) < TOLDBL*np.abs(exact_real[0])
    assert np.abs(e[0].imag - exact_imag[0]) < TOLDBL*np.abs(exact_imag[0])
    assert np.abs(e[1].real - exact_real[1]) < TOLDBL*np.abs(exact_real[1])
    assert np.abs(e[1].imag - exact_imag[1]) < TOLDBL*np.abs(exact_imag[1])


@pytest.mark.parametrize("precision", [30, 50, 80], indirect=True)
def test_system_eval_mp(precision):
    # tolerance scales with the working precision, replacing the old hard-coded 1e-27.
    tol = mpfr_float(10) ** (-(precision - 3))
    s = pb.parse.system('function f1, f2; variable_group x,y,z; f1 = x*y; f2 = x^2*y - z*x;')
    # variables are canonical-by-name and shared across tests, so they may carry a neighbor's
    # precision; set the system (and thus its variables) to this test's precision explicitly.
    s.precision(precision)
    exact_real = (mpfr_float('-32.841085'), mpfr_float('-62.9317230'))
    exact_imag = (mpfr_float('-26.66705'), mpfr_float('-196.39641065'))
    v = np.array((mpfr_complex('3.5', '2.89'), mpfr_complex('-9.32', '.0765'), mpfr_complex('5.4', '-2.13')))
    #
    e = s.eval(v)
    #
    assert mp.abs(e[0].real / exact_real[0] - 1) <= tol
    assert mp.abs(e[0].imag / exact_imag[0] - 1) <= tol
    assert mp.abs(e[1].real / exact_real[1] - 1) <= tol
    assert mp.abs(e[1].imag / exact_imag[1] - 1) <= tol


def test_system_Jac_double(variables, fg):
    x, y, z = variables
    f, g = fg
    exact_real = ((-9.32, 3.5, 0),
                  (-94.745870, 3.8979, -13.5848))
    exact_imag = ((.0765, 2.89, 0),
                  (-49.54549, 20.230, -18.45733))
    #
    s = System()

    vg = pb.VariableGroup()
    vg.append(x)
    vg.append(y)
    vg.append(z)
    s.add_variable_group(vg)

    s.add_function(f)
    s.add_function(g)
    #
    v = np.array([complex(0), complex(0), complex(0)])
    v[0] = complex(3.5, 2.89); v[1] = complex(-9.32, .0765); v[2] = complex(5.4, -2.13)
    #
    s.differentiate()
    e = s.eval_jacobian(v)
    #
    for r in range(2):
        for c in range(3):
            assert np.abs(e[r][c].real - exact_real[r][c]) <= TOLDBL*np.abs(exact_real[r][c])
            assert np.abs(e[r][c].imag - exact_imag[r][c]) <= TOLDBL*np.abs(exact_imag[r][c])


@pytest.mark.parametrize("precision", [30, 50, 80], indirect=True)
def test_system_Jac_mp(precision):
    tol = mpfr_float(10) ** (-(precision - 3))
    s = pb.parse.system('function f1, f2; variable_group x,y,z; f1 = x*y; f2 = x^2*y - z*x;')
    # shared canonical variables may carry a neighbor's precision; pin this test's precision.
    s.precision(precision)
    exact_real = ((mpfr_float('-9.32'), mpfr_float('3.5'), mpfr_float('0')),
                  (mpfr_float('-71.082170'), mpfr_float('3.8979'), mpfr_float('-3.5')))
    exact_imag = ((mpfr_float('.0765'), mpfr_float('2.89'), mpfr_float('0')),
                  (mpfr_float('-51.20410'), mpfr_float('20.230'), mpfr_float('-2.89')))
    v = np.array((mpfr_complex('3.5', '2.89'), mpfr_complex('-9.32', '.0765'), mpfr_complex('5.4', '-2.13')))
    #
    s.differentiate()
    e = s.eval_jacobian(v)
    #
    # the (0, 2) entry (d(x*y)/dz) is exactly zero -- compare absolutely there, relatively
    # elsewhere. (positions are hard-coded rather than testing `exact == 0`, which trips a
    # RecursionError in the mpfr_float comparison binding.)
    zero_positions = {(0, 2)}
    for r in range(2):
        for c in range(3):
            if (r, c) in zero_positions:
                assert mp.abs(e[r][c].real) <= tol
                assert mp.abs(e[r][c].imag) <= tol
            else:
                assert mp.abs(e[r][c].real / exact_real[r][c] - 1) <= tol
                assert mp.abs(e[r][c].imag / exact_imag[r][c] - 1) <= tol


def test_add_systems(variables):
    x, y, z = variables
    s1 = System(); s2 = System()
    #
    vars = VariableGroup()
    vars.append(x); vars.append(y)
    #
    s1.add_variable_group(vars)
    s1.add_function(y+1)
    s1.add_function(x*y)
    #
    s2.add_variable_group(vars)
    s2.add_function(-y-1)
    s2.add_function(-x*y)
    #
    s1 += s2
    values = np.array((2, 3))
    v = s1.eval(values)
    #
    assert v[0] == 0.0
    assert v[1] == 0.0
    #
    deg = s1.degrees()
    assert len(deg) == 2
    #
    assert deg[0] == 1
    assert deg[1] == 2


def test_homogenize_multiple_variable_groups():
    """Homogenizing a two-group system inserts one hom. variable per group.

    Mirrors system_class/system_homogenize_multiple_variable_groups.
    """
    x, y = Variable("x"), Variable("y")
    s = System()
    g1 = pb.VariableGroup(); g1.append(x)
    g2 = pb.VariableGroup(); g2.append(y)
    s.add_variable_group(g1)
    s.add_variable_group(g2)
    s.add_function(x + y - 1)
    s.add_function(x * y)
    assert not s.is_homogeneous()
    s.homogenize()
    assert s.is_homogeneous()
    # Each affine group gains a hom. variable → 2 original + 2 hom. = 4 variables.
    assert s.num_variables() == 4


def test_reorder_by_degree_decreasing(variables):
    """Reorder functions by decreasing degree.

    Mirrors system_class/system_reorder_by_degree_decreasing.
    """
    x, y, z = variables
    s = System()
    vg = pb.VariableGroup(); vg.append(x); vg.append(y)
    s.add_variable_group(vg)
    s.add_function(x + y)          # degree 1
    s.add_function(x**3 + y**3)    # degree 3
    s.add_function(x**2 - y)       # degree 2

    s.reorder_functions_by_degree_decreasing()
    degs = list(s.degrees())
    assert degs == sorted(degs, reverse=True), f"Not decreasing: {degs}"


def test_reorder_by_degree_increasing(variables):
    """Reorder functions by increasing degree.

    Mirrors system_class/system_reorder_by_degree_increasing.
    """
    x, y, z = variables
    s = System()
    vg = pb.VariableGroup(); vg.append(x); vg.append(y)
    s.add_variable_group(vg)
    s.add_function(x**3 + y**3)    # degree 3
    s.add_function(x + y)          # degree 1
    s.add_function(x**2 - y)       # degree 2

    s.reorder_functions_by_degree_increasing()
    degs = list(s.degrees())
    assert degs == sorted(degs), f"Not increasing: {degs}"


def test_eval_wrong_size_input_throws(variables):
    """Feeding the wrong number of variables raises RuntimeError.

    Mirrors system_class/eval_wrong_size_input_throws.
    """
    x, y, z = variables
    s = System()
    vg = pb.VariableGroup(); vg.append(x); vg.append(y)
    s.add_variable_group(vg)
    s.add_function(x + y)

    with pytest.raises(RuntimeError):
        s.eval(np.array([complex(1, 0)]))  # 1 value, need 2


def test_dehomogenize_one_affine_group():
    """Dehomogenizing a single-group homogenized system divides out the hom. var.

    Mirrors system_class/system_dehomogenize_FIFO_one_aff_group.
    After Homogenize(), the hom. variable is PREPENDED (FIFO), so the point
    vector is ordered [h, x, y, ...].  DehomogenizePoint returns [x/h, y/h, ...].
    """
    x, y = Variable("x"), Variable("y")
    s = System()
    vg = pb.VariableGroup(); vg.append(x); vg.append(y)
    s.add_variable_group(vg)
    s.add_function(x**2 + y - 1)
    s.homogenize()
    # Variable order after hom: [h, x, y]
    # v = [2+3i, 3+4i, 4+5i] → dehom = [v[1]/v[0], v[2]/v[0]]
    v = np.array([complex(2, 3), complex(3, 4), complex(4, 5)])
    dehom = s.dehomogenize_point(v)
    assert dehom.shape == (2,)
    assert abs(dehom[0] - v[1] / v[0]) < 1e-14
    assert abs(dehom[1] - v[2] / v[0]) < 1e-14


def test_dehomogenize_two_affine_groups():
    """Dehomogenizing a two-group system divides each group by its own hom. var.

    Mirrors system_class/system_dehomogenize_FIFO_two_aff_groups.
    After Homogenize(), each group gets its hom. var. prepended.
    Variable order: [h1, x, y, h2, z, w] for groups {x,y} and {z,w}.
    """
    x, y = Variable("x"), Variable("y")
    z, w = Variable("z"), Variable("w")
    s = System()
    g1 = pb.VariableGroup(); g1.append(x); g1.append(y)
    g2 = pb.VariableGroup(); g2.append(z); g2.append(w)
    s.add_variable_group(g1)
    s.add_variable_group(g2)
    s.add_function(x + y)
    s.add_function(z * w)
    s.add_function(x - z)
    s.add_function(y - w)
    s.homogenize()
    # Variable order: [h1, x, y, h2, z, w]
    v = np.array([complex(2,3), complex(3,4), complex(4,5),
                  complex(5,6), complex(6,7), complex(7,8)])
    dehom = s.dehomogenize_point(v)
    assert dehom.shape == (4,)
    assert abs(dehom[0] - v[1] / v[0]) < 1e-14
    assert abs(dehom[1] - v[2] / v[0]) < 1e-14
    assert abs(dehom[2] - v[4] / v[3]) < 1e-14
    assert abs(dehom[3] - v[5] / v[3]) < 1e-14


def test_mult_system_node():
    tol_d = TOLDBL
    sys = pb.parse.system('function f1, f2; variable_group x,y,z; f1 = x+2; f2 = y*y;')
    #
    sys *= Integer(2)
    #
    vals = np.array((complex(-2.43, .21), complex(4.84, -1.94), complex(-6.48, -.731)))
    sysEval = sys.eval(vals)
    #
    assert np.abs(sysEval[0].real / (-.86) - 1) <= tol_d
    assert np.abs(sysEval[0].imag / (0.42) - 1) <= tol_d
    assert np.abs(sysEval[1].real / (39.3240) - 1) <= tol_d
    assert np.abs(sysEval[1].imag / (-37.5584) - 1) <= tol_d
