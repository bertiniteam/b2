# This file is part of Bertini 2.
#
# python/test/function_tree_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/function_tree_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/function_tree_test.py.  If not, see <http://www.gnu.org/licenses/>.
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


from bertini import *
from bertini.function_tree.symbol import *
from bertini.function_tree.root import *
from bertini.function_tree import *
import numpy as np
import pytest

import bertini.multiprec as mp
from bertini.multiprec import Float as mpfr_float
from bertini.multiprec import Complex as mpfr_complex


# global default precision is reset to 30 before every test by the autouse
# _reset_precision fixture in python/test/conftest.py. These cases mix double and
# multiprecision evaluation with ~47-digit transcendental constants, so they stay at the
# baseline precision rather than parametrizing.


# ---------------------------------------------------------------------------- symbols

@pytest.fixture
def sym():
    """double and multiprecision scalar test values + tolerances (SymbolTest)."""
    x_d = complex(-2.43, .21)
    y_d = complex(4.84, -1.94)
    z_d = complex(-6.48, -.731)
    p_d = complex(-.321, -.72)
    tol_d = float(1e-15)
    #
    x_mp = mpfr_complex("-2.43", ".21")
    y_mp = mpfr_complex("4.84", "-1.94")
    z_mp = mpfr_complex("-6.48", "-.731")
    p_mp = mpfr_complex("-.321", "-.72")
    tol_mp = mpfr_float("1e-27")
    return x_d, y_d, z_d, p_d, tol_d, x_mp, y_mp, z_mp, p_mp, tol_mp


def test_Float_construct(sym):
    x_d, y_d, z_d, p_d, tol_d, x_mp, y_mp, z_mp, p_mp, tol_mp = sym
    x = Complex("4.3", "-9e-3")
    x = Complex(y_mp)
    x = Complex(mpfr_float("9.3"), mpfr_float("-3"))
    x = Complex("9.2", "-43.2e2")


def test_Float_funcs(sym):
    x_d, y_d, z_d, p_d, tol_d, x_mp, y_mp, z_mp, p_mp, tol_mp = sym
    x = Complex(x_mp); y = Complex(y_mp); z = Complex(z_mp)
    #
    assert x.degree() == 0
    y.differentiate()                 # a constant differentiates without error
    assert y.is_homogeneous()
    assert y.is_polynomial()


def test_Variable_construct():
    x = Variable("x")


def test_Variable_unicode_name():
    # Names are arbitrary UTF-8; a Unicode letter round-trips through str().
    x = Variable("Ω")  # Ω
    assert str(x) == "Ω"


def test_Variable_cjk_name():
    x = Variable("中")  # 中
    assert str(x) == "中"


def test_variables_count():
    v = variables('x', 3)
    assert len(v) == 3
    assert [str(z) for z in v] == ['x0', 'x1', 'x2']


def test_variables_iterable_indices():
    assert [str(z) for z in variables('x', range(2, 5))] == ['x2', 'x3', 'x4']
    assert [str(z) for z in variables('y', [0, 2])] == ['y0', 'y2']


def test_variables_custom_format():
    assert [str(z) for z in variables('x', 2, fmt='{base}_{index}')] == ['x_0', 'x_1']


def test_variables_compose():
    v = variables('x', 2)
    g = VariableGroup(v)
    assert len(g) == 2
    # generated entries are real Variable nodes usable in expressions
    assert (v[0]**2).degree() == 2


def test_Variable_funcs():
    x = Variable("x"); y = Variable("y")
    #
    assert x.degree() == 1
    assert x.degree(x) == 1
    assert x.degree(y) == 0
    d = y.differentiate()
    assert y.is_homogeneous()
    assert y.is_polynomial()


def test_special_constants_construct():
    # Construction smoke for the special-number nodes; their values are verified through the SLP
    # (see eval_expression_test.py, which evaluates Pi).
    Pi(); make_pi(); E(); make_e(); make_i()


# --------------------------------------------------------------------------- operators

@pytest.fixture
def op():
    """Variables and Complex constants for the structural operator tests (degree / homogeneity /
    polynomiality / homogenization).  Operator value-correctness is verified through the SLP in
    eval_expression_test.py, not by evaluating nodes here."""
    x = Variable("x")
    y = Variable("y")
    z = Variable("z")
    p = Variable("p")
    a = Complex(mpfr_complex("3.12", ".612"))
    b = Complex(mpfr_complex("-.823", "2.62"))
    tol_d = float(9e-14)
    tol_mp = mpfr_float("1e-27")
    return x, y, z, p, a, b, tol_d, tol_mp


def test_Operator_degree(op):
    x, y, z, p, a, b, tol_d, tol_mp = op
    #
    f = (x**3-9*y*y*z*pow(x, 4))
    w = VariableGroup(); w.append(x); w.append(y)
    assert f.degree() == 7
    assert f.degree(x) == 4
    assert f.degree(y) == 2
    assert f.degree(z) == 1
    assert f.degree(w) == 6
    w = VariableGroup(); w.append(y); w.append(z)
    assert f.degree(w) == 3


def test_Operator_ishom(op):
    x, y, z, p, a, b, tol_d, tol_mp = op
    #
    f = (y**2 - 9*y*z + 3*x*x - y*82*x - pow(z, 2))
    assert f.is_homogeneous()
    f = (y**2 - 9*y*z + 3*x*x - y*82*x - pow(z, 4))
    assert not f.is_homogeneous()
    f = (y**2 - 9*y*z + 3*x*x - 5)
    assert not f.is_homogeneous()


def test_Operator_ispoly(op):
    x, y, z, p, a, b, tol_d, tol_mp = op
    #
    f = (y**2 - 9*y*z + 3*x*x - y*82*x - pow(z, 2))
    assert f.is_polynomial()
    #
    f = (x**2*y - 9 + sin(z))
    assert not f.is_polynomial()


def test_Homogenize(op):
    x, y, z, p, a, b, tol_d, tol_mp = op
    #
    h = Variable("h")
    f = (x**2 + y**2 + z**2 - 1)
    #
    vars = VariableGroup()
    vars.append(y); vars.append(x); vars.append(z)
    #
    g = f.homogenized(vars, h)   # functional: returns a homogenized copy
    assert g.degree(h) == 2
    assert g.is_homogeneous()
    #
    assert not g.is_homogeneous(x)
    assert not g.is_homogeneous(y)
    assert not g.is_homogeneous(z)
    assert not g.is_homogeneous(h)
    #
    vars.append(h)
    assert g.is_homogeneous(vars)
    # the original is untouched (non-mutating homogenization)
    assert not f.is_homogeneous()


# ---------------------------------------------------------------------------- homogenization

def test_variable_is_homogeneous():
    """A single Variable node is degree-1 homogeneous.

    Mirrors homogenization_test.cpp/no_homogenization_needed_x.
    """
    x = Variable("x")
    assert x.is_homogeneous()


def test_constant_is_homogeneous():
    """A constant node (degree 0) is homogeneous.

    Mirrors homogenization_test.cpp/is_homogeneous_* cases.
    """
    assert Integer(0).is_homogeneous()
    assert Integer(1).is_homogeneous()
    # sin(0) = 0 is a constant
    assert sin(Integer(0)).is_homogeneous()
    assert cos(Integer(1) - Integer(1)).is_homogeneous()


def test_x_minus_1_is_not_homogeneous():
    """x - 1 has mixed degrees (degree 1 and degree 0), so not homogeneous.

    Mirrors homogenization_test.cpp/homogenization_needed_x_minus_1.
    """
    x = Variable("x")
    assert not (x - Integer(1)).is_homogeneous()


def test_transcendentals_are_not_homogeneous():
    """sin(x), cos(x), exp(x) are not homogeneous (non-polynomial).

    Mirrors homogenization_test.cpp/nothomogeneous_sin_x etc.
    """
    x = Variable("x")
    assert not sin(x).is_homogeneous()
    assert not cos(x).is_homogeneous()
    assert not tan(x).is_homogeneous()
    assert not exp(x).is_homogeneous()
    assert not log(x).is_homogeneous()
    assert not asin(x).is_homogeneous()
    assert not acos(x).is_homogeneous()
    assert not atan(x).is_homogeneous()


def test_homogenize_x_minus_1():
    """Homogenizing x - 1 with variable group [x] and hom var h gives x - h.

    Mirrors homogenization_test.cpp/homogenization_needed_x_minus_1.
    """
    x = Variable("x"); h = Variable("h")
    vg = VariableGroup(); vg.append(x)

    f = x - Integer(1)
    assert not f.is_homogeneous()
    g = f.homogenized(vg, h)        # x - h
    assert g.is_homogeneous()
    assert not f.is_homogeneous()   # original untouched

    # Evaluate g = x - h through the SLP: at (x=2, h=1) → 1, and at (x=3, h=2) → 1
    assert mp.abs(g.eval(x=2, h=1) - mpfr_complex("1")) < mpfr_float("1e-14")
    assert mp.abs(g.eval(x=3, h=2) - mpfr_complex("1")) < mpfr_float("1e-14")


def test_homogenize_leaves_already_homogeneous_unchanged():
    """Homogenizing a degree-1 expression should leave it homogeneous.

    Mirrors homogenization_test.cpp/no_homogenization_needed_x.
    """
    x = Variable("x"); h = Variable("h")
    vg = VariableGroup(); vg.append(x)

    f = x  # degree 1, already homogeneous
    assert f.is_homogeneous()
    g = f.homogenized(vg, h)  # nothing to pad; returns an equivalent (homogeneous) tree
    assert g.is_homogeneous()


def test_forbid_doubles(op):
    """
    make sure that we're correctly forbidding mixing in doubles to making symbolic expressions

    you should get exceptions when you try to do it.
    """
    x, y, z, p, a, b, tol_d, tol_mp = op

    with pytest.raises(TypeError):
        0.1 + x

    with pytest.raises(TypeError):
        0.1 - x

    with pytest.raises(TypeError):
        0.1 * x

    with pytest.raises(TypeError):
        0.1 / x

    with pytest.raises(TypeError):
        x / 0.1

    with pytest.raises(TypeError):
        0.1**x

    with pytest.raises(TypeError):
        x**0.1
