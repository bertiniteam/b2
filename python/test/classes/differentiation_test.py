# This file is part of Bertini 2.
#
# python/test/differentiation_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/differentiation_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/differentiation_test.py.  If not, see <http://www.gnu.org/licenses/>.
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
from bertini.symbolics import *
from bertini.symbolics import *
from bertini.symbolics import *
import numpy as np
import pytest

import bertini.multiprec as mp
from bertini.multiprec import real_mp as mpfr_float
from bertini.multiprec import complex_mp as mpfr_complex

from eval_helper import eval_at


# global default precision is reset to 30 before every test by the autouse
# _reset_precision fixture in python/test/conftest.py.  Partial derivatives are taken with
# node.differentiate(var) and evaluated through the SLP (eval_at) in multiple precision at the
# baseline precision; the transcendental expected constants are only ~47 digits, so this suite
# stays at the baseline precision rather than parametrizing.

TOL_D = float(1e-14)


@pytest.fixture
def tol_mp():
    return mpfr_float("1e-27")


@pytest.fixture
def diffvars():
    # variables, two Complex constants, and the evaluation point (a superset map; eval_at takes only
    # the variables each partial derivative actually contains).
    x = Variable("x")
    y = Variable("y")
    z = Variable('z')
    p = Variable('p')
    a = Complex(mpfr_complex("3.12", ".612"))
    b = Complex(mpfr_complex("-.823", "2.62"))
    pt = dict(x=mpfr_complex("-2.43", ".21"), y=mpfr_complex("4.84", "-1.94"),
              z=mpfr_complex("-6.48", "-.731"), p=mpfr_complex("-.321", "-.72"))
    return x, y, z, p, a, b, pt


def test_sum_rule(diffvars, tol_mp):
    x, y, z, p, a, b, pt = diffvars
    f = x + y + a
    #
    dfx = eval_at(f.differentiate(x), **pt)
    assert mp.abs(dfx.real / mpfr_float("1.0") - 1) <= tol_mp
    assert mp.abs(dfx.imag - mpfr_float("0")) <= tol_mp
    #
    dfy = eval_at(f.differentiate(y), **pt)
    assert mp.abs(dfy.real / mpfr_float("1.0") - 1) <= tol_mp
    assert mp.abs(dfy.imag - mpfr_float("0")) <= tol_mp
    #
    dfz = eval_at(f.differentiate(z), **pt)
    assert mp.abs(dfz.real - mpfr_float("0.0")) <= tol_mp
    assert mp.abs(dfz.imag - mpfr_float("0")) <= tol_mp


def test_power_rule(diffvars, tol_mp):
    x, y, z, p, a, b, pt = diffvars
    f = x**2 + y**3
    #
    dfx = eval_at(f.differentiate(x), **pt)
    assert mp.abs(dfx.real / mpfr_float("-4.86") - 1) <= tol_mp
    assert mp.abs(dfx.imag / mpfr_float("0.42") - 1) <= tol_mp
    #
    dfy = eval_at(f.differentiate(y), **pt)
    assert mp.abs(dfy.real / mpfr_float("58.9860") - 1) <= tol_mp
    assert mp.abs(dfy.imag / mpfr_float("-56.3376") - 1) <= tol_mp


def test_prod_rule(diffvars, tol_mp):
    x, y, z, p, a, b, pt = diffvars
    f = x**2*y**4 - a*x*y*z**2
    #
    dfx = eval_at(f.differentiate(x), **pt)
    assert mp.abs(dfx.real / mpfr_float("-559.28968169592") - 1) <= tol_mp
    assert mp.abs(dfx.imag / mpfr_float("3577.05276993648") - 1) <= tol_mp
    #
    dfy = eval_at(f.differentiate(y), **pt)
    assert mp.abs(dfy.real / mpfr_float("1161.85042980828") - 1) <= tol_mp
    assert mp.abs(dfy.imag / mpfr_float("-3157.24325320476") - 1) <= tol_mp
    #
    dfz = eval_at(f.differentiate(z), **pt)
    assert mp.abs(dfz.real / mpfr_float("-520.5265859088") - 1) <= tol_mp
    assert mp.abs(dfz.imag / mpfr_float("84.7479679056") - 1) <= tol_mp


def test_trancendental(diffvars, tol_mp):
    x, y, z, p, a, b, pt = diffvars
    f = sin(x*y) + exp(z*y) - log(x*x)
    #
    dfx = eval_at(f.differentiate(x), **pt)
    assert mp.abs(dfx.real / mpfr_float("-17.648420086229721902138620795382021306411662490") - 1) <= tol_mp
    assert mp.abs(dfx.imag / mpfr_float("-803.11883403426275105632833868183320319093878729") - 1) <= tol_mp
    #
    dfy = eval_at(f.differentiate(y), **pt)
    assert mp.abs(dfy.real / mpfr_float("-100.97157179433748763552280062599971478593963953") - 1) <= tol_mp
    assert mp.abs(dfy.imag / mpfr_float("361.98093991820979266721712882115615553425318528") - 1) <= tol_mp
    #
    dfz = eval_at(f.differentiate(z), **pt)
    assert mp.abs(dfz.real / mpfr_float("-2.1642907643013779167501866500194314960002972412e-14") - 1) <= tol_mp
    assert mp.abs(dfz.imag / mpfr_float("2.1105887207247540399884720817624768568595288922e-14") - 1) <= tol_mp


# --- derivatives come out already simplified ---
# differentiation builds trees through the Simplified* factories: literal
# zeros/ones never appear, exact Integer/Rational constants fold (no Floats),
# nested products flatten.  these exact-form assertions are the contract.

def test_derivative_trees_are_simplified():
    x = Variable('x')
    y = Variable('y')
    assert str((x * y).differentiate(x)) == 'y'
    assert str((x + y).differentiate(x)) == '1'
    assert str((x**3).differentiate(x)) == '3*x^2'
    assert str(sin(x).differentiate(x)) == 'cos(x)'
    assert str(cos(x).differentiate(x)) == '-sin(x)'
    assert str(log(x).differentiate(x)) == '1/x'
    assert str((x**3 * y).differentiate(x).differentiate(x)) == '6*x*y'
    assert str((x / y).differentiate(y)) == '-x/y^2'
