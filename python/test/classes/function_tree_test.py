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
    x = Float("4.3", "-9e-3")
    x = Float(y_mp)
    x = Float(mpfr_float("9.3"), mpfr_float("-3"))
    x = Float("9.2", "-43.2e2")


def test_Float_eval(sym):
    x_d, y_d, z_d, p_d, tol_d, x_mp, y_mp, z_mp, p_mp, tol_mp = sym
    x = Float(x_mp); y = Float(y_mp); z = Float(z_mp)
    #
    assert np.abs(x.eval_d().real/(-2.43)-1) <= tol_d
    assert np.abs(x.eval_d().imag/(.21)-1) <= tol_d
    #
    assert mp.abs(y.eval_mp().real/mpfr_float("4.84")-1) <= tol_mp
    assert mp.abs(y.eval_mp().imag/mpfr_float("-1.94")-1) <= tol_mp


def test_Float_funcs(sym):
    x_d, y_d, z_d, p_d, tol_d, x_mp, y_mp, z_mp, p_mp, tol_mp = sym
    x = Float(x_mp); y = Float(y_mp); z = Float(z_mp)
    #
    assert x.degree() == 0
    d = y.differentiate()
    assert mp.abs(d.eval_mp().real-mpfr_float("0")) <= tol_mp
    assert mp.abs(d.eval_mp().imag-mpfr_float("0")) <= tol_mp
    assert y.is_homogeneous()
    assert y.is_polynomial()


def test_Variable_construct():
    x = Variable("x")


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


def test_Variable_eval(sym):
    x_d, y_d, z_d, p_d, tol_d, x_mp, y_mp, z_mp, p_mp, tol_mp = sym
    x = Variable("x"); y = Variable("y")
    x.set_current_value(x_d); x.set_current_value(x_mp)
    #
    assert np.abs(x.eval_d().real/(-2.43)-1) <= tol_d
    assert np.abs(x.eval_d().imag/(.21)-1) <= tol_d
    #
    assert mp.abs(x.eval_mp().real/mpfr_float("-2.43")-1) <= tol_mp
    assert mp.abs(x.eval_mp().imag/mpfr_float(".21")-1) <= tol_mp


def test_Variable_funcs():
    x = Variable("x"); y = Variable("y")
    #
    assert x.degree() == 1
    assert x.degree(x) == 1
    assert x.degree(y) == 0
    d = y.differentiate()
    assert y.is_homogeneous()
    assert y.is_polynomial()


def test_Pi_construct(sym):
    x_d, y_d, z_d, p_d, tol_d, x_mp, y_mp, z_mp, p_mp, tol_mp = sym
    x = Pi()
    y = make_pi()
    #
    assert np.abs(x.eval_d().real/(3.1415926535897932384626433832795028841971693994)-1) <= tol_d
    assert np.abs(x.eval_d().imag - (0)) <= tol_d
    #
    assert mp.abs(x.eval_mp().real/mpfr_float("3.1415926535897932384626433832795028841971693994")-1) <= tol_mp
    assert mp.abs(x.eval_mp().imag - mpfr_float("0")) <= tol_mp
    #
    assert np.abs(y.eval_d().real/(3.1415926535897932384626433832795028841971693994)-1) <= tol_d
    assert np.abs(y.eval_d().imag - (0)) <= tol_d
    #
    assert mp.abs(y.eval_mp().real/mpfr_float("3.1415926535897932384626433832795028841971693994")-1) <= tol_mp
    assert mp.abs(y.eval_mp().imag - mpfr_float("0")) <= tol_mp


def test_E_construct(sym):
    x_d, y_d, z_d, p_d, tol_d, x_mp, y_mp, z_mp, p_mp, tol_mp = sym
    x = E()
    y = make_e()
    #
    assert np.abs(x.eval_d().real/(2.7182818284590452353602874713526624977572470937)-1) <= tol_d
    assert np.abs(x.eval_d().imag - (0)) <= tol_d
    #
    assert mp.abs(x.eval_mp().real - mpfr_float("2.7182818284590452353602874713526624977572470937")) <= tol_mp
    assert mp.abs(x.eval_mp().imag - mpfr_float("0")) <= tol_mp
    #
    assert np.abs(y.eval_d().real/(2.7182818284590452353602874713526624977572470937)-1) <= tol_d
    assert np.abs(y.eval_d().imag - (0)) <= tol_d
    #
    assert mp.abs(y.eval_mp().real/mpfr_float("2.7182818284590452353602874713526624977572470937")-1) <= tol_mp
    assert mp.abs(y.eval_mp().imag - mpfr_float("0")) <= tol_mp


def test_I_construct(sym):
    x_d, y_d, z_d, p_d, tol_d, x_mp, y_mp, z_mp, p_mp, tol_mp = sym
    y = make_i()
    #
    assert np.abs(y.eval_d().real - (0)) <= tol_d
    assert np.abs(y.eval_d().imag/(1.0)-1) <= tol_d
    #
    assert mp.abs(y.eval_mp().real - mpfr_float("0")) <= tol_mp
    assert mp.abs(y.eval_mp().imag/mpfr_float("1.0")-1) <= tol_mp


# --------------------------------------------------------------------------- operators

@pytest.fixture
def op():
    """Variables (with double + multiprecision current values), constants, tolerances."""
    x = Variable("x")
    x.set_current_value(complex(-2.43, .21))
    x.set_current_value(mpfr_complex("-2.43", ".21"))
    y = Variable("y")
    y.set_current_value(complex(4.84, -1.94))
    y.set_current_value(mpfr_complex("4.84", "-1.94"))
    z = Variable("z")
    z.set_current_value(complex(-6.48, -.731))
    z.set_current_value(mpfr_complex("-6.48", "-.731"))
    p = Variable("p")
    p.set_current_value(complex(-.321, -.72))
    p.set_current_value(mpfr_complex("-.321", "-.72"))
    a = Float(mpfr_complex("3.12", ".612"))
    b = Float(mpfr_complex("-.823", "2.62"))
    tol_d = float(9e-14)
    tol_mp = mpfr_float("1e-27")
    return x, y, z, p, a, b, tol_d, tol_mp


def test_plus(op):
    x, y, z, p, a, b, tol_d, tol_mp = op
    #
    f = (x+y+a)
    assert np.abs(f.eval_d().real/(5.53)-1) <= tol_d
    assert np.abs(f.eval_d().imag/(-1.118)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real/mpfr_float("5.53")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag/mpfr_float("-1.118")-1) <= tol_mp
    #
    f = (x+Float("3.87"))
    assert np.abs(f.eval_d().real/(1.44)-1) <= tol_d
    assert np.abs(f.eval_d().imag/(0.21)-1) <= tol_d
    #
    f = (x+mpfr_complex("3.87", "-2.1"))
    assert np.abs(f.eval_d().real/(1.44)-1) <= tol_d
    assert np.abs(f.eval_d().imag/(-1.89)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real/mpfr_float("1.44")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag/mpfr_float("-1.89")-1) <= tol_mp
    #
    f = (x+(-5))
    assert np.abs(f.eval_d().real/(-7.43)-1) <= tol_d
    assert np.abs(f.eval_d().imag/(0.21)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real/mpfr_float("-7.43")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag/mpfr_float("0.21")-1) <= tol_mp
    #
    f = (x); f += y; f += a
    assert np.abs(f.eval_d().real/(5.53)-1) <= tol_d
    assert np.abs(f.eval_d().imag/(-1.118)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real/mpfr_float("5.53")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag/mpfr_float("-1.118")-1) <= tol_mp
    #
    f = (x); f += Float("3.87")
    assert np.abs(f.eval_d().real/(1.44)-1) <= tol_d
    assert np.abs(f.eval_d().imag/(0.21)-1) <= tol_d


def test_sub(op):
    x, y, z, p, a, b, tol_d, tol_mp = op
    #
    f = (x-y-a)
    assert np.abs(f.eval_d().real/(-10.39)-1) <= tol_d
    assert np.abs(f.eval_d().imag/(1.538)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real/mpfr_float("-10.39")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag/mpfr_float("1.538")-1) <= tol_mp
    #
    f = (y-Float("3.87"))
    assert np.abs(f.eval_d().real/(0.97)-1) <= tol_d
    assert np.abs(f.eval_d().imag/(-1.94)-1) <= tol_d
    #
    f = (y-Float("3.87", "0"))
    assert np.abs(f.eval_d().real / (0.97)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-1.94)-1) <= tol_d
    #
    #
    f = (y-mpfr_complex("3.87", "-2.1"))
    assert np.abs(f.eval_d().real / (0.97)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (0.16)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("0.97")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("0.16")-1) <= tol_mp
    #
    f = (y-(-5))
    assert np.abs(f.eval_d().real / (9.84)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-1.94)-1) <= tol_d
    #
    f = (x); f -= y; f -= a
    assert np.abs(f.eval_d().real / (-10.39)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (1.538)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("-10.39")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("1.538")-1) <= tol_mp
    #
    f = (y); f -= Float("3.87")
    assert np.abs(f.eval_d().real / (0.97)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-1.94)-1) <= tol_d


def test_num_times_var(op):
    x, y, z, p, a, b, tol_d, tol_mp = op
    #
    f = (a*x*b*y)
    #
    assert np.abs(f.eval_d().real / (3.4011196056)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-110.9953448712)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("3.4011196056")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("-110.9953448712")-1) <= tol_mp


def test_var_times_var(op):
    x, y, z, p, a, b, tol_d, tol_mp = op
    #
    f = (x*y*z)
    assert np.abs(f.eval_d().real / (77.7616926)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-28.8346602)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("77.7616926")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("-28.8346602")-1) <= tol_mp


def test_var_div_var(op):
    x, y, z, p, a, b, tol_d, tol_mp = op
    #
    f = (x/y)
    assert np.abs(f.eval_d().real / (-.44755270475041560619657805305047592426404601827)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-.13600253041648890000441351713180233327939034617)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("-.44755270475041560619657805305047592426404601827")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("-.13600253041648890000441351713180233327939034617")-1) <= tol_mp


def test_trans_funcs(op):
    x, y, z, p, a, b, tol_d, tol_mp = op
    #
    f = (sin(x))
    assert np.abs(f.eval_d().real/(-.66749329633668695550441899166308616328986315948)-1) <= tol_d
    assert np.abs(f.eval_d().imag/(-.16020928942503633132090203927960650380076680938)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("-.66749329633668695550441899166308616328986315948")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("-.16020928942503633132090203927960650380076680938")-1) <= tol_mp
    #
    f = (cos(y))
    assert np.abs(f.eval_d().real / (0.45194679593300564730917329070452033759984813611)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-3.3798161097977088705360399142708324234265626016)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("0.45194679593300564730917329070452033759984813611")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("-3.3798161097977088705360399142708324234265626016")-1) <= tol_mp
    #
    f = (tan(z))
    assert np.abs(f.eval_d().real / (-.11998086808607765336591715593714295443402911227)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-.63859741450762243500349270264429166927927928889)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("-.11998086808607765336591715593714295443402911227")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("-.63859741450762243500349270264429166927927928889")-1) <= tol_mp
    #
    f = (asin(x))
    assert np.abs(f.eval_d().real / (-1.4763431474004472804452143435221887167393328861)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (1.5406263884278099750127157814537559611048741005)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("-1.4763431474004472804452143435221887167393328861")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("1.5406263884278099750127157814537559611048741005")-1) <= tol_mp
    #
    f = (acos(y))
    assert np.abs(f.eval_d().real / (0.38769800860408664087229892623614567735197135529)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (2.3379037587834289977359318611458042347281923566)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("0.38769800860408664087229892623614567735197135529")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("2.3379037587834289977359318611458042347281923566")-1) <= tol_mp
    #
    f = (atan(z))
    assert np.abs(f.eval_d().real / (-1.4195347801361539102032503530226060969949192059)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-0.16801358511827150554928904776095870747673962940e-1)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("-1.4195347801361539102032503530226060969949192059")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("-0.16801358511827150554928904776095870747673962940e-1")-1) <= tol_mp
    #
    f = (exp(y))
    assert np.abs(f.eval_d().real / (-45.639359208255772966298371983389382308765171859)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-117.94721623715960520658000231550940351946595854)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("-45.639359208255772966298371983389382308765171859")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("-117.94721623715960520658000231550940351946595854")-1) <= tol_mp
    #
    f = (log(z))
    assert np.abs(f.eval_d().real / (1.8750432590213669716781046977781508038070552297)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-3.0292589170775161726973168096174940982043177322)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("1.8750432590213669716781046977781508038070552297")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("-3.0292589170775161726973168096174940982043177322")-1) <= tol_mp


def test_power(op):
    x, y, z, p, a, b, tol_d, tol_mp = op
    #
    f = (y**3)
    assert np.abs(f.eval_d().real / (58.732432)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-129.035608)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("58.732432")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("-129.035608")-1) <= tol_mp
    #
    f = (pow(y, 3))
    assert np.abs(f.eval_d().real / (58.732432)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-129.035608)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("58.732432")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("-129.035608")-1) <= tol_mp
    #
    f = (x**p)
    assert np.abs(f.eval_d().real / (-.35190932545709434788093164550270669097948909024)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-6.7687858345625791466707575744042177964518271087)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("-.35190932545709434788093164550270669097948909024")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("-6.7687858345625791466707575744042177964518271087")-1) <= tol_mp
    #
    f = (pow(x, p))
    assert np.abs(f.eval_d().real / (-.35190932545709434788093164550270669097948909024)-1) <= tol_d
    assert np.abs(f.eval_d().imag / (-6.7687858345625791466707575744042177964518271087)-1) <= tol_d
    #
    assert mp.abs(f.eval_mp().real / mpfr_float("-.35190932545709434788093164550270669097948909024")-1) <= tol_mp
    assert mp.abs(f.eval_mp().imag / mpfr_float("-6.7687858345625791466707575744042177964518271087")-1) <= tol_mp


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

    # Evaluate g: at (x=2, h=1) → x - h = 2 - 1 = 1
    x.set_current_value(complex(2, 0))
    h.set_current_value(complex(1, 0))
    assert abs(g.eval_d() - complex(1, 0)) < 1e-14

    # At (x=3, h=2) → 3 - 2 = 1
    x.set_current_value(complex(3, 0))
    h.set_current_value(complex(2, 0))
    assert abs(g.eval_d() - complex(1, 0)) < 1e-14


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
