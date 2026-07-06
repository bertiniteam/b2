# This file is part of Bertini 2.
#
# python/test/mpfr_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/mpfr_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/mpfr_test.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
#
#   silviana amethyst
#   University of Wisconsin - Eau Claire
#   Fall 2017, Spring 2018
#
#   James Collins
#   West Texas A&M University
#   Spring 2016
#

import bertini as pb

from bertini import multiprec as mp

import pytest


# global default precision is reset to 30 before every test by the autouse
# _reset_precision fixture in python/test/conftest.py, AND restored afterward -- so the
# precision-changing tests below (test_change_prec, test_mp_complex_precision) can poke
# the global default freely without leaking into their neighbors. The old per-test "reset
# back to 30" lines are gone; the fixture owns that now.


def _prec(x):
    """Return the precision of an mpfr value.

    Boost.Python exposes ``precision`` as a read-write property on macOS/Linux,
    but as a method on Windows (clang-cl). This helper returns the integer
    precision regardless of which binding form is active.
    """
    p = x.precision
    return p() if callable(p) else p


dbltol = 1e-15


# ------------------------------------------------------------------------------- Float

@pytest.fixture
def fvals():
    x = mp.real_mp("4.23")
    y = mp.real_mp("-3.86")
    z = mp.real_mp("1.1495")
    p = mp.real_mp(".34")
    tol = mp.real_mp("1e-27")
    return x, y, z, p, tol


def test_arith_int(fvals):
    x, y, z, p, tol = fvals
    assert mp.abs((x+8) - mp.real_mp("12.23")) <= tol
    assert mp.abs((y-2) - mp.real_mp("-5.86")) <= tol
    assert mp.abs((8+x) - mp.real_mp("12.23")) <= tol
    assert mp.abs((2-y) - mp.real_mp("5.86")) <= tol
    assert mp.abs((z*6) - mp.real_mp("6.897")) <= tol
    assert mp.abs((6*z) - mp.real_mp("6.897")) <= tol
    assert mp.abs((y/3) - mp.real_mp("-1.2866666666666666666666666666666667")) <= tol
    assert mp.abs((3/y) - mp.real_mp("-.77720207253886010362694300518134714")) <= tol
    assert mp.abs((x**3) - mp.real_mp("75.686967")) <= tol
    #
    result = mp.real_mp(x)
    result += 8
    assert mp.abs(result - mp.real_mp("12.23")) <= tol
    result = mp.real_mp(y)
    result -= 2
    assert mp.abs(result - mp.real_mp("-5.86")) <= tol
    result = mp.real_mp(z)
    result *= 6
    assert mp.abs(result - mp.real_mp("6.897")) <= tol
    result = mp.real_mp(y)
    result /= 3
    assert mp.abs(result - mp.real_mp("-1.2866666666666666666666666666666667")) <= tol
    #
    assert mp.abs((-z) - mp.real_mp("-1.1495")) <= tol


def test_arith_mpfr(fvals):
    x, y, z, p, tol = fvals
    assert mp.abs((x+y) - mp.real_mp("0.37")) <= tol
    assert mp.abs((z-y) - mp.real_mp("5.0095")) <= tol
    assert mp.abs((z*y) - mp.real_mp("-4.437070")) <= tol
    assert mp.abs((y/x) - mp.real_mp("-.91252955082742316784869976359338061")) <= tol
    assert mp.abs((x**y) - mp.real_mp("0.0038223124228935822000384505727705508")) <= tol
    assert mp.abs((-z) - mp.real_mp("-1.1495")) <= tol
    #
    result = mp.real_mp(x)
    result += y
    assert mp.abs(result - mp.real_mp("0.37")) <= tol
    result = mp.real_mp(z)
    result -= y
    assert mp.abs(result - mp.real_mp("5.0095")) <= tol
    result = mp.real_mp(z)
    result *= y
    assert mp.abs(result - mp.real_mp("-4.437070")) <= tol
    result = mp.real_mp(y)
    result /= x
    assert mp.abs(result - mp.real_mp("-.91252955082742316784869976359338061")) <= tol


def test_trancendentals(fvals):
    x, y, z, p, tol = fvals
    assert mp.abs((mp.exp(x)) - mp.real_mp("68.717232173846461408252914213396109")) <= tol
    assert mp.abs((mp.log(z)) - mp.real_mp("0.13932706522109918666170810230684295")) <= tol
    assert mp.abs((mp.sqrt(z)) - mp.real_mp("1.0721473779289860297522254519889560")) <= tol
    #
    assert mp.abs((mp.sin(x)) - mp.real_mp("-.88588921129660245121088859729926237")) <= tol
    assert mp.abs((mp.cos(y)) - mp.real_mp("-.75285494656729525719980460936483635")) <= tol
    assert mp.abs((mp.tan(z)) - mp.real_mp("2.2315038042849919118711153687209483")) <= tol
    #
    assert mp.abs((mp.asin(p)) - mp.real_mp("0.34691689752716170922069696210451452")) <= tol
    assert mp.abs((mp.acos(p)) - mp.real_mp("1.2238794292677349100106247295352369")) <= tol
    assert mp.abs((mp.atan(z)) - mp.real_mp("0.85483739856328448882289109284144652")) <= tol
    #
    assert mp.abs((mp.sinh(x)) - mp.real_mp("34.351339891649022639414777866662100")) <= tol
    assert mp.abs((mp.cosh(y)) - mp.real_mp("23.743209684188284295743755381842167")) <= tol
    assert mp.abs((mp.tanh(z)) - mp.real_mp("0.81758837109637920976170104688035086")) <= tol


def test_change_prec(fvals):
    x, y, z, p, tol = fvals
    mp.default_precision(40)
    tol = mp.real_mp("1e-37")
    t = mp.real_mp("4.23")
    assert mp.abs(t**(-2) - mp.real_mp("0.055888089689206333238323580861682566828183245868")) <= tol


# ----------------------------------------------------------------------------- Complex

@pytest.fixture
def cvals():
    x = mp.complex_mp("-2.43", ".21")
    y = mp.complex_mp("4.84", "-1.94")
    z = mp.complex_mp("-6.48", "-.731")
    p = mp.complex_mp("-.321", "-.72")
    tol = mp.real_mp("1e-27")
    return x, y, z, p, tol


def test_construct():
    t = mp.complex_mp(3.452)
    t = mp.complex_mp(mp.real_mp("-5.6"))
    t = mp.complex_mp("3.89")
    t = mp.complex_mp(mp.real_mp("2.98"), mp.real_mp("-1e-4"))
    t = mp.complex_mp(3.4, 3.5)
    t = mp.complex_mp("6e2", mp.real_mp("4.32"))
    t = mp.complex_mp(mp.real_mp("4.32"), "6e2")


def test_arith_mp_float(cvals):
    x, y, z, p, tol = cvals
    a = mp.real_mp("3.12"); b = mp.real_mp("-5.92")
    res = mp.complex_mp(x+a)
    assert mp.abs(res.real - mp.real_mp("0.69")) <= tol
    assert mp.abs(res.imag - mp.real_mp("0.21")) <= tol
    res = mp.complex_mp(y-b)
    assert mp.abs(res.real - mp.real_mp("10.76")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-1.94")) <= tol
    res = mp.complex_mp(a+x)
    assert mp.abs(res.real - mp.real_mp("0.69")) <= tol
    assert mp.abs(res.imag - mp.real_mp("0.21")) <= tol
    res = mp.complex_mp(b-y)
    assert mp.abs(res.real - mp.real_mp("-10.76")) <= tol
    assert mp.abs(res.imag - mp.real_mp("1.94")) <= tol
    res = mp.complex_mp(z*a)
    assert mp.abs(res.real - mp.real_mp("-20.2176")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-2.28072")) <= tol
    res = mp.complex_mp(a*z)
    assert mp.abs(res.real - mp.real_mp("-20.2176")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-2.28072")) <= tol
    res = mp.complex_mp(y/b)
    assert mp.abs(res.real - mp.real_mp("-.81756756756756756756756756756756756756756756757")) <= tol
    assert mp.abs(res.imag - mp.real_mp(".3277027027027027027027027027027027027027027027")) <= tol
    res = mp.complex_mp(b/y)
    assert mp.abs(res.real - mp.real_mp("-1.0538301972842157915642975887484736586585850264")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-.42240301296102864372618539714298324334662292381")) <= tol
    res = mp.complex_mp(x**a)
    assert mp.abs(res.real - mp.real_mp("-16.054376621961088182387920766649714821973863952")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-1.7411284591111236754359685247799914985638458821")) <= tol
    #
    #
    #
    res = mp.complex_mp(x)
    res += a
    assert mp.abs(res.real - mp.real_mp("0.69")) <= tol
    assert mp.abs(res.imag - mp.real_mp("0.21")) <= tol
    res = mp.complex_mp(y)
    res -= b
    assert mp.abs(res.real - mp.real_mp("10.76")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-1.94")) <= tol
    res = mp.complex_mp(z)
    res *= a
    assert mp.abs(res.real - mp.real_mp("-20.2176")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-2.28072")) <= tol
    res = mp.complex_mp(y)
    res /= b
    assert mp.abs(res.real - mp.real_mp("-.81756756756756756756756756756756756756756756757")) <= tol
    assert mp.abs(res.imag - mp.real_mp(".3277027027027027027027027027027027027027027027")) <= tol
    #
    res = mp.complex_mp(x**4)
    assert mp.abs(res.real - mp.real_mp("33.30735228")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-11.96306496")) <= tol
    #
    res = mp.complex_mp(-z)
    assert mp.abs(res.real - mp.real_mp("6.48")) <= tol
    assert mp.abs(res.imag - mp.real_mp(".731")) <= tol


def test_mp_complex_precision():
    mp.default_precision(30)

    x = mp.complex_mp(1)

    mp.default_precision(40)

    y = x

    a = mp.complex_mp(4)
    b = mp.complex_mp(5)

    assert _prec(y) == 30

    mp.default_precision(50)
    z = mp.complex_mp(3)

    mp.default_precision(60)

    c = mp.complex_mp(6)
    d = mp.complex_mp(7)
    w = x+y
    assert _prec(w) == 60

    c = a
    assert _prec(c) == 40  # even though the source is 40, target is 60, and APPoT.
    d = a+b
    assert _prec(d) == 60  # even though the source is 30, target is 60, and APPoT.


def test_arith_mp_complex(cvals):
    x, y, z, p, tol = cvals
    #
    res = mp.complex_mp(x+y)
    assert mp.abs(res.real - mp.real_mp("2.41")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-1.73")) <= tol
    res = mp.complex_mp(y-z)
    assert mp.abs(res.real - mp.real_mp("11.32")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-1.209")) <= tol
    res = mp.complex_mp(y+x)
    assert mp.abs(res.real - mp.real_mp("2.41")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-1.73")) <= tol
    res = mp.complex_mp(z-y)
    assert mp.abs(res.real - mp.real_mp("-11.32")) <= tol
    assert mp.abs(res.imag - mp.real_mp("1.209")) <= tol
    res = mp.complex_mp(z*x)
    assert mp.abs(res.real - mp.real_mp("15.89991")) <= tol
    assert mp.abs(res.imag - mp.real_mp("0.41553")) <= tol
    res = mp.complex_mp(x*z)
    assert mp.abs(res.real - mp.real_mp("15.89991")) <= tol
    assert mp.abs(res.imag - mp.real_mp(".41553")) <= tol
    res = mp.complex_mp(y/x)
    assert mp.abs(res.real - mp.real_mp("-2.0454866364094805849722642460917801311144730207")) <= tol
    assert mp.abs(res.imag - mp.real_mp("0.62158345940494200706001008572869389813414019163")) <= tol
    res = mp.complex_mp(x/y)
    assert mp.abs(res.real - mp.real_mp("-.44755270475041560619657805305047592426404601827")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-.13600253041648890000441351713180233327939034617")) <= tol
    res = mp.complex_mp(y**z)
    assert mp.abs(res.real - mp.real_mp("0.0000051612634484879218649489640888954160904291899461")) <= tol
    assert mp.abs(res.imag - mp.real_mp("0.16242051733741136410199105656393042124100116889e-4")) <= tol
    #
    #
    #
    res = mp.complex_mp(x)
    res += y
    assert mp.abs(res.real - mp.real_mp("2.41")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-1.73")) <= tol
    res = mp.complex_mp(y)
    res -= z
    assert mp.abs(res.real - mp.real_mp("11.32")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-1.209")) <= tol
    res = mp.complex_mp(z)
    res *= x
    assert mp.abs(res.real - mp.real_mp("15.89991")) <= tol
    assert mp.abs(res.imag - mp.real_mp(".41553")) <= tol
    res = mp.complex_mp(y)
    res /= x
    assert mp.abs(res.real - mp.real_mp("-2.0454866364094805849722642460917801311144730207")) <= tol
    assert mp.abs(res.imag - mp.real_mp("0.62158345940494200706001008572869389813414019163")) <= tol


def test_trancendentals_complex(cvals):
    x, y, z, p, tol = cvals
    #
    res = mp.exp(x)
    assert mp.abs(res.real - mp.real_mp("0.086102743899954532232498058731947255067424332219")) <= tol
    assert mp.abs(res.imag - mp.real_mp("0.018352149302889219131202317785160051395400089327")) <= tol
    res = mp.log(y)
    assert mp.abs(res.real - mp.real_mp("1.6514099178148475691128039277241118340531698491")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-0.38121862770417378405072154507774424569831993182")) <= tol
    res = mp.sqrt(z)
    assert mp.abs(res.real - mp.real_mp("0.14335482322754515813189359093523204816185445814")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-2.5496177371015053245565185485769652617478797292")) <= tol
    res = mp.sin(x)
    assert mp.abs(res.real - mp.real_mp("-0.66749329633668695550441899166308616328986315948")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-0.16020928942503633132090203927960650380076680938")) <= tol
    res = mp.cos(y)
    assert mp.abs(res.real - mp.real_mp("0.45194679593300564730917329070452033759984813611")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-3.3798161097977088705360399142708324234265626016")) <= tol
    res = mp.tan(z)
    assert mp.abs(res.real - mp.real_mp("-0.11998086808607765336591715593714295443402911227")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-0.63859741450762243500349270264429166927927928889")) <= tol
    res = mp.asin(x)
    assert mp.abs(res.real - mp.real_mp("-1.4763431474004472804452143435221887167393328861")) <= tol
    assert mp.abs(res.imag - mp.real_mp("1.5406263884278099750127157814537559611048741005")) <= tol
    res = mp.acos(y)
    assert mp.abs(res.real - mp.real_mp("0.38769800860408664087229892623614567735197135529")) <= tol
    assert mp.abs(res.imag - mp.real_mp("2.3379037587834289977359318611458042347281923566")) <= tol
    res = mp.atan(z)
    assert mp.abs(res.real - mp.real_mp("-1.4195347801361539102032503530226060969949192059")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-0.016801358511827150554928904776095870747673962940")) <= tol
    res = mp.sinh(x)
    assert mp.abs(res.real - mp.real_mp("-5.5116175435238027338707341903303682792175349461")) <= tol
    assert mp.abs(res.imag - mp.real_mp("1.1931117850318239903301857156967336540973461581")) <= tol
    res = mp.cosh(y)
    assert mp.abs(res.real - mp.real_mp("-22.821106324812396153305984541517740047129741299")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-58.969920999917202046449299163531439999586349377")) <= tol
    res = mp.tanh(z)
    assert mp.abs(res.real - mp.real_mp("-.99999948909538256503828034023523935671055767287")) <= tol
    assert mp.abs(res.imag - mp.real_mp("-0.0000046773288796255165542679839050497228616710086412")) <= tol


def test_complex_abs(cvals):
    x, y, z, p, tol = cvals
    #
    # test absolute value computation
    res = mp.abs(x)
    assert mp.abs(res - mp.real_mp("2.4390571949013413818264861306502")) <= tol


def test_complex_conj(cvals):
    x, y, z, p, tol = cvals
    # placeholder kept from the original suite (conj asserts were commented out)


def test_complex_construct_from_polar(cvals):
    x, y, z, p, tol = cvals
    # test construction of complex from polar coordinates
    res = mp.polar(mp.real_mp("3.21"), mp.real_mp("-5.62"))

    assert mp.abs(res.real - mp.real_mp("2.5295931897050156212406076422629449344206513531")) <= tol
    assert mp.abs(res.imag - mp.real_mp("1.9761726378527774897831544771425943545375239972")) <= tol


def test_complex_arg(cvals):
    x, y, z, p, tol = cvals
    # compute the argument of a complex number
    res = mp.arg(y)
    assert mp.abs(res - mp.real_mp("-.38121862770417378405072154507774424569831993182")) <= tol


def test_change_prec_complex(cvals):
    x, y, z, p, tol = cvals
    mp.default_precision(45)
    tol = mp.real_mp("1e-37")
    t = mp.complex_mp("-2.43", ".21")
    t = t**(-2)
    assert mp.abs(t.real - mp.real_mp("0.16560329111110602501494676510297183141930819429")) <= tol
    assert mp.abs(t.imag - mp.real_mp("0.028838165251841866149715852522538399390278791818")) <= tol


# --- conversion to python builtins ---
# regression: without __float__/__complex__ on the bound scalar types, CPython's
# conversion fell into the numpy user-dtype dispatch and recursed until the C
# stack overflowed (SIGSEGV).  values chosen exactly representable in binary.

def test_float_of_Float():
    assert float(mp.real_mp("2.5")) == 2.5


def test_complex_of_Float():
    assert complex(mp.real_mp("-0.25")) == -0.25 + 0j


def test_complex_of_Complex():
    assert complex(mp.complex_mp("2.5", "-0.25")) == complex(2.5, -0.25)


def test_float_of_Complex_raises():
    with pytest.raises(TypeError):
        float(mp.complex_mp("2.5", "-0.25"))


def test_complex_of_numpy_element():
    import numpy as np
    a = np.zeros((2,), dtype=mp.complex_mp)
    assert complex(a[0]) == 0j
