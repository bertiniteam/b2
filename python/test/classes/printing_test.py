# This file is part of Bertini 2.
#
# python/test/classes/printing_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/classes/printing_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/classes/printing_test.py.  If not, see <http://www.gnu.org/licenses/>.
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

"""Tests for precedence-aware printing of function trees.

Printing wraps a child in parentheses only when its precedence is too low for
the position it occupies, instead of wrapping every operator unconditionally
(x^2+2*x*y-1, not (((x^2)+((2*x)*y))-1)).  The printed form must remain
re-parseable with unchanged meaning; the round-trip test at the bottom is the
load-bearing check.
"""

import numpy as np
import pytest

import bertini as pb
import bertini.parse as parse
import bertini.multiprec as mp
from bertini.multiprec import complex_mp as mpfr_complex
from bertini.symbolics import Variable, Integer, Rational
from bertini.symbolics import sin

from eval_helper import eval_at


@pytest.fixture
def xyz():
    return Variable('x'), Variable('y'), Variable('z')


def test_sum_prints_flat(xyz):
    x, y, z = xyz
    assert str(x + y + z) == 'x+y+z'
    assert str(x - y) == 'x-y'
    assert str(x + y - z) == 'x+y-z'


def test_subtraction_groups_sums(xyz):
    x, y, z = xyz
    assert str(x - (y + z)) == 'x-(y+z)'
    assert str(x - (y - z)) == 'x-(y-z)'


def test_mult_binds_tighter_than_sum(xyz):
    x, y, z = xyz
    assert str(x * y + z) == 'x*y+z'
    assert str((x + y) * z) == '(x+y)*z'
    # canonical ordering sorts the graded factor (y+z) ahead of the degree-1 x
    assert str(x * (y + z)) == '(y+z)*x'


def test_division_groups(xyz):
    x, y, z = xyz
    assert str(x / (y * z)) == 'x/(y*z)'
    assert str(x / y / z) == 'x/y/z'
    assert str(x / (y + z)) == 'x/(y+z)'


def test_power_printing(xyz):
    x, y, z = xyz
    assert str(x ** 2) == 'x^2'
    assert str((x + y) ** 2) == '(x+y)^2'
    assert str(x ** 2 * y) == 'x^2*y'
    assert str(x ** -2) == 'x^(-2)'
    assert str(x ** y) == 'x^y'
    assert str((x ** y) ** z) == '(x^y)^z'


def test_negation_printing(xyz):
    x, y, z = xyz
    assert str(-x) == '-x'
    assert str(-(x + y)) == '-(x+y)'
    assert str(x - -y) == 'x-(-y)'
    assert str(-x * y) == '(-x)*y'


def test_function_call_operators_self_delimit(xyz):
    x, y, _ = xyz
    assert str(sin(x * y)) == 'sin(x*y)'
    assert str(sin(x) ** 2) == 'sin(x)^2'


def test_real_constants_print_bare(xyz):
    x, _, _ = xyz
    assert str(Rational('1/3') * x) == '1/3*x'
    assert str(x / Rational('1/3')) == 'x/(1/3)'  # a divisor that prints with '/' must group
    assert str(Rational('1/3', '1/2') * x) == '(1/3,1/2)*x'  # genuinely complex: pair form


def test_the_motivating_example(xyz):
    x, y, _ = xyz
    assert str(x**2 + 2 * x * y - Integer(1)) == 'x^2+2*x*y-1'


def test_printed_form_reparses_to_same_values(xyz):
    """The load-bearing property: dropping parentheses must not change meaning."""
    x, y, z = xyz
    expr = (x + y) * z - x / (y + z) + 3 * x**2 * y - Rational('1/3') * (x - (y - z))

    text = f'function f; variable_group x,y,z; f = {expr};'
    reparsed = parse.system(text)

    vals = np.array([complex(-2.43, .21), complex(4.84, -1.94), complex(-6.48, -.731)])

    got = reparsed.eval(vals)[0]                                          # System eval (double)
    want = eval_at(expr, x=vals[0], y=vals[1], z=vals[2])                 # node expr through the SLP
    got_mp = mpfr_complex(str(got.real), str(got.imag))
    assert mp.abs(got_mp - want) / mp.abs(want) <= 1e-13
