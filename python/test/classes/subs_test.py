# This file is part of Bertini 2.
#
# python/test/classes/subs_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# Bertini 2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Bertini 2.
# If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license, as well as COPYING.

"""node.subs(...): symbolic substitution (variable -> node), returning a NEW expression.

Substitution is simultaneous and single-pass; results are simplified (constants fold).  It is
kept distinct from eval (which is numeric).  Structural assertions use str() (the classic '^'
form); numeric cross-checks go through eval.  Non-identity values (x=3, y=5) are used so terms
stay separable.
"""

import pytest
import bertini as pb
from bertini import Variable
from bertini.symbolics import Integer, Complex
from fractions import Fraction
import bertini.multiprec as mp
from bertini.multiprec import real_mp as mpfr_float
from bertini.multiprec import complex_mp as mpfr_complex

TOL = mpfr_float("1e-25")


def test_freeze_forms_agree():
    x, y = Variable('x'), Variable('y')
    f = x * x * y
    # the three call forms agree, and the constant power folds (x**2 -> 3**2 -> 9)
    assert str(f.subs(x, 3)) == '9*y'
    assert str(f.subs({x: 3})) == '9*y'
    assert str(f.subs(x=3)) == '9*y'


def test_dict_keys_may_be_variables_or_names():
    x, y = Variable('x'), Variable('y')
    f = x * x * y
    assert str(f.subs({x: 3})) == str(f.subs({'x': 3}))


def test_rename_and_compose():
    x, y, z = Variable('x'), Variable('y'), Variable('z')
    f = x * x * y
    assert str(f.subs(x, z)) == 'z^2*y'          # rename
    assert str(f.subs(x, y + 1)) == '(y+1)^2*y'  # compose with an expression


def test_simultaneous_swap_does_not_cascade():
    x, y = Variable('x'), Variable('y')
    f = x * x * y
    g = f.subs({x: y, y: x})                      # -> y^2 x, NOT a cascade
    assert str(g) == 'y^2*x'
    # numerically: at x=3, y=5, y^2 x = 25*3 = 75
    assert mp.abs(g.eval(x=3, y=5) - mpfr_complex("75")) < TOL


def test_absent_variable_is_noop():
    x, y, z = Variable('x'), Variable('y'), Variable('z')
    f = x * x * y
    assert str(f.subs({z: 9})) == 'x^2*y'


def test_fraction_folds_to_rational():
    x, y = Variable('x'), Variable('y')
    f = x * x * y
    assert str(f.subs(x, Fraction(1, 2))) == '1/4*y'   # (1/2)^2 = 1/4, exact


def test_substitutes_into_exponent():
    x, y = Variable('x'), Variable('y')
    assert str((x ** y).subs(y, 2)) == 'x^2'


def test_subs_then_eval_equals_joint_eval():
    x, y = Variable('x'), Variable('y')
    f = x * x * y
    assert mp.abs(f.subs(x, 3).eval(y=5) - f.eval(x=3, y=5)) < TOL


def test_imaginary_unit_squared_is_minus_one():
    x = Variable('x')
    i = Complex("0", "1")
    # i^2 = -1, delivered by constant-power folding wherever simplification runs
    assert str((x ** 2).subs(x, i)) == '-1'
    assert str((i ** 2).simplify()) == '-1'


def test_simplify_folds_constant_powers():
    assert str((Integer(3) ** 2).simplify()) == '9'
    assert str((Integer(2) ** 10).simplify()) == '1024'


def test_repr_uses_python_power_operator():
    x, y = Variable('x'), Variable('y')
    # str keeps the classic '^' (round-trips with the parser); repr is copy-pasteable Python '**'
    assert str(x ** 2) == 'x^2'
    assert repr(x ** 2) == 'x**2'
    assert str((y + 1) ** 2) == '(y+1)^2'
    assert repr((y + 1) ** 2) == '(y+1)**2'
