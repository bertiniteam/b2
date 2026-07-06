# This file is part of Bertini 2.
#
# python/test/classes/differentiation_safety_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/classes/differentiation_safety_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/classes/differentiation_safety_test.py.  If not, see <http://www.gnu.org/licenses/>.
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

"""The user-visible contract of differentiation: no surprises.

If you hold a function f (or any subexpression of it), differentiating -- in any flavor:
explicit-variable, no-argument Jacobian form, or through a System with default settings -- must
never change what you hold (its printed form, its degree, or its value through the SLP), and the
derivative must not sprout variables that were not in f.
"""

import pytest

import bertini as pb
from bertini.symbolics import Variable, Rational
from bertini.symbolics import sin, gather_variables

from eval_helper import eval_at


@pytest.fixture
def held():
    """f, plus a subexpression g the user also holds (shared inside f), and an evaluation point."""
    x, y = Variable('x'), Variable('y')
    g = x**2 + y
    f = sin(g) + g * y - Rational('1/3') * x
    pt = dict(x=complex(1.25, -0.3), y=complex(-0.5, 0.75))
    return x, y, g, f, pt


def test_explicit_diff_leaves_f_alone(held):
    x, y, g, f, pt = held
    s_f, s_g = str(f), str(g)
    v_f = eval_at(f, **pt)
    deg_f = f.degree()

    d = f.differentiate(x)
    d.differentiate(y)  # second derivative too

    assert str(f) == s_f
    assert str(g) == s_g
    assert f.degree() == deg_f
    assert eval_at(f, **pt) == v_f


def test_jacobian_form_diff_leaves_f_alone(held):
    x, y, g, f, pt = held
    s_f, s_g = str(f), str(g)
    v_f = eval_at(f, **pt)

    f.differentiate()  # no-arg form: builds Differential-leaf tree

    assert str(f) == s_f
    assert str(g) == s_g
    assert eval_at(f, **pt) == v_f


def test_system_differentiate_leaves_functions_alone(held):
    """Regression for the auto-simplify mutation channel: System::Differentiate
    used to run an in-place Simplify over derivative trees that share subtrees
    with f, restructuring parts of f itself."""
    x, y, g, f, pt = held
    sys = pb.System()
    sys.add_function(f)
    sys.add_variable_group(pb.VariableGroup([x, y]))

    s_f, s_g = str(f), str(g)
    s_fn = str(sys.function(0))

    sys.differentiate()

    assert str(f) == s_f
    assert str(g) == s_g
    assert str(sys.function(0)) == s_fn


def test_no_new_variables(held):
    x, y, g, f, pt = held
    vars_f = {str(v) for v in gather_variables(f)}

    d = f.differentiate(x)
    assert {str(v) for v in gather_variables(d)} <= vars_f

    dd = d.differentiate(y)
    assert {str(v) for v in gather_variables(dd)} <= vars_f

    # the Jacobian form introduces Differential leaves -- they are not
    # Variables and must not appear as variables
    dj = f.differentiate()
    assert {str(v) for v in gather_variables(dj)} <= vars_f


def test_derivative_wrt_foreign_variable_is_variable_free(held):
    x, y, g, f, pt = held
    z = Variable('z')
    d = f.differentiate(z)
    assert len(gather_variables(d)) == 0
    assert eval_at(d) == 0


def test_repeat_differentiation_consistent(held):
    x, y, g, f, pt = held
    d1 = f.differentiate(x)
    d2 = f.differentiate(x)
    assert str(d1) == str(d2)
    assert eval_at(d1, **pt) == eval_at(d2, **pt)

    # taking a second derivative does not perturb the first
    s = str(d1)
    d1.differentiate(y)
    assert str(d1) == s
