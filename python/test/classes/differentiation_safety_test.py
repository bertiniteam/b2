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

If you hold a function f (or any subexpression of it), differentiating --
in any flavor: explicit-variable, no-argument Jacobian form, or through a
System with default settings -- must never change what you hold, and the
derivative must not sprout variables that were not in f.  Derivative trees
deliberately share the SAME Variable objects as f, so setting a value once
drives both.
"""

import pytest

import bertini as pb
from bertini.function_tree.symbol import Variable, Rational
from bertini.function_tree import sin, gather_variables


@pytest.fixture
def held():
    """f, plus a subexpression g the user also holds, shared inside f."""
    x, y = Variable('x'), Variable('y')
    g = x**2 + y
    f = sin(g) + g * y - Rational('1/3') * x
    x.set_current_value(complex(1.25, -0.3))
    y.set_current_value(complex(-0.5, 0.75))
    return x, y, g, f


def test_explicit_diff_leaves_f_alone(held):
    x, y, g, f = held
    s_f, s_g = str(f), str(g)
    v_f = f.eval_d()
    deg_f = f.degree()

    d = f.differentiate(x)
    d.differentiate(y)  # second derivative too

    assert str(f) == s_f
    assert str(g) == s_g
    assert f.degree() == deg_f
    f.reset()
    assert f.eval_d() == v_f


def test_jacobian_form_diff_leaves_f_alone(held):
    x, y, g, f = held
    s_f, s_g = str(f), str(g)
    v_f = f.eval_d()

    f.differentiate()  # no-arg form: builds Differential-leaf tree

    assert str(f) == s_f
    assert str(g) == s_g
    f.reset()
    assert f.eval_d() == v_f


def test_system_differentiate_leaves_functions_alone(held):
    """Regression for the auto-simplify mutation channel: System::Differentiate
    used to run an in-place Simplify over derivative trees that share subtrees
    with f, restructuring parts of f itself."""
    x, y, g, f = held
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
    x, y, g, f = held
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
    x, y, g, f = held
    z = Variable('z')
    d = f.differentiate(z)
    assert len(gather_variables(d)) == 0
    assert d.eval_d() == 0


def test_shared_variables_are_live(held):
    """f' references the SAME Variable objects as f -- by design."""
    x, y, g, f = held
    d = f.differentiate(x)
    v1 = d.eval_d()
    x.set_current_value(complex(2.0, 0.125))
    d.reset()
    v2 = d.eval_d()
    assert v1 != v2


def test_repeat_differentiation_consistent(held):
    x, y, g, f = held
    d1 = f.differentiate(x)
    d2 = f.differentiate(x)
    assert str(d1) == str(d2)
    assert d1.eval_d() == d2.eval_d()

    # taking a second derivative does not perturb the first
    s = str(d1)
    d1.differentiate(y)
    assert str(d1) == s


def test_eval_cache_hygiene(held):
    x, y, g, f = held
    v = f.eval_d()
    f.differentiate(x)
    f.differentiate()
    assert f.eval_d() == v
