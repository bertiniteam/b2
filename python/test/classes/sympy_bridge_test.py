# This file is part of Bertini 2.
#
# python/test/classes/sympy_bridge_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/classes/sympy_bridge_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/classes/sympy_bridge_test.py.  If not, see <http://www.gnu.org/licenses/>.
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

"""Tests for the sympy bridge: exact two-way conversion and the solve round trip.

The whole suite skips cleanly when sympy is not installed -- it is an optional
dependency.
"""

import numpy as np
import pytest

sp = pytest.importorskip('sympy')

import bertini as pb
import bertini.multiprec as mp
from bertini.function_tree import sin, gather_variables
from bertini.function_tree.symbol import Rational
from bertini.sympy_bridge import from_sympy, to_sympy, system_from_sympy

from eval_helper import eval_at


@pytest.fixture
def sxy():
    return sp.symbols('x y')


# --- forward: sympy -> bertini ---

def test_forward_values_match(sxy):
    sx, sy = sxy
    expr = sx**3 * sy + 2 * sx * sy - sp.Rational(1, 3) + sp.sin(sx * sy) * sp.pi

    x, y = pb.Variable('x'), pb.Variable('y')
    tree = from_sympy(expr, [x, y])

    x0, y0 = complex(1.25, -0.3), complex(-0.5, 0.75)

    want_py = complex(expr.subs({sx: x0, sy: y0}).evalf())   # sympy oracle
    want = mp.Complex(str(want_py.real), str(want_py.imag))
    got = eval_at(tree, x=x0, y=y0)                          # bertini tree through the SLP
    assert mp.abs(got - want) / mp.abs(want) < 1e-13


def test_forward_reuses_supplied_variables(sxy):
    sx, sy = sxy
    x = pb.Variable('x')
    tree = from_sympy(sx**2, [x])
    # the supplied Variable is reused (not a fresh one), and the tree evaluates correctly
    assert {str(v) for v in gather_variables(tree)} == {'x'}
    assert mp.abs(eval_at(tree, x=complex(3.0, 0)) - mp.Complex('9')) < mp.Float('1e-14')
    assert mp.abs(eval_at(tree, x=complex(2.0, 0)) - mp.Complex('4')) < mp.Float('1e-14')


def test_forward_exact_rational(sxy):
    sx, _ = sxy
    tree = from_sympy(sp.Rational(1, 3) * sx)
    # the coefficient leaf is an exact bertini Rational, not a float64 dump
    leaf = next(tree.operand(i) for i in range(tree.num_operands())
                if isinstance(tree.operand(i), Rational))
    assert leaf.value_real() == mp.Rational('1/3')


def test_forward_big_integer():
    big = sp.Integer(10)**40 + 1
    tree = from_sympy(big)
    assert tree.value() == mp.Int('1' + '0' * 39 + '1')


def test_forward_hyperbolic_raises(sxy):
    sx, _ = sxy
    with pytest.raises(NotImplementedError, match='rewrite'):
        from_sympy(sp.sinh(sx))
    # ...and the suggested rewrite makes it convertible
    tree = from_sympy(sp.sinh(sx).rewrite(sp.exp), [pb.Variable('x')])
    assert tree is not None


def test_forward_unknown_head_raises(sxy):
    sx, _ = sxy
    with pytest.raises(NotImplementedError):
        from_sympy(sp.gamma(sx))


# --- reverse: bertini -> sympy ---

def test_reverse_symbolic_equality(sxy):
    sx, sy = sxy
    x, y = pb.Variable('x'), pb.Variable('y')
    tree = x**3 * y + 2 * x * y - Rational('1/3') * sin(x * y)
    expr = to_sympy(tree)
    want = sx**3 * sy + 2 * sx * sy - sp.Rational(1, 3) * sp.sin(sx * sy)
    assert sp.simplify(expr - want) == 0


def test_reverse_of_derivative(sxy):
    """the original motivation: symbolic access to bertini's derivative trees."""
    sx, sy = sxy
    x, y = pb.Variable('x'), pb.Variable('y')
    f = x**3 * y + sin(x * y)
    fx = to_sympy(f.differentiate(x))
    want = sp.diff(sx**3 * sy + sp.sin(sx * sy), sx)
    assert sp.simplify(fx - want) == 0


def test_reverse_jacobian_form_raises():
    x = pb.Variable('x')
    f = x**2
    with pytest.raises(NotImplementedError, match='explicit'):
        to_sympy(f.differentiate())  # no-arg form has Differential leaves


def test_reverse_float_precision_carry():
    mp.default_precision(50)
    fifty = '1.' + '3' * 49
    from bertini.function_tree.symbol import Complex
    expr = to_sympy(Complex(fifty))
    # the mpfr precision rides into sympy's Complex; allow the last digits to round
    assert str(expr)[:48] == fifty[:48]


# --- round trips ---

def test_round_trip_sympy_to_bertini_to_sympy(sxy):
    sx, sy = sxy
    cases = [
        sx**2 + 2 * sx * sy - 1,
        sp.Rational(3, 7) * sx - sp.exp(sy) + sp.sqrt(sx),
        sp.cos(sx) * sp.pi + sp.E,
    ]
    for e in cases:
        back = to_sympy(from_sympy(e))
        assert sp.simplify(back - e) == 0


def test_round_trip_bertini_to_sympy_to_bertini():
    x, y = pb.Variable('x'), pb.Variable('y')
    tree = x**2 * y - Rational('2/5') + sin(x)
    rebuilt = from_sympy(to_sympy(tree), [x, y])
    pt = dict(x=complex(0.6, -0.2), y=complex(-1.1, 0.4))
    assert mp.abs(eval_at(rebuilt, **pt) - eval_at(tree, **pt)) <= mp.Float('1e-25')


# --- the acceptance test: define in sympy, solve with bertini ---

def test_end_to_end_solve_matches_sympy(sxy):
    from bertini.nag_algorithm import ZeroDimCauchyAdaptivePrecisionTotalDegree

    sx, sy = sxy
    eqs = [sx**2 + sy**2 - 1, sx + sy]

    sys = system_from_sympy(eqs, [sx, sy])
    solver = ZeroDimCauchyAdaptivePrecisionTotalDegree(sys)
    solver.solve()
    got = [np.array([complex(s[i]) for i in range(len(s))]) for s in solver.all_solutions()]

    exact = sp.solve(eqs, [sx, sy])
    want = [np.array([complex(r[0]), complex(r[1])]) for r in exact]

    assert len(got) == len(want) == 2
    for g in got:
        assert min(np.linalg.norm(g - w) for w in want) < 1e-8
