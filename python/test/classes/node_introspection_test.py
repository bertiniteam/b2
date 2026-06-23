# This file is part of Bertini 2.
#
# python/test/classes/node_introspection_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/classes/node_introspection_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/classes/node_introspection_test.py.  If not, see <http://www.gnu.org/licenses/>.
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

"""Tests for walking/reconstructing function trees from Python.

These accessors exist so exporters (e.g. a future sympy bridge) can traverse a
tree: indexed operands with signs/flags on the n-ary operators, base/exponent
on the power operators, exact literal values on the number leaves, and names
on named symbols.
"""

import pytest

import bertini as pb
import bertini.function_tree as ft
from bertini.function_tree import operator as op
from bertini.function_tree import symbol as sym
from bertini import multiprec as mp

from eval_helper import eval_at


@pytest.fixture
def xy():
    return pb.Variable('x'), pb.Variable('y')


# --- n-ary operands and signs/flags ---

def test_sum_operands_and_signs(xy):
    x, y = xy
    s = op.Sum(x, True, y, False)  # x - y
    assert s.num_operands() == 2
    assert str(s.operand(0)) == 'x'
    assert str(s.operand(1)) == 'y'
    assert s.sign(0) is True
    assert s.sign(1) is False


def test_mult_operands_and_flags(xy):
    x, y = xy
    m = op.Mult(x, True, y, False)  # x / y
    assert m.num_operands() == 2
    assert str(m.operand(0)) == 'x'
    assert str(m.operand(1)) == 'y'
    assert m.mult_or_div(0) is True
    assert m.mult_or_div(1) is False


def test_operand_index_out_of_range_raises(xy):
    x, y = xy
    s = op.Sum(x, y)
    with pytest.raises(IndexError):
        s.operand(2)
    with pytest.raises(IndexError):
        s.sign(2)
    m = op.Mult(x, y)
    with pytest.raises(IndexError):
        m.mult_or_div(2)


# --- power operators ---

def test_power_base_and_exponent(xy):
    x, _ = xy
    p = x ** sym.Rational('1/2')
    assert isinstance(p, op.Power)
    assert str(p.get_base()) == 'x'
    e = p.get_exponent()
    assert isinstance(e, sym.Rational)
    assert e.value_real() == mp.Rational('1/2')


def test_integer_power_exponent(xy):
    x, _ = xy
    p = x ** 3
    assert isinstance(p, op.IntegerPower)
    assert p.exponent == 3
    assert str(p.operand()) == 'x'
    p.exponent = 5
    assert p.exponent == 5


# --- literal values of number leaves, exactly ---

def test_integer_value_exact():
    n = sym.Integer('123456789012345678901234567890')
    assert n.value() == mp.Int('123456789012345678901234567890')


def test_rational_value_exact():
    n = sym.Rational('1/3')
    assert n.value_real() == mp.Rational('1/3')
    assert n.value_imag() == mp.Rational(0)


def test_float_value():
    n = sym.Float('2.5')  # exactly representable in binary
    v = n.value()
    assert v.real == mp.Float('2.5')
    assert v.imag == mp.Float(0)


# --- names ---

def test_variable_name(xy):
    x, _ = xy
    assert x.name == 'x'
    x.name = 'renamed'
    assert x.name == 'renamed'
    assert str(x) == 'renamed'


def test_differential_name(xy):
    x, _ = xy
    d = x.differentiate()
    assert d.get_variable().name == 'x'


def test_function_root_walkable(xy):
    x, y = xy
    sys = pb.System()
    sys.add_function(x + y)
    sys.add_variable_group(pb.VariableGroup([x, y]))
    # function(i) now returns the bare expression root itself (functions are no longer wrapped).
    root = sys.function(0)
    assert isinstance(root, op.Sum)
    assert root.num_operands() == 2


# --- the gap-is-closed proof: rebuild a tree from introspection alone ---

def _rebuild(n):
    if isinstance(n, op.Sum):
        r = None
        for i in range(n.num_operands()):
            child = _rebuild(n.operand(i))
            if r is None:
                r = child if n.sign(i) else -child
            else:
                r = (r + child) if n.sign(i) else (r - child)
        return r
    if isinstance(n, op.Mult):
        r = None
        for i in range(n.num_operands()):
            child = _rebuild(n.operand(i))
            if r is None:
                r = child if n.mult_or_div(i) else (sym.Integer(1) / child)
            else:
                r = (r * child) if n.mult_or_div(i) else (r / child)
        return r
    if isinstance(n, op.Power):
        return _rebuild(n.get_base()) ** _rebuild(n.get_exponent())
    if isinstance(n, op.IntegerPower):
        return _rebuild(n.operand()) ** n.exponent
    if isinstance(n, op.Negate):
        return -_rebuild(n.operand())
    if isinstance(n, op.Sin):
        return ft.sin(_rebuild(n.operand()))
    if isinstance(n, sym.Variable):
        return n  # reuse: rebuilt tree shares the variables
    if isinstance(n, sym.Integer):
        return sym.Integer(n.value())
    if isinstance(n, sym.Rational):
        return sym.Rational(n.value_real(), n.value_imag())
    if isinstance(n, sym.Float):
        return sym.Float(n.value())
    raise NotImplementedError(type(n).__name__)


def test_rebuild_tree_from_introspection(xy):
    x, y = xy
    expr = x**2 + 2 * x * y - sym.Rational('1/3') * ft.sin(x)
    rebuilt = _rebuild(expr)

    assert str(rebuilt) == str(expr)

    # the rebuilt tree evaluates identically to the original (through the SLP)
    pt = dict(x=complex(1.25, 0), y=complex(-0.5, 0.75))
    assert mp.abs(eval_at(rebuilt, **pt) - eval_at(expr, **pt)) < 1e-15
