# This file is part of Bertini 2.
#
# python/test/classes/canonicalization_test.py is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
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

"""Canonical operand ordering and how it shows up in printed expressions.

With canonicalization on (the default), a Sum/Mult normalizes its operand order, so
structurally-equal expressions print identically (and intern to one node).  The order is a
selectable monomial order, and it can be turned off to preserve the authored operand order.
"""

import pytest
import bertini as pb
from bertini import Variable, MonomialOrder
from bertini.function_tree.symbol import Integer


@pytest.fixture(autouse=True)
def _restore_canonicalization():
    """Canonicalization is session-global; save and restore it around every test here."""
    on, order = pb.canonicalize(), pb.monomial_order()
    yield
    pb.canonicalize(on)
    pb.monomial_order(order)


def test_canonicalization_is_on_by_default():
    assert pb.canonicalize() is True


def test_commutative_sum_prints_the_same_either_way():
    x, y = Variable('x'), Variable('y')
    assert str(x + y) == str(y + x) == 'x+y'


def test_commutative_product_prints_the_same_either_way():
    x, y = Variable('x'), Variable('y')
    assert str(x * y) == str(y * x) == 'x*y'


def test_coefficient_prints_first_within_a_monomial():
    x = Variable('x')
    # the degree-0 constant sorts ahead of the variables: 3*x^2, not x^2*3
    assert str(3 * x**2) == '3*x^2'
    assert str(x**2 * 3) == '3*x^2'


def test_constant_term_of_a_sum_stays_last():
    x = Variable('x')
    # within a sum, terms order by degree (the constant term is last), so this is conventional
    assert str(x**2 + 2 * x - Integer(1)) == 'x^2+2*x-1'


def test_subtracted_sum_leads_with_a_positive_term():
    t = Variable('t')
    # the natural homotopy factor reads "1-t", never "-t+1"
    assert str(Integer(1) - t) == '1-t'


def test_turning_canonicalization_off_preserves_authored_order():
    pb.canonicalize(False)
    x, y = Variable('x'), Variable('y')
    assert str(x + y) == 'x+y'
    assert str(y + x) == 'y+x'      # authored order kept (the opt-out)
    assert str(y * x) == 'y*x'


def test_monomial_order_is_selectable():
    x, y = Variable('x'), Variable('y')
    pb.monomial_order(MonomialOrder.Lex)
    lex = str(x**2 + y**3)          # lexicographic: x before y
    pb.monomial_order(MonomialOrder.GrevLex)
    grev = str(x**2 + y**3)         # graded: higher total degree (y^3) leads
    assert lex == 'x^2+y^3'
    assert grev == 'y^3+x^2'
    assert lex != grev
