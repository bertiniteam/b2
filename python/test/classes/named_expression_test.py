# This file is part of Bertini 2.
#
# python/test/classes/named_expression_test.py is free software: you can redistribute it
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

"""Interface checks for Named(expr, name).  Numeric correctness lives in the C++
named_expression suite; here we only exercise the binding: construction, that it prints
as its name, and that it is usable as an ordinary node in further expressions.
"""

import bertini as pb


def test_named_constructs_and_prints_as_its_name():
    x, y = pb.Variable('x'), pb.Variable('y')
    a = pb.Named(x * x + y * y, "a")
    assert str(a) == "a"


def test_named_is_usable_as_a_subexpression():
    x = pb.Variable('x')
    a = pb.Named(x * x, "a")
    f = a * a + a            # a node built from the named expression
    # it still prints by name within the larger expression
    assert "a" in str(f)
