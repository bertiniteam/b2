# This file is part of Bertini 2.
#
# python/test/classes/eval_expression_test.py is free software: you can redistribute it
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

"""f.eval(x=2, y=5): evaluating a bare expression at a point with no System.

Values are bound to variables by name (variables are canonical by name across the whole
session, so the name is unambiguous).  Evaluation runs through the one evaluation engine
(the compiled SLP) via a throwaway adapter System built internally.  Every variable of the
expression must be supplied a value, and every keyword must name a variable of it.
"""

import pytest
import bertini as pb
from bertini import Variable
from bertini.function_tree.symbol import Integer
import bertini.multiprec as mp
from bertini.multiprec import Float as mpfr_float
from bertini.multiprec import Complex as mpfr_complex

TOL = mpfr_float("1e-25")


def test_eval_bare_polynomial_by_name():
    x, y = Variable('x'), Variable('y')
    f = x * x + y
    assert mp.abs(f.eval(x=2, y=5) - mpfr_complex("9")) < TOL


def test_eval_binds_by_name_not_order():
    a, b = Variable('a'), Variable('b')
    f = a - b
    # 'b' is passed first, but binding is by name: 7 - 3 == 4
    assert mp.abs(f.eval(b=3, a=7) - mpfr_complex("4")) < TOL


def test_eval_complex_valued_input():
    x = Variable('x')
    f = x * x
    # i^2 == -1
    assert mp.abs(f.eval(x=mpfr_complex("0", "1")) - mpfr_complex("-1")) < TOL


def test_eval_accepts_native_python_float():
    x = Variable('x')
    f = x * x
    # 1.5^2 == 2.25; native float carries only float64 of information
    assert mp.abs(f.eval(x=1.5) - mpfr_complex("2.25")) < mpfr_float("1e-12")


def test_eval_constant_with_no_variables():
    f = Integer(3) * Integer(4)
    assert mp.abs(f.eval() - mpfr_complex("12")) < TOL


def test_eval_missing_variable_is_an_error():
    x, y = Variable('x'), Variable('y')
    f = x + y
    with pytest.raises(Exception):
        f.eval(x=1)  # y omitted


def test_eval_unknown_variable_is_an_error():
    x = Variable('x')
    f = x * x
    with pytest.raises(Exception):
        f.eval(x=2, z=9)  # z is not in the expression (typo guard)


def test_eval_positional_argument_is_an_error():
    x = Variable('x')
    f = x * x
    with pytest.raises(Exception):
        f.eval(2)  # values must be passed by keyword
