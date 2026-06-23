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
from bertini.function_tree.symbol import Integer, Pi
from bertini.function_tree import sin, cos, tan, asin, acos, atan, exp, log, sqrt
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


def test_every_operator_evaluates_through_the_slp():
    """Each operator node compiles and evaluates through the SLP (f.eval).  Focused successor to the
    node-level operator matrix removed from function_tree_test.py."""
    x, y = Variable('x'), Variable('y')
    tol = mpfr_float("1e-12")
    def close(got, want):
        return mp.abs(got - mpfr_complex(want)) < tol
    assert close((x + y + Integer(3)).eval(x=2, y=5), "10")   # Sum
    assert close((x - y - Integer(1)).eval(x=5, y=2), "2")    # Subtract
    assert close((x * y * Integer(2)).eval(x=3, y=4), "24")   # Multiply
    assert close((x / y).eval(x=6, y=2), "3")                 # Divide
    assert close((-x).eval(x=3), "-3")                        # Negate
    assert close((x ** 3).eval(x=2), "8")                     # IntPower
    assert close((x ** y).eval(x=2, y=3), "8")                # Power (variable exponent)
    assert close(sqrt(x).eval(x=4), "2")                      # Sqrt
    assert close(exp(x).eval(x=0), "1")                       # Exp
    assert close(log(x).eval(x=1), "0")                       # Log
    assert close(sin(x).eval(x=0), "0")                       # Sin
    assert close(cos(x).eval(x=0), "1")                       # Cos
    assert close(tan(x).eval(x=0), "0")                       # Tan
    assert close(asin(x).eval(x=0), "0")                      # Asin
    assert close(acos(x).eval(x=1), "0")                      # Acos
    assert close(atan(x).eval(x=0), "0")                      # Atan


def test_pi_constant_evaluates_through_the_slp():
    x = Variable('x')
    assert mp.abs((Pi() * x).eval(x=1) - mpfr_complex("3.14159265358979323846")) < mpfr_float("1e-15")
