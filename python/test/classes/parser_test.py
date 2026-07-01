# This file is part of Bertini 2.
#
# python/test/parser_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/parser_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/parser_test.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
#
#   James Collins
#   West Texas A&M University
#   Spring 2016
#
#  silviana amethyst
#  UWEC
#  Spring 2018
#


__author__ = 'jcollins'

from bertini import *
import numpy as np

import bertini.parse as pp


def test_create_system():
    tol_d = 1e-14
    input = 'function f, g; variable_group x,y,z; f = 3*x*y*z; g = x^2 + y^2 + z^2 - 1;'
    sys = pp.system(input)
    #
    vals = np.array((complex(-2.43, .21), complex(4.84, -1.94), complex(-6.48, -.731)))
    sysEval = sys.eval(vals)
    #
    assert np.abs(sysEval[0].real / (233.2850778)-1) <= tol_d
    assert np.abs(sysEval[0].imag / (-86.5039806)-1) <= tol_d
    assert np.abs(sysEval[1].real / (65.978839)-1) <= tol_d
    assert np.abs(sysEval[1].imag / (-10.32604)-1) <= tol_d
    #
    sys.differentiate()
    sysJac = sys.eval_jacobian(vals)


def test_parse_unicode_variable():
    # The classic parser accepts Unicode letters (Ω, α) as variable names.
    sys = pp.system('variable_group Ω, α; function f; f = Ω^2 + α^2 - 1;')
    assert sys.num_variables() == 2
    vals = np.array((complex(0.3, 0.0), complex(0.4, 0.0)))  # 0.09 + 0.16 - 1 = -0.75
    assert abs(sys.eval(vals)[0] - (-0.75)) < 1e-12


def test_parse_utf8_bom():
    # A leading UTF-8 BOM is stripped before parsing.
    sys = pp.system('﻿variable_group x, y; function f; f = x^2 + y^2 - 1;')
    assert sys.num_variables() == 2


def _f_eval(expr, vals):
    """Parse 'f = <expr>' over x,y,z and evaluate at vals."""
    return pp.system(f'function f; variable_group x,y,z; f = {expr};').eval(vals)[0]


# Regression: a leading unary '-' used to negate the ENTIRE following expression ("-y+x"
# parsed as "-(y+x)") instead of just its operand.  This surfaced once canonical
# ordering let a printed sum lead with a subtracted term, breaking print->parse round-trips.
def test_leading_unary_minus_negates_only_its_operand():
    vals = np.array((complex(2, 0), complex(5, 0), complex(3, 0)))  # x=2, y=5, z=3
    assert abs(_f_eval('-y+x', vals) - _f_eval('x-y', vals)) < 1e-12
    assert abs(_f_eval('-y+x', vals) - (-3)) < 1e-12               # (-y)+x, not -(y+x)=-7


def test_leading_minus_on_parenthesized_sum_round_trips():
    vals = np.array((complex(2, 0), complex(5, 0), complex(3, 0)))
    assert abs(_f_eval('-(y-z)+x', vals) - _f_eval('x-(y-z)', vals)) < 1e-12
    assert abs(_f_eval('-(y-z)+x', vals) - 0) < 1e-12             # x-(y-z) == 0, not -((y-z)+x)


def test_unary_minus_binds_looser_than_power():
    vals = np.array((complex(2, 0), complex(5, 0), complex(3, 0)))
    assert abs(_f_eval('-x^2', vals) - (-4)) < 1e-12              # -(x^2), not (-x)^2 == 4


def test_unary_minus_then_product():
    vals = np.array((complex(2, 0), complex(5, 0), complex(3, 0)))
    assert abs(_f_eval('-x*y', vals) - (-10)) < 1e-12            # -(x*y)
