# This file is part of Bertini 2.
#
# python/test/mpfr_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/mpfr_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/mpfr_test.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
#
#   silviana amethyst
#   Max Planck Institute of Molecular Cell Biology and Genetics
#   Fall 2024, Spring 2025

# the purpose of this test suite is to make ensure numpy functionality with the custom types defined by Bertini (via EigenPy).

import numpy as np
import bertini as pb

from bertini import multiprec as mp


# global default precision is reset to 30 before every test by the autouse
# _reset_precision fixture in python/test/conftest.py.

SHAPE = (5, 10)


# --------------------------------------------------------------------------------- Float

def test_float_make_array_empty():
    """check that we can call `np.empty` using Float"""
    A = np.empty(SHAPE, dtype=mp.Float)
    A[0, 0] = mp.Float(0)


def test_float_make_array_zeros():
    """check that we can make variable precision real's directly using np.zeros"""
    A = np.zeros(SHAPE, dtype=mp.Float)
    A[0, 0] = mp.Float(1)


def test_float_make_array_zeros_with_conversion():
    """check that we can make variable precision real's by converting from a previously
    constructed array of zeros, without specifying the type in the converted-from array"""
    A = np.array(np.zeros(SHAPE), dtype=mp.Float)


def test_float_make_array_zeros_with_conversion_and_astype_int64():
    """check that we can make variable precision real's by converting from a previously
    constructed array of zeros, by first passing through int64"""
    A = np.array(np.zeros(SHAPE).astype(np.int64), dtype=mp.Float)


def test_float_make_array_ones():
    """check that we can make variable precision real's directly using np.ones"""
    A = np.ones(SHAPE, dtype=mp.Float)
    A[0, 0] = mp.Float(2)


def test_float_make_array_ones_with_conversion():
    """check that we can make variable precision real's by converting from a previously
    constructed array of ones, without specifying the type in the converted-from array"""
    A = np.array(np.ones(SHAPE), dtype=mp.Float)


def test_float_make_array_ones_with_conversion_and_astype_int64():
    """check that we can make variable precision real's by converting from a previously
    constructed array of ones, by first passing through int64"""
    A = np.array(np.ones(SHAPE).astype(np.int64), dtype=mp.Float)


def test_float_make_array_point_one():
    """check that if we make an array of variable precision reals from double 0.1, we
    don't get the same thing as if we constructed the high-precision type from a string."""
    A = np.array(0.1 * np.ones(SHAPE), dtype=mp.Float)
    assert np.all(A != mp.Float('0.1'))


# ------------------------------------------------------------------------------- Complex

def test_complex_make_array_empty():
    """check that we can call `np.empty` using Complex"""
    A = np.empty(SHAPE, dtype=mp.Complex)
    A[0, 0] = mp.Complex(0)


def test_complex_make_array_zeros():
    """check that we can make variable precision complex's directly using np.zeros"""
    A = np.zeros(SHAPE, dtype=mp.Complex)
    A[0, 0] = mp.Complex(1)


def test_complex_make_array_zeros_with_conversion():
    """check that we can make variable precision complex's by converting from a previously
    constructed array of zeros"""
    intermediary = np.zeros(SHAPE)
    assert isinstance(intermediary[0, 0], np.float64)
    A = np.array(intermediary, dtype=mp.Complex)


def test_complex_make_array_zeros_with_conversion_and_astype_int64():
    """check that we can make variable precision complex's by converting from a previously
    constructed array of zeros, by first passing through int64"""
    A = np.array(np.zeros(SHAPE).astype(np.int64), dtype=mp.Complex)


def test_complex_make_array_ones():
    """check that we can make variable precision complex's directly using np.ones"""
    A = np.ones(SHAPE, dtype=mp.Complex)
    A[0, 0] = mp.Complex(2)


def test_complex_make_array_ones_with_conversion():
    """check that we can make variable precision complex's by converting from a previously
    constructed array of ones"""
    A = np.array(np.ones(SHAPE), dtype=mp.Complex)


def test_complex_make_array_ones_with_conversion_and_astype_int64():
    """check that we can make variable precision complex's by converting from a previously
    constructed array of ones, by first passing through int64"""
    A = np.array(np.ones(SHAPE).astype(np.int64), dtype=mp.Complex)
