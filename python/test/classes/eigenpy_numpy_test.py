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

import warnings

import numpy as np
import bertini as pb

from bertini import multiprec as mp


# global default precision is reset to 30 before every test by the autouse
# _reset_precision fixture in python/test/conftest.py.

SHAPE = (5, 10)


# --------------------------------------------------------------------------------- Float

def test_float_make_array_empty():
    """check that we can call `np.empty` using Float"""
    A = np.empty(SHAPE, dtype=mp.real_mp)
    A[0, 0] = mp.real_mp(0)


def test_float_make_array_zeros():
    """check that we can make variable precision real's directly using np.zeros"""
    A = np.zeros(SHAPE, dtype=mp.real_mp)
    A[0, 0] = mp.real_mp(1)


def test_float_make_array_zeros_with_conversion():
    """check that we can make variable precision real's by converting from a previously
    constructed array of zeros, without specifying the type in the converted-from array"""
    A = np.array(np.zeros(SHAPE), dtype=mp.real_mp)


def test_float_make_array_zeros_with_conversion_and_astype_int64():
    """check that we can make variable precision real's by converting from a previously
    constructed array of zeros, by first passing through int64"""
    A = np.array(np.zeros(SHAPE).astype(np.int64), dtype=mp.real_mp)


def test_float_make_array_ones():
    """check that we can make variable precision real's directly using np.ones"""
    A = np.ones(SHAPE, dtype=mp.real_mp)
    A[0, 0] = mp.real_mp(2)


def test_float_make_array_ones_with_conversion():
    """check that we can make variable precision real's by converting from a previously
    constructed array of ones, without specifying the type in the converted-from array"""
    A = np.array(np.ones(SHAPE), dtype=mp.real_mp)


def test_float_make_array_ones_with_conversion_and_astype_int64():
    """check that we can make variable precision real's by converting from a previously
    constructed array of ones, by first passing through int64"""
    A = np.array(np.ones(SHAPE).astype(np.int64), dtype=mp.real_mp)


def test_float_make_array_point_one():
    """check that if we make an array of variable precision reals from double 0.1, we
    don't get the same thing as if we constructed the high-precision type from a string."""
    A = np.array(0.1 * np.ones(SHAPE), dtype=mp.real_mp)
    assert np.all(A != mp.real_mp('0.1'))


# ------------------------------------------------------------------------------- Complex

def test_complex_make_array_empty():
    """check that we can call `np.empty` using Complex"""
    A = np.empty(SHAPE, dtype=mp.complex_mp)
    A[0, 0] = mp.complex_mp(0)


def test_complex_make_array_zeros():
    """check that we can make variable precision complex's directly using np.zeros"""
    A = np.zeros(SHAPE, dtype=mp.complex_mp)
    A[0, 0] = mp.complex_mp(1)


def test_complex_make_array_zeros_with_conversion():
    """check that we can make variable precision complex's by converting from a previously
    constructed array of zeros"""
    intermediary = np.zeros(SHAPE)
    assert isinstance(intermediary[0, 0], np.float64)
    A = np.array(intermediary, dtype=mp.complex_mp)


def test_complex_make_array_zeros_with_conversion_and_astype_int64():
    """check that we can make variable precision complex's by converting from a previously
    constructed array of zeros, by first passing through int64"""
    A = np.array(np.zeros(SHAPE).astype(np.int64), dtype=mp.complex_mp)


def test_complex_make_array_ones():
    """check that we can make variable precision complex's directly using np.ones"""
    A = np.ones(SHAPE, dtype=mp.complex_mp)
    A[0, 0] = mp.complex_mp(2)


def test_complex_make_array_ones_with_conversion():
    """check that we can make variable precision complex's by converting from a previously
    constructed array of ones"""
    A = np.array(np.ones(SHAPE), dtype=mp.complex_mp)


def test_complex_make_array_ones_with_conversion_and_astype_int64():
    """check that we can make variable precision complex's by converting from a previously
    constructed array of ones, by first passing through int64"""
    A = np.array(np.ones(SHAPE).astype(np.int64), dtype=mp.complex_mp)


# ------------------------------------------------ mp -> narrower dtype casts (no SIGABRT)
#
# Regression guard for the crash where storing a complex_mp carrying a nonzero imaginary
# part into a REAL (float64) numpy array hard-crashed the interpreter (SIGABRT), instead of
# behaving like numpy's builtin complex->real cast.  Root cause: the registered
# complex_mp->double cast did static_cast<double>(complex_mp), which routes through
# boost.multiprecision's complex->scalar conversion and THROWS
# "Could not convert imaginary number to scalar." for any nonzero imaginary part -- and that
# C++ throw, escaping numpy's C cast loop, called std::terminate().  The fix converts the
# real component, mirroring numpy (which discards the imaginary part).  It bit naturally, e.g.
# `M = np.zeros(...); M[i,j] = solver.real_solutions()[k][0]` (solution coords carry ~1e-13
# imaginary noise).  See python_bindings/include/eigenpy_interaction.hpp (cast<complex_mp,To>).

def test_complex_scalar_with_imag_into_real_array_element():
    """Assigning a complex_mp with nonzero imaginary part into a float64 array element must
    NOT crash; it discards the imaginary part (like numpy's complex128 -> float64)."""
    z = mp.complex_mp('0.25', '1e-13')          # deliberately nonzero imaginary part
    M = np.zeros((2, 2))
    M[0, 1] = z                                  # the exact crashing operation, pre-fix
    assert M[0, 1] == 0.25


def test_complex_scalar_zero_imag_into_real_array_element():
    """The zero-imaginary case (which never crashed) still assigns the real value."""
    M = np.zeros((2, 2))
    M[1, 0] = mp.complex_mp('-3.5', '0')
    assert M[1, 0] == -3.5


def test_complex_array_with_imag_astype_float64_discards_imag():
    """The vectorized cast loop (arr.astype) over complex_mp with nonzero imaginary parts
    also discards the imaginary part rather than throwing through numpy."""
    A = np.array([mp.complex_mp('1.5', '2e-13'), mp.complex_mp('-3.25', '9e-14')],
                 dtype=mp.complex_mp)
    B = A.astype(np.float64)
    assert B.dtype == np.float64
    assert np.array_equal(B, np.array([1.5, -3.25]))


def test_complex_vector_with_imag_broadcast_into_real_row():
    """Broadcast assignment of a complex_mp vector (nonzero imag) into a float64 row keeps
    the real parts and does not crash."""
    A = np.array([mp.complex_mp('1.5', '2e-13'), mp.complex_mp('-3.25', '9e-14')],
                 dtype=mp.complex_mp)
    R = np.zeros(2)
    R[:] = A
    assert np.array_equal(R, np.array([1.5, -3.25]))


def test_complex_to_real_matches_numpy_builtin():
    """The mp complex -> real cast produces the same values as numpy's builtin
    complex128 -> float64 cast (both discard the imaginary part)."""
    A = np.array([mp.complex_mp('1.5', '2e-13'), mp.complex_mp('-3.25', '9e-14')],
                 dtype=mp.complex_mp)
    mp_real = A.astype(np.float64)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')          # numpy emits ComplexWarning on this cast
        builtin_real = np.array([1.5 + 2e-13j, -3.25 + 9e-14j]).astype(np.float64)
    assert np.array_equal(mp_real, builtin_real)


def test_complex_astype_complex128_preserves_imag():
    """Regression the other way: complex_mp -> complex128 must still KEEP the imaginary part
    (its own cast specialization), so the real-target fix did not clobber it."""
    A = np.array([mp.complex_mp('1.5', '0.5'), mp.complex_mp('-3.25', '-0.75')],
                 dtype=mp.complex_mp)
    C = A.astype(complex)
    assert C.dtype == np.complex128
    assert np.array_equal(C, np.array([1.5 + 0.5j, -3.25 - 0.75j]))


def test_real_mp_astype_float64():
    """real_mp -> float64 (the always-real path) is unaffected and correct."""
    A = np.array([mp.real_mp('2.5'), mp.real_mp('-1.25')], dtype=mp.real_mp)
    B = A.astype(np.float64)
    assert B.dtype == np.float64
    assert np.array_equal(B, np.array([2.5, -1.25]))
