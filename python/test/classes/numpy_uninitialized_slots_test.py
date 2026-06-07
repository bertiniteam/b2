# This file is part of Bertini 2.
#
# python/test/classes/numpy_uninitialized_slots_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/classes/numpy_uninitialized_slots_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/classes/numpy_uninitialized_slots_test.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""
Regression tests for SIGSEGV/SIGABRT on numpy arrays with never-written slots.

numpy zero-fills fresh buffers for the mpfr dtypes (NPY_NEEDS_INIT), but an
all-zero mpfr_t/mpc_t is Boost.Multiprecision's "uninitialized" sentinel, not
a valid value.  eigenpy's stock ufunc loops and cast loops read such slots
unguarded and crash inside libmpfr/libmpc.  The guarded loops in
python_bindings/include/eigenpy_interaction.hpp substitute an exact zero on
the read side; these tests pin that behavior.  Before the guards, the
arithmetic/equality/matmul tests below crashed the interpreter outright.
"""

import numpy as np
import pytest

from bertini.multiprec import Complex as mpfr_complex, Float as mpfr_float


@pytest.fixture(params=[mpfr_float, mpfr_complex], ids=["mpfr_float", "mpfr_complex"])
def dtype(request):
    return request.param


class TestUnwrittenSlotsAreExactZero:
    """Never-written slots of np.zeros/np.empty must act as exact zero, not crash."""

    def test_binary_arithmetic_on_unwritten_zeros(self, dtype):
        a = np.zeros(4, dtype=dtype)
        b = np.zeros(4, dtype=dtype)
        assert (a + b)[0] == dtype(0)
        assert (a - b)[1] == dtype(0)
        assert (a * b)[2] == dtype(0)

    def test_arithmetic_mixing_unwritten_and_written(self, dtype):
        a = np.zeros(3, dtype=dtype)
        b = np.array([dtype(1), dtype(2), dtype(3)])
        c = a + b
        assert c[0] == dtype(1)
        assert c[2] == dtype(3)
        d = a * b
        assert d[1] == dtype(0)

    def test_equality_on_unwritten_arrays(self, dtype):
        a = np.zeros(3, dtype=dtype)
        b = np.zeros(3, dtype=dtype)
        assert (a == b).all()
        assert not (a != b).any()

    def test_matmul_on_unwritten_arrays(self, dtype):
        a = np.zeros((2, 2), dtype=dtype)
        p = a @ a
        assert p[0, 0] == dtype(0)

    def test_matmul_mixing_unwritten_and_written(self, dtype):
        eye = np.array([[dtype(1), dtype(0)], [dtype(0), dtype(1)]])
        z = np.zeros((2, 2), dtype=dtype)
        assert (eye @ z)[0, 0] == dtype(0)
        assert (eye @ eye)[0, 0] == dtype(1)

    def test_unary_ops_on_unwritten_array(self, dtype):
        a = np.zeros(3, dtype=dtype)
        assert (-a)[0] == dtype(0)
        assert np.square(a)[1] == dtype(0)
        assert np.sqrt(a)[2] == dtype(0)

    def test_setitem_getitem_roundtrip_through_empty(self, dtype):
        a = np.empty(5, dtype=dtype)
        for i in range(5):
            a[i] = dtype(i)
        assert a[3] == dtype(3)

    def test_getitem_from_unwritten_empty_slot(self, dtype):
        a = np.empty(4, dtype=dtype)
        assert a[2] == dtype(0)  # getitem heals the slot to zero

    def test_copy_of_partially_written_array(self, dtype):
        a = np.empty(4, dtype=dtype)
        a[0] = dtype(7)
        b = a.copy()
        assert b[0] == dtype(7)

    def test_repr_of_unwritten_array(self, dtype):
        repr(np.zeros(3, dtype=dtype))


class TestOrderingComparitorsRealOnly:
    """mpfr_float arrays get guarded ordering ufuncs; complex types have none."""

    def test_orderings_on_unwritten_float_arrays(self):
        a = np.zeros(3, dtype=mpfr_float)
        b = np.array([mpfr_float(-1), mpfr_float(0), mpfr_float(1)])
        assert (a > b)[0]
        assert (a <= b)[1]
        assert (a < b)[2]
        assert (a >= b)[0]


class TestCastsFromUnwrittenSlots:
    """eigenpy cast loops read From-slots; unwritten slots must cast as zero."""

    def test_cast_unwritten_float_to_double(self):
        a = np.zeros(3, dtype=mpfr_float)
        d = np.array(a, dtype=np.float64)
        assert d[0] == 0.0

    def test_cast_float_to_complex_with_unwritten_slots(self):
        a = np.empty(3, dtype=mpfr_float)
        a[0] = mpfr_float(2)
        c = np.array(a, dtype=mpfr_complex)
        assert c[0] == mpfr_complex(2)
        assert c[1] == mpfr_complex(0)

    def test_cast_int64_zeros_to_mpfr_complex(self):
        # the exact pattern that used to SIGABRT in CI (pre-4b7f1689 tests)
        a = np.array(np.zeros(4, dtype=np.int64), dtype=mpfr_complex)
        a[0] = mpfr_complex(1)
        assert a[0] == mpfr_complex(1)
        assert a[1] == mpfr_complex(0)
