# This file is part of Bertini 2.
#
# python/test/classes/random_matrix_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/classes/random_matrix_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this file.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""bertini.random_matrix: random linear-coefficient / projection / slice / patch building block."""

import numpy as np
import pytest

import bertini as pb
from bertini._pybertini.function_tree import AbstractNode


def _as_complex(M):
    return np.array([[complex(e) for e in row] for row in M])


def test_shape_is_always_two_dimensional():
    assert pb.random_matrix(1, 3).shape == (1, 3)          # eigenpy would collapse a 1-row matrix to 1-D
    assert pb.random_matrix(3, 1).shape == (3, 1)
    assert pb.random_matrix(2, 4).shape == (2, 4)


def test_orthonormal_complex_rows_are_conjugate_orthonormal():
    pb.random.set_random_seed(7)
    M = _as_complex(pb.random_matrix(3, 3))                 # orthonormal=True is the default
    assert np.allclose(M @ M.conj().T, np.eye(3), atol=1e-10)


def test_orthonormal_real_is_real_orthogonal():
    pb.random.set_random_seed(11)
    M = _as_complex(pb.random_matrix(4, 4, real=True))
    assert np.allclose(M.imag, 0)
    assert np.allclose(M @ M.T, np.eye(4), atol=1e-10)


def test_bounded_modulus_entries():
    pb.random.set_random_seed(3)
    M = pb.random_matrix(2, 3, orthonormal=False)          # complex, bounded modulus
    assert M.shape == (2, 3)
    assert all(abs(complex(e)) <= np.sqrt(2) + 1e-9 for e in M.ravel())


def test_units_have_modulus_one():
    pb.random.set_random_seed(5)
    M = pb.random_matrix(3, 2, orthonormal=False, units=True)
    assert all(abs(abs(complex(e)) - 1) < 1e-12 for e in M.ravel())


def test_real_units_are_plus_or_minus_one():
    pb.random.set_random_seed(9)
    M = pb.random_matrix(2, 3, real=True, units=True, orthonormal=False)
    for e in M.ravel():
        c = complex(e)
        assert c.imag == 0 and abs(abs(c.real) - 1) < 1e-12


def test_real_bounded_modulus_has_zero_imaginary_part():
    pb.random.set_random_seed(13)
    M = pb.random_matrix(2, 2, real=True, orthonormal=False)
    assert all(complex(e).imag == 0 for e in M.ravel())


def test_symbolic_returns_coefficient_nodes():
    pb.random.set_random_seed(17)
    M = pb.random_matrix(1, 3, symbolic=True)
    assert M.shape == (1, 3)
    assert all(isinstance(e, AbstractNode) for e in M.ravel())
    # the nodes combine with a variable to build an expression
    x = pb.Variable('x')
    expr = M[0, 0] * x
    assert isinstance(expr, AbstractNode)


def test_reproducible_under_seed():
    pb.random.set_random_seed(123)
    a = pb.random_matrix(2, 2)
    pb.random.set_random_seed(123)
    b = pb.random_matrix(2, 2)
    assert all(complex(p) == complex(q) for p, q in zip(a.ravel(), b.ravel()))


def test_projection_row_stacks_onto_a_jacobian():
    """A projection is a linear functional (zero constant); its gradient row is random_matrix(1, n)."""
    pb.random.set_random_seed(21)
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    J = pb.jacobian([x * y, x - z], [x, y, z])             # (2, 3)
    pi = pb.random_matrix(1, 3, symbolic=True)             # (1, 3) coefficient nodes
    M = np.vstack([J, pi])
    assert M.shape == (3, 3)
    assert all(isinstance(e, AbstractNode) for e in M.ravel())
