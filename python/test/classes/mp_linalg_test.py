"""Dense linear algebra for the multiprecision types: bertini.linalg.solve / lu.

numpy's np.linalg cannot touch complex_mp/real_mp arrays (LAPACK is float/complex128 only).
bertini.linalg fills the gap at full multiprecision, backed by eigenpy's own Eigen decomposition
wrappers instantiated on the mp scalars.  These tests pin the everyday surface: solve (both scalar
types), the reusable LU object (determinant/inverse/solve), the dtype-dispatching lu() factory,
and that the arithmetic really stays multiprecision (not silently truncated to double).
"""

import numpy as np
import pytest

import bertini as pb
from bertini.multiprec import complex_mp, real_mp


def _cmat(M):
    """A double-precision numpy copy of an mp matrix/vector (for np.allclose comparisons)."""
    return np.array([[complex(e) for e in np.atleast_1d(row)] for row in np.atleast_2d(M)])


# A = [[2,1],[1,3]], b = [3,5]  ->  x = [0.8, 1.4], det(A) = 5
def _complex_system():
    A = np.array([[complex_mp('2'), complex_mp('1')],
                  [complex_mp('1'), complex_mp('3')]])
    b = np.array([complex_mp('3'), complex_mp('5')])
    return A, b


def _real_system():
    A = np.array([[real_mp('2'), real_mp('1')],
                  [real_mp('1'), real_mp('3')]])
    b = np.array([real_mp('3'), real_mp('5')])
    return A, b


def test_solve_complex_mp():
    A, b = _complex_system()
    x = pb.linalg.solve(A, b)
    assert x.dtype == np.dtype(complex_mp)                 # stays mp, not cast to double
    assert [complex(v) for v in x] == [0.8, 1.4]
    resid = max(abs(complex(r)) for r in (A @ x - b))
    assert resid < 1e-25                                   # to mp precision (0.8, 1.4 not binary-exact)


def test_solve_real_mp():
    A, b = _real_system()
    x = pb.linalg.solve(A, b)
    assert x.dtype == np.dtype(real_mp)
    assert [float(v) for v in x] == [0.8, 1.4]
    assert max(abs(float(r)) for r in (A @ x - b)) < 1e-25


def test_lu_object_determinant_inverse_solve():
    A, b = _complex_system()
    lu = pb.linalg.PartialPivLU(A)
    assert complex(lu.determinant()) == 5.0
    # a factorization can be reused to solve
    x = lu.solve(b)
    assert max(abs(complex(r)) for r in (A @ x - b)) < 1e-25
    # inverse round-trips
    assert np.allclose(_cmat(A @ lu.inverse()), np.eye(2))


def test_lu_factory_dispatches_on_dtype():
    Ac, _ = _complex_system()
    Ar, _ = _real_system()
    assert isinstance(pb.linalg.lu(Ac), pb.linalg.PartialPivLU)
    assert isinstance(pb.linalg.lu(Ar), pb.linalg.PartialPivLUReal)


def test_lu_factory_rejects_double():
    with pytest.raises(TypeError):
        pb.linalg.lu(np.array([[1.0, 0.0], [0.0, 1.0]]))   # double array -> use numpy.linalg


@pytest.mark.parametrize("precision", [30, 60, 100], indirect=True)
def test_solve_is_genuinely_multiprecision(precision):
    # A = [[3,1],[1,3]], b = [1,0]  ->  x = [3/8, -1/8] exactly, at whatever precision is set.
    A = np.array([[complex_mp('3'), complex_mp('1')],
                  [complex_mp('1'), complex_mp('3')]])
    b = np.array([complex_mp('1'), complex_mp('0')])
    x = pb.linalg.solve(A, b)
    assert x[0].precision == precision                      # result carries the working precision
    assert str(x[0]) == '0.375' and str(x[1]) == '-0.125'  # exact rational, no double rounding


def test_singular_matrix_has_zero_determinant():
    # partial-pivot LU cannot flag singularity, but the determinant is exactly zero.
    A = np.array([[complex_mp('1'), complex_mp('2')],
                  [complex_mp('2'), complex_mp('4')]])
    assert complex(pb.linalg.lu(A).determinant()) == 0.0
