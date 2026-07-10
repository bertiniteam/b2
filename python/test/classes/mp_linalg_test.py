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


# ---- QR ------------------------------------------------------------------------------------

def test_qr_solve_and_rank():
    A, b = _complex_system()
    qr = pb.linalg.qr(A)
    assert qr.rank() == 2
    x = qr.solve(b)
    assert max(abs(complex(r)) for r in (A @ x - b)) < 1e-18


def test_lstsq_overdetermined_exact_fit():
    # A is 3x2 of rank 2; b lies in the column space, so the least-squares fit is exact: x = [1, 2].
    A = np.array([[complex_mp('1'), complex_mp('0')],
                  [complex_mp('0'), complex_mp('1')],
                  [complex_mp('1'), complex_mp('1')]])
    b = np.array([complex_mp('1'), complex_mp('2'), complex_mp('3')])
    x = pb.linalg.lstsq(A, b)
    assert [complex(v) for v in x] == [1.0, 2.0]
    assert max(abs(complex(r)) for r in (A @ x - b)) < 1e-18


def test_qr_rank_deficient():
    A = np.array([[complex_mp('1'), complex_mp('2')],
                  [complex_mp('2'), complex_mp('4')]])
    assert pb.linalg.qr(A).rank() == 1


# ---- SVD -----------------------------------------------------------------------------------

def test_svd_singular_values_and_reconstruction():
    A, _ = _complex_system()                       # [[2,1],[1,3]], det 5
    s = pb.linalg.svd(A)
    sv = s.singularValues()
    assert sv.dtype == np.dtype(real_mp)           # singular values are real
    # singular values are sorted descending and positive
    assert float(sv[0]) >= float(sv[1]) > 0
    # product of singular values == |det|
    assert abs(float(sv[0]) * float(sv[1]) - 5.0) < 1e-12
    # U S V* reconstructs A
    U, V = _cmat(s.matrixU()), _cmat(s.matrixV())
    S = np.diag([float(v) for v in sv])
    assert np.allclose(U @ S @ V.conj().T, _cmat(A))


def test_svd_least_squares_solve():
    A, b = _complex_system()
    x = pb.linalg.svd(A).solve(b)
    assert max(abs(complex(r)) for r in (A @ x - b)) < 1e-18


def test_svd_real():
    A, _ = _real_system()
    sv = pb.linalg.svd(A).singularValues()
    assert sv.dtype == np.dtype(real_mp)
    assert abs(float(sv[0]) * float(sv[1]) - 5.0) < 1e-12


# ---- factories dispatch on dtype (mp -> eigenpy classes) ------------------------------------

def test_qr_svd_factories_dispatch_real_vs_complex():
    Ac, _ = _complex_system()
    Ar, _ = _real_system()
    assert isinstance(pb.linalg.qr(Ac), pb.linalg.ColPivHouseholderQR)
    assert isinstance(pb.linalg.qr(Ar), pb.linalg.ColPivHouseholderQRReal)
    assert isinstance(pb.linalg.svd(Ac), pb.linalg.JacobiSVD)
    assert isinstance(pb.linalg.svd(Ar), pb.linalg.JacobiSVDReal)


# ---- dtype-agnostic: the same entry points also handle float64 / complex128 ----------------

def _double_system(dtype):
    A = np.array([[2, 1], [1, 3]], dtype=dtype)
    b = np.array([3, 5], dtype=dtype)
    return A, b


@pytest.mark.parametrize("dtype", [float, complex])
def test_solve_and_lstsq_on_double(dtype):
    A, b = _double_system(dtype)
    x = pb.linalg.solve(A, b)                       # routes to numpy for double
    assert np.allclose(A @ x - b, 0)
    assert np.allclose(pb.linalg.lstsq(A, b), x)


@pytest.mark.parametrize("dtype", [float, complex])
def test_lu_qr_svd_objects_on_double(dtype):
    A, b = _double_system(dtype)
    # lu: solve / determinant / inverse
    lu = pb.linalg.lu(A)
    assert np.allclose(A @ lu.solve(b) - b, 0)
    assert abs(lu.determinant() - 5) < 1e-9
    assert np.allclose(A @ lu.inverse(), np.eye(2))
    # qr: solve / rank
    qr = pb.linalg.qr(A)
    assert qr.rank() == 2
    assert np.allclose(A @ qr.solve(b) - b, 0)
    # svd: singular values / U,V reconstruction / solve -- same method names as the mp object
    s = pb.linalg.svd(A)
    sv = s.singularValues()
    assert abs(sv[0] * sv[1] - 5) < 1e-9
    assert np.allclose(s.matrixU() @ np.diag(sv) @ s.matrixV().conj().T, A)
    assert np.allclose(A @ s.solve(b) - b, 0)


def test_double_result_matches_mp_result():
    # the polymorphic entry point agrees, mp vs double, on the same system
    Amp, bmp = _complex_system()
    Ad = np.array([[2, 1], [1, 3]], dtype=complex)
    bd = np.array([3, 5], dtype=complex)
    xmp = [complex(v) for v in pb.linalg.solve(Amp, bmp)]
    xd = list(pb.linalg.solve(Ad, bd))
    assert np.allclose(xmp, xd)


def test_no_type_iffing_one_call_site():
    # one function, every dtype -- the whole point
    def solve_anything(A, b):
        return pb.linalg.solve(A, b)
    Amp, bmp = _complex_system()
    Ad, bd = _double_system(complex)
    assert max(abs(complex(r)) for r in (Amp @ solve_anything(Amp, bmp) - bmp)) < 1e-18
    assert np.allclose(Ad @ solve_anything(Ad, bd) - bd, 0)
