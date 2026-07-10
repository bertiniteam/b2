# This file is part of Bertini 2.
#
# python/bertini/linalg.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/linalg.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/bertini/linalg.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team

"""Dense linear algebra that works across all four scalar types, so you never type-if on dtype.

For the multiprecision types (``complex_mp`` / ``real_mp``) numpy's ``np.linalg`` is unavailable
(LAPACK is float/complex128 only), so this module drives eigenpy's own Eigen decomposition wrappers
instantiated on the mp scalars inside the native module -- the decompositions a stock
``import eigenpy`` cannot do on these custom types.  For the native double types
(``float64`` / ``complex128``) the same entry points transparently route to ``numpy.linalg`` (the
right, fast tool there).  So one call site handles every dtype::

    x  = bertini.linalg.solve(A, b)    # square system A x = b
    x  = bertini.linalg.lstsq(A, b)    # least squares (A may be rectangular)
    lu = bertini.linalg.lu(A)          # .solve / .determinant / .inverse
    qr = bertini.linalg.qr(A)          # .solve / .rank
    s  = bertini.linalg.svd(A)         # .singularValues / .matrixU / .matrixV / .solve / .rank

The factories return objects sharing that common method surface for every dtype.  For mp arrays the
object is the full eigenpy decomposition (extra methods: ``matrixLU``, ``permutationP``, ``rcond``,
``matrixQR``, ...); for double arrays it is a thin numpy-backed adapter exposing the common methods.
"""

import numpy as _np

from bertini.multiprec import complex_mp as _complex_mp, real_mp as _real_mp

# the native mp surface: solve() (overloaded for complex_mp / real_mp) and the decomposition classes.
from bertini._pybertini.linalg import solve as _native_solve
from bertini._pybertini.linalg import (
    PartialPivLU, PartialPivLUReal,
    HouseholderQR, HouseholderQRReal,
    ColPivHouseholderQR, ColPivHouseholderQRReal,
    JacobiSVD, JacobiSVDReal,
)

# Eigen decomposition options (Eigen/src/Core/util/Constants.h), for JacobiSVD's U/V computation.
_COMPUTE_THIN_U = 0x08
_COMPUTE_THIN_V = 0x20
_COMPUTE_FULL_U = 0x04
_COMPUTE_FULL_V = 0x10


def _classify(A):
    """Return ('complex_mp' | 'real_mp' | 'double', array) for A, casting integers to double."""
    A = _np.asarray(A)
    dt = A.dtype
    if dt == _np.dtype(_complex_mp):
        return 'complex_mp', A
    if dt == _np.dtype(_real_mp):
        return 'real_mp', A
    if _np.issubdtype(dt, _np.complexfloating) or _np.issubdtype(dt, _np.floating):
        return 'double', A
    if _np.issubdtype(dt, _np.integer):
        return 'double', A.astype(float)
    raise TypeError(f"bertini.linalg does not handle dtype {dt!r}")


# --- numpy-backed adapters for the double types, mirroring the eigenpy decomposition API ---------

class _DoubleLU:
    """LU-decomposition adapter over numpy for float64 / complex128 (mirrors eigenpy PartialPivLU)."""
    def __init__(self, A):
        self._A = A
    def solve(self, b):
        return _np.linalg.solve(self._A, _np.asarray(b))
    def determinant(self):
        return _np.linalg.det(self._A)
    def inverse(self):
        return _np.linalg.inv(self._A)


class _DoubleQR:
    """Column-pivoting-QR adapter over numpy (mirrors eigenpy ColPivHouseholderQR)."""
    def __init__(self, A):
        self._A = A
    def solve(self, b):
        return _np.linalg.lstsq(self._A, _np.asarray(b), rcond=None)[0]
    def rank(self):
        return int(_np.linalg.matrix_rank(self._A))


class _DoubleSVD:
    """SVD adapter over numpy (mirrors eigenpy JacobiSVD: singularValues / matrixU / matrixV / solve).

    numpy returns ``V**H`` (``Vh``); ``matrixV()`` returns ``V`` (``= Vh.conj().T``), matching the
    eigenpy convention ``A = U S V**H``.
    """
    def __init__(self, A, full_matrices=False):
        self._A = A
        self._U, self._s, self._Vh = _np.linalg.svd(A, full_matrices=full_matrices)
    def singularValues(self):
        return self._s
    def matrixU(self):
        return self._U
    def matrixV(self):
        return self._Vh.conj().T
    def solve(self, b):
        return _np.linalg.lstsq(self._A, _np.asarray(b), rcond=None)[0]
    def rank(self):
        return int(_np.linalg.matrix_rank(self._A))


# --- the dtype-agnostic entry points ------------------------------------------------------------

def solve(A, b):
    """Solve the square system ``A x = b`` (partial-pivot LU), for mp or double A/b."""
    kind, A = _classify(A)
    if kind == 'double':
        return _np.linalg.solve(A, _np.asarray(b))
    return _native_solve(A, _np.asarray(b))


def lstsq(A, b):
    """Least-squares solution of ``A x = b`` (A may be rectangular), for mp or double A/b."""
    kind, A = _classify(A)
    if kind == 'double':
        return _np.linalg.lstsq(A, _np.asarray(b), rcond=None)[0]
    return qr(A).solve(_np.asarray(b))


def lu(A):
    """Partial-pivot LU of a square matrix; returns an object with ``.solve/.determinant/.inverse``."""
    kind, A = _classify(A)
    if kind == 'complex_mp':
        return PartialPivLU(A)
    if kind == 'real_mp':
        return PartialPivLUReal(A)
    return _DoubleLU(A)


def qr(A):
    """Column-pivoting (rank-revealing) QR; returns an object with ``.solve`` (least squares) / ``.rank``.

    For mp the object is eigenpy's ColPivHouseholderQR (also ``.matrixQR``, ``.absDeterminant``, ...);
    for double it is a numpy adapter.  The plain full-rank QR is the ``HouseholderQR`` class.
    """
    kind, A = _classify(A)
    if kind == 'complex_mp':
        return ColPivHouseholderQR(A)
    if kind == 'real_mp':
        return ColPivHouseholderQRReal(A)
    return _DoubleQR(A)


def svd(A, full_matrices=False):
    """SVD; returns an object with ``.singularValues/.matrixU/.matrixV/.solve/.rank``.

    Computes (thin, by default) U and V so least-squares ``.solve`` works.  For mp the object is
    eigenpy's JacobiSVD; for double it is a numpy adapter.  Pass ``full_matrices=True`` for full U/V.
    """
    kind, A = _classify(A)
    if kind in ('complex_mp', 'real_mp'):
        options = (_COMPUTE_FULL_U | _COMPUTE_FULL_V) if full_matrices else (_COMPUTE_THIN_U | _COMPUTE_THIN_V)
        cls = JacobiSVD if kind == 'complex_mp' else JacobiSVDReal
        return cls(A, options)
    return _DoubleSVD(A, full_matrices)


__all__ = ['solve', 'lstsq', 'lu', 'qr', 'svd',
           'PartialPivLU', 'PartialPivLUReal',
           'HouseholderQR', 'HouseholderQRReal',
           'ColPivHouseholderQR', 'ColPivHouseholderQRReal',
           'JacobiSVD', 'JacobiSVDReal']
