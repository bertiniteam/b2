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

"""Dense linear algebra for bertini's multiprecision types (``real_mp`` / ``complex_mp``).

numpy's ``np.linalg`` routes into LAPACK, which only knows ``float``/``complex128``, so it
cannot solve or factor arrays of the multiprecision dtypes.  This module fills that gap **at full
multiprecision**, backed by eigenpy's own Eigen decomposition wrappers instantiated on the mp
scalars inside the native module -- the decompositions that a stock ``import eigenpy`` cannot do on
these custom types (its compiled module only baked in the standard scalars).

Everyday entry points::

    x  = bertini.linalg.solve(A, b)    # square system A x = b (partial-pivot LU)
    x  = bertini.linalg.lstsq(A, b)    # least-squares (column-pivoting QR), possibly rectangular
    lu = bertini.linalg.lu(A)          # reusable LU  (.solve/.determinant/.inverse)
    qr = bertini.linalg.qr(A)          # reusable QR  (.solve/.rank/.matrixQR), rank-revealing
    s  = bertini.linalg.svd(A)         # SVD          (.singularValues/.matrixU/.matrixV/.solve)

Each factory dispatches on the array dtype (complex_mp / real_mp).  For a specific variant, the
underlying eigenpy classes are exposed directly: ``PartialPivLU``, ``HouseholderQR``,
``ColPivHouseholderQR``, ``JacobiSVD`` (and their ``...Real`` counterparts).

For a start point that only needs to seed path tracking, casting to double and using
``numpy.linalg`` is faster and enough; reach for this module when you need the answer in mp.
"""

import numpy as _np

from bertini.multiprec import complex_mp as _complex_mp, real_mp as _real_mp

# the native surface: solve() (overloaded for complex_mp / real_mp) and the decomposition classes.
from bertini._pybertini.linalg import (
    solve,
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


def _dispatch(A, complex_cls, real_cls, what):
    """Pick the complex_mp or real_mp class for array A, or raise for anything else."""
    A = _np.asarray(A)
    if A.dtype == _np.dtype(_complex_mp):
        return complex_cls, A
    if A.dtype == _np.dtype(_real_mp):
        return real_cls, A
    raise TypeError(
        f"bertini.linalg.{what} needs a complex_mp or real_mp matrix, not dtype {A.dtype!r}; "
        "for double-precision matrices use numpy.linalg")


def lu(A):
    """Partial-pivot LU factorization of a square multiprecision matrix.

    Returns the eigenpy LU object (``.solve(b)``, ``.determinant()``, ``.inverse()``,
    ``.matrixLU()``, ``.permutationP()``), dispatched on dtype.
    """
    cls, A = _dispatch(A, PartialPivLU, PartialPivLUReal, "lu")
    return cls(A)


def qr(A):
    """Column-pivoting (rank-revealing) Householder QR of a multiprecision matrix.

    Handles rectangular and rank-deficient matrices; exposes ``.solve(b)`` (least squares),
    ``.rank()``, ``.matrixQR()``, ``.absDeterminant()``.  For the plain full-rank QR use the
    ``HouseholderQR`` class directly.
    """
    cls, A = _dispatch(A, ColPivHouseholderQR, ColPivHouseholderQRReal, "qr")
    return cls(A)


def svd(A, full_matrices=False):
    """Two-sided Jacobi SVD of a multiprecision matrix.

    Computes the singular values and (by default, thin) U and V, so ``.singularValues()``,
    ``.matrixU()``, ``.matrixV()`` and least-squares ``.solve(b)`` all work.  Pass
    ``full_matrices=True`` for full U and V.
    """
    cls, A = _dispatch(A, JacobiSVD, JacobiSVDReal, "svd")
    if full_matrices:
        options = _COMPUTE_FULL_U | _COMPUTE_FULL_V
    else:
        options = _COMPUTE_THIN_U | _COMPUTE_THIN_V
    return cls(A, options)


def lstsq(A, b):
    """Least-squares solution of ``A x = b`` (A may be rectangular), via column-pivoting QR."""
    return qr(A).solve(b)


__all__ = ['solve', 'lstsq', 'lu', 'qr', 'svd',
           'PartialPivLU', 'PartialPivLUReal',
           'HouseholderQR', 'HouseholderQRReal',
           'ColPivHouseholderQR', 'ColPivHouseholderQRReal',
           'JacobiSVD', 'JacobiSVDReal']
