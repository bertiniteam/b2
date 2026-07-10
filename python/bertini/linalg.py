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
scalars inside the native module -- the LU that a stock ``import eigenpy`` cannot do on these
custom types (its compiled module only baked in the standard scalars).

Everyday entry points::

    x  = bertini.linalg.solve(A, b)   # solve the square system A x = b
    lu = bertini.linalg.lu(A)         # a reusable LU factorization (.solve/.determinant/.inverse)

For a start point that only needs to seed path tracking, casting to double and using
``numpy.linalg`` is faster and enough; reach for this module when you need the answer in mp.
"""

import numpy as _np

from bertini.multiprec import complex_mp as _complex_mp, real_mp as _real_mp

# the native surface: solve() (overloaded for complex_mp / real_mp) and the LU classes.
from bertini._pybertini.linalg import solve, PartialPivLU, PartialPivLUReal


def lu(A):
    """Partial-pivot LU factorization of a square multiprecision matrix.

    Dispatches on the array dtype and returns the matching eigenpy LU object (a
    :class:`PartialPivLU` for ``complex_mp``, a :class:`PartialPivLUReal` for ``real_mp``).  The
    result exposes ``solve(b)``, ``determinant()``, ``inverse()``, ``matrixLU()`` and
    ``permutationP()``.
    """
    A = _np.asarray(A)
    if A.dtype == _np.dtype(_complex_mp):
        return PartialPivLU(A)
    if A.dtype == _np.dtype(_real_mp):
        return PartialPivLUReal(A)
    raise TypeError(
        f"bertini.linalg.lu needs a complex_mp or real_mp matrix, not dtype {A.dtype!r}; "
        "for double-precision matrices use numpy.linalg")


__all__ = ['solve', 'lu', 'PartialPivLU', 'PartialPivLUReal']
