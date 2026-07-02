# This file is part of Bertini 2.
#
# python/bertini/_slice_ops.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/_slice_ops.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/bertini/_slice_ops.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""Friendly ``Slice.from_coefficients`` classmethod (accepts exact values + a variable list).

Replaces the old ``bertini.linalg.slice_from_coefficients``.  ``install`` captures the native
static ``from_coefficients`` (which takes an already-coerced ``(VariableGroup, mpfr_matrix)``) so
the friendly override -- which takes exact values first, coerces them, and accepts a plain iterable
of variables -- can still reach it.
"""

from bertini._coefficients import _coerce_mpfr_matrix
from bertini._pybertini.container import VariableGroup as _VariableGroup

_native = {}


def from_coefficients(cls, coefficients, variables, homogeneous=False):
    """Build a :class:`Slice` from an exact augmented coefficient matrix.

    ``coefficients`` is an ``(m x n+1)`` array/list of EXACT values (see
    :func:`bertini.coefficient`; Python floats are refused) -- one row per linear form, the trailing
    column being each form's constant term (give ``0`` there for a homogeneous slice).
    ``variables`` is the length-``n`` vector of variables the slice is over (a plain list/iterable of
    :class:`~bertini.Variable`, or a :class:`~bertini.VariableGroup`).

    Returns a ``bertini.Slice``.  Its rows are ready-made factors for a products-of-linears block
    (see :meth:`bertini.System.add_slices_as_products`)::

        s = bertini.Slice.from_coefficients([[2, 1, -1]], [x, y])   # 2x + y - 1 = 0
    """
    vg = variables if isinstance(variables, _VariableGroup) else _VariableGroup(list(variables))
    M, _, ncol = _coerce_mpfr_matrix(coefficients, "slice coefficients")
    if ncol != len(vg) + 1:
        raise ValueError(
            f"slice coefficient matrix has {ncol} columns but needs num_variables+1 = {len(vg) + 1} "
            "(one column per variable plus a trailing constant-term column)"
        )
    return _native['from_coefficients'](vg, M, homogeneous)


def install(Slice):
    """Replace ``Slice.from_coefficients`` with the friendly, coercing classmethod (idempotent)."""
    if getattr(Slice, "_b2_slice_ops_installed", False):
        return
    _native['from_coefficients'] = Slice.from_coefficients   # native static: (variables, mpfr_matrix)
    Slice.from_coefficients = classmethod(from_coefficients)
    Slice._b2_slice_ops_installed = True
