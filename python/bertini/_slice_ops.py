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

from bertini._coefficients import _coerce_mpfr_matrix, _coerce_mpfr_vector
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


def through_point(cls, variables, point, dim=1, coefficients=None, real=False,
                  orthogonal=True, homogeneous=False):
    """Build a :class:`Slice` through a given point -- the one place to do this.

    ``variables`` is the vector of variables the slice is over (a :class:`~bertini.VariableGroup`
    or a flat list of :class:`~bertini.Variable`).  ``point`` is the point every linear form must
    vanish at -- a length-``num_variables`` vector of EXACT values (see :func:`bertini.coefficient`;
    Python floats are refused).

    With ``coefficients=None`` (the default) a random block of ``dim`` linear forms is generated
    (complex, or real if ``real=True``; orthonormalized when ``orthogonal``) and anchored so it
    passes through ``point``::

        s = bertini.Slice.through_point([x, y], pt)          # a random line through pt

    Pass ``coefficients`` (a bare ``dim x num_variables`` block of EXACT values, NOT augmented) to
    use exactly those directional coefficients; ``dim``/``real``/``orthogonal`` are then ignored.

    ``homogeneous=True`` builds a projective slice through ``point`` (each form ``a.x = 0`` with
    ``a.point = 0``); it is random-only, so it cannot be combined with an explicit ``coefficients``.

    Returns a ``bertini.Slice``.
    """
    vg = _coerce_slice_variables(variables, "through_point")
    p, plen = _coerce_mpfr_vector(point, "through_point point")
    if plen != len(vg):
        raise ValueError(
            f"through_point point has {plen} entries but needs num_variables = {len(vg)}")

    if coefficients is not None:
        if homogeneous:
            raise ValueError(
                "through_point: homogeneous=True is random-only for now -- do not also pass "
                "coefficients (projecting supplied rows onto the point's complement would "
                "silently rewrite them)")
        M, _, ncol = _coerce_mpfr_matrix(coefficients, "through_point coefficients")
        if ncol != len(vg):
            raise ValueError(
                f"through_point coefficients has {ncol} columns but needs num_variables = "
                f"{len(vg)} (a bare block, not augmented)")
        return _native['through_point'](vg, M, p)

    native = _native['random_real'] if real else _native['random_complex']
    return native(vg, dim, homogeneous, orthogonal, p)


def _coerce_slice_variables(variables, method_name):
    """Coerce the `variables` argument of a random Slice factory to a single VariableGroup (issue #293).

    Accepts a VariableGroup, a one-element list of groups (unwrapped), or a flat list of Variables
    (wrapped).  A list of *several* groups -- the common mistake of passing ``sys.variable_groups()``
    -- raises a precise message instead of the opaque converter TypeError.
    """
    from bertini._pybertini.container import VariableGroup as _VariableGroup
    from bertini._pybertini.function_tree.symbol import Variable as _Variable

    if isinstance(variables, _VariableGroup):
        return variables
    # Any other iterable: a list/tuple/ndarray of Variables, or a sequence of VariableGroups (e.g. the
    # ListOfVariableGroup that system.variable_groups() returns -- NOT a plain Python list).
    try:
        seq = list(variables)
    except TypeError:
        seq = None
    if seq is not None:
        if len(seq) == 1 and isinstance(seq[0], _VariableGroup):
            return seq[0]                                   # a one-element sequence of groups -- unwrap
        if seq and all(isinstance(g, _VariableGroup) for g in seq):
            raise TypeError(
                "Slice.{m} wants ONE variable group, but was given a sequence of {n} of them (this is "
                "what system.variable_groups() returns).  Pass a single group, e.g. "
                "Slice.{m}(sys.variable_groups()[0], ...), or a flat list of Variables."
                .format(m=method_name, n=len(seq)))
        if seq and all(isinstance(v, _Variable) for v in seq):
            return _VariableGroup(seq)                      # a flat list of Variables -- wrap
    raise TypeError(
        "Slice.{m}: `variables` must be a VariableGroup or a list of Variables (got a {t})"
        .format(m=method_name, t=type(variables).__name__))


def _make_random_slice_factory(native, method_name):
    def factory(cls, variables, dim, homogeneous=False, orthogonal=True):
        vg = _coerce_slice_variables(variables, method_name)
        return native(vg, dim, homogeneous, orthogonal)
    factory.__name__ = method_name
    factory.__doc__ = (
        "Make a random {kind} slice of `dim` linear forms over `variables` (a VariableGroup or a flat "
        "list of Variables).  homogeneous=True zeroes the constant column; orthogonal=True (default) "
        "orthonormalizes the coefficient block.".format(
            kind='real' if method_name == 'random_real' else 'complex'))
    return classmethod(factory)


def install(Slice):
    """Replace ``Slice.from_coefficients`` with the friendly, coercing classmethod (idempotent);
    give the random factories clear variable-group coercion + errors (issue #293)."""
    if getattr(Slice, "_b2_slice_ops_installed", False):
        return
    _native['from_coefficients'] = Slice.from_coefficients   # native static: (variables, mpfr_matrix)
    Slice.from_coefficients = classmethod(from_coefficients)
    for _name in ('random_complex', 'random_real'):
        if hasattr(Slice, _name):
            _native[_name] = getattr(Slice, _name)
            setattr(Slice, _name, _make_random_slice_factory(_native[_name], _name))
    # through_point dispatches to the native through_point (exact) or random factories (random).
    _native['through_point'] = Slice.through_point           # native static: (variables, coeffs, point)
    Slice.through_point = classmethod(through_point)
    Slice._b2_slice_ops_installed = True
