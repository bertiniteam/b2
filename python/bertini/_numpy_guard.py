# This file is part of Bertini 2.
#
# python/bertini/_numpy_guard.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/_numpy_guard.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this file.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team

"""Make numpy's component accessors fail LOUDLY on multiprecision complex arrays.

numpy's ``ndarray.real`` / ``.imag`` are hardwired to its three built-in complex types
(``PyArray_ISCOMPLEX`` in ``getset.c``): on any other dtype, ``.real`` returns the array
itself and ``.imag`` returns zeros -- **silently wrong** for ``complex_mp``, and there is
no user-dtype hook to fix or even detect it at the numpy level.  ``np.real`` / ``np.imag``
/ ``np.angle`` are thin wrappers over those attributes and inherit the lie.

The attributes themselves are C-level and untouchable, but the module *functions* are
plain Python -- so importing bertini wraps them: called on a plain ndarray of the
multiprecision complex dtype they raise ``TypeError`` (wrong-by-construction is converted
to loud), and every other input passes straight through to the original numpy functions.
A crash is better than incorrect values.

Correct spellings, always:

* ``bertini.real(x)`` / ``bertini.imag(x)`` / ``bertini.multiprec.arg(x)``
* solution points (:class:`bertini.records.Solution`) override ``.real`` / ``.imag``
  at the subclass level and are simply correct -- the guard lets them through.

The remaining untouchable spelling is the raw ``.real`` / ``.imag`` attribute on a plain
ndarray you built yourself; see the "Multiprecision numbers and NumPy" docs page.
"""

import functools as _functools

import numpy as _np

from bertini.multiprec import complex_mp as _complex_mp

_CPLX_MP = _np.dtype(_complex_mp)

_MSG = (
    "numpy.{name}() cannot work on an array of the multiprecision complex dtype: "
    "numpy's component access is hardwired to its built-in complex types and would "
    "silently return WRONG values for user dtypes (there is no hook for bertini to fix "
    "it).  Use bertini.real(x) / bertini.imag(x) / bertini.multiprec.arg(x) instead -- "
    "or hold a bertini Solution, whose .real/.imag are correct."
)


def _is_plain_complex_mp_array(x):
    # exact-type check: subclasses (e.g. bertini's Solution) override .real/.imag
    # correctly and must pass through
    return type(x) is _np.ndarray and x.dtype == _CPLX_MP


def _would_lie(val):
    if _is_plain_complex_mp_array(val):
        return True
    # a list/tuple of complex_mp converts to a plain mp-dtype array INSIDE numpy's
    # real()/imag(), landing on the same wrong attribute path
    if isinstance(val, (list, tuple)):
        try:
            return _is_plain_complex_mp_array(_np.asanyarray(val))
        except Exception:
            return False
    return False


def _is_complex_mp_anything(val):
    # np.angle never consults .real/.imag attributes: it branches on the DTYPE and
    # would die in arctan2 with a cryptic error for every mp-complex input --
    # Solutions included, since the subclass property cannot help it.  Catch them
    # all and say what to use instead.
    if isinstance(val, _complex_mp):
        return True
    try:
        return _np.asanyarray(val).dtype == _CPLX_MP
    except Exception:
        return False


def _guarded(orig, name, applies):
    @_functools.wraps(orig)
    def wrapper(val, *args, **kwargs):
        if applies(val):
            raise TypeError(_MSG.format(name=name))
        return orig(val, *args, **kwargs)

    wrapper._bertini_guarded = True
    wrapper._bertini_original = orig
    return wrapper


def install():
    """Wrap ``np.real`` / ``np.imag`` / ``np.angle``.  Idempotent."""
    for name, applies in (("real", _would_lie),
                          ("imag", _would_lie),
                          ("angle", _is_complex_mp_anything)):
        orig = getattr(_np, name)
        if getattr(orig, "_bertini_guarded", False):
            continue
        setattr(_np, name, _guarded(orig, name, applies))
