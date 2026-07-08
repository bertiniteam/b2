# This file is part of Bertini 2.
#
# python/bertini/_numpy_helpers.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/_numpy_helpers.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this file.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team

"""Vectorized helpers for the multiprecision types (issues #298, #301).

The mp dtypes have full native numpy support (ufuncs, reductions, sorting -- see the
"Multiprecision numbers and NumPy" docs page); these helpers are the everyday sugar on
top.  They accept a scalar, a list, or a numpy array, and return **mp-native** results,
so values stay arbitrary-precision instead of collapsing to float64.  On an mp-dtype
ndarray they ride the native numpy loops; on lists/tuples and mixed input they work
element-wise.

    bertini.real(pt)       # real parts, as real_mp
    bertini.abs(pt)        # magnitudes, as real_mp
    bertini.is_real(pt)    # is every coordinate real (imag within tol)?

They also paper over the one true numpy limitation for user dtypes: the ndarray
``.real`` / ``.imag`` attributes return silently wrong values on complex_mp arrays
(numpy cannot know a user dtype is complex-like), so component access must go through
these helpers or ``bertini.multiprec.real/imag/arg``.
"""

import numpy as _np
from decimal import Decimal as _Decimal, ROUND_HALF_EVEN as _ROUND_HALF_EVEN

import bertini.multiprec as _mp
from bertini.multiprec import real_mp as _real_mp, complex_mp as _complex_mp
from bertini.multiprec import abs as _mp_abs, conj as _mp_conj

_MP_DTYPES = (_np.dtype(_real_mp), _np.dtype(_complex_mp))


def _is_container(x):
    return isinstance(x, (list, tuple)) or (isinstance(x, _np.ndarray) and x.ndim > 0)


def _is_mp_array(x):
    """A numpy array of one of the mp dtypes -- eligible for the native ufunc loops."""
    return isinstance(x, _np.ndarray) and x.ndim > 0 and x.dtype in _MP_DTYPES


def _elementwise(fn, x):
    """Apply scalar ``fn`` over a scalar / list / numpy array, returning mp-native results
    (a numpy object array, preserving shape, for container input)."""
    if not _is_container(x):
        return fn(x)
    arr = _np.asarray(x, dtype=object)
    out = _np.empty(arr.shape, dtype=object)
    a_flat = arr.reshape(-1)
    o_flat = out.reshape(-1)
    for i in range(a_flat.size):
        o_flat[i] = fn(a_flat[i])
    return out


# --- scalar operations (robust across complex_mp / real_mp / python numbers) ------------------

def _real_scalar(v):
    if isinstance(v, _complex_mp):
        return v.real
    if isinstance(v, _real_mp):
        return v
    return getattr(v, 'real', v)


def _imag_scalar(v):
    if isinstance(v, _complex_mp):
        return v.imag
    if isinstance(v, _real_mp):
        return _real_mp(0)
    return getattr(v, 'imag', 0.0)


def _abs_scalar(v):
    if isinstance(v, (_complex_mp, _real_mp)):
        return _mp_abs(v)
    return abs(v)


def _conj_scalar(v):
    if isinstance(v, _complex_mp):
        return _mp_conj(v)
    if isinstance(v, _real_mp):
        return v
    return v.conjugate() if hasattr(v, 'conjugate') else v


def _round_real_mp(r, decimals):
    """Round a real_mp to ``decimals`` places, staying arbitrary-precision (via Decimal)."""
    q = _Decimal(1).scaleb(-decimals)
    return _real_mp(str(_Decimal(repr(r)).quantize(q, rounding=_ROUND_HALF_EVEN)))


def _round_scalar(v, decimals):
    if isinstance(v, _complex_mp):
        return _complex_mp(_round_real_mp(v.real, decimals), _round_real_mp(v.imag, decimals))
    if isinstance(v, _real_mp):
        return _round_real_mp(v, decimals)
    return round(v, decimals)


# --- the public helpers -----------------------------------------------------------------------
# each takes the native numpy path when handed an mp-dtype array (the C++ loops --
# fast, and the result is a proper mp-dtype array that keeps working with sort,
# reductions, and friends), and falls back to element-wise work for lists and
# mixed input.

def real(x):
    """Real part(s), as ``real_mp`` -- over a scalar / list / array (replaces numpy ``.real``,
    which returns silently wrong values on complex_mp arrays)."""
    if _is_mp_array(x) and x.ndim == 1:
        if x.dtype == _np.dtype(_complex_mp):
            return _mp.real(x)
        return x.copy()
    return _elementwise(_real_scalar, x)


def imag(x):
    """Imaginary part(s), as ``real_mp`` -- over a scalar / list / array (replaces numpy ``.imag``,
    which returns silently wrong values on complex_mp arrays)."""
    if _is_mp_array(x) and x.ndim == 1:
        if x.dtype == _np.dtype(_complex_mp):
            return _mp.imag(x)
        return _np.zeros(x.shape, dtype=_real_mp)
    return _elementwise(_imag_scalar, x)


def abs(x):
    """Magnitude(s), as ``real_mp`` -- over a scalar / list / array (same as ``np.abs``)."""
    if _is_mp_array(x):
        return _np.abs(x)
    return _elementwise(_abs_scalar, x)


def conj(x):
    """Complex conjugate(s) -- over a scalar / list / array (same as ``np.conj``)."""
    if _is_mp_array(x):
        return _np.conj(x)
    return _elementwise(_conj_scalar, x)


def round(x, decimals=0):
    """Round to ``decimals`` places, staying mp-native -- over a scalar / list / array (``np.round``)."""
    return _elementwise(lambda v: _round_scalar(v, decimals), x)


def is_real(point, tol=1e-10):
    """Is every coordinate of ``point`` real -- i.e. is each ``|imag| < tol``?  Returns a bool.

    The one-liner behind the notebook's real-solution filter::

        just_real = [pt for pt in solutions if bertini.is_real(pt)]
    """
    if _is_mp_array(point) and point.ndim == 1:
        if point.dtype == _np.dtype(_real_mp):
            return True
        # mp-vs-float tolerance orderings are exact and native
        return bool(_np.all(_np.abs(_mp.imag(point)) < tol))
    flat = _np.asarray(point, dtype=object).reshape(-1)
    return all(float(_abs_scalar(_imag_scalar(c))) < tol for c in flat)


def sum(x):
    """Sum of a collection, staying mp-native."""
    if _is_mp_array(x) and x.size > 0:
        return _np.sum(x)
    flat = list(_np.asarray(x, dtype=object).reshape(-1))
    if not flat:
        return 0
    total = flat[0]
    for v in flat[1:]:
        total = total + v
    return total


def norm(x):
    """Euclidean (2-)norm of a 1-D collection, as ``real_mp`` -- ``sqrt(sum |x_i|^2)``."""
    if _is_mp_array(x) and x.ndim == 1 and x.size > 0:
        a = _np.abs(x)                     # magnitudes, as real_mp
        return _mp.sqrt(_np.dot(a, a))     # dot slot + scalar sqrt, all mp
    flat = _np.asarray(x, dtype=object).reshape(-1)
    acc = None
    for v in flat:
        a = _abs_scalar(v)
        term = a * a
        acc = term if acc is None else acc + term
    if acc is None:
        return _real_mp(0)
    # mp square root via a half-power (stays arbitrary-precision)
    return acc ** _real_mp("0.5")
