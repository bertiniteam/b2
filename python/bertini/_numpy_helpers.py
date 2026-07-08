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

"""Vectorized element-wise helpers for the multiprecision types (issues #298, #301).

NumPy's ``.real`` / ``np.abs`` / ``np.round`` / identity-seeded reductions are unreliable on the
custom ``real_mp`` / ``complex_mp`` dtypes -- that boundary cannot be patched in the bindings (see
``docs/source/known_gotchas.rst``).  These helpers do the elementwise work themselves, over a scalar,
a list, or a numpy array, and return **mp-native** results (a numpy object array for array input), so
your values stay arbitrary-precision instead of collapsing to float64.

    bertini.real(pt)       # real parts, as real_mp
    bertini.abs(pt)        # magnitudes, as real_mp
    bertini.is_real(pt)    # is every coordinate real (imag within tol)?
"""

import numpy as _np
from decimal import Decimal as _Decimal, ROUND_HALF_EVEN as _ROUND_HALF_EVEN

from bertini.multiprec import real_mp as _real_mp, complex_mp as _complex_mp
from bertini.multiprec import abs as _mp_abs, conj as _mp_conj


def _is_container(x):
    return isinstance(x, (list, tuple)) or (isinstance(x, _np.ndarray) and x.ndim > 0)


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

def real(x):
    """Real part(s), as ``real_mp`` -- over a scalar / list / array (replaces numpy ``.real``)."""
    return _elementwise(_real_scalar, x)


def imag(x):
    """Imaginary part(s), as ``real_mp`` -- over a scalar / list / array (replaces numpy ``.imag``)."""
    return _elementwise(_imag_scalar, x)


def abs(x):
    """Magnitude(s), as ``real_mp`` -- over a scalar / list / array (replaces ``np.abs``)."""
    return _elementwise(_abs_scalar, x)


def conj(x):
    """Complex conjugate(s) -- over a scalar / list / array."""
    return _elementwise(_conj_scalar, x)


def round(x, decimals=0):
    """Round to ``decimals`` places, staying mp-native -- over a scalar / list / array (``np.round``)."""
    return _elementwise(lambda v: _round_scalar(v, decimals), x)


def is_real(point, tol=1e-10):
    """Is every coordinate of ``point`` real -- i.e. is each ``|imag| < tol``?  Returns a bool.

    The one-liner behind the notebook's real-solution filter::

        just_real = [pt for pt in solutions if bertini.is_real(pt)]
    """
    flat = _np.asarray(point, dtype=object).reshape(-1)
    return all(float(_abs_scalar(_imag_scalar(c))) < tol for c in flat)


def sum(x):
    """Sum of a 1-D collection, staying mp-native (sidesteps numpy's identity-reduction gotcha)."""
    flat = list(_np.asarray(x, dtype=object).reshape(-1))
    if not flat:
        return 0
    total = flat[0]
    for v in flat[1:]:
        total = total + v
    return total


def norm(x):
    """Euclidean (2-)norm of a 1-D collection, as ``real_mp`` -- ``sqrt(sum |x_i|^2)``."""
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
