# This file is part of Bertini 2.
#
# python/bertini/_points.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/_points.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/bertini/_points.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""Private: the one place a point vector handed in from Python is put into a shape the native
overloads accept.

The native functions that take a point take an ``Eigen::Matrix`` of ``complex<double>`` or of
``complex_mp``, and Boost.Python matches those by the array's dtype.  A caller has several
faithful ways to write the same vector that do *not* match:

* a plain ``list`` or ``tuple`` of numbers;
* an OBJECT-dtype array of multiprecision values.  ``np.array`` over ``complex_mp`` values infers
  the registered ``complex_mp`` dtype, which does match -- but an array built with
  ``dtype=object``, or grown by assignment, or sliced out of a larger object array, holds the same
  values behind a dtype no overload matches.

Unconverted, each produces an argument error quoting C++ Eigen types, which is a bug by the rule
that no C++ converter error reaches a Python caller raw (issues #348, #367).

This module holds the conversion once, rather than in each module that owns such a seam.
"""

import numpy as _np


def coerce_point_vector(value):
    """A point vector in a shape the native overloads accept.

    Returns the value unchanged unless it is a list/tuple, or a one-dimensional object-dtype
    numpy array whose every entry is a number.  Converted vectors come back as a numpy array:
    ``complex_mp`` when any entry is already multiprecision, double otherwise -- so the same call
    reaches the multiprecision overload or the double one according to what the caller supplied.

    Anything else passes through untouched, including an object array holding something other than
    numbers: widening what works must not silently reinterpret what did not.

    Parameters
    ----------
    value : object
        The candidate point vector, or any other argument.

    Returns
    -------
    object
        The value, converted if it was a convertible vector and unchanged otherwise.
    """
    from bertini._pybertini.multiprec import complex_mp as _complex_mp, real_mp as _real_mp

    if isinstance(value, _np.ndarray):
        if value.dtype != object or value.ndim != 1:
            return value
        if not all(isinstance(e, (_complex_mp, _real_mp, int, float, complex)) for e in value):
            return value
        value = list(value)
    elif not isinstance(value, (list, tuple)):
        return value

    if any(isinstance(e, (_complex_mp, _real_mp)) for e in value):
        return _np.array([e if isinstance(e, _complex_mp) else _complex_mp(e) for e in value])
    return _np.asarray(value, dtype=complex)
