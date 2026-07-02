# This file is part of Bertini 2.
#
# python/bertini/operators.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/operators.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/bertini/operators.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""The symbolic math vocabulary, gathered for ``from bertini.operators import *``.

These are the elementary functions (``sin``, ``cos``, ...) and constants (``E``, ``Pi``, ``I``) used
to *build* symbolic systems on :class:`~bertini.Variable`\\ s.  They also live at the top level
(``bertini.sin``, ``bertini.Pi``, ...); this module exists only so you can pull the math vocabulary
into your namespace without importing the rest of ``bertini``::

    from bertini.operators import *
    f = sin(x) + Pi*y - E

These are the *symbolic* operators; the numeric elementary functions (acting on multiprecision
numbers rather than expression nodes) live in :mod:`bertini.multiprec`.
"""

# Re-exported from the top-level bertini package (defined there before this module is imported).
from . import sin, cos, tan, asin, acos, atan, exp, log, sqrt, E, Pi, I  # noqa: F401

__all__ = ['sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'exp', 'log', 'sqrt', 'E', 'Pi', 'I']
