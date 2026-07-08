# This file is part of Bertini 2.
# 
# python/bertini/multiprec/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/bertini/multiprec/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/bertini/multiprec/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Copyright(C) Bertini2 Development Team
# 
#  See <http://www.gnu.org/licenses/> for a copy of the license, 
#  as well as COPYING.  Bertini2 is provided with permitted 
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
# 
#  silviana amethyst
#  UWEC
#  Spring 2018, summer 2023, winter/spring 2025
# 





"""
Multiprecision types, and functions that operate on them.

Numeric types exposed are

* complex_mp (Boost.Multiprecision mpc)
* real_mp (Boost.Multiprecision mpfr)
* int_mp (Boost.Multiprecision mpz)
* rational_mp (Boost.Multiprecision.mpq)

This namespace also includes the mathematical operators, like `cos`, etc.
"""

from bertini._pybertini import multiprec as _pybmp

from bertini._pybertini.multiprec import *

# (no Vector helper: eigenpy makes the mp number types work as numpy dtypes directly, so a plain
# numpy array -- e.g. np.zeros(n, dtype=bertini.complex_mp) -- is the vector.  numpy ufuncs
# (np.abs, np.exp, np.sum, ...) work on such arrays; see the "Multiprecision numbers and NumPy"
# docs page.  For the real/imaginary parts or argument of a COMPLEX ARRAY use this module's
# real()/imag()/arg() -- the ndarray .real/.imag attributes and np.angle return silently wrong
# values for user-defined dtypes, a numpy limitation.)

__all__ = dir(_pybmp)


