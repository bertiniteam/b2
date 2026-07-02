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

* Complex (Boost.Multiprecision mpc)
* Float (Boost.Multiprecision mpfr)
* Int (Boost.Multiprecision mpz)
* Rational (Boost.Multiprecision.mpq)

This namespace also includes the mathematical operators, like `cos`, etc.
"""

from bertini._pybertini import multiprec as _pybmp
import numpy as np

from bertini._pybertini.multiprec import *

def Vector(n=0):
	"""DEPRECATED: a numpy array of multiprecision complexes.  Prefer a plain numpy array."""
	return np.zeros((n,), dtype = complex_mp)

__all__ = dir(_pybmp)


