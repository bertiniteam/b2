# This file is part of Bertini 2.
# 
# python/pybertini/minieigen/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/pybertini/minieigen/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/pybertini/minieigen/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Copyright(C) 2018 by Bertini2 Development Team
# 
#  See <http://www.gnu.org/licenses/> for a copy of the license, 
#  as well as COPYING.  Bertini2 is provided with permitted 
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
# 
#  Danielle Brake
#  UWEC
#  Spring 2018
# 




"""
The `MiniEigen <https://github.com/eudoxos/minieigen>`_ submodule, for PyBertini

Bertini2 uses `Eigen <https://eigen.tuxfamily.org/>`_ for linear algebra, enabling expression templates on linear algebraic objects and operations for arbitrary types -- most importantly the multiprecision complex numbers needed for working near singularities in algebraic geometry.

This module exposes some functionality.  Not near all of Eigen is exposed.  We could use help on this.  If you know some Boost.Python and some Eigen, please consider helping expose more of Eigen through MiniEigen.  Bertini2 is not the host of MiniEigen -- user Eudoxus on GitHub provided MiniEigen.

If you encounter any problems with this functionality, please ask for help via a `GitHub issue <https://github.com/bertiniteam/b2/issues>`_.

"""


import _pybertini
import _pybertini.minieigen

from _pybertini.minieigen import *

__all__ = dir(_pybertini.minieigen)

