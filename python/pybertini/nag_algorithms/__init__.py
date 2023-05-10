# This file is part of Bertini 2.
# 
# python/pybertini/nag_algorithms/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/pybertini/nag_algorithms/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/pybertini/nag_algorithms/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Copyright(C) 2023 by Bertini2 Development Team
# 
#  See <http://www.gnu.org/licenses/> for a copy of the license, 
#  as well as COPYING.  Bertini2 is provided with permitted 
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
# 
#  silviana amethyst
#  UWEC
#  Spring 2023
# 





"""
nag_algorithms
"""


import _pybertini.nag_algorithms
from _pybertini.nag_algorithms import *

__all__ = dir(_pybertini.nag_algorithms)


DoublePrecisionTotalDegree = _pybertini.nag_algorithms.ZeroDimDoublePrecisionTotalDegree