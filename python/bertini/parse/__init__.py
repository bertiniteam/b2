# This file is part of Bertini 2.
# 
# python/bertini/parse/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/bertini/parse/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/bertini/parse/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Copyright(C) Bertini2 Development Team
# 
#  See <http://www.gnu.org/licenses/> for a copy of the license, 
#  as well as COPYING.  Bertini2 is provided with permitted 
#  additional terms in the b2/licenses/ directory.

"""
Parsing functions, taking strings and producing various other things
"""

from bertini._pybertini import parse as _pybparse
from bertini._pybertini.parse import *

# The native parser accepts a full Bertini 1 classic file -- comments, a leading CONFIG
# section and the INPUT/END; wrappers are all handled in C++ (#407, #396) -- and the classic
# writer spells complex constants the one way Bertini 1 and this parser read, (re+im*I), so
# parse.system(sys.to_classic_input()) round-trips with nothing left for Python to patch.
system = _pybparse.system


__all__ = dir(_pybparse)

