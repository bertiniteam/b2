# This file is part of Bertini 2.
# 
# python/pybertini/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/pybertini/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/pybertini/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
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
PyBertini -- Python bindings for Bertini2.

This code is licensed under the GNU Public License, Version 3, with
additional clauses under section 7 as permitted, to protect the 
Bertini name.  See b2/licenses/ for a complete copy of the license,
and the licenses of software upon which Bertini depends.

See the source at https://github.com/bertiniteam/b2
"""


### this __init__.py is strongly inspired by that for GalSim 
### https://github.com/GalSim-developers/GalSim

from ._version import __version__, __version_info__
version = __version__


# put stuff in the pybertini namespace
import pybertini.system as system
import pybertini.function_tree as function_tree
import pybertini.minieigen as minieigen
import pybertini.doubleprec as doubleprec
import pybertini.multiprec as multiprec
import pybertini.tracking as tracking
import pybertini.endgame as endgame
import pybertini.parse as parse




# some convenience assignments
Variable = function_tree.symbol.Variable
VariableGroup = function_tree.VariableGroup
System = system.System
default_precision = multiprec.default_precision

__all__ = ['tracking','endgame','multiprec','function_tree','system','default_precision','System','Variable','VariableGroup']