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

# from .system import System
import pybertini.system as system
# from .system import start_system


import pybertini.function_tree as function_tree

# #bring some things to the front namespace, because they are nested several deep
# from .function_tree import Variable
# from .function_tree import VariableGroup

from .multiprec import multiprec

import pybertini.tracking as tracking
import pybertini.endgame as endgame


__all__ = ['tracking','endgame','multiprec','function_tree','system','start_system']

# some convenience assignments
Variable = function_tree.symbol.Variable
VariableGroup = function_tree.VariableGroup
System = system.System
