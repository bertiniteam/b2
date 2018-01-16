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




# some convenience assignments
Variable = function_tree.symbol.Variable
VariableGroup = function_tree.VariableGroup
System = system.System
default_precision = multiprec.default_precision

__all__ = ['tracking','endgame','multiprec','function_tree','system','default_precision','System','Variable','VariableGroup']