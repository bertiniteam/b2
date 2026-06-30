# This file is part of Bertini 2.
# 
# python/bertini/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/bertini/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/bertini/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
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
#  Spring 2018
#
#  MPI-CBG 
#  2025
# 




"""
bertini -- Python bindings for Bertini 2.

This code is licensed under the GNU Public License, Version 3, with
additional clauses under section 7 as permitted, to protect the 
Bertini name.  See b2/licenses/ for a complete copy of the license,
and the licenses of software upon which Bertini depends.

See the source at https://github.com/bertiniteam/b2
"""


### this __init__.py is strongly inspired by that for GalSim 
### https://github.com/GalSim-developers/GalSim

from importlib.metadata import version
__version__ = version("bertini2")

import os
import sys

if sys.platform == "win32":
    from .windows_dll_manager import get_dll_paths, build_directory_manager
    _dll_manager = build_directory_manager()
    _dll_manager.__enter__()
    for p in get_dll_paths():
        _dll_manager.add_dll_directory(p)

# put stuff in the bertini namespace

from . import function_tree
from .function_tree import sin, cos, tan, asin, acos, atan, exp, log, sqrt
from .function_tree import canonicalize, monomial_order, MonomialOrder

from . import system
from . import tracking
from . import endgame
from . import parse
from . import container
from . import logging
from . import nag_algorithm
from . import random
from . import parallel

from ._pybertini import info

from . import multiprec

from . import linalg

# some convenience assignments
Variable = function_tree.symbol.Variable
variables = function_tree.variables
gather_variables = function_tree.gather_variables
VariableGroup = function_tree.VariableGroup
Named = function_tree.root.NamedExpression   # Named(expr, "a"): a user-named subexpression
System = system.System
default_precision = multiprec.default_precision

from ._calculus import jacobian



# https://stackoverflow.com/questions/44834/what-does-all-mean-in-python
# "a list of strings defining what symbols in a module will be exported when from <module> import * is used on the module"
__all__ = ['Variable','variables','gather_variables','VariableGroup','Named','system','System',
           'jacobian',
           'nag_algorithm','container','default_precision',
           'tracking','endgame','logging','function_tree','parse','multiprec','random','parallel',
           'linalg',
           'sin','cos','tan','asin','acos','atan','exp','log','sqrt',
           'canonicalize','monomial_order','MonomialOrder']




