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

import sys

if sys.platform == "win32":
    from .windows_dll_manager import get_dll_paths, build_directory_manager
    _dll_manager = build_directory_manager()
    _dll_manager.__enter__()
    for p in get_dll_paths():
        _dll_manager.add_dll_directory(p)

del sys  # used only for the platform check above; don't leak it into the bertini.* namespace

# put stuff in the bertini namespace

from . import symbolics
from .symbolics import sin, cos, tan, asin, acos, atan, exp, log, sqrt
from .symbolics import canonicalize, monomial_order, MonomialOrder

from . import system
from . import tracking
from . import endgame
from . import parse
from . import logging
from . import nag_algorithm
from . import random
from . import parallel

from . import multiprec

# some convenience assignments
Variable = symbolics.Variable
variables = symbolics.variables
gather_variables = symbolics.gather_variables
VariableGroup = symbolics.VariableGroup
Named = symbolics.NamedExpression            # Named(expr, "a"): a user-named subexpression
System = system.System
default_precision = multiprec.default_precision

# the multiprecision number types, hoisted to the top level (they also live in bertini.multiprec)
from .multiprec import complex_mp, real_mp, int_mp, rational_mp

# symbolic constants, ready to drop straight into expressions (bertini.E, bertini.Pi, bertini.I)
E = symbolics.E()
Pi = symbolics.Pi()
I = symbolics.Complex(0, 1)                     # imaginary unit -- no dedicated node, a complex leaf

# the everyday classes, hoisted to the top level for tab-completion.  They still live in their
# submodules (nag_algorithm.*, tracking.*); this just spares users the deep path.
from .nag_algorithm import ZeroDimSolver, HomotopySolver, SolutionPathCollector, Slice, StartSystemType
from .tracking import (AMPTracker, DoublePrecisionTracker, MultiplePrecisionTracker,
                       SuccessCode, Predictor)

from ._calculus import jacobian
from .random import random_matrix

# exact-coefficient coercion at the top level (was bertini.linalg.coefficient / as_coefficients)
from ._coefficients import coefficient, coefficients

# the casual records surface: solve / save / load over the structured output directory
from .records import solve, save, load, records_dir, Solution, SolveResult
from . import records

# attach the friendly system-building methods and Slice.from_coefficients (was bertini.linalg.*)
from . import _system_ops as _system_ops
_system_ops.install(system.System)
from . import _slice_ops as _slice_ops
_slice_ops.install(nag_algorithm.Slice)

from . import operators                          # `from bertini.operators import *` -> just the math ops



# https://stackoverflow.com/questions/44834/what-does-all-mean-in-python
# "a list of strings defining what symbols in a module will be exported when from <module> import * is used on the module"
__all__ = ['solve','save','load','records_dir','Solution','SolveResult','records',
           'Variable','variables','gather_variables','VariableGroup','Named','system','System',
           'jacobian','random_matrix','coefficient','coefficients',
           'complex_mp','real_mp','int_mp','rational_mp',
           'nag_algorithm','default_precision',
           'tracking','endgame','logging','symbolics','parse','multiprec','random','parallel',
           'operators',
           # everyday classes hoisted to the top level
           'ZeroDimSolver','HomotopySolver','SolutionPathCollector','Slice',
           'AMPTracker','DoublePrecisionTracker','MultiplePrecisionTracker',
           # enums at the root
           'SuccessCode','Predictor','MonomialOrder','StartSystemType',
           # symbolic constants
           'E','Pi','I',
           'sin','cos','tan','asin','acos','atan','exp','log','sqrt',
           'canonicalize','monomial_order']




