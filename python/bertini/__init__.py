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
bertini -- Python bindings for Bertini 2
========================================

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
is_distinct_up_to = multiprec.is_distinct_up_to   # tolerance point-inequality, infinity norm (#304)

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
from ._randomize import randomize
# dense linear algebra for the mp types (solve / LU), backed by eigenpy's decompositions
from . import linalg
from .random import random_matrix
from .random import random_vector, random_real, random_complex

# exact-coefficient coercion at the top level (was bertini.linalg.coefficient / as_coefficients)
from ._coefficients import coefficient, coefficients

# the casual records surface: solve / save / load over the structured output directory
from .records import (solve, save, load, annotate, solutions_of, provenance,
                      recording, records_dir, runs, tracks, provenance_graph,
                      plot_chain, Solution, SolveResult)
from . import records

# attach the friendly system-building methods and Slice.from_coefficients (was bertini.linalg.*)
from . import _system_ops as _system_ops
_system_ops.install(system.System)
from . import _slice_ops as _slice_ops
_slice_ops.install(nag_algorithm.Slice)

# --- sympy interop (#295) ----------------------------------------------------------------------
# Make every function-tree node auto-convert to sympy (the `_sympy_` protocol), so sympy.sympify(node),
# sympy.Matrix(array_of_nodes), and sympy.det(J) work directly.  sympy is an optional dependency; the
# import only happens when the method is actually called (so importing bertini never needs sympy).
def _node_to_sympy(self):
    from bertini.sympy_bridge import to_sympy
    return to_sympy(self)


symbolics.AbstractNode._sympy_ = _node_to_sympy


# `bertini.sympy_bridge` is exposed lazily: accessing it imports the module (which raises a helpful
# ImportError if sympy is missing) without making sympy a hard dependency of `import bertini`.
def __getattr__(name):
    if name == 'sympy_bridge':
        import importlib
        return importlib.import_module('bertini.sympy_bridge')   # import_module avoids re-entering __getattr__
    raise AttributeError("module {!r} has no attribute {!r}".format(__name__, name))


# --- numpy-friendly elementwise helpers for the mp types (#298, #301) --------------------------
# Vectorized real/imag/abs/conj/round/sum/norm/is_real that stay mp-native, riding the native numpy
# ufunc loops on mp-dtype arrays and working element-wise on lists/mixed input (see the
# "Multiprecision numbers and NumPy" docs page, docs/source/numpy.rst).  These are attributes of the
# top-level module (bertini.abs, bertini.real, ...).  The builtin-shadowing names (abs, round, sum) are
# deliberately kept OUT of __all__, so `from bertini import *` never clobbers the Python builtins.
from . import _numpy_helpers as _numpy_helpers

# make np.real/np.imag/np.angle raise (instead of silently returning wrong values)
# on plain mp-complex arrays -- numpy has no user-dtype hook for component access,
# and a crash is better than incorrect values.  See _numpy_guard for the story.
from . import _numpy_guard as _numpy_guard
_numpy_guard.install()

real = _numpy_helpers.real
imag = _numpy_helpers.imag
conj = _numpy_helpers.conj
arg = _numpy_helpers.arg
norm = _numpy_helpers.norm
is_real = _numpy_helpers.is_real
abs = _numpy_helpers.abs        # noqa: A001  (bertini.abs; not exported via *)
round = _numpy_helpers.round    # noqa: A001
sum = _numpy_helpers.sum        # noqa: A001

# --- the one-stop math vocabulary ---------------------------------------------------------------
# `from bertini.operators import *` gives sin/cos/.../abs/arg/real/imag/... that work on symbolic
# expressions AND numbers AND numpy containers alike, dispatching per argument.  The top-level
# elementary functions are rebound to the polymorphic versions (a strict superset of the symbolic
# ones bound above): bertini.sin(x) works for a Variable, a real_mp, or an array.
from . import operators
from .operators import (sin, cos, tan, asin, acos, atan, exp, log, sqrt,   # noqa: F811
                        sinh, cosh, tanh, asinh, acosh, atanh)



# https://stackoverflow.com/questions/44834/what-does-all-mean-in-python
# "a list of strings defining what symbols in a module will be exported when from <module> import * is used on the module"
__all__ = ['solve','save','load','annotate','solutions_of','provenance','recording','records_dir','runs','tracks','provenance_graph','plot_chain','Solution','SolveResult','records',
           'Variable','variables','gather_variables','VariableGroup','Named','system','System',
           'jacobian','randomize','linalg','random_matrix','random_vector','random_real','random_complex','coefficient','coefficients',
           'complex_mp','real_mp','int_mp','rational_mp',
           'nag_algorithm','default_precision','is_distinct_up_to',
           'real','imag','conj','arg','norm','is_real',
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
           'sinh','cosh','tanh','asinh','acosh','atanh',
           'canonicalize','monomial_order']




