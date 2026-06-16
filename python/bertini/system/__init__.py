# coding : utf-8
#
# This file is part of Bertini 2.
# 
# python/bertini/system/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/bertini/system/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/bertini/system/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
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





"""
Provides utilities for working with systems of functions -- polynomials are intended, although you can work with functions involving things like trig functions, arbitrary powers, etc.

Making a new `System` is the starting point you want, probably some of these things:

::

	sys = bertini.system.System()
	sys.add_function(...)
	sys.add_variable_group(...)
	
	x = sys.dehomogenize_point(z)

"""

from bertini._pybertini import system as _pybsys

from bertini._pybertini.system import *

# --- unified builder: System.add(*objects) ---
# One fluent verb instead of remembering add_function / add_variable_group: dispatch each
# argument by type.  Returns self for chaining.  Additive -- the explicit methods still work.
import numpy as _np
from bertini._pybertini.function_tree import AbstractNode as _AbstractNode
from bertini._pybertini.container import VariableGroup as _VariableGroup


def _system_add(self, *objects):
    """Add functions and/or variable groups to the System, dispatched by type.

    Each argument may be:
      * a function-tree expression (e.g. ``x**2 + y - 1``, or a lone ``Variable``) -> added
        as a function;
      * a :class:`~bertini.container.VariableGroup` -> added as an affine variable group;
      * a numpy array / list / tuple of the above -> each element is added (so
        ``sys.add(A @ x - lam*x)``, ``sys.add([f, g])``, and ``sys.add(grp, f, g)`` work).

    Projective groups still use :meth:`add_hom_variable_group`; structured blocks use
    :func:`bertini.linalg.add_linear`.  Returns ``self`` for chaining.
    """
    for obj in objects:
        if isinstance(obj, _AbstractNode):
            self.add_function(obj)
        elif isinstance(obj, _VariableGroup):
            self.add_variable_group(obj)
        elif isinstance(obj, (_np.ndarray, list, tuple)):
            for elt in (obj.ravel() if isinstance(obj, _np.ndarray) else obj):
                self.add(elt)
        else:
            raise TypeError(
                f"System.add does not know how to add a {type(obj).__name__}; pass a "
                "function-tree expression, a VariableGroup, or an array/list of those"
            )
    return self


System.add = _system_add

# Override C++ submodule reference with the Python wrapper (which has AbstractStartSystem removed).
# Can't use 'from . import start_system': the star import already set that name to the C++ submodule.
import importlib as _importlib
start_system = _importlib.import_module('bertini.system.start_system')
del _importlib

__all__ = dir(_pybsys)
__all__.extend(['start_system'])
