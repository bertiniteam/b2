# This file is part of Bertini 2.
# 
# python/bertini/function_tree/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/bertini/function_tree/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/bertini/function_tree/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
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



from bertini._pybertini import function_tree as _pybft
from bertini._pybertini.container import VariableGroup

from bertini._pybertini.function_tree import *

del AbstractNode

# Override C++ submodule references with their Python wrappers (which have abstracts removed).
# Can't use 'from . import symbol, root' here: the star import already put the C++ submodules
# under those names, and Python skips the subpackage import when the name is already defined.
import importlib as _importlib
symbol = _importlib.import_module('bertini.function_tree.symbol')
root = _importlib.import_module('bertini.function_tree.root')
del _importlib

from bertini._pybertini.function_tree.operator import Sqrt as _Sqrt

def sqrt(x):
    """Symbolic square-root operator."""
    return _Sqrt(x)


VariableGroup.__str__ = lambda vg: '[{}]'.format( ','.join([str(v) for v in vg]) )


def variables(base, indices, fmt='{base}{index}'):
    """Make a list of integer-indexed Variables.

    base    -- name prefix, e.g. 'x'
    indices -- an int n (shorthand for range(n)) or any iterable of ints
    fmt     -- str.format template using {base} and {index};
               default '{base}{index}' gives x0, x1, x2, ...

    Returns a list[Variable].  Wrap in a VariableGroup if desired:
        pb.VariableGroup(pb.variables('x', 5))
    """
    from bertini._pybertini.function_tree.symbol import Variable
    if isinstance(indices, int):
        indices = range(indices)
    return [Variable(fmt.format(base=base, index=i)) for i in indices]


_ABSTRACT = {'AbstractNode'}
__all__ = [n for n in dir(_pybft) if n not in _ABSTRACT] + ['variables', 'sqrt']

