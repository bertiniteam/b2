# This file is part of Bertini 2.
#
# python/bertini/symbolics/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/symbolics/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/bertini/symbolics/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""The symbolic expression system for building polynomial systems -- a FLAT namespace.

Renamed from ``function_tree`` (the underlying C++ module keeps that name).  Everything lives
directly here, with no ``.symbol`` / ``.operator`` / ``.root`` sub-levels:

* symbols -- :class:`Variable`, :class:`Complex`, :class:`Integer`, :class:`Rational`, ``E``, ``Pi``
* operator node types -- :class:`Sum`, :class:`Mult`, :class:`Power`, :class:`Sqrt`, ... (for
  inspecting/walking expression trees)
* roots -- :class:`NamedExpression`

plus :func:`variables`, :func:`sqrt`, and :class:`VariableGroup`.
"""

from bertini._pybertini import function_tree as _pybft
from bertini._pybertini.container import VariableGroup

# Flatten: pull the top-level utilities and every node type from symbol/operator/root up to here.
from bertini._pybertini.function_tree import *          # noqa: F401,F403  (sin, cos, canonicalize, ...)
from bertini._pybertini.function_tree.symbol import *   # noqa: F401,F403  (Variable, Complex, E, Pi, ...)
from bertini._pybertini.function_tree.operator import * # noqa: F401,F403  (Sum, Mult, Power, ...)
from bertini._pybertini.function_tree.root import *     # noqa: F401,F403  (NamedExpression, ...)

# The top-level `import *` dragged the C++ sub-submodule objects in as names (symbol/operator/root);
# drop them so this namespace stays flat.  Also drop the abstract bases -- they exist in C++ only to
# make isinstance-based tree walking possible, and are not part of the public surface.
_HIDDEN = {'symbol', 'operator', 'root',
           'AbstractNode', 'AbstractSymbol', 'AbstractNamedSymbol', 'AbstractNumber', 'AbstractOp',
           # constants factories superseded by the bertini.E / bertini.Pi / bertini.I constants
           'make_e', 'make_i', 'make_pi',
           # Differential is an internal differentiation artifact, not a user-facing symbol
           'Differential'}
for _n in _HIDDEN:
    globals().pop(_n, None)
del _n

from bertini._pybertini.function_tree.operator import Sqrt as _Sqrt


def sqrt(x):
    """Symbolic square-root operator."""
    return _Sqrt(x)


VariableGroup.__str__ = lambda vg: '[{}]'.format(','.join([str(v) for v in vg]))


def variables(base, indices, fmt='{base}{index}'):
    """Make a list of integer-indexed Variables.

    ::

        base    -- name prefix, e.g. 'x'
        indices -- an int n (shorthand for range(n)) or any iterable of ints
        fmt     -- str.format template using {base} and {index};
                   default '{base}{index}' gives x0, x1, x2, ...

    Returns a ``list[Variable]``.  Wrap in a VariableGroup if desired::

        pb.VariableGroup(pb.variables('x', 5))
    """
    from bertini._pybertini.function_tree.symbol import Variable
    if isinstance(indices, int):
        indices = range(indices)
    return [Variable(fmt.format(base=base, index=i)) for i in indices]


__all__ = list(dict.fromkeys(
    [n for n in dir(_pybft) if n not in _HIDDEN]
    + [n for n in dir(_pybft.symbol) if n not in _HIDDEN]
    + [n for n in dir(_pybft.operator) if n not in _HIDDEN]
    + [n for n in dir(_pybft.root) if n not in _HIDDEN]
    + ['VariableGroup', 'variables', 'sqrt']))
