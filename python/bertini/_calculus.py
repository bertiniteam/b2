# This file is part of Bertini 2.
#
# bertini/_calculus.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# bertini/_calculus.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with bertini/_calculus.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""Symbolic calculus over function trees (top-level ``bertini`` helpers)."""

import numpy as _np

from bertini._pybertini import system as _pybsys
from bertini._pybertini.function_tree import AbstractNode as _AbstractNode

__all__ = ['jacobian']


def jacobian(functions, variables):
    """The symbolic Jacobian of ``functions`` with respect to ``variables``.

    ``J[i, j]`` is the partial derivative of ``functions[i]`` with respect to ``variables[j]``,
    returned as a 2-D numpy object array of function-tree nodes -- always 2-D, even for a single
    function, so it is ready to ``numpy.vstack`` onto a coefficient row and ``@`` a vector of
    variables.  This is purely symbolic (it builds expression trees via differentiation); it does
    not evaluate.  For a whole :class:`~bertini.System` use :meth:`bertini.System.jacobian`, which
    is homogenization-aware.

    Parameters
    ----------
    functions : node or iterable of nodes
        A single expression node or a sequence of them (a list or numpy object array).
    variables : iterable of Variable
        The variables to differentiate with respect to (a list, a numpy object array such as
        :func:`bertini.linalg.variable_vector`, or a :class:`~bertini.container.VariableGroup`).

    Examples
    --------
    >>> import bertini as pb                                   # doctest: +SKIP
    >>> x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')   # doctest: +SKIP
    >>> J = pb.jacobian([x*y, x - z], [x, y, z])              # doctest: +SKIP
    >>> J.shape                                                # doctest: +SKIP
    (2, 3)
    """
    if isinstance(functions, _AbstractNode):
        functions = [functions]
    else:
        functions = list(functions)
    variables = list(variables)

    rows = _pybsys.symbolic_jacobian(functions, variables)
    arr = _np.array(rows, dtype=object)
    if arr.ndim != 2:   # empty function list collapses to 1-D; restore the (m, n) shape
        arr = arr.reshape(len(functions), len(variables))
    return arr
