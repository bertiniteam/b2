# This file is part of Bertini 2.
#
# bertini/linalg.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# bertini/linalg.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with bertini/linalg.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""DEPRECATED module -- kept as a thin shim while call sites migrate.

``bertini.linalg`` was never linear algebra; it dissolved into the everyday API:

* ``linalg.add_functions/add_linear/add_linear_forms/add_products_of_linears/randomize/
  add_slices_as_products(sys, ...)``  ->  ``sys.<same>(...)``  (methods on :class:`~bertini.System`)
* ``linalg.coefficient`` / ``linalg.as_coefficients``  ->  ``bertini.coefficient`` /
  ``bertini.coefficients``
* ``linalg.slice_from_coefficients(coeffs, vars)``  ->  ``bertini.Slice.from_coefficients(coeffs, vars)``
* ``linalg.variable_vector('v', n)``  ->  ``bertini.variables('v', n)`` (list) or a
  :class:`~bertini.VariableGroup`

The functions below simply forward to those homes, so old code keeps working.  This module will be
removed once the tutorials and tests are swept.
"""

import numpy as np

from bertini.symbolics import Variable
# re-exported so the exact-coefficient helpers keep their old import path during the transition
from bertini._coefficients import (coefficient, coefficients as as_coefficients,
                                    _exact_to_mpfr, _coerce_mpfr_matrix, _MP_VALUE_TYPES)

__all__ = ['variable_vector', 'variable_matrix', 'coefficient', 'as_coefficients',
           'add_functions', 'add_linear_forms', 'add_linear', 'add_products_of_linears',
           'randomize', 'slice_from_coefficients', 'add_slices_as_products']


def variable_vector(name, n, start=0):
    """DEPRECATED: a length-``n`` numpy object array of variables ``name{start}..``.

    Prefer :func:`bertini.variables` (returns a list) or a :class:`~bertini.VariableGroup`.
    """
    return np.array([Variable(f'{name}{i}') for i in range(start, start + n)], dtype=object)


def variable_matrix(name, rows, cols):
    """DEPRECATED: a ``rows`` x ``cols`` numpy object array of variables ``name_{i}_{j}``."""
    return np.array(
        [[Variable(f'{name}_{i}_{j}') for j in range(cols)] for i in range(rows)],
        dtype=object,
    )


def add_functions(system, expressions):
    """DEPRECATED: use ``system.add_functions(expressions)``."""
    return system.add_functions(expressions)


def add_linear_forms(system, coefficients):
    """DEPRECATED: use ``system.add_linear_forms(coefficients)``."""
    return system.add_linear_forms(coefficients)


def add_linear(system, A, x, b=None):
    """DEPRECATED: use ``system.add_linear(A, x, b)``."""
    return system.add_linear(A, x, b)


def add_products_of_linears(system, factors):
    """DEPRECATED: use ``system.add_products_of_linears(factors)``."""
    return system.add_products_of_linears(factors)


def randomize(system, matrix=None):
    """DEPRECATED: use ``system.randomize(matrix)``."""
    return system.randomize(matrix)


def add_slices_as_products(system, slices):
    """DEPRECATED: use ``system.add_slices_as_products(slices)``."""
    return system.add_slices_as_products(slices)


def slice_from_coefficients(coefficients, variables, homogeneous=False):
    """DEPRECATED: use ``bertini.Slice.from_coefficients(coefficients, variables, homogeneous)``."""
    from bertini.nag_algorithm import Slice
    return Slice.from_coefficients(coefficients, variables, homogeneous)
