# This file is part of Bertini 2.
#
# python/bertini/_system_ops.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/_system_ops.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/bertini/_system_ops.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""Friendly System-building methods, attached to the bound ``System`` class.

These are the everyday ways to add blocks of conditions to a :class:`~bertini.System`: pass numpy
arrays / lists of expressions or exact coefficient matrices and they do the flattening, exact-value
coercion, and homogenization-aware block construction, delegating to the low-level C++ methods
(``add_function``, ``add_linear_forms_block``, ``add_products_of_linears_block``, ``randomize``).

They replace the old free functions in :mod:`bertini.linalg` (now a deprecated shim).  ``install``
captures the native ``randomize`` before overriding it, so the friendly ``randomize`` (which accepts
an optional coefficient matrix) can still reach the C++ implementation without recursing.
"""

import numpy as np

from bertini import multiprec as _mp
from bertini._coefficients import _exact_to_mpfr, _coerce_mpfr_matrix
from bertini._pybertini.function_tree import AbstractNode as _AbstractNode

_native = {}


def add_functions(self, expressions):
    """Add a vector/array of expressions as individual functions; returns the count added.

    ``expressions`` may be a numpy object array, a (possibly nested) list, or a single
    expression node::

        sys.add_functions(A @ x - lam * x)
    """
    if isinstance(expressions, _AbstractNode):
        expressions = [expressions]
    flat = np.array(expressions, dtype=object).reshape(-1)
    for i in range(flat.size):
        e = flat[i]
        if not isinstance(e, _AbstractNode):
            raise TypeError(
                f"expression {i} is a {type(e).__name__}, not a function-tree node; "
                "did a coefficient fail to combine with a variable?"
            )
        self.add_function(e)
    return int(flat.size)


def add_linear_forms(self, coefficients):
    """Add an affine-linear-forms block f(x) = M [x;1] (one LinearFormsBlock).

    ``coefficients`` is an (m x num_vars+1) array/list of EXACT values -- one row per function,
    the trailing column being each form's constant term.  Python floats are refused.  Returns self.
    """
    rows = [list(r) for r in coefficients]
    if not rows:
        raise ValueError("coefficients must have at least one row")
    ncol = len(rows[0])
    M = np.empty((len(rows), ncol), dtype=_mp.Complex)
    for i, row in enumerate(rows):
        if len(row) != ncol:
            raise ValueError("coefficient matrix is ragged (rows of differing length)")
        for j, entry in enumerate(row):
            M[i, j] = _exact_to_mpfr(entry)
    self.add_linear_forms_block(ncol - 1, M)
    return self


def add_linear(self, A, x, b=None):
    """Add the linear conditions ``A @ x + b == 0`` as one LinearFormsBlock.

    ``A`` is an exact (m x n) coefficient matrix; ``x`` a length-n vector of the system's variables
    (each already in a variable group); ``b`` an optional length-m exact constant vector (default
    zero).  Coefficients must be exact.  Returns self.  Homogenization-aware.
    """
    rows = [list(r) for r in A]
    if not rows:
        raise ValueError("A must have at least one row")
    m = len(rows)
    n = len(rows[0])
    xs = list(x)
    if len(xs) != n:
        raise ValueError(f"A has {n} columns but x has {len(xs)} variables")
    if b is not None:
        b = list(b)
        if len(b) != m:
            raise ValueError(f"b has length {len(b)} but A has {m} rows")

    ordering = list(self.variable_ordering())
    col_of = {v.name: i for i, v in enumerate(ordering)}
    num_vars = len(ordering)

    M = np.zeros((m, num_vars + 1), dtype=_mp.Complex)
    for i in range(m):
        if len(rows[i]) != n:
            raise ValueError("A is ragged (rows of differing length)")
        for j in range(n):
            name = xs[j].name
            if name not in col_of:
                raise ValueError(
                    f"variable {name!r} is not in the system's variable ordering; "
                    "add its variable group before calling add_linear"
                )
            M[i, col_of[name]] = _exact_to_mpfr(rows[i][j])
        if b is not None:
            M[i, num_vars] = _exact_to_mpfr(b[i])

    self.add_linear_forms_block(num_vars, M)
    return self


def add_products_of_linears(self, factors):
    """Add a products-of-linear-forms block: f_i(x) = prod_r ( c_{i,r} . [x;1] ).

    ``factors`` is a list with one entry per function; entry i is an exact (k_i x num_vars+1)
    matrix (one row per linear factor, trailing column the factor's constant).  A product of k
    factors is one function of degree k.  Coefficients must be exact.  Returns self.
    """
    mats = []
    ncol = None
    for fi, factor_matrix in enumerate(factors):
        rows = [list(r) for r in factor_matrix]
        if not rows:
            raise ValueError(f"function {fi} has no linear factors")
        if ncol is None:
            ncol = len(rows[0])
        M = np.empty((len(rows), ncol), dtype=_mp.Complex)
        for i, row in enumerate(rows):
            if len(row) != ncol:
                raise ValueError("products-of-linears coefficient matrices are ragged: every factor "
                                 "row across every function must have num_vars+1 columns")
            for j, entry in enumerate(row):
                M[i, j] = _exact_to_mpfr(entry)
        mats.append(M)
    if ncol is None:
        raise ValueError("factors must contain at least one function")
    self.add_products_of_linears_block(ncol - 1, mats)
    return self


def randomize(self, matrix=None):
    """Randomize an overdetermined System down to a square one, returning a NEW system.

    ``matrix`` is an optional exact coefficient matrix R (one row per randomized function, one
    column per natural function); when ``None`` bertini builds a generic R.  The original is not
    mutated.  Coefficients must be exact.
    """
    if matrix is None:
        return _native['randomize'](self)
    rows = [list(r) for r in matrix]
    if not rows:
        raise ValueError("randomization matrix must have at least one row")
    ncol = len(rows[0])
    M = np.empty((len(rows), ncol), dtype=_mp.Complex)
    for i, row in enumerate(rows):
        if len(row) != ncol:
            raise ValueError("randomization matrix is ragged (rows of differing length)")
        for j, entry in enumerate(row):
            M[i, j] = _exact_to_mpfr(entry)
    return _native['randomize'](self, M)


def add_slices_as_products(self, slices):
    """Add one products-of-linears function per slice (the regeneration bridge).

    Each :class:`Slice` becomes a single function that is the product of its linear forms.  Slice
    coefficient columns must follow the system's variable ordering.  Returns self.
    """
    factors = [np.asarray(s.coefficients()) for s in slices]
    return add_products_of_linears(self, factors)


def install(System):
    """Attach the friendly methods to the bound ``System`` class (idempotent).

    Captures the native ``randomize`` first so the friendly override can still reach it.
    """
    if getattr(System, "_b2_system_ops_installed", False):
        return
    _native['randomize'] = System.randomize
    System.add_functions = add_functions
    System.add_linear_forms = add_linear_forms
    System.add_linear = add_linear
    System.add_products_of_linears = add_products_of_linears
    System.randomize = randomize
    System.add_slices_as_products = add_slices_as_products
    System._b2_system_ops_installed = True
