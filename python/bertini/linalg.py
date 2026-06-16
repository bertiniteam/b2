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

"""
Express conditions over vectors and matrices of variables.

Vectors and matrices of variables are ordinary numpy object arrays of
:class:`~bertini.function_tree.symbol.Variable` nodes, so the usual numpy operators build
function-tree expressions::

    import numpy as np
    import bertini
    from bertini import linalg

    x   = linalg.variable_vector('x', 3)        # array([x0, x1, x2], dtype=object)
    lam = bertini.Variable('lam')
    A   = np.array([[2, 1, 0], [1, 3, 1], [0, 1, 4]])

    equations = A @ x - lam * x                  # an object array of 3 expressions

    sys = bertini.System()
    linalg.add_functions(sys, equations)         # add them all at once

**Coefficients must be exact.**  Python ints (including numpy integers) and bertini nodes
work directly inside numpy arithmetic, so an *integer* matrix needs nothing special.  For
any other coefficient use :func:`as_coefficients`, which accepts ``fractions.Fraction``,
exact strings (``'2.5'`` decimal or ``'3/4'`` rational), and bertini multiprecision values
(:class:`bertini.multiprec.Complex` / :class:`bertini.multiprec.Float`).  Python ``float``
and ``complex`` are **refused**: a 64-bit float literal would inject only ~16 correct
digits into the arbitrary-precision function tree, silently capping the precision of every
downstream computation.  Pass the value exactly instead.
"""

import numpy as np
from fractions import Fraction

from bertini.function_tree.symbol import Variable, Integer, Rational, Float
from bertini import multiprec as _mp

try:
    from bertini.function_tree import AbstractNode as _AbstractNode
except ImportError:  # pragma: no cover - fall back to the native module
    from bertini._pybertini.function_tree import AbstractNode as _AbstractNode

# bertini multiprecision value types (not function-tree nodes) that we accept as exact
# coefficients by wrapping them in a Float node.
_MP_VALUE_TYPES = tuple(
    t for t in (getattr(_mp, n, None) for n in ('Complex', 'Float', 'Int', 'Rational'))
    if isinstance(t, type)
)

__all__ = ['variable_vector', 'variable_matrix', 'coefficient', 'as_coefficients',
           'add_functions', 'add_linear_forms', 'add_linear', 'add_products_of_linears']


def variable_vector(name, n, start=0):
    """A length-``n`` numpy object array of variables ``name{start} .. name{start+n-1}``.

    >>> x = variable_vector('x', 3)            # doctest: +SKIP
    >>> [v.name for v in x]                    # doctest: +SKIP
    ['x0', 'x1', 'x2']
    """
    return np.array([Variable(f'{name}{i}') for i in range(start, start + n)], dtype=object)


def variable_matrix(name, rows, cols):
    """A ``rows`` x ``cols`` numpy object array of variables ``name_{i}_{j}``."""
    return np.array(
        [[Variable(f'{name}_{i}_{j}') for j in range(cols)] for i in range(rows)],
        dtype=object,
    )


def coefficient(value):
    """Coerce a single exact value into a function-tree coefficient node.

    Accepts bertini nodes (returned unchanged), ``int`` / numpy integers,
    ``fractions.Fraction``, exact strings (``'2.5'`` or ``'3/4'``), and bertini
    multiprecision values.  Raises ``TypeError`` for Python ``float`` / ``complex`` (and
    numpy floating types) -- see the module docstring for why.
    """
    if isinstance(value, _AbstractNode):
        return value
    if isinstance(value, bool):
        # bool is a subclass of int; refuse it explicitly rather than coercing True -> 1.
        raise TypeError("a bool is not a valid coefficient")
    if isinstance(value, (int, np.integer)):
        return Integer(int(value))
    if isinstance(value, Fraction):
        return Rational(str(value))            # 'p/q', or 'p' when the denominator is 1
    if isinstance(value, str):
        return Rational(value) if '/' in value else Float(value)
    if _MP_VALUE_TYPES and isinstance(value, _MP_VALUE_TYPES):
        return Float(value)
    if isinstance(value, (float, complex, np.floating, np.complexfloating)):
        raise TypeError(
            f"refusing to use the Python {type(value).__name__} {value!r} as a coefficient: "
            "a floating-point literal carries only ~16 digits and would cap the precision of "
            "the arbitrary-precision function tree.  Pass an exact value instead -- an int, a "
            "fractions.Fraction, an exact string such as '2.5' or '3/4', or a "
            "bertini.multiprec value."
        )
    raise TypeError(f"cannot use {type(value).__name__} as a coefficient")


def as_coefficients(array_like):
    """Coerce an array (or nested list) of exact values to an object array of coefficient nodes.

    Use this to bring a non-integer matrix or vector into the function tree::

        A = as_coefficients([['5/2', '1'], ['0', '3']])     # exact rationals
        A = as_coefficients(fraction_matrix)                # fractions.Fraction entries
        A = as_coefficients(mpfr_complex_matrix)            # bertini multiprecision values

    Integer arrays do not need this: numpy ``int * Variable`` already builds ``Integer``
    coefficients.  Refuses Python floats (see :func:`coefficient`).
    """
    arr = np.array(array_like, dtype=object)
    out = np.empty(arr.shape, dtype=object)
    flat_in = arr.reshape(-1)
    flat_out = out.reshape(-1)
    for i in range(flat_in.size):
        flat_out[i] = coefficient(flat_in[i])
    return out


def add_functions(system, expressions, basename=None):
    """Add a vector/array of expressions to a ``System`` as individual functions.

    ``expressions`` may be a numpy object array, a (possibly nested) list, or a single
    expression node.  With ``basename`` the functions are named ``basename0, basename1,
    ...``.  Returns the number of functions added.

        equations = A @ x - lam * x
        add_functions(sys, equations, basename='eig')
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
        if basename is None:
            system.add_function(e)
        else:
            system.add_function(e, f'{basename}{i}')
    return int(flat.size)


def _exact_to_mpfr(value):
    """Convert a single exact value to a multiprec.Complex (refusing Python floats).

    Mirrors :func:`coefficient`'s exact-only rule, but produces a multiprecision *value*
    (for a coefficient matrix) rather than a function-tree node.
    """
    if _MP_VALUE_TYPES and isinstance(value, _MP_VALUE_TYPES):
        return _mp.Complex(value)
    if isinstance(value, bool):
        raise TypeError("a bool is not a valid coefficient")
    if isinstance(value, (int, np.integer)):
        return _mp.Complex(str(int(value)))
    if isinstance(value, Fraction):
        return _mp.Complex(str(value.numerator)) / _mp.Complex(str(value.denominator))
    if isinstance(value, str):
        if '/' in value:
            num, den = value.split('/')
            return _mp.Complex(num) / _mp.Complex(den)
        return _mp.Complex(value)
    if isinstance(value, (float, complex, np.floating, np.complexfloating)):
        raise TypeError(
            f"refusing to use the Python {type(value).__name__} {value!r} as a coefficient: "
            "a floating-point literal would cap the precision of the block.  Pass an exact "
            "value -- an int, a fractions.Fraction, an exact string, or a bertini.multiprec "
            "value."
        )
    raise TypeError(f"cannot use {type(value).__name__} as a coefficient")


def add_linear_forms(system, coefficients):
    """Add an affine-linear-forms block f(x) = M [x;1] to ``system``.

    This is the efficient, first-class form of a stack of linear forms: instead of expanding
    each row into a scalar function-tree expression, the whole block is evaluated as one
    matrix-vector product (the C++ LinearFormsBlock).

    ``coefficients`` is an (m x num_vars+1) array/list of EXACT values -- one row per
    function, the trailing column being each form's constant term.  Entries are converted to
    an mpfr_complex coefficient matrix at the current default precision (set
    :func:`bertini.default_precision` higher beforehand if you want the block's master to
    carry more digits for adaptive-precision tracking).  Python floats are refused.

    Returns ``system`` for chaining.
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
    system.add_linear_forms_block(ncol - 1, M)
    return system


def add_linear(system, A, x, b=None):
    """Add the linear conditions ``A @ x + b == 0`` to ``system`` as a LinearFormsBlock.

    This is the "auto-target" for constant-coefficient linear forms: instead of expanding
    each row of ``A @ x`` into a scalar function-tree expression, the whole stack is added
    as one block evaluated by a single matrix-vector product.

    Parameters
    ----------
    system : the System to add to.  Its variable groups must already be set -- the block is
        built over the system's current variable ordering.
    A : an exact (m x n) coefficient matrix (array/list of lists).
    x : a length-n vector of the system's variables (numpy object array of Variable, e.g.
        from :func:`variable_vector`).  Each x[j] must already belong to a variable group of
        ``system``.
    b : optional length-m exact constant vector (default all zero).

    Coefficients must be exact (see :func:`coefficient`); Python floats are refused.

    The resulting linear-forms block is homogenization-aware, so it survives the
    homogenization the zero-dim solver performs: a system mixing polynomial functions with an
    ``add_linear`` block solves end to end (currently for a single affine variable group).

    Returns ``system`` for chaining.
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

    ordering = list(system.variable_ordering())
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

    system.add_linear_forms_block(num_vars, M)
    return system


def add_products_of_linears(system, factors):
    """Add a products-of-linear-forms block to ``system``: f_i(x) = prod_r ( c_{i,r} . [x;1] ).

    This is the first-class C++ ``ProductsOfLinearsBlock`` -- the same evaluation block the
    multihomogeneous start system uses -- evaluated as matrix-multiplies-then-row-products rather
    than as expanded scalar function-tree expressions.  It is the natural way to author your own
    start system as a product of linears: each factor ``c.[x;1] = 0`` is a hyperplane, so the start
    solutions are exact intersections of one hyperplane per function.

    Distinct from :func:`add_linear_forms`: a single ``c.[x;1]`` is one degree-1 linear form (a
    ``LinearFormsBlock``); a *product* of k such factors is one function of **degree k**.

    Parameters
    ----------
    system : the System to add to.
    factors : a list with one entry per function.  Entry i is an exact (k_i x (num_vars+1)) matrix
        (array/list of lists): one row per linear factor, the trailing column being that factor's
        constant term.  Different functions may have different numbers of factors, but every matrix
        must have the same number of columns (num_vars+1).

    Coefficients must be exact (see :func:`coefficient`); Python floats are refused.

    Works across multiple variable groups: the coefficient columns follow the system's variable
    ordering and the trailing column is the affine constant.  For a *projective* (homogeneous)
    variable group, write homogeneous factors -- give the constant column as ``0`` -- since the
    block is evaluated as authored (it does not re-homogenize).

    Returns ``system`` for chaining.
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
    system.add_products_of_linears_block(ncol - 1, mats)
    return system
