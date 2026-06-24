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

from bertini.function_tree import VariableGroup as _VariableGroup

__all__ = ['variable_vector', 'variable_matrix', 'coefficient', 'as_coefficients',
           'add_functions', 'add_linear_forms', 'add_linear', 'add_products_of_linears',
           'randomize', 'slice_from_coefficients', 'add_slices_as_products']


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

    Examples
    --------
    There are four ways to spell an *exact* non-integer coefficient (a Python ``float`` is not one
    of them -- it would silently cap precision)::

        >>> from fractions import Fraction
        >>> from bertini import linalg, multiprec
        >>> _ = linalg.coefficient('2.5')                       # exact decimal string
        >>> _ = linalg.coefficient('3/4')                       # exact rational string
        >>> _ = linalg.coefficient(Fraction(3, 4))              # a fractions.Fraction
        >>> _ = linalg.coefficient(multiprec.Complex('0.1'))    # a full-precision multiprec value

    A Python ``float`` is refused: ``linalg.coefficient(0.1)`` raises ``TypeError`` (0.1 is not
    exactly 1/10 in binary, so it would inject only ~16 correct digits).  An exact *complex*
    coefficient (e.g. an off-axis gamma) is a ``multiprec.Complex`` with real and imaginary parts:
    ``linalg.coefficient(multiprec.Complex('0.6', '0.8'))``.
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


def add_functions(system, expressions):
    """Add a vector/array of expressions to a ``System`` as individual functions.

    ``expressions`` may be a numpy object array, a (possibly nested) list, or a single
    expression node.  Returns the number of functions added.

        equations = A @ x - lam * x
        add_functions(sys, equations)
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
        system.add_function(e)
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
    system
        the System to add to.  Its variable groups must already be set -- the block is
        built over the system's current variable ordering.
    A
        an exact (m x n) coefficient matrix (array/list of lists).
    x
        a length-n vector of the system's variables (numpy object array of Variable, e.g.
        from :func:`variable_vector`).  Each x[j] must already belong to a variable group of
        ``system``.
    b
        optional length-m exact constant vector (default all zero).

    Returns
    -------
    System
        ``system``, for chaining.

    Notes
    -----
    Coefficients must be exact (see :func:`coefficient`); Python floats are refused.

    The resulting linear-forms block is homogenization-aware, so it survives the
    homogenization the zero-dim solver performs: a system mixing polynomial functions with an
    ``add_linear`` block solves end to end (currently for a single affine variable group).

    Examples
    --------
    The unit circle meeting the line ``2x + y = 1`` -- a polynomial function plus one
    constant-coefficient linear condition (a degree-1 block)::

        >>> import numpy as np, bertini
        >>> from bertini import linalg
        >>> x, y = bertini.Variable('x'), bertini.Variable('y')
        >>> S = bertini.System()
        >>> S.add_variable_group(bertini.VariableGroup([x, y]))
        >>> S.add_function(x*x + y*y - 1)
        >>> _ = linalg.add_linear(S, np.array([[2, 1]]), np.array([x, y]), [-1])  # 2x + y - 1
        >>> list(S.degrees())
        [2, 1]
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
    system
        the System to add to.
    factors
        a list with one entry per function.  Entry i is an exact (k_i x (num_vars+1)) matrix
        (array/list of lists): one row per linear factor, the trailing column being that factor's
        constant term.  Different functions may have different numbers of factors, but every matrix
        must have the same number of columns (num_vars+1).

    Returns
    -------
    System
        ``system``, for chaining.

    Notes
    -----
    Coefficients must be exact (see :func:`coefficient`); Python floats are refused.

    Works across multiple variable groups: the coefficient columns follow the system's variable
    ordering and the trailing column is the affine constant.  For a *projective* (homogeneous)
    variable group, write homogeneous factors -- give the constant column as ``0`` -- since the
    block is evaluated as authored (it does not re-homogenize).

    Examples
    --------
    Two functions, each a product of two linear forms.  A product's degree is its number of
    factors, so these are degree 2::

        >>> import bertini
        >>> from bertini import linalg
        >>> x, y = bertini.Variable('x'), bertini.Variable('y')
        >>> S = bertini.System()
        >>> S.add_variable_group(bertini.VariableGroup([x, y]))
        >>> _ = linalg.add_products_of_linears(S, [
        ...     [[1, 0, -1], [1, 0, 1]],     # (x - 1)(x + 1)
        ...     [[0, 1, -1], [0, 1, -2]],    # (y - 1)(y - 2)
        ... ])
        >>> list(S.degrees())
        [2, 2]

    The start solutions are then the obvious hyperplane intersections -- here x in {1, -1} times
    y in {1, 2}.

    A non-integer coefficient must be given *exactly* -- a Python ``float`` like ``1.5`` is refused
    (its ~16 digits would cap the arbitrary-precision tree).  Give it instead as an exact decimal or
    rational string (``'1.5'``, ``'3/2'``), a :class:`fractions.Fraction` (``Fraction(3, 2)``), or a
    :mod:`bertini.multiprec` value (``bertini.multiprec.Complex('1.5')`` /
    ``bertini.multiprec.Float('0.1')``, which carry full working precision).  See :func:`coefficient`.
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


def randomize(system, matrix=None):
    """Randomize an overdetermined ``System`` down to a square one, returning a NEW system.

    An overdetermined system (N functions, n variables, N > n) cannot be fed to a total-degree
    start system, which needs a square system.  Randomization replaces the N functions with n
    generic combinations whose isolated solutions still contain the original's; you then solve
    the square result and discard the extraneous solutions by re-evaluating the original system.

    The original ``system`` is **not** mutated.

    Parameters
    ----------
    system : the overdetermined System to randomize.
    matrix : optional exact coefficient matrix R (array/list of lists), one row per randomized
        function and one column per natural function of ``system``.  When ``None`` (the default)
        bertini builds a generic R for you -- for a single affine variable group it sorts the
        functions by descending degree and uses ``R = [I | C]`` (random ``C``), giving the optimal
        total-degree path count (the product of the n largest degrees); for several variable groups
        it uses a dense random R.  When supplied, R is used verbatim and the functions are kept in
        their current order.  Coefficients must be exact (see :func:`coefficient`); Python floats
        are refused.

    Returns
    -------
    A new square ``System`` carrying a single randomization block.

    Examples
    --------
    Three quadrics in two variables (overdetermined), squared up to two functions::

        >>> import bertini
        >>> from bertini import linalg
        >>> x, y = bertini.Variable('x'), bertini.Variable('y')
        >>> S = bertini.System()
        >>> S.add_variable_group(bertini.VariableGroup([x, y]))
        >>> S.add_function(x*x + y*y - 1)
        >>> S.add_function(x*y)
        >>> S.add_function(x*x + y*y - x - y)
        >>> R = linalg.randomize(S)
        >>> R.num_functions(), S.num_functions()
        (2, 3)
    """
    if matrix is None:
        return system.randomize()

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
    return system.randomize(M)


def _coerce_mpfr_matrix(coefficients, context):
    """Coerce a (possibly nested list / numpy) array of exact values to an mpfr_complex matrix.

    Shared by the slice helpers; mirrors the per-entry coercion of :func:`add_linear_forms`.
    Returns ``(numpy_object_matrix, num_rows, num_cols)``.  Refuses Python floats.
    """
    rows = [list(r) for r in coefficients]
    if not rows:
        raise ValueError(f"{context} must have at least one row")
    ncol = len(rows[0])
    M = np.empty((len(rows), ncol), dtype=_mp.Complex)
    for i, row in enumerate(rows):
        if len(row) != ncol:
            raise ValueError(f"{context} is ragged (rows of differing length)")
        for j, entry in enumerate(row):
            M[i, j] = _exact_to_mpfr(entry)
    return M, len(rows), ncol


def slice_from_coefficients(coefficients, variables, homogeneous=False):
    """Build a :class:`Slice` from an exact augmented coefficient matrix.

    The linear part of a witness set.  ``coefficients`` is an ``(m x n+1)`` array/list of EXACT
    values (see :func:`coefficient`; Python floats are refused) -- one row per linear form, the
    trailing column being each form's constant term (give ``0`` there for a homogeneous slice).
    ``variables`` is the length-``n`` vector of variables the slice is a function of (a numpy object
    array of :class:`~bertini.function_tree.symbol.Variable`, e.g. from :func:`variable_vector`, or a
    :class:`~bertini.container.VariableGroup`).

    Returns a ``bertini.nag_algorithm.Slice``.  Its rows are also ready-made factors for a
    products-of-linears block -- see :func:`add_slices_as_products`.

    Examples
    --------
    The line ``2x + y - 1 = 0`` as a one-form slice on (x, y)::

        >>> import bertini                                            # doctest: +SKIP
        >>> from bertini import linalg                                # doctest: +SKIP
        >>> x, y = bertini.Variable('x'), bertini.Variable('y')      # doctest: +SKIP
        >>> s = linalg.slice_from_coefficients([[2, 1, -1]], [x, y]) # doctest: +SKIP
        >>> s.dimension(), s.num_variables()                         # doctest: +SKIP
        (1, 2)
    """
    from bertini.nag_algorithm import Slice  # local import: linalg loads before this is needed

    vg = variables if isinstance(variables, _VariableGroup) else _VariableGroup(list(variables))
    M, _, ncol = _coerce_mpfr_matrix(coefficients, "slice coefficients")
    if ncol != len(vg) + 1:
        raise ValueError(
            f"slice coefficient matrix has {ncol} columns but needs num_variables+1 = {len(vg) + 1} "
            "(one column per variable plus a trailing constant-term column)"
        )
    return Slice.from_coefficients(vg, M, homogeneous)


def add_slices_as_products(system, slices):
    """Add one products-of-linears function per slice to ``system`` (the regen bridge).

    Each :class:`Slice` is a stack of linear forms; this turns each slice into a single function that
    is the *product* of its forms -- a product of hyperplanes -- using the same first-class
    ``ProductsOfLinearsBlock`` as :func:`add_products_of_linears`.  Slice rows and product-of-linears
    factor rows share the augmented ``num_vars+1`` layout, so this is a direct hand-off.

    The slices' coefficient columns must follow ``system``'s variable ordering (build each slice over
    the system's variables).  Returns ``system`` for chaining.
    """
    # atleast_2d: a one-form slice's coefficients come back from eigenpy as a 1-D array; keep it a
    # (1 x num_vars+1) factor matrix so a single-hyperplane slice is still one degree-1 product.
    factors = [np.atleast_2d(np.asarray(s.coefficients())) for s in slices]
    return add_products_of_linears(system, factors)
