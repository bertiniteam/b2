# This file is part of Bertini 2.
#
# python/bertini/_coefficients.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/_coefficients.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/bertini/_coefficients.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""Exact-coefficient coercion for building systems.

The single home for turning exact values into function-tree coefficient *nodes*
(:func:`coefficient` / :func:`coefficients`) and into multiprecision coefficient *matrices* (the
private ``_exact_to_mpfr`` / ``_coerce_mpfr_matrix`` used by the ``System`` block methods).

**Coefficients must be exact.**  Python ints (including numpy integers) and bertini nodes work
directly inside numpy arithmetic, so an *integer* matrix needs nothing special.  For any other
coefficient use :func:`coefficient` / :func:`coefficients`, which accept ``fractions.Fraction``,
exact strings (``'2.5'`` decimal or ``'3/4'`` rational), and bertini multiprecision values.  Python
``float`` and ``complex`` are **refused**: a 64-bit float literal would inject only ~16 correct
digits into the arbitrary-precision function tree, silently capping downstream precision.
"""

import numpy as np
from fractions import Fraction

from bertini.symbolics import Integer, Rational, Complex
from bertini import multiprec as _mp
from bertini._pybertini.function_tree import AbstractNode as _AbstractNode

# bertini multiprecision value types (not function-tree nodes) that we accept as exact
# coefficients by wrapping them in a Complex node.
_MP_VALUE_TYPES = tuple(
    t for t in (getattr(_mp, n, None) for n in ('real_mp', 'complex_mp', 'int_mp', 'rational_mp'))
    if isinstance(t, type)
)


def coefficient(value):
    """Coerce a single exact value into a function-tree coefficient node.

    Accepts bertini nodes (returned unchanged), ``int`` / numpy integers,
    ``fractions.Fraction``, exact strings (``'2.5'`` or ``'3/4'``), and bertini
    multiprecision values.  Raises ``TypeError`` for Python ``float`` / ``complex`` (and
    numpy floating types) -- a floating-point literal would cap precision.

    There are four ways to spell an *exact* non-integer coefficient::

        >>> from fractions import Fraction
        >>> import bertini
        >>> from bertini import multiprec
        >>> _ = bertini.coefficient('2.5')                      # exact decimal string
        >>> _ = bertini.coefficient('3/4')                      # exact rational string
        >>> _ = bertini.coefficient(Fraction(3, 4))             # a fractions.Fraction
        >>> _ = bertini.coefficient(multiprec.complex_mp('0.1'))   # a full-precision multiprec value

    A Python ``float`` is refused: ``bertini.coefficient(0.1)`` raises ``TypeError``.  An exact
    *complex* coefficient is a ``multiprec.complex_mp`` with real and imaginary parts:
    ``bertini.coefficient(multiprec.complex_mp('0.6', '0.8'))``.
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
        return Rational(value) if '/' in value else Complex(value)
    if _MP_VALUE_TYPES and isinstance(value, _MP_VALUE_TYPES):
        return Complex(value)
    if isinstance(value, (float, complex, np.floating, np.complexfloating)):
        raise TypeError(
            f"refusing to use the Python {type(value).__name__} {value!r} as a coefficient: "
            "a floating-point literal carries only ~16 digits and would cap the precision of "
            "the arbitrary-precision function tree.  Pass an exact value instead -- an int, a "
            "fractions.Fraction, an exact string such as '2.5' or '3/4', or a "
            "bertini.multiprec value."
        )
    raise TypeError(f"cannot use {type(value).__name__} as a coefficient")


def coefficients(array_like):
    """Coerce an array (or nested list) of exact values to an object array of coefficient nodes.

    Use this to bring a non-integer matrix or vector into the function tree::

        A = bertini.coefficients([['5/2', '1'], ['0', '3']])     # exact rationals
        A = bertini.coefficients(fraction_matrix)                # fractions.Fraction entries

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


def _exact_to_mpfr(value):
    """Convert a single exact value to a multiprec.complex_mp (refusing Python floats).

    Mirrors :func:`coefficient`'s exact-only rule, but produces a multiprecision *value*
    (for a coefficient matrix) rather than a function-tree node.  Accepts the same exact
    spellings as :func:`coefficient`, including a bertini constant *node* -- e.g. a
    ``symbolics.Complex`` -- so a value good enough for :func:`coefficient` is good enough
    here too (issue #326).
    """
    # A function-tree constant node (as coefficient() / Complex() / Integer() / Rational()
    # produce): accept it, the symmetric partner to the complex_mp case below.  A Complex node
    # carries a full-precision complex_mp directly; any other constant node is evaluated (it
    # holds no variables).  This closes the gap where a complex_mp was accepted but the Complex
    # *node* wrapping the same value was not (issue #326).
    if isinstance(value, _AbstractNode):
        if isinstance(value, Complex):
            return _mp.complex_mp(value.value())
        try:
            return _mp.complex_mp(value.eval())
        except Exception as e:
            raise TypeError(
                "cannot use a non-constant node as a coefficient (it still contains "
                f"variables): {value!r}"
            ) from e
    if _MP_VALUE_TYPES and isinstance(value, _MP_VALUE_TYPES):
        return _mp.complex_mp(value)
    if isinstance(value, bool):
        raise TypeError("a bool is not a valid coefficient")
    if isinstance(value, (int, np.integer)):
        return _mp.complex_mp(str(int(value)))
    if isinstance(value, Fraction):
        return _mp.complex_mp(str(value.numerator)) / _mp.complex_mp(str(value.denominator))
    if isinstance(value, str):
        if '/' in value:
            num, den = value.split('/')
            return _mp.complex_mp(num) / _mp.complex_mp(den)
        return _mp.complex_mp(value)
    if isinstance(value, (float, complex, np.floating, np.complexfloating)):
        raise TypeError(
            f"refusing to use the Python {type(value).__name__} {value!r} as a coefficient: "
            "a floating-point literal would cap the precision of the block.  Pass an exact "
            "value -- an int, a fractions.Fraction, an exact string, or a bertini.multiprec "
            "value."
        )
    raise TypeError(f"cannot use {type(value).__name__} as a coefficient")


def _coerce_mpfr_matrix(coefficients, context):
    """Coerce a (possibly nested list / numpy) array of exact values to an mpfr_complex matrix.

    Returns ``(numpy_object_matrix, num_rows, num_cols)``.  Refuses Python floats.
    """
    rows = [list(r) for r in coefficients]
    if not rows:
        raise ValueError(f"{context} must have at least one row")
    ncol = len(rows[0])
    M = np.empty((len(rows), ncol), dtype=_mp.complex_mp)
    for i, row in enumerate(rows):
        if len(row) != ncol:
            raise ValueError(f"{context} is ragged (rows of differing length)")
        for j, entry in enumerate(row):
            M[i, j] = _exact_to_mpfr(entry)
    return M, len(rows), ncol
