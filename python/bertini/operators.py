# This file is part of Bertini 2.
#
# python/bertini/operators.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/operators.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/bertini/operators.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""The whole math vocabulary in one namespace -- symbols and numbers alike.

One import gives you functions that work on *everything*: a symbolic
:class:`~bertini.Variable`/expression, a multiprecision number, a numpy array of them,
or a plain python number::

    from bertini.operators import *

    f = sin(x) + Pi*y - E          # symbolic       (x, y Variables -> an expression)
    v = sin(real_mp('0.5'))        # numeric        (full precision)
    m = abs(solutions[0])          # numpy arrays   (multiprecision dtypes included)
    t = arg(complex_mp(1, 1))      # components     (arg/real/imag/conj)

Dispatch is by argument: a function-tree node builds a symbolic node; everything else
takes the numeric path (multiprecision scalars and mp-dtype numpy arrays go through the
native precision-preserving loops; python numbers and float arrays are plain numpy).
You never have to remember whether a name lives in ``bertini.multiprec`` or at the top
level -- it is here.

The functions with no symbolic counterpart (``abs``, ``arg``, ``real``, ``imag``,
``conj``, ``round``, ``sum``, ``norm``, ``is_real``, and the hyperbolics) raise a clear
``TypeError`` when handed a symbolic expression.

``abs``, ``round``, and ``sum`` shadow the python builtins **within your namespace**
when you star-import this module -- that is the point (they fall back to builtin
behavior on plain python input), but it is opt-in: ``from bertini import *`` never
shadows builtins.
"""

import numpy as _np

from bertini._pybertini.function_tree import AbstractNode as _AbstractNode

from . import symbolics as _sym
from . import _numpy_helpers as _nh

# the symbolic constants, ready to drop into expressions
from . import E, Pi, I  # noqa: F401


def _polymorphic(name, sym_fn, num_fn, doc):
    def f(x):
        if isinstance(x, _AbstractNode):
            return sym_fn(x)
        return num_fn(x)
    f.__name__ = name
    f.__qualname__ = name
    f.__doc__ = doc + ("\n\nPolymorphic: builds a symbolic node for a function-tree "
                       "argument, computes numerically (precision-preserving, numpy "
                       "containers included) for everything else.")
    return f


def _numeric_only(name, num_fn, doc):
    def f(x, *args, **kwargs):
        if isinstance(x, _AbstractNode):
            raise TypeError(
                f"{name}() is not defined for symbolic expressions -- it is a numeric "
                "operation.  Evaluate the expression first, or use the symbolic "
                "functions (sin, cos, exp, ...) to build systems.")
        return num_fn(x, *args, **kwargs)
    f.__name__ = name
    f.__qualname__ = name
    f.__doc__ = doc + ("\n\nNumeric: multiprecision scalars, numpy arrays (mp dtypes "
                       "included), lists, and plain python numbers.")
    return f


# --- the elementary functions with symbolic twins: full polymorphic dispatch -------------------

sin  = _polymorphic('sin',  _sym.sin,  _np.sin,    "Sine.")
cos  = _polymorphic('cos',  _sym.cos,  _np.cos,    "Cosine.")
tan  = _polymorphic('tan',  _sym.tan,  _np.tan,    "Tangent.")
asin = _polymorphic('asin', _sym.asin, _np.arcsin, "Arcsine.")
acos = _polymorphic('acos', _sym.acos, _np.arccos, "Arccosine.")
atan = _polymorphic('atan', _sym.atan, _np.arctan, "Arctangent.")
exp  = _polymorphic('exp',  _sym.exp,  _np.exp,    "Exponential, base e.")
log  = _polymorphic('log',  _sym.log,  _np.log,    "Natural logarithm.")
sqrt = _polymorphic('sqrt', _sym.sqrt, _np.sqrt,   "Square root.")

# --- numeric-only elementary functions (no symbolic node exists) -------------------------------

sinh  = _numeric_only('sinh',  _np.sinh,    "Hyperbolic sine.")
cosh  = _numeric_only('cosh',  _np.cosh,    "Hyperbolic cosine.")
tanh  = _numeric_only('tanh',  _np.tanh,    "Hyperbolic tangent.")
asinh = _numeric_only('asinh', _np.arcsinh, "Hyperbolic arcsine.")
acosh = _numeric_only('acosh', _np.arccosh, "Hyperbolic arccosine.")
atanh = _numeric_only('atanh', _np.arctanh, "Hyperbolic arctangent.")

# --- components, magnitudes, and friends (numeric-only) ----------------------------------------

abs     = _numeric_only('abs',     _nh.abs,     "Magnitude(s), as real_mp for mp input.")        # noqa: A001
arg     = _numeric_only('arg',     _nh.arg,     "Argument(s) (angle from 0), as real_mp.  Beware the branch cut.")
real    = _numeric_only('real',    _nh.real,    "Real part(s), as real_mp for mp input.")
imag    = _numeric_only('imag',    _nh.imag,    "Imaginary part(s), as real_mp for mp input.")
conj    = _numeric_only('conj',    _nh.conj,    "Complex conjugate(s).")
round   = _numeric_only('round',   _nh.round,   "Round to N DECIMAL digits, staying mp-native.") # noqa: A001
sum     = _numeric_only('sum',     _nh.sum,     "Sum of a collection, staying mp-native.")       # noqa: A001
norm    = _numeric_only('norm',    _nh.norm,    "Euclidean (2-)norm, as real_mp.")
is_real = _numeric_only('is_real', _nh.is_real, "Is every coordinate real (abs(imag) < tol)?")


__all__ = [
    'sin', 'cos', 'tan', 'asin', 'acos', 'atan',
    'sinh', 'cosh', 'tanh', 'asinh', 'acosh', 'atanh',
    'exp', 'log', 'sqrt',
    'abs', 'arg', 'real', 'imag', 'conj', 'round', 'sum', 'norm', 'is_real',
    'E', 'Pi', 'I',
]
