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
      * a :class:`~bertini.VariableGroup` -> added as an affine variable group;
      * a :class:`~bertini.Slice` -> its linear forms are added as a linear-forms block;
      * a numpy array / list / tuple of the above -> each element is added (so
        ``sys.add(A @ x - lam*x)``, ``sys.add([f, g])``, and ``sys.add(grp, f, g)`` work).

    Projective groups still use :meth:`add_hom_variable_group`; structured blocks use
    :meth:`~bertini.System.add_linear`.  Returns ``self`` for chaining.
    """
    from bertini._pybertini.nag_algorithms import Slice as _Slice
    for obj in objects:
        if isinstance(obj, _AbstractNode):
            self.add_function(obj)
        elif isinstance(obj, _VariableGroup):
            self.add_variable_group(obj)
        elif isinstance(obj, _Slice):
            obj.add_to(self)                     # its linear forms, as a linear-forms block (#372)
        elif isinstance(obj, (_np.ndarray, list, tuple)):
            for elt in (obj.ravel() if isinstance(obj, _np.ndarray) else obj):
                self.add(elt)
        else:
            raise TypeError(
                f"System.add does not know how to add a {type(obj).__name__}; pass a "
                "function-tree expression, a VariableGroup, a Slice, or an array/list of those"
            )
    return self


System.add = _system_add


_native_system_eval = System.eval


def _eval_argument(value):
    """A point given as a list/tuple (or a numpy array of Python numbers) as the vector the
    native overloads accept: double when every entry is an ordinary number, multiprecision when
    any entry already is one.  Other arguments (numpy vectors, scalars, times) pass through."""
    if not isinstance(value, (list, tuple)):
        return value
    from bertini._pybertini.multiprec import complex_mp as _complex_mp, real_mp as _real_mp
    if any(isinstance(e, (_complex_mp, _real_mp)) for e in value):
        return _np.array([e if isinstance(e, _complex_mp) else _complex_mp(e) for e in value])
    return _np.asarray(value, dtype=complex)


def _system_eval(self, *args):
    """Evaluate the system: ``eval()`` at the values already set, ``eval(point)``, or
    ``eval(point, time)`` for a system with a path variable.

    ``point`` may be a numpy array (double or multiprecision) or a plain list/tuple of numbers;
    a list is converted here, so ``sys.eval([0.3, 0.7])`` works instead of failing inside the
    native overload resolution with an argument error (issue #367).  A list holding any
    multiprecision entry evaluates in multiple precision; a list of ordinary numbers in double.
    """
    return _native_system_eval(self, *(_eval_argument(a) for a in args))


_system_eval.__doc__ = (_native_system_eval.__doc__ or "").rstrip() + "\n\n" + _system_eval.__doc__
System.eval = _system_eval


def _system_jacobian(self, usercoordinates=True):
    """The symbolic Jacobian of the system, as a 2-D numpy object array of expression nodes.

    ``J[i, j]`` is the partial derivative of function ``i`` with respect to variable ``j`` -- an
    expression tree, not a number (contrast :meth:`eval_jacobian`, which is numeric).  Ready to
    ``numpy.vstack`` onto a coefficient row and ``@`` a vector of variables.

    Parameters
    ----------
    usercoordinates : bool, default True
        When True, differentiate the functions *as authored* (the natural, pre-homogenization
        functions) with respect to the user-declared affine/projective variable groups: the
        solver-added homogenizing variables never appear and patches are omitted.  When False,
        differentiate the functions *as currently stored* (possibly homogenized) with respect to
        the full internal variable ordering (homogenizing variables included), with the patch's
        rows appended when the system is patched.

    Notes
    -----
    For a system that has *already* been homogenized, the user-coordinate Jacobian relies on the
    natural functions snapshotted at homogenization time.  Build the system affinely and call
    ``jacobian`` before homogenizing/solving for the cleanest result.
    """
    rows = self.symbolic_jacobian(usercoordinates)
    return _np.array(rows, dtype=object)


System.jacobian = _system_jacobian

# Override C++ submodule reference with the Python wrapper (which has AbstractStartSystem removed).
# Can't use 'from . import start_system': the star import already set that name to the C++ submodule.
import importlib as _importlib
start_system = _importlib.import_module('bertini.system.start_system')
del _importlib

__all__ = dir(_pybsys)
__all__.extend(['start_system'])
