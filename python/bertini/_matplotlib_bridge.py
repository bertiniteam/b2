# This file is part of Bertini 2.
#
# python/bertini/_matplotlib_bridge.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/_matplotlib_bridge.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this file.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team

"""Teach matplotlib to draw multiprecision numbers, and to refuse to draw them wrongly.

Most of pyplot already works on ``real_mp``, because the type converts to a float and
matplotlib asks each element for one.  Everything that computes with the data before it
draws -- ``hist``, ``bar`` with a plain float width -- instead fails with a dtype
promotion error, because there is no common dtype between ``real_mp`` and ``float64``
and (deliberately) never will be: an automatic promotion would silently drop every digit
past the 16th in ordinary arithmetic.

matplotlib's own extension point for a number type it does not know is the units
registry, so that is what this module fills in.  A converter is consulted only for data
on its way to an axis, which is the one place where dropping digits is right: a screen
resolves about three.  The values themselves are never touched.

The complex converter *refuses* instead of converting, and that is the point of it.
numpy's ``.real`` / ``.imag`` on an array of a user complex dtype return the array
itself and an array of zeros -- see :mod:`bertini._numpy_guard` -- so
``scatter(pts.real, pts.imag)`` draws every point on the x axis, with no error and no
warning.  Because ``.real`` of a ``complex_mp`` array is still ``complex_mp``, refusing
complex data at the axis catches exactly that mistake and says what to write instead.

Registration happens when :func:`install` runs, whether matplotlib is imported before
bertini or long after it: if matplotlib is already loaded the converters go in
immediately, and if it is not, an import hook waits for it.
"""

import sys as _sys

import numpy as _np

from bertini.multiprec import real_mp as _real_mp, complex_mp as _complex_mp

_REAL_MP = _np.dtype(_real_mp)
_CPLX_MP = _np.dtype(_complex_mp)

_COMPLEX_MSG = (
    "a multiprecision complex value cannot be an axis coordinate.  Drawing it would show "
    "its real parts and call them the data -- and numpy's .imag on a multiprecision "
    "complex array is zeros rather than the imaginary parts, so the usual "
    "scatter(z.real, z.imag) would put every point on the x axis.  Plot the parts "
    "explicitly instead: scatter(bertini.real(z), bertini.imag(z))."
)


def _as_float(value):
    """One value, as a float, whatever multiprecision type it arrived as."""
    if isinstance(value, _complex_mp):
        raise TypeError(_COMPLEX_MSG)
    return float(value)


def _floats(value):
    """A float, or an array of floats, matching the shape of what came in.

    Complex data is refused here as well as in the complex converter, because an axis
    keeps the first converter it is given: once real data has set one, later complex
    data on that same axis arrives here rather than at the refusal.
    """
    if isinstance(value, _np.ndarray):
        if value.dtype == _CPLX_MP:
            raise TypeError(_COMPLEX_MSG)
        if value.dtype == _REAL_MP:
            return value.astype(float)
        if value.dtype == object:
            flat = [_as_float(v) for v in value.ravel()]
            return _np.array(flat, dtype=float).reshape(value.shape)
        return value.astype(float)

    if isinstance(value, (list, tuple)):
        return _np.array([_as_float(v) for v in value], dtype=float)

    return _as_float(value)


def _converters():
    """The two converter instances, built against the matplotlib that is now loaded."""
    import matplotlib.units as munits

    class RealConverter(munits.ConversionInterface):
        """Hand matplotlib float64 for multiprecision real data bound for an axis."""

        @staticmethod
        def convert(value, unit, axis):
            """Convert one value, or a whole array of them, to float64."""
            return _floats(value)

        @staticmethod
        def axisinfo(unit, axis):
            """No special ticks, labels or limits: these are ordinary numbers."""
            return munits.AxisInfo()

        @staticmethod
        def default_units(x, axis):
            """These values carry no unit."""
            return None

    class ComplexRefusal(munits.ConversionInterface):
        """Refuse multiprecision complex data at an axis, saying what to write instead."""

        @staticmethod
        def convert(value, unit, axis):
            """Always raises ``TypeError``."""
            raise TypeError(_COMPLEX_MSG)

        @staticmethod
        def axisinfo(unit, axis):
            """No axis information; the data never gets this far."""
            return munits.AxisInfo()

        @staticmethod
        def default_units(x, axis):
            """Always raises ``TypeError``; this runs before ``convert``."""
            raise TypeError(_COMPLEX_MSG)

    return RealConverter(), ComplexRefusal()


def _register(units_module):
    """Put the converters in a loaded ``matplotlib.units`` registry.  Idempotent."""
    registry = units_module.registry
    if getattr(registry.get(_real_mp), "_bertini_converter", False):
        return

    real_converter, complex_refusal = _converters()
    real_converter._bertini_converter = True
    complex_refusal._bertini_converter = True

    registry[_real_mp] = real_converter
    registry[_complex_mp] = complex_refusal


class _RegisterWhenMatplotlibArrives:
    """A meta path finder that registers the converters as ``matplotlib.units`` loads.

    It finds nothing itself.  It lets the ordinary finders produce the spec, wraps that
    spec's loader so registration runs the moment the module body finishes, and then
    takes itself back off ``sys.meta_path``.
    """

    WATCHED = "matplotlib.units"

    def find_spec(self, fullname, path=None, target=None):
        """Wrap the loader for ``matplotlib.units``; return ``None`` for everything else."""
        if fullname != self.WATCHED:
            return None

        spec = None
        for finder in _sys.meta_path:
            if finder is self:
                continue
            find = getattr(finder, "find_spec", None)
            if find is None:
                continue
            spec = find(fullname, path, target)
            if spec is not None:
                break

        if spec is None or spec.loader is None:
            return None

        loader = spec.loader
        original = loader.exec_module

        def exec_module(module, _original=original):
            _original(module)
            try:
                _register(module)
            except Exception:
                pass   # a plotting convenience must never break someone's import
            self.uninstall()

        loader.exec_module = exec_module
        return spec

    def uninstall(self):
        """Take this finder off ``sys.meta_path``, if it is still on it."""
        if self in _sys.meta_path:
            _sys.meta_path.remove(self)


def install():
    """Register the converters now, or arrange for it when matplotlib is imported.

    Idempotent, and safe whether or not matplotlib is installed at all.
    """
    units_module = _sys.modules.get("matplotlib.units")
    if units_module is not None:
        _register(units_module)
        return

    for finder in _sys.meta_path:
        if isinstance(finder, _RegisterWhenMatplotlibArrives):
            return

    _sys.meta_path.insert(0, _RegisterWhenMatplotlibArrives())
