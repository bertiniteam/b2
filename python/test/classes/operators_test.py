# This file is part of Bertini 2.
#
# python/test/classes/operators_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/classes/operators_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this file.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team

"""bertini.operators: one star-import, functions that work on symbols AND numbers AND
numpy containers alike, dispatched per argument."""

import builtins

import numpy as np
import pytest

import bertini as pb
import bertini.multiprec as mp
import bertini.operators as ops
from bertini.multiprec import complex_mp, real_mp
from bertini._pybertini.function_tree import AbstractNode


class TestPolymorphicDispatch:
    """sin & friends: symbolic on expressions, numeric on everything else."""

    def test_symbolic_on_variables(self):
        x = pb.Variable('x')
        f = ops.sin(x) + ops.Pi * x - ops.E
        assert isinstance(ops.sin(x), AbstractNode)
        assert isinstance(f, AbstractNode)

    def test_numeric_on_mp_scalars(self):
        v = real_mp('0.5')
        assert ops.sin(v) == mp.sin(v)
        assert ops.exp(v) == mp.exp(v)
        z = complex_mp('0.5', '0.25')
        assert ops.sqrt(z) == mp.sqrt(z)

    def test_numeric_on_mp_arrays(self):
        v = np.array([real_mp('0.25'), real_mp('0.5')])
        out = ops.cos(v)
        assert out.dtype == np.dtype(real_mp)
        assert out[0] == mp.cos(v[0])

    def test_numeric_on_lists_and_python_numbers(self):
        assert ops.sin(0.0) == 0.0
        out = ops.tan([real_mp('0.25'), real_mp('0.5')])
        assert out[1] == mp.tan(real_mp('0.5'))

    def test_asin_maps_to_arcsin(self):
        v = real_mp('0.5')
        assert ops.asin(v) == mp.asin(v)
        assert ops.acos(v) == mp.acos(v)
        assert ops.atan(v) == mp.atan(v)

    def test_hyperbolics_numeric(self):
        v = real_mp('0.5')
        assert ops.sinh(v) == mp.sinh(v)
        assert ops.atanh(v) == mp.atanh(v)

    def test_hyperbolics_reject_symbols(self):
        x = pb.Variable('x')
        with pytest.raises(TypeError, match="symbolic"):
            ops.sinh(x)


class TestComponentsAndFriends:
    """abs/arg/real/imag/conj/round/sum/norm/is_real, all in the same namespace."""

    def test_components_on_arrays(self):
        w = np.array([complex_mp(1, 2), complex_mp(3, 4)])
        assert [str(t) for t in ops.real(w)] == ['1', '3']
        assert [str(t) for t in ops.imag(w)] == ['2', '4']
        assert str(ops.abs(np.array([complex_mp(3, 4)]))[0]) == '5'
        assert ops.arg(w)[0] == mp.arg(w[0])
        assert complex(ops.conj(w)[1]) == complex(3, -4)

    def test_arg_on_scalars_and_reals(self):
        assert ops.arg(complex_mp(0, 1)) == mp.arg(complex_mp(0, 1))
        # arg of a negative real is pi
        assert mp.abs(ops.arg(np.array([real_mp(-2)]))[0] - mp.arg(complex_mp(-2))) == 0
        # plain python numbers give floats
        assert ops.arg(1j) == pytest.approx(np.pi / 2)

    def test_sum_norm_is_real(self):
        v = np.array([real_mp(3), real_mp(4)])
        assert ops.sum(v) == real_mp(7)
        assert ops.norm(v) == real_mp(5)
        assert ops.is_real(np.array([complex_mp(1)])) is True

    def test_round_stays_decimal(self):
        assert str(ops.round(real_mp('2.34567'), 2)) == '2.35'

    def test_numeric_only_reject_symbols(self):
        x = pb.Variable('x')
        for fn in (ops.abs, ops.arg, ops.real, ops.imag, ops.conj,
                   ops.round, ops.sum, ops.norm, ops.is_real):
            with pytest.raises(TypeError, match="symbolic"):
                fn(x)

    def test_builtin_fallback_on_plain_python(self):
        # the shadowing names still behave sanely on plain python input
        assert ops.abs(-3) == 3
        assert ops.sum([1, 2, 3]) == 6
        assert ops.round(2.345, 1) == builtins.round(2.345, 1)


class TestStarImportSurface:
    def test_star_import_gives_the_whole_vocabulary(self):
        ns = {}
        exec("from bertini.operators import *", ns)
        for name in ('sin', 'cos', 'tan', 'asin', 'acos', 'atan',
                     'sinh', 'cosh', 'tanh', 'asinh', 'acosh', 'atanh',
                     'exp', 'log', 'sqrt',
                     'abs', 'arg', 'real', 'imag', 'conj',
                     'round', 'sum', 'norm', 'is_real',
                     'E', 'Pi', 'I'):
            assert name in ns, name

    def test_one_import_covers_symbols_and_numbers(self):
        # the point of the module, as a single flow
        ns = {}
        exec("from bertini.operators import *", ns)
        x = pb.Variable('x')
        f = ns['sin'](x)                        # symbolic
        assert isinstance(f, AbstractNode)
        val = ns['sin'](real_mp('0.5'))         # mp scalar
        assert val == mp.sin(real_mp('0.5'))
        w = np.array([complex_mp(1, 2)])        # numpy container
        assert ns['imag'](w)[0] == real_mp(2)

    def test_top_level_functions_are_polymorphic_too(self):
        x = pb.Variable('x')
        assert isinstance(pb.sin(x), AbstractNode)
        assert pb.sin(real_mp('0.5')) == mp.sin(real_mp('0.5'))
        assert pb.arg(complex_mp(1, 1)) == mp.arg(complex_mp(1, 1))

    def test_from_bertini_star_still_never_shadows_builtins(self):
        ns = {}
        exec("from bertini import *", ns)
        assert 'abs' not in ns
        assert 'round' not in ns
        assert 'sum' not in ns
