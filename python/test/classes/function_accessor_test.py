# This file is part of Bertini 2.
#
# python/test/classes/function_accessor_test.py is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/classes/function_accessor_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this file.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""``system.function(i)`` for every kind of evaluation block.

A System is block-composed: its rows come from a PolynomialBlock and/or structured blocks
(linear forms, products of linears, randomization, blend).  ``function(i)`` must return the
i-th natural function as a function-tree node for *every* block type -- structured blocks are
expanded to nodes on demand.  Regression for issue #263, where ``function(0)`` segfaulted on a
system built from a slice (a pure LinearFormsBlock, so no PolynomialBlock to dereference).

``system.slices()`` is the companion accessor: it backs out the linear-form slices embedded in
a system -- the inverse of ``Slice.as_system()``.
"""

import numpy as np
import pytest

import bertini as pb
from bertini.symbolics import AbstractNode
from bertini import nag_algorithm as na
from bertini import Slice


def _vg(*vs):
    return pb.VariableGroup(list(vs))


def _all_functions_are_nodes(system):
    """function(i) returns a function-tree node for every i in range, and overshoot raises."""
    n = system.num_functions()
    assert n > 0
    for i in range(n):
        f = system.function(i)
        assert isinstance(f, AbstractNode)
        assert str(f)                         # has a printable form
    with pytest.raises(IndexError):
        system.function(n)                    # past the end raises, no longer segfaults
    return [system.function(i) for i in range(n)]


# ---- function(i) across every block type ----------------------------------------------------

def test_function_polynomial_block():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System(); s.add_variable_group(_vg(x, y))
    s.add_function(x * x + y * y - 1); s.add_function(x - y)
    fns = _all_functions_are_nodes(s)
    assert 'x' in str(fns[0]) and 'y' in str(fns[0])


def test_function_linear_forms_block_from_slice():
    # issue #263: a system built from a slice has only a LinearFormsBlock (no PolynomialBlock).
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    s = Slice.random_complex(_vg(x, y, z), 1)
    sys = s.as_system()
    fns = _all_functions_are_nodes(sys)       # used to segfault here
    text = str(fns[0])
    assert 'x' in text and 'y' in text and 'z' in text   # the expanded linear form


def test_function_products_of_linears_block():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System(); s.add_variable_group(_vg(x, y))
    s.add_products_of_linears([[[1, 0, -1], [1, 0, 1]]])   # (x-1)(x+1)
    fns = _all_functions_are_nodes(s)
    assert 'x' in str(fns[0])


def test_function_randomization_block():
    x, y = pb.Variable('x'), pb.Variable('y')
    o = pb.System(); o.add_variable_group(_vg(x, y))
    o.add_function(x * x + y * y - 1); o.add_function(x * y); o.add_function(x * x - y)
    r = o.randomize()                   # only a RandomizationBlock, no PolynomialBlock
    _all_functions_are_nodes(r)


def test_function_blend_block():
    x, y = pb.Variable('x'), pb.Variable('y')
    fx = pb.System(); fx.add_variable_group(_vg(x, y)); fx.add_function(x * x + y * y - 1)
    sm = pb.System(); sm.add_variable_group(_vg(x, y)); sm.add_function(y)
    em = pb.System(); em.add_variable_group(_vg(x, y)); em.add_function(y - x)
    H = na.moving_homotopy(fx, sm, em,
                           gamma=pb.coefficient(pb.multiprec.complex_mp('0.6', '0.8')))
    _all_functions_are_nodes(H)               # a BlendBlock row expands recursively


def test_function_matches_eval_for_slice_system():
    # the node returned by function(i) must evaluate to the same value as the block.
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.Slice.from_coefficients([[2, 3, 1], [1, -1, 4]], [x, y]).as_system()
    one = pb.multiprec.complex_mp('1')
    pt = np.array([one, one], dtype=pb.multiprec.complex_mp)
    from_block = sys.eval(pt)                  # System eval: a vector in variable order
    for i in range(sys.num_functions()):
        node_val = sys.function(i).eval(x=one, y=one)   # node eval: variable values by keyword
        assert abs(complex(node_val) - complex(from_block[i])) < 1e-12


# ---- slices(): back out the slice structure -------------------------------------------------

def test_slices_empty_for_polynomial_system():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System(); s.add_variable_group(_vg(x, y)); s.add_function(x * x - y)
    assert s.slices() == []


def test_slices_roundtrips_coefficients():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.Slice.from_coefficients([[2, 3, 1], [1, -1, 4]], [x, y])
    recovered = s.as_system().slices()
    assert len(recovered) == 1
    orig = np.asarray(s.coefficients())
    back = np.asarray(recovered[0].coefficients())
    assert back.shape == orig.shape
    assert all(abs(complex(orig[i, j]) - complex(back[i, j])) < 1e-25
               for i in range(orig.shape[0]) for j in range(orig.shape[1]))


def test_slices_one_per_linear_forms_block():
    # a system that mixes a polynomial row and a linear-forms block exposes exactly one slice.
    x, y = pb.Variable('x'), pb.Variable('y')
    m = pb.System(); m.add_variable_group(_vg(x, y))
    m.add_function(x * x + y * y - 1)
    m.add_linear(np.array([[2, 1]]), np.array([x, y]), [-1])   # 2x + y - 1
    slices = m.slices()
    assert len(slices) == 1
    assert slices[0].dimension() == 1
