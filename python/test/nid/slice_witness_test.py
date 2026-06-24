# This file is part of Bertini 2.
#
# python/test/nid/slice_witness_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/nid/slice_witness_test.py is distributed in the hope that it will be useful,
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

"""The block-backed Slice and the WitnessSet exposure: building, evaluating, row-subsetting,
adding a slice to a System, the regen products-of-linears bridge, and constructing a witness
set both incrementally and all-at-once."""

import copy
import pickle

import numpy as np

import bertini as pb
from bertini import linalg
from bertini import multiprec as mp
from bertini.nag_algorithm import Slice, WitnessSetMultiplePrecision


def _mpvec(*entries):
    return np.array([mp.Complex(str(e)) for e in entries], dtype=mp.Complex)


def _known_slice():
    # f0 = 2x + 3y + 1,  f1 = x - y + 4   on (x, y)
    x, y = pb.Variable('x'), pb.Variable('y')
    return linalg.slice_from_coefficients([[2, 3, 1], [1, -1, 4]], [x, y])


# ---- Slice ---------------------------------------------------------------------------------

def test_slice_from_coefficients_shape():
    s = _known_slice()
    assert s.dimension() == 2
    assert s.num_variables() == 2
    assert not s.is_homogeneous()
    assert np.asarray(s.coefficients()).shape == (2, 3)


def test_slice_eval_mp():
    s = _known_slice()
    v = s.eval(_mpvec(1, 1))            # f0 = 6, f1 = 4
    assert abs(complex(v[0]) - 6) < 1e-12
    assert abs(complex(v[1]) - 4) < 1e-12


def test_slice_eval_double():
    s = _known_slice()
    v = s.eval(np.array([1 + 0j, 1 + 0j]))
    assert abs(v[0] - 6) < 1e-12
    assert abs(v[1] - 4) < 1e-12


def test_slice_coefficients_match_first_row():
    s = _known_slice()
    C = np.asarray(s.coefficients())
    assert [complex(C[0, j]).real for j in range(3)] == [2.0, 3.0, 1.0]


def test_random_complex_slice():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    vg = pb.VariableGroup([x, y, z])
    s = Slice.random_complex(vg, 2)
    assert s.dimension() == 2
    assert s.num_variables() == 3
    assert not s.is_homogeneous()


def test_homogeneous_slice_has_zero_constant_column():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    vg = pb.VariableGroup([x, y, z])
    s = Slice.random_complex(vg, 2, True)
    assert s.is_homogeneous()
    C = np.asarray(s.coefficients())
    assert all(abs(complex(C[i, 3])) < 1e-25 for i in range(2))


def test_random_real_slice():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    s = Slice.random_real(pb.VariableGroup([x, y, z]), 2)
    assert s.dimension() == 2


# ---- row subsetting / composition ----------------------------------------------------------

def test_slice_len_and_getitem_slice():
    s = _known_slice()
    assert len(s) == 2
    head = s[:1]
    assert head.dimension() == 1
    # a one-form slice's coefficients come back from eigenpy as a 1-D array (numpy convention).
    assert np.atleast_2d(np.asarray(head.coefficients())).shape == (1, 3)
    tail = s[-1:]
    assert tail.dimension() == 1


def test_slice_getitem_list_and_int():
    s = _known_slice()
    reordered = s[[1, 0]]
    assert reordered.dimension() == 2
    single = s[0]
    assert single.dimension() == 1


def test_head_tail_rows_methods():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    s = Slice.random_complex(pb.VariableGroup([x, y, z]), 4)
    assert s.head(2).dimension() == 2
    assert s.tail(1).dimension() == 1
    assert s.rows([0, 3]).dimension() == 2


def test_subslice_rows_match_parent():
    s = _known_slice()
    C = np.asarray(s.coefficients())
    sub = np.atleast_2d(np.asarray(s[[1]].coefficients()))   # one-form slice -> 1-D from eigenpy
    assert all(abs(complex(C[1, j]) - complex(sub[0, j])) < 1e-25 for j in range(3))


# ---- a slice as a System block -------------------------------------------------------------

def test_slice_add_to_system_agrees_with_slice_eval():
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    s = linalg.slice_from_coefficients([[2, 3, 1], [1, -1, 4]], [x, y])
    s.add_to(sys)
    assert sys.num_functions() == 2

    pt = _mpvec(1, 1)
    from_sys = sys.eval(pt)
    from_slice = s.eval(pt)
    assert abs(complex(from_sys[0]) - complex(from_slice[0])) < 1e-12
    assert abs(complex(from_sys[1]) - complex(from_slice[1])) < 1e-12


def test_add_slices_as_products_is_the_regen_bridge():
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    s1 = linalg.slice_from_coefficients([[1, 0, -1], [1, 0, 1]], [x, y])    # (x - 1)(x + 1)
    s2 = linalg.slice_from_coefficients([[0, 1, -1], [0, 1, -2]], [x, y])   # (y - 1)(y - 2)
    linalg.add_slices_as_products(sys, [s1, s2])
    assert list(sys.degrees()) == [2, 2]
    # at (x, y) = (1, 1):  f0 = (0)(2) = 0,  f1 = (0)(-1) = 0
    v = sys.eval(_mpvec(1, 1))
    assert abs(complex(v[0])) < 1e-12
    assert abs(complex(v[1])) < 1e-12


# ---- WitnessSet ----------------------------------------------------------------------------

def test_witness_set_built_incrementally():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    vg = pb.VariableGroup([x, y, z])

    w = WitnessSetMultiplePrecision()
    assert w.degree() == 0

    w.add_point(_mpvec(1, 0, 0))
    w.add_point(_mpvec(0, 1, 0))
    assert w.degree() == 2

    w.set_slice(Slice.random_complex(vg, 1))
    assert w.dimension() == 1

    assert len(w.get_points()) == 2
    assert abs(complex(w.get_point(0)[0]) - 1) < 1e-12
    assert w.get_slice().dimension() == 1


def test_witness_set_from_parts_and_consistency():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    vg = pb.VariableGroup([x, y, z])

    sys = pb.System()
    sys.add_variable_group(vg)
    sys.add_function(x * x + y * y + z * z - 1)   # the sphere: a dim-2 component in 3-space

    s = Slice.random_complex(vg, 2)
    pts = [_mpvec(1, 0, 0), _mpvec(0, 1, 0)]
    w = WitnessSetMultiplePrecision(pts, s, sys)

    assert w.degree() == 2
    assert w.dimension() == 2
    # 3 variables - 1 natural function == slice dimension 2
    assert w.is_consistent()
    assert w.get_system().num_functions() == 1


# ---- serialization (pickle / deepcopy) -----------------------------------------------------

def test_slice_pickle_roundtrip():
    s = _known_slice()
    s2 = pickle.loads(pickle.dumps(s))
    assert s2.dimension() == 2
    assert s2.num_variables() == 2
    C, C2 = np.asarray(s.coefficients()), np.asarray(s2.coefficients())
    assert all(abs(complex(C[i, j]) - complex(C2[i, j])) < 1e-25
               for i in range(2) for j in range(3))
    v = s2.eval(_mpvec(1, 1))
    assert abs(complex(v[0]) - 6) < 1e-12 and abs(complex(v[1]) - 4) < 1e-12


def test_slice_deepcopy():
    s = _known_slice()
    s2 = copy.deepcopy(s)
    assert s2.dimension() == 2
    assert np.asarray(s2.coefficients()).shape == (2, 3)


def test_witness_set_pickle_roundtrip():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    vg = pb.VariableGroup([x, y, z])

    sys = pb.System()
    sys.add_variable_group(vg)
    sys.add_function(x * x + y * y + z * z - 1)

    s = Slice.random_complex(vg, 2)
    w = WitnessSetMultiplePrecision([_mpvec(1, 0, 0), _mpvec(0, 1, 0)], s, sys)

    w2 = pickle.loads(pickle.dumps(w))
    assert w2.degree() == 2
    assert w2.dimension() == 2
    assert w2.is_consistent()
    assert w2.get_slice().dimension() == 2
    assert len(w2.get_points()) == 2

    # the deserialized system is usable: its load() re-differentiated it.
    pt = _mpvec(1, 1, 1)
    v = w2.get_system().eval(pt)
    assert abs(complex(v[0]) - 2) < 1e-12     # 1 + 1 + 1 - 1 = 2
    # a 1-function Jacobian comes back from eigenpy as a 1-D array (see other tests); the point is
    # that differentiating the deserialized system works at all -- load() restored it.
    assert np.atleast_2d(np.asarray(w2.get_system().eval_jacobian(pt))).shape == (1, 3)
