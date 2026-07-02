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
import pytest

import bertini as pb
from bertini import multiprec as mp
from bertini import Slice
from bertini.nag_algorithm import WitnessSetMultiplePrecision


def _mpvec(*entries):
    return np.array([mp.Complex(str(e)) for e in entries], dtype=mp.Complex)


def _known_slice():
    # f0 = 2x + 3y + 1,  f1 = x - y + 4   on (x, y)
    x, y = pb.Variable('x'), pb.Variable('y')
    return pb.Slice.from_coefficients([[2, 3, 1], [1, -1, 4]], [x, y])


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

def test_slice_getitem_slice_is_a_subcollection():
    # s[i:j] selects a SUB-COLLECTION -> a (sub-)Slice (Python list semantics).
    s = _known_slice()
    assert len(s) == 2
    head = s[:1]
    assert head.dimension() == 1
    assert head.coefficients().shape == (1, 3)   # coefficients() is always 2-D
    tail = s[-1:]
    assert tail.dimension() == 1


def test_slice_getitem_list_is_a_subcollection():
    s = _known_slice()
    reordered = s[[1, 0]]
    assert reordered.dimension() == 2
    rows = np.asarray(reordered.coefficients())
    assert [complex(rows[0, j]).real for j in range(3)] == [1.0, -1.0, 4.0]   # form 1 came first


def test_slice_int_index_is_a_form_vector():
    # s[i] selects an ELEMENT -> the i-th form's coefficient VECTOR (1-D), "a line is a vector".
    s = _known_slice()
    v0 = np.asarray(s[0])
    assert v0.shape == (3,)
    assert [complex(c).real for c in v0] == [2.0, 3.0, 1.0]
    v_last = np.asarray(s[-1])
    assert [complex(c).real for c in v_last] == [1.0, -1.0, 4.0]
    with pytest.raises(IndexError):
        _ = s[5]


def test_iterating_slice_yields_form_vectors():
    s = _known_slice()
    forms = list(s)               # iteration yields elements = form vectors
    assert len(forms) == 2
    assert all(np.asarray(f).shape == (3,) for f in forms)


def test_coefficients_always_2d_even_for_one_form():
    x, y = pb.Variable('x'), pb.Variable('y')
    one = pb.Slice.from_coefficients([[2, 3, 1]], [x, y])
    assert one.coefficients().shape == (1, 3)     # NOT (3,): coefficients() never collapses


def test_slice_concatenate_and_add_operator():
    x, y = pb.Variable('x'), pb.Variable('y')
    sa = pb.Slice.from_coefficients([[2, 3, 1]], [x, y])    # 2x + 3y + 1
    sb = pb.Slice.from_coefficients([[1, -1, 4]], [x, y])   # x - y + 4

    combined = sa.concatenate(sb)
    assert combined.dimension() == 2
    v = combined.eval(_mpvec(1, 1))
    assert abs(complex(v[0]) - 6) < 1e-12 and abs(complex(v[1]) - 4) < 1e-12

    # the + operator is concatenation
    plus = sa + sb
    assert plus.dimension() == 2


def test_slice_as_system_matches_eval():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.Slice.from_coefficients([[2, 3, 1], [1, -1, 4]], [x, y])
    sys = s.as_system()
    assert sys.num_functions() == 2
    pt = _mpvec(1, 1)
    from_sys, from_slice = sys.eval(pt), s.eval(pt)
    assert all(abs(complex(from_sys[i]) - complex(from_slice[i])) < 1e-12 for i in range(2))


def test_slice_repr_is_readable():
    s = _known_slice()
    text = repr(s)
    assert 'slice' in text.lower()
    assert 'coefficient' in text.lower()


def test_head_tail_rows_methods():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    s = Slice.random_complex(pb.VariableGroup([x, y, z]), 4)
    assert s.head(2).dimension() == 2
    assert s.tail(1).dimension() == 1
    assert s.rows([0, 3]).dimension() == 2


def test_subslice_rows_match_parent():
    s = _known_slice()
    C = np.asarray(s.coefficients())
    sub = np.asarray(s[[1]].coefficients())   # a sub-Slice; coefficients() is always 2-D (1, 3)
    assert sub.shape == (1, 3)
    assert all(abs(complex(C[1, j]) - complex(sub[0, j])) < 1e-25 for j in range(3))


# ---- a slice as a System block -------------------------------------------------------------

def test_slice_add_to_system_agrees_with_slice_eval():
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    s = pb.Slice.from_coefficients([[2, 3, 1], [1, -1, 4]], [x, y])
    s.add_to(sys)
    assert sys.num_functions() == 2

    pt = _mpvec(1, 1)
    from_sys = sys.eval(pt)
    from_slice = s.eval(pt)
    assert abs(complex(from_sys[0]) - complex(from_slice[0])) < 1e-12
    assert abs(complex(from_sys[1]) - complex(from_slice[1])) < 1e-12


def test_add_to_rejects_variable_count_mismatch():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    s = Slice.random_complex(pb.VariableGroup([x, y, z]), 1)   # three variables
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))           # two variables
    with pytest.raises(RuntimeError):
        s.add_to(sys)


def test_add_slices_as_products_is_the_regen_bridge():
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    s1 = pb.Slice.from_coefficients([[1, 0, -1], [1, 0, 1]], [x, y])    # (x - 1)(x + 1)
    s2 = pb.Slice.from_coefficients([[0, 1, -1], [0, 1, -2]], [x, y])   # (y - 1)(y - 2)
    sys.add_slices_as_products([s1, s2])
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


# ---- Bertini 1 / classic emission ----------------------------------------------------------

def test_slice_to_classic_input():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.Slice.from_coefficients([[2, 3, 1]], [x, y])     # 2x + 3y + 1
    text = s.to_classic_input()
    assert 'CONFIG' in text and 'INPUT' in text
    assert text.count('END;') == 2
    assert 'x' in text and 'y' in text
    assert 'f0' in text


def test_concatenate_accepts_a_slice_system():
    # concatenate must append a slice's structured linear-forms block, not just polynomial rows.
    from bertini.system import concatenate
    x, y = pb.Variable('x'), pb.Variable('y')
    poly = pb.System()
    poly.add_variable_group(pb.VariableGroup([x, y]))
    poly.add_function(x * x + y * y - 1)

    slice_sys = pb.Slice.from_coefficients([[2, 3, 1], [1, -1, 4]], [x, y]).as_system()
    combined = concatenate(poly, slice_sys)
    assert combined.num_functions() == 3
    v = combined.eval(_mpvec(1, 1))     # f0=1, f1=6, f2=4 at (1,1)
    assert abs(complex(v[0]) - 1) < 1e-12
    assert abs(complex(v[1]) - 6) < 1e-12
    assert abs(complex(v[2]) - 4) < 1e-12


def test_witness_system_preserves_patch_on_a_projective_system():
    # A slice rides along with a system; it owns no homogenization.  For an already-projective,
    # patched system the witness (square) system must keep the system's patch and stay evaluable.
    x0, x1, x2 = pb.Variable('x0'), pb.Variable('x1'), pb.Variable('x2')
    vg = pb.VariableGroup([x0, x1, x2])

    sys = pb.System()
    sys.add_hom_variable_group(vg)                 # projective P^2
    sys.add_function(x0 * x0 + x1 * x1 - x2 * x2)  # a homogeneous conic
    sys.auto_patch()
    assert sys.is_patched()

    s = Slice.random_complex(vg, 1, True)          # one homogeneous hyperplane (no constant term)
    w = WitnessSetMultiplePrecision([], s, sys)

    wsys = w.witness_system()
    assert wsys.is_patched()                       # the patch carried through concatenate
    assert wsys.eval(_mpvec(1, 1, 1)) is not None  # the combined projective system still evaluates


def test_witness_system_homogenizes_an_affine_slice_for_a_homogenized_system():
    # Adding a slice to an already-homogenized system folds the slice's constant onto the
    # homogenizing variable (Slice.add_to is homogenization-aware), so witness_system stays
    # consistent even when the witness set's system was homogenized and patched.
    x, y = pb.Variable('x'), pb.Variable('y')
    vg = pb.VariableGroup([x, y])

    sys = pb.System()
    sys.add_variable_group(vg)
    sys.add_function(x * x + y * y - 1)
    sys.homogenize()                 # mints a homogenizing variable -> 3 variables
    sys.auto_patch()
    assert sys.is_patched()
    assert sys.num_variables() == 3

    s = pb.Slice.from_coefficients([[2, 3, 1]], [x, y])   # an affine slice over the 2 natural vars
    w = WitnessSetMultiplePrecision([], s, sys)

    wsys = w.witness_system()        # clone + add_to: the slice's constant folds onto the hom var
    assert wsys.is_patched()         # the system's patch carried through
    assert wsys.is_homogeneous()     # the appended slice form is homogeneous too -- it was folded
    assert wsys.eval(_mpvec(1, 1, 1)) is not None


def test_affine_slice_constant_rides_on_hom_var_after_homogenize():
    # The slice does not homogenize itself; the system does.  When a system carrying an affine slice
    # form is homogenized, the form's constant term becomes the coefficient on the homogenizing
    # variable (LinearFormsBlock::Homogenize) -- so the slice stays unaware of homogenization.
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.Slice.from_coefficients([[2, 3, 1]], [x, y]).as_system()   # 2x + 3y + 1

    # affine: at the origin the form is its constant term, 1.
    assert abs(complex(sys.eval(_mpvec(0, 0))[0]) - 1) < 1e-12

    sys.homogenize()
    assert sys.num_variables() == 3                # a homogenizing variable was added

    # now homogeneous of degree 1: a*x + b*y + c*h with c the old constant (1).  Across the three
    # coordinate axes the form takes the values {2, 3, 1}: the hom-var axis yields the old constant.
    axis_vals = [complex(sys.eval(_mpvec(*([1 if i == k else 0 for i in range(3)])))[0])
                 for k in range(3)]
    assert any(abs(v - 1) < 1e-12 for v in axis_vals)     # the constant rode onto the hom var
    assert any(abs(v - 2) < 1e-12 for v in axis_vals)     # x coefficient
    assert any(abs(v - 3) < 1e-12 for v in axis_vals)     # y coefficient
    # and the standalone constant is gone: the homogeneous form vanishes at the origin.
    assert abs(complex(sys.eval(_mpvec(0, 0, 0))[0])) < 1e-12


def test_witness_set_classic_emission():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    vg = pb.VariableGroup([x, y, z])
    sys = pb.System()
    sys.add_variable_group(vg)
    sys.add_function(x * x + y * y + z * z - 1)

    s = Slice.random_complex(vg, 2)
    w = WitnessSetMultiplePrecision([_mpvec(1, 0, 0), _mpvec(0, 1, 0)], s, sys)

    # the witness (square) system: 1 sphere function + 2 slice forms in 3 variables.
    wsys = w.witness_system()
    assert wsys.num_functions() == 3

    text = w.to_classic_input()
    assert 'CONFIG' in text and 'INPUT' in text
    assert text.count('END;') == 2


def test_witness_set_repr_is_readable():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    vg = pb.VariableGroup([x, y, z])
    sys = pb.System()
    sys.add_variable_group(vg)
    sys.add_function(x * x + y * y + z * z - 1)
    s = Slice.random_complex(vg, 2)
    w = WitnessSetMultiplePrecision([_mpvec(1, 0, 0), _mpvec(0, 1, 0)], s, sys)

    text = repr(w)
    assert 'WitnessSet' in text
    assert 'dimension 2' in text
    assert 'degree 2' in text
    assert 'consistent' in text
    assert 'slice' in text and 'system' in text
    # single vs plural: 1 function reads "1 function", not "1 functions"
    assert '1 function ' in text and '2 linear forms' in text


def test_empty_witness_set_repr():
    w = WitnessSetMultiplePrecision()
    text = repr(w)
    assert 'WitnessSet' in text
    assert 'degree 0' in text


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
