"""UI ergonomics regressions: variable/group construction (#293), System accessors (#296/#297),
and node.eval accepting a point/dict (#300)."""

import numpy as np
import pytest

import bertini as pb
from bertini.multiprec import complex_mp


# --- #293: variables([...]) list form ---------------------------------------------------------

def test_variables_explicit_names_list():
    x, y, z = pb.variables(['x', 'y', 'z'])
    assert [str(v) for v in (x, y, z)] == ['x', 'y', 'z']


def test_variables_integer_indexed_still_works():
    vs = pb.variables('v', 3)
    assert [str(v) for v in vs] == ['v0', 'v1', 'v2']


def test_variables_prefix_without_count_is_a_clear_error():
    with pytest.raises(TypeError):
        pb.variables('x')


# --- #293: add_variable_group accepts several forms -------------------------------------------

@pytest.mark.parametrize("build", [
    lambda x, y, z: (lambda s: s.add_variable_group(x, y, z)),          # loose variadic
    lambda x, y, z: (lambda s: s.add_variable_group([x, y, z])),        # a list
    lambda x, y, z: (lambda s: s.add_variable_group(pb.VariableGroup([x, y, z]))),  # explicit
])
def test_add_variable_group_forms(build):
    x, y, z = pb.variables(['x', 'y', 'z'])
    sys = pb.System()
    build(x, y, z)(sys)
    assert sys.num_variable_groups() == 1
    assert sys.num_variables() == 3


def test_add_variable_group_single_variable():
    (x,) = pb.variables(['x'])
    sys = pb.System()
    sys.add_variable_group(x)
    assert sys.num_variable_groups() == 1


def test_add_variable_group_rejects_mixed_variadic():
    x, y, z = pb.variables(['x', 'y', 'z'])
    sys = pb.System()
    with pytest.raises(TypeError):
        sys.add_variable_group(x, [y, z])


# --- #293: Slice.random_* clear error on a list of groups -------------------------------------

def test_slice_random_unwraps_single_group_sequence():
    x, y, z = pb.variables(['x', 'y', 'z'])
    sys = pb.System()
    sys.add_variable_group([x, y, z])
    sys.add_functions([x * x + y * y + z * z - 1])
    s = pb.Slice.random_complex(sys.variable_groups(), 1)   # a one-group sequence unwraps
    assert s.dimension() == 1


def test_slice_random_multi_group_sequence_raises_clear_error():
    x, y, u, v = pb.variables(['x', 'y', 'u', 'v'])
    sys = pb.System()
    sys.add_variable_group([x, y])
    sys.add_variable_group([u, v])
    with pytest.raises(TypeError, match="ONE variable group"):
        pb.Slice.random_complex(sys.variable_groups(), 1)


# --- #297 / #296: System.functions(), copy_functions(), clone() ------------------------------

def _sphere_system():
    x, y, z = pb.variables(['x', 'y', 'z'])
    sys = pb.System()
    sys.add_variable_group([x, y, z])
    sys.add_functions([x * x + y * y + z * z - 1, x - y])
    return sys, (x, y, z)


def test_functions_returns_all_functions():
    sys, _ = _sphere_system()
    fns = sys.functions()
    assert len(fns) == sys.num_functions() == 2


def test_copy_functions_builds_from_another_system():
    sys, (x, y, z) = _sphere_system()
    other = pb.System()
    other.copy_functions(sys)
    assert len(other.functions()) == 2
    # add_functions(functions()) is the equivalent long form
    other2 = pb.System()
    other2.add_functions(sys.functions())
    assert len(other2.functions()) == 2


def test_clone_returns_extendable_system():
    sys, (x, y, z) = _sphere_system()
    c = sys.clone()
    assert len(c.functions()) == 2
    c.add_function(x + y + z)          # clone is independently extendable
    assert len(c.functions()) == 3
    assert len(sys.functions()) == 2   # original untouched


# --- #300: node.eval accepts a point array and a dict ----------------------------------------

def test_eval_accepts_array_dict_and_kwargs_consistently():
    x, y, z = pb.variables(['x', 'y', 'z'])
    f = x * x + y - z
    pt = [complex_mp(2), complex_mp(3), complex_mp(4)]   # variables() order is alphabetical: x, y, z
    by_kwargs = complex(f.eval(x=pt[0], y=pt[1], z=pt[2]))
    by_array = complex(f.eval(pt))
    by_nparray = complex(f.eval(np.array(pt, dtype=object)))
    by_dict = complex(f.eval({x: pt[0], y: pt[1], z: pt[2]}))
    by_dict_names = complex(f.eval({'x': pt[0], 'y': pt[1], 'z': pt[2]}))
    assert by_kwargs == by_array == by_nparray == by_dict == by_dict_names == (2 * 2 + 3 - 4)


def test_eval_returns_mp_native_not_float():
    x = pb.Variable('x')
    f = x * x
    out = f.eval([complex_mp(3)])
    assert isinstance(out, complex_mp)           # NOT a python float -- the projection through-line
    assert isinstance(out.real, type(out.real))  # .real is real_mp


def test_eval_array_wrong_length_raises():
    x, y = pb.variables(['x', 'y'])
    f = x + y
    with pytest.raises(RuntimeError):
        f.eval([complex_mp(1)])                   # 1 value, 2 variables
