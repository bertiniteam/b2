"""The unified System.add(*objects) builder.

One fluent verb that dispatches by type: function-tree expressions become functions,
VariableGroups become variable groups, arrays/lists of those are added element-wise.
It is sugar over add_function / add_variable_group -- a System built with .add must be
identical to one built the explicit way.
"""

import numpy as np
import pytest

import bertini as pb
from bertini import linalg


def test_add_expression_and_group_with_chaining():
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    out = sys.add(pb.VariableGroup([x, y])).add(x**2 + y**2 - 1).add(x + y)
    assert out is sys                       # returns self for chaining
    assert sys.num_functions() == 2
    assert sys.num_variable_groups() == 1


def test_add_array_of_expressions():
    xs = linalg.variable_vector('x', 2)
    A = np.array([[2, 1], [0, 3]])
    sys = pb.System()
    sys.add(pb.VariableGroup(list(xs)))
    sys.add(A @ xs - 1)                      # numpy object array of expressions
    assert sys.num_functions() == 2


def test_add_mixed_varargs():
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add(pb.VariableGroup([x, y]), x * y - 1, x + y)   # group + two exprs in one call
    assert sys.num_functions() == 2
    assert sys.num_variable_groups() == 1


def test_add_list_of_expressions():
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add(pb.VariableGroup([x, y]))
    assert sys.add([x * y - 1, x + y]) is sys
    assert sys.num_functions() == 2


def test_add_matches_explicit_construction():
    # the same system, built with .add vs the explicit methods, evaluates identically.
    def via_add():
        x, y = pb.Variable('x'), pb.Variable('y')
        s = pb.System()
        s.add(pb.VariableGroup([x, y]), x**2 + y**2 - 1, x + y)
        return s

    def via_explicit():
        x, y = pb.Variable('x'), pb.Variable('y')
        s = pb.System()
        s.add_variable_group(pb.VariableGroup([x, y]))
        s.add_function(x**2 + y**2 - 1)
        s.add_function(x + y)
        return s

    pt = np.array([1 + 0j, 1 + 0j])
    va = via_add().eval(pt)
    ve = via_explicit().eval(pt)
    assert len(va) == len(ve) == 2
    for i in range(2):
        assert abs(complex(va[i]) - complex(ve[i])) < 1e-12


@pytest.mark.parametrize("bad", [3.5, "x", 7, None])
def test_add_rejects_unsupported_types(bad):
    with pytest.raises(TypeError):
        pb.System().add(bad)
