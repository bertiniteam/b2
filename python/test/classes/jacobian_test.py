# This file is part of Bertini 2.
#
# python/test/classes/jacobian_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/classes/jacobian_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/classes/jacobian_test.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""The symbolic Jacobian: bertini.jacobian (free) and System.jacobian (method)."""

import numpy as np
import pytest

import bertini as pb
from bertini._pybertini.function_tree import AbstractNode

from eval_helper import eval_at


def test_free_jacobian_shape_and_entries():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    J = pb.jacobian([x * y, x - z], [x, y, z])
    assert J.shape == (2, 3)
    assert all(isinstance(e, AbstractNode) for e in J.ravel())
    # d(xy)/dx = y, d(xy)/dy = x, d(xy)/dz = 0; d(x-z)/dx = 1, /dy = 0, /dz = -1
    assert str(J[0, 0]) == 'y'
    assert str(J[0, 1]) == 'x'
    assert eval_at(J[0, 2], x=2, y=3, z=5) == 0
    assert eval_at(J[1, 0], x=2, y=3, z=5) == 1
    assert eval_at(J[1, 2], x=2, y=3, z=5) == -1


def test_single_function_is_two_dimensional():
    x, y = pb.Variable('x'), pb.Variable('y')
    J = pb.jacobian(x * y, [x, y])      # a lone node, not a list
    assert J.shape == (1, 2)
    assert str(J[0, 0]) == 'y'
    assert str(J[0, 1]) == 'x'


def test_jacobian_accepts_variable_vector_and_numpy_functions():
    from bertini import linalg
    x = linalg.variable_vector('x', 3)            # numpy object array of Variables
    F = np.array([x[0] * x[1], x[2] ** 2], dtype=object)
    J = pb.jacobian(F, x)
    assert J.shape == (2, 3)
    assert eval_at(J[1, 2], x0=0, x1=0, x2=4) == 8     # d(x2^2)/dx2 = 2*x2 = 8


def test_system_jacobian_matches_free_function():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    f = x * (x ** 2 + y ** 2 - 1)
    g = z * ((y - 1) ** 2 + z ** 2 - 1)
    S = pb.System()
    S.add(pb.VariableGroup([x, y, z]), f, g)

    Js = S.jacobian()                              # usercoordinates=True default
    Jf = pb.jacobian([f, g], [x, y, z])
    assert Js.shape == Jf.shape == (2, 3)
    for a, b in zip(Js.ravel(), Jf.ravel()):
        assert str(a) == str(b)


def test_system_jacobian_matmul_nullvector():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    f = x * (x ** 2 + y ** 2 - 1)                  # df/dx = 3x^2 + y^2 - 1, df/dy = 2xy, df/dz = 0
    g = z * ((y - 1) ** 2 + z ** 2 - 1)
    S = pb.System()
    S.add(pb.VariableGroup([x, y, z]), f, g)
    J = S.jacobian()
    v = np.array([pb.Variable('v0'), pb.Variable('v1'), pb.Variable('v2')], dtype=object)
    prod = J @ v
    assert prod.shape == (2,)
    # row 0 . v at (x,y,z)=(2,3,5), v=(1,0,0) -> df/dx = 3*4 + 9 - 1 = 20
    assert eval_at(prod[0], x=2, y=3, z=5, v0=1, v1=0, v2=0) == 20


def test_usercoordinates_after_homogenize_has_no_homvar():
    x, y = pb.Variable('x'), pb.Variable('y')
    f = x ** 2 * y
    g = x * y - y ** 2
    S = pb.System()
    S.add(pb.VariableGroup([x, y]), f, g)
    S.homogenize()

    J = S.jacobian(usercoordinates=True)
    assert J.shape == (2, 2)                        # only x, y -- not the homogenizing variable
    assert all('HOM_VAR' not in str(e) for e in J.ravel())
    # entries equal the affine partials: d(x^2 y)/dx = 2xy = 12 at (2,3)
    assert eval_at(J[0, 0], x=2, y=3) == 12
    assert eval_at(J[0, 1], x=2, y=3) == 4          # d(x^2 y)/dy = x^2


def test_internal_coordinates_includes_homvar_and_patch_rows():
    x, y = pb.Variable('x'), pb.Variable('y')
    f = x ** 2 * y
    g = x * y - y ** 2
    S = pb.System()
    S.add(pb.VariableGroup([x, y]), f, g)
    S.homogenize()
    S.auto_patch()

    Juser = S.jacobian(usercoordinates=True)
    assert Juser.shape == (2, 2)                     # two functions, two user variables, no patch row

    Jint = S.jacobian(usercoordinates=False)
    assert Jint.shape[1] == S.num_variables()       # homogenizing variable is a column (3)
    assert Jint.shape[0] == Juser.shape[0] + 1      # the two functions plus one patch row (one variable group)


# ---- the symbolic Jacobian must work for every block type, not just polynomials -------------
#
# SymbolicJacobian differentiates NaturalFunctionsAsNodes(), which expands every block variant
# (polynomial, linear-forms / slice, products-of-linears, randomization, blend) to nodes.  The
# oracle below evaluates the symbolic Jacobian and checks it against the numeric eval_jacobian,
# so each block type is covered end to end.

def _assert_symbolic_matches_numeric(S, **point):
    """The symbolic user-coordinate Jacobian, evaluated, equals the numeric eval_jacobian."""
    order = [v.name for g in S.variable_groups() for v in g]
    vec = np.array([complex(point[name]) for name in order], dtype=complex)
    J = S.jacobian(usercoordinates=True)
    # numeric eval_jacobian collapses to 1-D for a single-function system; our symbolic Jacobian is
    # always 2-D, so reshape the numeric one to compare element-for-element.
    Jn = np.array(S.eval_jacobian(vec)).reshape(J.shape)
    for i in range(J.shape[0]):
        for j in range(J.shape[1]):
            got = complex(eval_at(J[i, j], **point))
            assert abs(got - complex(Jn[i, j])) < 1e-9, (i, j, got, Jn[i, j])


def test_jacobian_linear_forms_block():
    from bertini import linalg
    x, y = pb.Variable('x'), pb.Variable('y')
    S = pb.System()
    S.add_variable_group(pb.VariableGroup([x, y]))
    linalg.add_linear(S, np.array([[2, 1]]), np.array([x, y]), [-1])   # 2x + y - 1
    J = S.jacobian()
    assert J.shape == (1, 2)
    assert eval_at(J[0, 0], x=9, y=9) == 2 and eval_at(J[0, 1], x=9, y=9) == 1   # constant coefficients
    _assert_symbolic_matches_numeric(S, x=1.3, y=-0.7)


def test_jacobian_slice_block():
    from bertini import linalg
    x, y = pb.Variable('x'), pb.Variable('y')
    S = pb.System()
    S.add_variable_group(pb.VariableGroup([x, y]))
    sl = linalg.slice_from_coefficients([[2, 1, -1]], [x, y])          # a slice is a linear-forms block
    sl.add_to(S)
    J = S.jacobian()
    assert J.shape == (1, 2)
    _assert_symbolic_matches_numeric(S, x=0.4, y=2.1)


def test_jacobian_products_of_linears_block():
    from bertini import linalg
    x, y = pb.Variable('x'), pb.Variable('y')
    S = pb.System()
    S.add_variable_group(pb.VariableGroup([x, y]))
    linalg.add_products_of_linears(S, [
        [[1, 0, -1], [1, 0, 1]],     # (x - 1)(x + 1)
        [[0, 1, -1], [0, 1, -2]],    # (y - 1)(y - 2)
    ])
    J = S.jacobian()
    assert J.shape == (2, 2)
    assert eval_at(J[0, 0], x=3, y=5) == 6     # d/dx (x^2 - 1) = 2x = 6
    assert eval_at(J[0, 1], x=3, y=5) == 0
    _assert_symbolic_matches_numeric(S, x=3.0, y=5.0)


def test_jacobian_randomization_block():
    from bertini import linalg
    pb.random.set_random_seed(91)
    x, y = pb.Variable('x'), pb.Variable('y')
    S = pb.System()
    S.add_variable_group(pb.VariableGroup([x, y]))
    S.add_function(x * x + y * y - 1)
    S.add_function(x * y)
    S.add_function(x * x + y * y - x - y)       # overdetermined: 3 functions, 2 variables
    R = linalg.randomize(S)                       # a randomization block
    J = R.jacobian()
    assert J.shape == (2, 2)
    _assert_symbolic_matches_numeric(R, x=0.7, y=-1.2)


def test_jacobian_blend_block_homotopy():
    import bertini.system as bsys
    x, y = pb.Variable('x'), pb.Variable('y')
    target = pb.System()
    target.add_variable_group(pb.VariableGroup([x, y]))
    target.add_function(x * x - 1)
    target.add_function(y - 2)
    start = pb.System()
    start.add_variable_group(pb.VariableGroup([x, y]))
    start.add_function(x - 1)
    start.add_function(y - 1)
    H = bsys.make_homotopy(target, start)         # a blend block, with a path variable t

    J = H.jacobian(usercoordinates=True)          # space Jacobian: differentiates w.r.t. x, y only (not t)
    assert J.shape == (2, 2)
    # compare against the numeric space Jacobian at a fixed (space, time)
    order = [v.name for g in H.variable_groups() for v in g]
    pt = dict(x=0.5, y=1.5, t=0.3)
    vec = np.array([complex(pt[name]) for name in order], dtype=complex)
    Jn = np.array(H.eval_jacobian(vec, complex(pt['t'])))
    for i in range(2):
        for j in range(2):
            assert abs(complex(eval_at(J[i, j], **pt)) - complex(Jn[i, j])) < 1e-9
