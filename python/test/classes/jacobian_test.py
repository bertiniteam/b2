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
