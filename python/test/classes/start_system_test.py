# This file is part of Bertini 2.
#
# python/test/classes/start_system_test.py is free software: you can
# redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This file is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this file.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

# individual authors of this file include:
#   silviana amethyst

"""Tests for the TotalDegree start system.

Mirrors C++ tests in core/test/classes/start_system_test.cpp:
  - system_class/make_total_degree_system_linear
  - system_class/make_total_degree_system_quadratic
  - system_class/linear_total_degree_start_system
  - system_class/quadratic_cubic_quartic_total_degree_start_system
  - system_class/quadratic_cubic_quartic_start_points
  - system_class/quadratic_cubic_quartic_start_points_homogenized_patched
  - system_class/total_degree_start_system_precision_16
  - system_class/total_degree_start_system_homogenized_patched_precision_16
"""

import numpy as np
import pytest

import bertini
from bertini import System, VariableGroup
from bertini.function_tree.symbol import Variable
from bertini.function_tree import *
import bertini.system.start_system as ss

import bertini.multiprec as mp
from bertini.multiprec import Complex as mpfr_complex


@pytest.fixture
def linear_system():
    """Single-variable linear system: 2x - 3."""
    x = Variable("x")
    s = System()
    vg = VariableGroup(); vg.append(x)
    s.add_variable_group(vg)
    s.add_function(2 * x - 3)
    return s


@pytest.fixture
def quad_cubic_quartic_system():
    """3-variable system with degrees (2, 3, 4): 24 start points."""
    x = Variable("x"); y = Variable("y"); z = Variable("z")
    s = System()
    vg = VariableGroup(); vg.append(x); vg.append(y); vg.append(z)
    s.add_variable_group(vg)
    s.add_function(x**2 + y**2 - 1)
    s.add_function(x**3 - z)
    s.add_function(y**4 + z**2 - 2)
    return s, x, y, z


# ---------------------------------------------------------------------------
# num_start_points and degrees
# ---------------------------------------------------------------------------

def test_total_degree_linear(linear_system):
    """Linear system has exactly 1 start point.

    Mirrors start_system_test.cpp/make_total_degree_system_linear.
    """
    td = ss.TotalDegree(linear_system)
    assert td.num_start_points() == 1
    assert list(td.degrees()) == [1]


def test_total_degree_quad_cubic_quartic(quad_cubic_quartic_system):
    """Degrees (2,3,4) yield 24 start points.

    Mirrors start_system_test.cpp/make_total_degree_system_quadratic
    (generalised to the 3-variable case).
    """
    s, *_ = quad_cubic_quartic_system
    td = ss.TotalDegree(s)
    assert list(td.degrees()) == [2, 3, 4]
    assert td.num_start_points() == 24


# ---------------------------------------------------------------------------
# start_point_d — double precision start points
# ---------------------------------------------------------------------------

def test_start_points_evaluate_to_zero_double(quad_cubic_quartic_system):
    """All 24 double start points evaluate to ~0 on the TD system.

    Mirrors start_system_test.cpp/quadratic_cubic_quartic_start_points.
    """
    s, *_ = quad_cubic_quartic_system
    td = ss.TotalDegree(s)
    n = td.num_start_points()
    assert n == 24
    for i in range(n):
        p = td.start_point_d(i)
        ev = td.eval(p)
        assert np.linalg.norm(ev) < 1e-5, f"start_point_d({i}) residual {np.linalg.norm(ev)}"


def test_total_degree_jacobian_at_111(quad_cubic_quartic_system):
    """Jacobian diagonal at (1,1,1) equals the degree vector (2,3,4).

    Mirrors start_system_test.cpp/quadratic_cubic_quartic_total_degree_start_system.
    The TD Jacobian at (1,1,1) should be diag(2,3,4) since the TD functions
    are d_i * x_i^{d_i - 1}.
    """
    s, *_ = quad_cubic_quartic_system
    td = ss.TotalDegree(s)
    td.differentiate()
    pt = np.array([complex(1, 0)] * 3)
    J = td.eval_jacobian(pt)
    assert abs(J[0, 0] - complex(2, 0)) < 1e-14
    assert abs(J[1, 1] - complex(3, 0)) < 1e-14
    assert abs(J[2, 2] - complex(4, 0)) < 1e-14


# ---------------------------------------------------------------------------
# start_point_mp — multiprecision start points
# ---------------------------------------------------------------------------

def test_start_points_evaluate_to_zero_mp(quad_cubic_quartic_system):
    """All 24 mp start points evaluate to ~0 on the TD system.

    Mirrors start_system_test.cpp/quadratic_cubic_quartic_start_points (mp variant).
    """
    s, *_ = quad_cubic_quartic_system
    td = ss.TotalDegree(s)
    n = td.num_start_points()
    for i in range(n):
        p = td.start_point_mp(i)
        ev = td.eval(p)
        norm = max(mp.abs(v) for v in ev)
        assert norm < 1e-5, f"start_point_mp({i}) residual {norm}"


# ---------------------------------------------------------------------------
# Homogenized + patched total degree system
# ---------------------------------------------------------------------------

def test_total_degree_homogenized_patched(quad_cubic_quartic_system):
    """After homogenize + auto_patch the TD system is hom. and patched.

    Mirrors start_system_test.cpp/quadratic_cubic_quartic_start_points_homogenized_patched.
    """
    s, *_ = quad_cubic_quartic_system
    s.homogenize()
    s.auto_patch()
    assert s.is_homogeneous()
    assert s.is_patched()

    td = ss.TotalDegree(s)
    td.homogenize()
    assert td.is_homogeneous()
    assert td.is_patched()

    n = td.num_start_points()
    assert n == 24
    for i in range(n):
        p = td.start_point_d(i)
        ev = td.eval(p)
        assert np.linalg.norm(ev) < 1e-5, f"start_point_d({i}) residual {np.linalg.norm(ev)}"


# ---------------------------------------------------------------------------
# Precision of start points
# ---------------------------------------------------------------------------

def test_start_point_precision_16(quad_cubic_quartic_system):
    """All start points inherit the current default precision (16).

    Mirrors start_system_test.cpp/total_degree_start_system_precision_16.
    """
    bertini.default_precision(16)

    s, *_ = quad_cubic_quartic_system
    td = ss.TotalDegree(s)
    n = td.num_start_points()
    for i in range(n):
        p = td.start_point_mp(i)
        assert all(c.precision == 16 for c in p), (
            f"start_point_mp({i}) has wrong precision: {[c.precision for c in p]}"
        )
