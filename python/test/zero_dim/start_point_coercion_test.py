"""Start points arrive as whatever numbers the user has; the solver seam coerces them (issue #347).

A HomotopySolver's start points are TRANSPORTED values -- the tracker refines them -- so any
faithful numeric representation must be accepted: multiprecision scalars, Python/numpy numbers,
and constant symbolic nodes.  Before the coercion seam, an ndarray of ``symbolics.Complex`` nodes
(the natural thing to build when your coefficients are already exact nodes) crashed with a raw
eigenpy converter TypeError naming a C++ Eigen type; a complex128 array died the same way.

Non-numbers are refused with errors that say what to do: an expression still containing
variables is a math error (it is not a number), anything else names the accepted kinds.
"""

import numpy as np
import pytest

import bertini as pb
from bertini import nag_algorithm as na

C = pb.multiprec.complex_mp
GAMMA = pb.symbolics.Complex('0.6', '0.8')   # exact, off the real axis -> reproducible path


def _quadratic_homotopy():
    """H deforming x^2-4 (roots +-2) into x^2-1 (roots +-1); returns (homotopy, target)."""
    x, t = pb.Variable('x'), pb.Variable('t')
    target = pb.System()
    target.add_variable_group(pb.VariableGroup([x]))
    target.add_function(x * x - 1)
    hom = pb.System()
    hom.add_variable_group(pb.VariableGroup([x]))
    one = pb.symbolics.Integer(1)
    hom.add_function((one - t) * (x * x - 1) + GAMMA * t * (x * x - 4))
    hom.add_path_variable(t)
    return hom, target


def _solved_roots(start_points):
    hom, target = _quadratic_homotopy()
    solver = na.HomotopySolver(hom, start_points, target)
    solver.solve()
    return sorted(round(complex(s[0]).real, 6) for s in solver.all_solutions())


def test_start_points_as_mp_arrays():
    # the previously-working spelling keeps working
    assert _solved_roots([np.array([C('2')]), np.array([C('-2')])]) == [-1.0, 1.0]


def test_start_points_as_constant_symbolic_nodes():
    # issue #347: ndarrays of symbolics.Complex crashed with a raw eigenpy converter TypeError
    sp = [np.array([pb.symbolics.Complex('2', '0')]),
          np.array([pb.symbolics.Complex('-2', '0')])]
    assert _solved_roots(sp) == [-1.0, 1.0]


def test_start_points_as_complex128_arrays():
    # doubles are a faithful transport format for start points; the tracker refines them
    sp = [np.array([2.0 + 0.0j]), np.array([-2.0 + 0.0j])]
    assert _solved_roots(sp) == [-1.0, 1.0]


def test_start_points_as_plain_python_lists():
    assert _solved_roots([[2.0], [-2]]) == [-1.0, 1.0]


def test_start_points_of_mixed_kinds():
    # one point per representation, and mixed coordinates within a point elsewhere: all one seam
    sp = [np.array([pb.symbolics.Complex('2', '0')]), [C('-2')]]
    assert _solved_roots(sp) == [-1.0, 1.0]


def test_start_point_with_variables_is_a_math_error():
    x = pb.Variable('x')
    with pytest.raises(ValueError, match=r"symbolic expression.*variable.*x"):
        _solved_roots([np.array([x + 2])])


def test_start_point_junk_names_the_accepted_kinds():
    with pytest.raises(TypeError, match=r"complex_mp.*numpy numbers.*constant symbolic"):
        _solved_roots([["two"]])


def test_start_point_must_be_a_vector():
    with pytest.raises(TypeError, match=r"vector of coordinates"):
        _solved_roots([2.0, -2.0])
