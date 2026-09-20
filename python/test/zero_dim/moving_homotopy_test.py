"""Moving-row homotopies: move only the slices/regenerated rows, evaluate the fixed system once.

A regeneration / moving-slice homotopy keeps the polynomial system and any static linear slices
FIXED -- they are their own evaluation blocks, evaluated once per point and contributing zero to
dH/dt -- while only the moving rows (a sliding linear slice, or a products-of-linears deforming into
a polynomial) carry the path variable.  nag_algorithm.moving_homotopy builds it; this verifies the
solutions land where hand computation says, and that the fixed rows really are left out of dH/dt.
"""

import numpy as np
import pytest

import bertini as pb
from bertini import nag_algorithm as na

C = pb.multiprec.complex_mp
GAMMA = C('0.6', '0.8')   # exact, off the real axis -> reproducible path


def _vg(*vs):
    return pb.VariableGroup(list(vs))


def _roots(solutions, n):
    return sorted(tuple(round(complex(s[i]).real, 4) for i in range(n)) for s in solutions)


def _solve(fixed, start_moving, end_moving, start_points, nvars):
    # affine on purpose: these cases are about the BLEND structure -- which rows move, which stay
    # out of dH/dt -- and projectivizing (the builder's default, b2#382) would add a patch row and
    # a homogenizing coordinate to every count below without changing what is under test
    b = na.straight_line_homotopy(end_moving, start_moving, fixed=fixed,
                                  gamma=pb.coefficient(GAMMA), projectivize=False)
    solver = na.user_homotopy(b.homotopy, start_points, b.target)
    solver.solve()
    return b.homotopy, _roots(solver.all_solutions(), nvars)


def _fixed_rows_have_zero_dHdt(H, num_fixed_rows, num_vars):
    """dH/dt must be exactly zero on every fixed (leading) row -- the fixed system is left out."""
    p = np.array([C(str(0.3 + 0.1 * k)) for k in range(num_vars)])
    dhdt = H.eval_time_derivative(p, C('0.5'))
    return all(abs(complex(dhdt[i])) == 0.0 for i in range(num_fixed_rows))


def test_move_one_slice_2var():
    # fixed unit circle; slice moves from the x-axis (y=0) to the diagonal (y-x=0).
    x, y = pb.Variable('x'), pb.Variable('y')
    fixed = pb.System(); fixed.add_variable_group(_vg(x, y)); fixed.add_function(x*x + y*y - 1)
    start_moving = pb.System(); start_moving.add_variable_group(_vg(x, y)); start_moving.add_function(y)
    end_moving = pb.System(); end_moving.add_variable_group(_vg(x, y)); end_moving.add_function(y - x)

    sp = [np.array([C('1'), C('0')]), np.array([C('-1'), C('0')])]   # circle ∩ x-axis
    H, got = _solve(fixed, start_moving, end_moving, sp, 2)

    r = round(1 / np.sqrt(2), 4)
    assert got == sorted([(r, r), (-r, -r)])           # circle ∩ diagonal
    assert H.num_functions() == 2
    assert _fixed_rows_have_zero_dHdt(H, num_fixed_rows=1, num_vars=2)   # circle row out of dH/dt


def test_static_and_moving_slice_3var():
    # fixed unit sphere + a STATIC slice z=0; a second slice moves y=0 -> y-x=0.  Two fixed blocks.
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    fixed = pb.System(); fixed.add_variable_group(_vg(x, y, z))
    fixed.add_function(x*x + y*y + z*z - 1)                                  # sphere (PolynomialBlock)
    fixed.add_linear(np.array([[0, 0, 1]]), np.array([x, y, z]))    # static slice z=0 (LinearFormsBlock)
    start_moving = pb.System(); start_moving.add_variable_group(_vg(x, y, z)); start_moving.add_function(y)
    end_moving = pb.System(); end_moving.add_variable_group(_vg(x, y, z)); end_moving.add_function(y - x)

    sp = [np.array([C('1'), C('0'), C('0')]), np.array([C('-1'), C('0'), C('0')])]
    H, got = _solve(fixed, start_moving, end_moving, sp, 3)

    r = round(1 / np.sqrt(2), 4)
    assert got == sorted([(r, r, 0.0), (-r, -r, 0.0)])
    assert H.num_functions() == 3
    # BOTH fixed rows (sphere AND static slice) are left out of dH/dt
    assert _fixed_rows_have_zero_dHdt(H, num_fixed_rows=2, num_vars=3)


def test_deform_products_of_linears_into_polynomial():
    # the regeneration "add a degree" step: deform the product (x-1)(x+1) into the circle, while a
    # static slice y=1/2 stays put.  start points are the product's roots on the slice: (±1, 1/2).
    x, y = pb.Variable('x'), pb.Variable('y')
    fixed = pb.System(); fixed.add_variable_group(_vg(x, y))
    fixed.add_linear(np.array([[0, 1]]), np.array([x, y]), ['-1/2'])     # static slice y=1/2
    start_moving = pb.System(); start_moving.add_variable_group(_vg(x, y))
    start_moving.add_products_of_linears([[[1, 0, -1], [1, 0, 1]]])      # (x-1)(x+1) as a block
    end_moving = pb.System(); end_moving.add_variable_group(_vg(x, y)); end_moving.add_function(x*x + y*y - 1)

    sp = [np.array([C('1'), C('0.5')]), np.array([C('-1'), C('0.5')])]
    H, got = _solve(fixed, start_moving, end_moving, sp, 2)

    s3 = round(np.sqrt(3) / 2, 4)
    assert got == sorted([(s3, 0.5), (-s3, 0.5)])      # circle ∩ {y=1/2}
    assert H.num_functions() == 2
    assert _fixed_rows_have_zero_dHdt(H, num_fixed_rows=1, num_vars=2)   # static slice row out of dH/dt


def test_the_builder_projectivizes_by_default():
    # b2#382.  A homotopy is made projective where it is BUILT -- once the blend exists its
    # operands are fixed and cannot be homogenized -- so the option lives here rather than on the
    # solver.  The default is on: infinity should be an ordinary place without anybody asking.
    x, y = pb.Variable('x'), pb.Variable('y')
    fixed = pb.System(); fixed.add_variable_group(_vg(x, y)); fixed.add_function(x*x + y*y - 1)
    start_moving = pb.System(); start_moving.add_variable_group(_vg(x, y)); start_moving.add_function(y)
    end_moving = pb.System(); end_moving.add_variable_group(_vg(x, y)); end_moving.add_function(y - x)

    b = na.straight_line_homotopy(end_moving, start_moving, fixed=fixed,
                                  gamma=pb.coefficient(GAMMA))

    assert b.homotopy.is_homogeneous()
    assert b.homotopy.num_variables() == 3         # x, y, and the homogenizing coordinate
    assert b.homotopy.num_functions() == 3         # the two rows, plus the patch
    # the CALLER's systems are untouched: a System's content is its identity, and records key a
    # solve on it, so converting one in place would change what their own system is
    assert not fixed.is_homogeneous() and not end_moving.is_homogeneous()

    # and it still solves, to the same roots as the affine spelling -- the builder hands back the
    # end systems so nothing is concatenated by hand, and the solver lifts the affine start points
    sp = [np.array([C('1'), C('0')]), np.array([C('-1'), C('0')])]
    solver = na.user_homotopy(b.homotopy, sp, b.target)
    solver.solve()
    r = round(1 / np.sqrt(2), 4)
    assert _roots(solver.finite_solutions(), 2) == sorted([(r, r), (-r, -r)])


def test_projectivize_false_builds_the_homotopy_as_written():
    x, y = pb.Variable('x'), pb.Variable('y')
    fixed = pb.System(); fixed.add_variable_group(_vg(x, y)); fixed.add_function(x*x + y*y - 1)
    sm = pb.System(); sm.add_variable_group(_vg(x, y)); sm.add_function(y)
    em = pb.System(); em.add_variable_group(_vg(x, y)); em.add_function(y - x)

    b = na.straight_line_homotopy(em, sm, fixed=fixed, gamma=pb.coefficient(GAMMA),
                                  projectivize=False)

    assert b.homotopy.num_variables() == 2
    assert b.homotopy.num_functions() == 2
    assert not fixed.is_homogeneous()              # the caller's systems are left alone


def test_mixing_projective_and_affine_systems_is_refused():
    # A homotopy between a projective system and an affine one deforms between points that do not
    # correspond, and nothing downstream notices -- the shapes can agree and the tracking runs.
    x, y = pb.Variable('x'), pb.Variable('y')
    fixed = pb.System(); fixed.add_variable_group(_vg(x, y)); fixed.add_function(x*x + y*y - 1)
    sm = pb.System(); sm.add_variable_group(_vg(x, y)); sm.add_function(y)
    em = pb.System(); em.add_variable_group(_vg(x, y)); em.add_function(y - x)

    sm.homogenize()                                # one of the three, in other coordinates
    with pytest.raises(RuntimeError, match=r"projective and the other is affine"):
        na.straight_line_homotopy(em, sm, fixed=fixed)
    with pytest.raises(RuntimeError, match=r"projective and the other is affine"):
        na.straight_line_homotopy(em, sm, fixed=fixed, projectivize=False)   # refused either way

    # and with no held rows at all
    T = pb.System(); T.add_variable_group(_vg(x, y)); T.add_functions([x*y - 1, x - 1])
    S = pb.System(); S.add_variable_group(_vg(x, y)); S.add_functions([x*x - 1, y - 1])
    S.homogenize()
    with pytest.raises(RuntimeError, match=r"projective and the other is affine"):
        na.straight_line_homotopy(T, S)


def test_held_and_moving_rows_must_agree_in_count():
    x, y = pb.Variable('x'), pb.Variable('y')
    fixed = pb.System(); fixed.add_variable_group(_vg(x, y)); fixed.add_function(x*x + y*y - 1)
    sm = pb.System(); sm.add_variable_group(_vg(x, y)); sm.add_function(y)
    em = pb.System(); em.add_variable_group(_vg(x, y)); em.add_function(y - x); em.add_function(x)  # 2 != 1
    with pytest.raises(RuntimeError):
        na.straight_line_homotopy(em, sm, fixed=fixed)


def test_rows_may_be_given_as_lists_of_functions():
    # b2#382: with fixed= supplying the variable structure, the moving rows need not be dressed up
    # as Systems -- the common case is one or two functions, and building a System for them is
    # ceremony.
    x, y = pb.Variable('x'), pb.Variable('y')
    fixed = pb.System(); fixed.add_variable_group(_vg(x, y)); fixed.add_function(x*x + y*y - 1)

    b = na.straight_line_homotopy([y - x], [y], fixed=fixed, gamma=pb.coefficient(GAMMA))

    assert b.target.num_functions() == b.homotopy.num_functions()
    assert b.start.num_functions() == b.homotopy.num_functions()
    assert b.fixed is not None

    sp = [np.array([C('1'), C('0')]), np.array([C('-1'), C('0')])]
    solver = na.user_homotopy(b.homotopy, sp, b.target)
    solver.solve()
    r = round(1 / np.sqrt(2), 4)
    assert _roots(solver.finite_solutions(), 2) == sorted([(r, r), (-r, -r)])


def test_rows_as_lists_need_a_variable_structure():
    x, y = pb.Variable('x'), pb.Variable('y')
    with pytest.raises(ValueError, match=r"variable structure"):
        na.straight_line_homotopy([y - x], [y])


def test_the_start_system_is_what_you_solve_for_start_points():
    # the builder hands back both ends, so neither concatenation is the caller's to get right
    x, y = pb.Variable('x'), pb.Variable('y')
    fixed = pb.System(); fixed.add_variable_group(_vg(x, y)); fixed.add_function(x*x + y*y - 1)

    b = na.straight_line_homotopy([y - x], [y], fixed=fixed, gamma=pb.coefficient(GAMMA))

    pts = na.ZeroDimSolver(b.start).solve().solutions          # circle meets y = 0
    assert len(pts) == 2
    solver = na.user_homotopy(b.homotopy, pts, b.target)
    solver.solve()
    r = round(1 / np.sqrt(2), 4)
    assert _roots(solver.finite_solutions(), 2) == sorted([(r, r), (-r, -r)])


def test_the_record_remembers_the_gamma_it_drew():
    # .gamma is the coefficient ACTUALLY used, including the random default -- without it a path
    # drawn by the builder could not be reproduced
    x, y = pb.Variable('x'), pb.Variable('y')
    fixed = pb.System(); fixed.add_variable_group(_vg(x, y)); fixed.add_function(x*x + y*y - 1)

    b = na.straight_line_homotopy([y - x], [y], fixed=fixed)
    assert b.gamma is not None
    assert b.path_variable == 't'

    again = na.straight_line_homotopy([y - x], [y], fixed=fixed, gamma=b.gamma)
    assert str(again.homotopy) == str(b.homotopy)
