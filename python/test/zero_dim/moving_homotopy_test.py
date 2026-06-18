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
from bertini import linalg, nag_algorithm as na

C = pb.multiprec.Complex
GAMMA = C('0.6', '0.8')   # exact, off the real axis -> reproducible path


def _vg(*vs):
    return pb.VariableGroup(list(vs))


def _roots(solutions, n):
    return sorted(tuple(round(complex(s[i]).real, 4) for i in range(n)) for s in solutions)


def _solve(fixed, start_moving, end_moving, start_points, nvars):
    H = na.moving_homotopy(fixed, start_moving, end_moving, gamma=linalg.coefficient(GAMMA))
    target = pb.system.concatenate(fixed, end_moving)
    solver = na.user_homotopy(H, start_points, target)
    solver.solve()
    return H, _roots(solver.solutions(), nvars)


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
    linalg.add_linear(fixed, np.array([[0, 0, 1]]), np.array([x, y, z]))    # static slice z=0 (LinearFormsBlock)
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
    linalg.add_linear(fixed, np.array([[0, 1]]), np.array([x, y]), ['-1/2'])     # static slice y=1/2
    start_moving = pb.System(); start_moving.add_variable_group(_vg(x, y))
    linalg.add_products_of_linears(start_moving, [[[1, 0, -1], [1, 0, 1]]])      # (x-1)(x+1) as a block
    end_moving = pb.System(); end_moving.add_variable_group(_vg(x, y)); end_moving.add_function(x*x + y*y - 1)

    sp = [np.array([C('1'), C('0.5')]), np.array([C('-1'), C('0.5')])]
    H, got = _solve(fixed, start_moving, end_moving, sp, 2)

    s3 = round(np.sqrt(3) / 2, 4)
    assert got == sorted([(s3, 0.5), (-s3, 0.5)])      # circle ∩ {y=1/2}
    assert H.num_functions() == 2
    assert _fixed_rows_have_zero_dHdt(H, num_fixed_rows=1, num_vars=2)   # static slice row out of dH/dt


def test_moving_homotopy_rejects_mismatched_endpoints():
    x, y = pb.Variable('x'), pb.Variable('y')
    fixed = pb.System(); fixed.add_variable_group(_vg(x, y)); fixed.add_function(x*x + y*y - 1)
    sm = pb.System(); sm.add_variable_group(_vg(x, y)); sm.add_function(y)
    em = pb.System(); em.add_variable_group(_vg(x, y)); em.add_function(y - x); em.add_function(x)  # 2 != 1
    with pytest.raises(RuntimeError):
        na.moving_homotopy(fixed, sm, em)
