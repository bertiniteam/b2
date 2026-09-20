"""The Python boundary refuses what the native layer cannot survive, and accepts what it should.

Every case here used to end in one of three bad ways: the process aborting from inside the
tracker (a start point of the wrong length, an overdetermined homotopy -- issues #369, #383),
heap corruption inside Eigen (a non-square matrix handed to the multiprecision LU -- #390), or a
raw Boost.Python argument error for an input that is obviously fine (a list as an evaluation
point -- #367; a Slice where a System is expected -- #372, #381; a config NAME where
get_config wanted the class -- #392).  These are interface tests; the arithmetic is covered in C++.
"""

import numpy as np
import pytest

import bertini as pb
from bertini import nag_algorithm as na
from bertini.nag_algorithm import TolerancesConfig


def _circle_and_line():
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_function(x * x + y * y - 1)
    sys.add_function(x - y)
    return x, y, sys


def _straight_homotopy():
    """H(x, t) = x^2 - (9 - 5 t): roots +/-2 at t = 1, +/-3 at t = 0."""
    x, t = pb.Variable('x'), pb.Variable('t')
    H = pb.System()
    H.add_variable_group(pb.VariableGroup([x]))
    H.add_function(x * x - (9 - 5 * t))
    H.add_path_variable(t)
    target = pb.System()
    target.add_variable_group(pb.VariableGroup([x]))
    target.add_function(x * x - 9)
    return H, target


# --- #369 / #383: shape mismatches are Python exceptions, never an abort ------------------------

def test_start_point_of_the_wrong_length_is_a_value_error():
    H, target = _straight_homotopy()
    with pytest.raises(ValueError, match="coordinate"):
        na.HomotopySolver(H, [[2, 0]], target)          # two coordinates for a one-variable homotopy


def test_overdetermined_homotopy_is_a_value_error():
    x, t = pb.Variable('x'), pb.Variable('t')
    H = pb.System()
    H.add_variable_group(pb.VariableGroup([x]))
    H.add_function(x * x - (9 - 5 * t))
    H.add_function(x - 3 * (1 - t) - 2 * t)             # a second equation in one variable
    H.add_path_variable(t)
    target = pb.System()
    target.add_variable_group(pb.VariableGroup([x]))
    target.add_function(x * x - 9)
    with pytest.raises(ValueError, match="SQUARE"):
        na.HomotopySolver(H, [[2]], target)


def test_well_shaped_homotopy_still_solves():
    H, target = _straight_homotopy()
    solver = na.HomotopySolver(H, [[2], [-2]], target, mptype='double')
    solver.solve()
    roots = sorted(round(complex(s[0]).real, 6) for s in solver.all_solutions())
    assert roots == [-3.0, 3.0]


# --- #390: a non-square system to linalg.solve raises, for every dtype alike --------------------

@pytest.mark.parametrize("make", [
    lambda: (np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]), np.array([1.0, 2.0])),
    lambda: (np.array([[pb.multiprec.complex_mp(1), pb.multiprec.complex_mp(2), pb.multiprec.complex_mp(3)],
                       [pb.multiprec.complex_mp(4), pb.multiprec.complex_mp(5), pb.multiprec.complex_mp(6)]]),
             np.array([pb.multiprec.complex_mp(1), pb.multiprec.complex_mp(2)])),
])
def test_linalg_solve_refuses_a_non_square_matrix(make):
    A, b = make()
    with pytest.raises(np.linalg.LinAlgError, match="square"):
        pb.linalg.solve(A, b)


def test_linalg_solve_refuses_a_right_hand_side_that_does_not_fit():
    A = np.array([[pb.multiprec.complex_mp(2), pb.multiprec.complex_mp(0)],
                  [pb.multiprec.complex_mp(0), pb.multiprec.complex_mp(2)]])
    b = np.array([pb.multiprec.complex_mp(1)])
    with pytest.raises(np.linalg.LinAlgError, match="right-hand side"):
        pb.linalg.solve(A, b)


def test_native_solve_guard_is_a_value_error_not_a_crash():
    """The native layer refuses too (Boost.Python maps invalid_argument to ValueError), so the
    Python check is a convenience, not the only thing between a wrong shape and the heap."""
    from bertini._pybertini.linalg import solve as native_solve
    A = np.array([[pb.multiprec.complex_mp(1), pb.multiprec.complex_mp(2), pb.multiprec.complex_mp(3)],
                  [pb.multiprec.complex_mp(4), pb.multiprec.complex_mp(5), pb.multiprec.complex_mp(6)]])
    b = np.array([pb.multiprec.complex_mp(1), pb.multiprec.complex_mp(2)])
    with pytest.raises(ValueError, match="square"):
        native_solve(A, b)


# --- #367: a plain list is an evaluation point -------------------------------------------------

def test_system_eval_accepts_a_list():
    x, y, sys = _circle_and_line()
    from_list = sys.eval([0.3, 0.7])
    from_array = sys.eval(np.array([0.3 + 0j, 0.7 + 0j]))
    assert np.allclose(from_list, from_array)
    assert from_list.dtype == np.complex128

    mp_from_list = sys.eval([pb.multiprec.complex_mp('0.3'), pb.multiprec.complex_mp('0.7')])
    assert str(mp_from_list.dtype) == 'complex_mp'
    assert abs(complex(mp_from_list[0]) - complex(from_list[0])) < 1e-12


# --- #372 / #381: a Slice is accepted where a System is expected --------------------------------

def test_system_add_accepts_a_slice():
    x, y, sys = _circle_and_line()
    s = pb.Slice.from_coefficients([[2, 1, -1]], [x, y])     # 2x + y - 1 = 0
    n = sys.num_functions()
    sys.add(s)
    assert sys.num_functions() == n + 1
    pt = np.array([0.25 + 0j, 0.5 + 0j])
    assert np.isclose(sys.eval(pt)[-1], 2 * 0.25 + 0.5 - 1)


def test_the_builder_accepts_slices_for_the_moving_rows():
    x, y = pb.Variable('x'), pb.Variable('y')
    fixed = pb.System()
    fixed.add_variable_group(pb.VariableGroup([x, y]))
    fixed.add_function(x * x + y * y - 1)
    start = pb.Slice.from_coefficients([[0, 1, 0]], [x, y])   # y = 0
    end = pb.Slice.from_coefficients([[-1, 1, 0]], [x, y])    # y - x = 0
    # affine on purpose: this checks that a Slice is accepted where a System is expected, and the
    # builder's projectivize default (b2#382) would change the row count the assertions pin
    b = na.straight_line_homotopy(end, start, fixed=fixed, projectivize=False,
                                  gamma=pb.coefficient(pb.multiprec.complex_mp('0.6', '0.8')))
    assert b.homotopy.num_functions() == 2
    solver = na.HomotopySolver(b, [[1, 0], [-1, 0]], mptype='double')   # target comes with it
    solver.solve()
    r = round(1 / np.sqrt(2), 4)
    roots = sorted((round(complex(s[0]).real, 4), round(complex(s[1]).real, 4)) for s in solver.all_solutions())
    assert roots == sorted([(r, r), (-r, -r)])


# --- #392: config names compose with get_config; final_tolerance has a deterministic owner -------

def test_get_config_accepts_the_names_config_names_lists():
    _, _, sys = _circle_and_line()
    zd = na.ZeroDimSolver(sys, mptype='double')
    for name in zd.config_names():
        cfg = zd.get_config(name)
        assert type(cfg) is type(zd.get_config(type(cfg)))
    assert zd.get_config('tolerances').final_tolerance == zd.get_config(TolerancesConfig).final_tolerance


def test_final_tolerance_set_on_the_endgame_survives_solve_and_solver_setting_wins():
    _, _, sys = _circle_and_line()
    zd = na.ZeroDimSolver(sys, mptype='double')
    eg = zd.get_endgame()
    cfg = eg.get_endgame_settings()
    cfg.update(final_tolerance='1e-8')
    eg.set_endgame_settings(cfg)
    zd.solve()
    assert float(zd.get_endgame().get_endgame_settings().final_tolerance) == pytest.approx(1e-8)   # was reverted

    zd.set(final_tolerance='1e-9')                       # the solver-level setting, set later, wins
    zd.solve()
    assert float(zd.get_endgame().get_endgame_settings().final_tolerance) == pytest.approx(1e-9)
