import pytest

import bertini as pb
from bertini import ZeroDimSolver
from bertini.nag_algorithm import TolerancesConfig


@pytest.fixture
def circle_intersection_solver():
    """x^2 + y^2 - 1 = 0 and x + y = 0: two solutions (±1/√2, ∓1/√2)."""
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_function(x**2 + y**2 - 1)
    sys.add_function(x + y)
    sys.add_variable_group(pb.VariableGroup([x, y]))
    return ZeroDimSolver(sys, endgame='cauchy', mptype='adaptive', startsystem='binomial')


@pytest.fixture
def solver(circle_intersection_solver):
    """Alias kept for backward-compatibility with the original test."""
    return circle_intersection_solver


def test_can_solve_multiple_times(solver):
    """Multiple solver.solve() calls must not crash.

    Regression: calling solve() twice triggered an MPFR precision assertion
    in patch rescaling.  Mirrors zero_dim C++ can_run_griewank_osborn (basic
    run test).
    """
    solver.solve()
    solver.solve()
    solver.solve()


def test_solution_count(circle_intersection_solver):
    """The circle-intersection system has exactly 2 solutions.

    Mirrors zero_dim/can_run_griewank_osborn — checks the result, not just
    that the solver runs.
    """
    circle_intersection_solver.solve()
    solns = circle_intersection_solver.all_solutions()
    assert len(solns) == 2


def test_custom_tolerances(circle_intersection_solver):
    """Changing tolerances before solving should still yield correct solutions.

    Mirrors zero_dim/can_run_change_some_settings.
    """
    solver = circle_intersection_solver
    tols = solver.get_config(TolerancesConfig)
    tols.newton_before_endgame = 1e-4
    tols.newton_during_endgame = 1e-4
    solver.set_config(tols)

    solver.solve()
    solns = solver.all_solutions()
    assert len(solns) == 2


def test_power_series_endgame_variant():
    """ZeroDimSolver solves a simple system.

    Mirrors zero_dim/can_run_griewank_osborn (PowerSeriesEndgame variant).
    """
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_function(x**2 + y**2 - 1)
    sys.add_function(x + y)
    sys.add_variable_group(pb.VariableGroup([x, y]))

    solver = ZeroDimSolver(sys, endgame='powerseries', mptype='adaptive', startsystem='binomial')
    solver.solve()
    solns = solver.all_solutions()
    assert len(solns) == 2


# --- coordinate representations of solutions ---
# solutions() returns USER coordinates by default (the variables you wrote);
# solutions(user_coords=False) is the explicit opt-out giving the solver's
# internal representation: homogenized, lying on the target system's patch.
# see ADR-0013.

import numpy as np
from bertini import multiprec as mp

INV_SQRT2 = 1 / np.sqrt(2)
KNOWN_SOLUTIONS = (
    np.array([complex(INV_SQRT2), complex(-INV_SQRT2)]),
    np.array([complex(-INV_SQRT2), complex(INV_SQRT2)]),
)


def _as_complex(pt):
    return np.array([complex(pt[i]) for i in range(len(pt))])


def _distance_to_known(pt):
    return min(np.linalg.norm(_as_complex(pt) - k) for k in KNOWN_SOLUTIONS)


@pytest.fixture
def solved():
    """the circle/line system, its variables, and a solved solver."""
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_function(x**2 + y**2 - 1)
    sys.add_function(x + y)
    sys.add_variable_group(pb.VariableGroup([x, y]))
    solver = ZeroDimSolver(sys, endgame='cauchy', mptype='adaptive', startsystem='binomial')
    solver.solve()
    return sys, solver


def test_solutions_are_in_user_coordinates_by_default(solved):
    _, solver = solved
    sols = solver.all_solutions()
    assert len(sols) == 2
    for s in sols:
        assert len(s) == 2  # the user's variables, not [h, x, y]
        assert _distance_to_known(s) < 1e-8


def test_solutions_internal_coords_are_explicit_optout(solved):
    _, solver = solved
    internal = solver.all_solutions(user_coords=False)
    assert len(internal) == 2
    ts = solver.target_system()
    user = solver.all_solutions()
    for i in range(2):
        assert len(internal[i]) == 3  # [hom_var, x, y]
        dehomed = ts.dehomogenize_point(internal[i])
        # compare in mpfr (a float64 round-trip puts an ~1e-16 ulp floor on
        # the difference; x86_64 CI caught this).  the AMP solver may emit
        # precision-16 values, so the identity holds at the value's own
        # working precision -- hence 1e-12, not the ambient 30 digits.
        for j in range(2):
            assert mp.abs(dehomed[j] - user[i][j]) < mp.Float('1e-12')


def test_homogenize_point_reenters_internal_coordinates(solved):
    """the lift: user coords -> homogenized, on the target system's patch."""
    _, solver = solved
    ts = solver.target_system()
    user = solver.all_solutions()
    internal = solver.all_solutions(user_coords=False)
    for i in range(2):
        lifted = ts.homogenize_point(user[i])
        assert len(lifted) == 3
        # projectively the same point, on the same patch -> numerically equal
        assert np.linalg.norm(_as_complex(lifted) - _as_complex(internal[i])) < 1e-8
        # already on the patch: rescaling is the identity -- compare in mpfr
        # at the value's own working precision (see note in the previous test)
        rescaled = ts.rescale_point_to_fit_patch(lifted)
        for j in range(3):
            assert mp.abs(rescaled[j] - lifted[j]) < mp.Float('1e-12')


def test_variable_orderings_label_the_representations(solved):
    sys, solver = solved
    user_names = [v.name for v in sys.variable_ordering()]
    assert user_names == ['x', 'y']

    internal_names = [v.name for v in solver.target_system().variable_ordering()]
    assert len(internal_names) == 3
    assert internal_names[1:] == ['x', 'y']  # leading entry is the homogenizing variable
