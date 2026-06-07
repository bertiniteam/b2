import pytest

import bertini as pb
from bertini.nag_algorithm import (
    ZeroDimCauchyAdaptivePrecisionTotalDegree,
    ZeroDimPowerSeriesAdaptivePrecisionTotalDegree,
    TolerancesConfig,
)


@pytest.fixture
def circle_intersection_solver():
    """x^2 + y^2 - 1 = 0 and x + y = 0: two solutions (±1/√2, ∓1/√2)."""
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_function(x**2 + y**2 - 1)
    sys.add_function(x + y)
    sys.add_variable_group(pb.VariableGroup([x, y]))
    return ZeroDimCauchyAdaptivePrecisionTotalDegree(sys)


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
    solns = circle_intersection_solver.solutions()
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
    solns = solver.solutions()
    assert len(solns) == 2


def test_power_series_endgame_variant():
    """ZeroDimPowerSeriesAdaptivePrecisionTotalDegree solves a simple system.

    Mirrors zero_dim/can_run_griewank_osborn (PSEG variant).
    """
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_function(x**2 + y**2 - 1)
    sys.add_function(x + y)
    sys.add_variable_group(pb.VariableGroup([x, y]))

    solver = ZeroDimPowerSeriesAdaptivePrecisionTotalDegree(sys)
    solver.solve()
    solns = solver.solutions()
    assert len(solns) == 2
