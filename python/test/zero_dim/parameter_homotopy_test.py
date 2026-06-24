"""Parameter homotopy: solve once, re-solve many — via nag_algorithm.user_homotopy.

The user-homotopy entry runs the full zero-dim pipeline (pre-endgame tracking, midpath check,
endgame, post-processing) on a homotopy you constructed, starting from a list of start points you
already have.  The flagship use is the parameter homotopy: solve a system once at a generic
parameter, then reuse those solutions to move the parameter wherever you like -- cheaply.
"""

import numpy as np

import bertini as pb


def _roots_real(solutions):
    return sorted(round(complex(s[0]).real, 4) for s in solutions)


def test_solve_once_then_sweep_a_parameter():
    # Solve x^2 - 4 = 0 once (generic parameter p = 4): the two roots are +/- 2.
    generic = pb.System()
    xg = pb.Variable('x')
    generic.add_variable_group(pb.VariableGroup([xg]))
    generic.add_function(xg * xg - 4)
    zd0 = pb.nag_algorithm.ZeroDimCauchyAdaptivePrecisionTotalDegree(generic)
    zd0.solve()
    start_points = zd0.all_solutions()
    assert _roots_real(start_points) == [-2.0, 2.0]

    # Now reuse those two solutions as start points to move the parameter p: 4 -> p_target,
    # tracking the homotopy H(x,t) = x^2 - ((1-t) p_target + t * 4).  At t=1 it is x^2 - 4 (the
    # start points); at t=0 it is x^2 - p_target.  No new ab-initio solve per parameter.
    for p_target, expected in [(9, [-3.0, 3.0]), (16, [-4.0, 4.0]), (25, [-5.0, 5.0])]:
        x, t = pb.Variable('x'), pb.Variable('t')
        H = pb.System()
        H.add_variable_group(pb.VariableGroup([x]))
        H.add_function(x * x - ((1 - t) * p_target + t * 4))
        H.add_path_variable(t)

        target = pb.System()
        target.add_variable_group(pb.VariableGroup([x]))
        target.add_function(x * x - p_target)

        solver = pb.nag_algorithm.user_homotopy(H, start_points, target)
        solver.solve()
        assert _roots_real(solver.all_solutions()) == expected


def test_user_homotopy_rejects_bad_precision():
    x, t = pb.Variable('x'), pb.Variable('t')
    H = pb.System(); H.add_variable_group(pb.VariableGroup([x])); H.add_function(x * x - (4 - 3 * t)); H.add_path_variable(t)
    target = pb.System(); target.add_variable_group(pb.VariableGroup([x])); target.add_function(x * x - 1)
    import pytest
    with pytest.raises(ValueError):
        pb.nag_algorithm.user_homotopy(H, [], target, precision='quadruple')


def test_user_homotopy_rejects_solver_as_start_points():
    # Regression for issue #258: passing the start-point *solver* (not its .all_solutions())
    # used to fail with a cryptic "object is not iterable"; now it explains the mistake.
    x, t = pb.Variable('x'), pb.Variable('t')
    H = pb.System(); H.add_variable_group(pb.VariableGroup([x])); H.add_function(x * x - (4 - 3 * t)); H.add_path_variable(t)
    target = pb.System(); target.add_variable_group(pb.VariableGroup([x])); target.add_function(x * x - 1)

    generic = pb.System(); generic.add_variable_group(pb.VariableGroup([x])); generic.add_function(x * x - 4)
    solver = pb.nag_algorithm.ZeroDim(generic, mptype='adaptive')   # a solver, not solutions

    import pytest
    with pytest.raises(TypeError, match="start_points must be"):
        pb.nag_algorithm.user_homotopy(H, solver, target)


def test_coefficient_parameter_homotopy_helper():
    # the coefficient_parameter_homotopy helper builds (1-t)*target + t*generic for you.
    x = pb.Variable('x')
    generic = pb.System(); generic.add_variable_group(pb.VariableGroup([x])); generic.add_function(x * x - 4)
    target = pb.System(); target.add_variable_group(pb.VariableGroup([x])); target.add_function(x * x - 9)

    gen_solver = pb.nag_algorithm.ZeroDimCauchyAdaptivePrecisionTotalDegree(generic)
    gen_solver.solve()

    H = pb.nag_algorithm.coefficient_parameter_homotopy(target, generic)
    solver = pb.nag_algorithm.user_homotopy(H, gen_solver.all_solutions(), target)
    solver.solve()
    assert _roots_real(solver.all_solutions()) == [-3.0, 3.0]
