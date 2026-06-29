"""ZeroDimSolver feasibility + square-up behaviors through the Python interface.

The correctness gate is the C++ suite (test/nag_algorithms/zero_dim.cpp); this pins the
Python-facing contract: an over-determined system is squared up and its extraneous solutions
filtered out, and an under-determined system raises a helpful error.
"""

import pytest

import bertini as pb
from bertini.nag_algorithm import ZeroDimSolver


def _xy():
    return pb.Variable('x'), pb.Variable('y')


def test_square_system_is_not_randomized():
    x, y = _xy()
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_function(x * x + y * y - 1)
    sys.add_function(x - y)
    solver = ZeroDimSolver(sys, mptype='adaptive')
    assert solver.was_randomized() is False
    solver.solve()
    assert len(solver.finite_solutions()) == 2


def test_overdetermined_system_is_squared_and_filtered():
    x, y = _xy()
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_function(x * x + y * y - 2)   # three equations, two variables
    sys.add_function(x - y)
    sys.add_function(x * y - 1)
    # genuine common solutions: (1, 1) and (-1, -1)

    solver = ZeroDimSolver(sys, mptype='adaptive')
    assert solver.was_randomized() is True
    solver.solve()

    # the randomized (square) system has MORE endpoints than genuine solutions ...
    n_all = len(solver.all_solutions())
    assert n_all > 2
    # ... but solutions() / finite_solutions returns only the two that satisfy the ORIGINAL system.
    assert len(solver.solutions()) == 2
    assert len(solver.finite_solutions()) == 2
    # how many of the squaring's extra roots stay finite (nonsolutions) vs. diverge is RNG-dependent,
    # so assert the ROBUST partition rather than an exact nonsolution count: every endpoint is a genuine
    # finite solution, a finite nonsolution, or at infinity.
    assert (len(solver.solutions()) + len(solver.nonsolutions())
            + len(solver.infinite_solutions())) == n_all
    # opting the nonsolutions back in returns everything except the at-infinity endpoints
    assert len(solver.solutions(nonsolution=True)) == len(solver.solutions()) + len(solver.nonsolutions())
    # nonsolution endpoints are flagged is_nonsolution but remain geometrically finite
    md = solver.solution_metadata()
    assert all(m.is_finite for m in md if m.is_nonsolution)
    assert sum(1 for m in md if m.is_nonsolution) == len(solver.nonsolutions())


def test_solutions_filter_by_realness():
    # circle meeting the line x = y: two REAL solutions
    x, y = _xy()
    real_sys = pb.System()
    real_sys.add_variable_group(pb.VariableGroup([x, y]))
    real_sys.add_function(x * x + y * y - 1)
    real_sys.add_function(x - y)
    rl = ZeroDimSolver(real_sys, mptype='adaptive')
    rl.solve()
    assert len(rl.solutions()) == 2                       # both, by default
    assert len(rl.solutions(nonreal=False)) == 2          # real only
    assert len(rl.solutions(real=False)) == 0             # complex only

    # circle meeting the hyperbola xy = 1: four COMPLEX (non-real) solutions
    x, y = _xy()
    cx_sys = pb.System()
    cx_sys.add_variable_group(pb.VariableGroup([x, y]))
    cx_sys.add_function(x * x + y * y - 1)
    cx_sys.add_function(x * y - 1)
    cx = ZeroDimSolver(cx_sys, mptype='adaptive')
    cx.solve()
    assert len(cx.solutions()) == 4                       # all four
    assert len(cx.solutions(nonreal=False)) == 0          # none are real
    assert len(cx.solutions(real=False)) == 4             # all complex


def test_underdetermined_system_raises_helpful_error():
    x, y = _xy()
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_function(x * x + y * y - 1)   # one equation, two variables

    with pytest.raises(RuntimeError, match='under-determined'):
        ZeroDimSolver(sys, mptype='adaptive')
