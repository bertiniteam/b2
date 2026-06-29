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
    assert len(solver.all_solutions()) > 2
    # ... but finite_solutions returns only the two that satisfy the ORIGINAL system.
    finite = solver.finite_solutions()
    assert len(finite) == 2


def test_underdetermined_system_raises_helpful_error():
    x, y = _xy()
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_function(x * x + y * y - 1)   # one equation, two variables

    with pytest.raises(RuntimeError, match='under-determined'):
        ZeroDimSolver(sys, mptype='adaptive')
