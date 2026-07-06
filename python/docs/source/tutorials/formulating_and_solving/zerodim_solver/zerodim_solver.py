"""Solving a system with ZeroDimSolver -- Bertini 2 tutorial.

Demonstrates ZeroDimSolver on a square, an over-determined, and an
under-determined polynomial system.

Run:  python zerodim_solver.py
"""

import bertini
from bertini import ZeroDimSolver


def square_system():
    """A pleasant square system: the unit circle meeting the line x = y."""
    x, y = bertini.Variable('x'), bertini.Variable('y')

    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([x, y]))
    sys.add_function(x**2 + y**2 - 1)     # the unit circle
    sys.add_function(x - y)               # the line x = y

    solver = ZeroDimSolver(sys, mptype='adaptive')
    solver.solve()

    assert solver.was_randomized() is False    # already square -- nothing to do
    assert len(solver.finite_solutions()) == 2

    return solver


def over_determined_system():
    """An over-determined system: extra equations, extraneous solutions filtered out."""
    x, y = bertini.Variable('x'), bertini.Variable('y')
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([x, y]))
    sys.add_function(x**2 + y**2 - 2)
    sys.add_function(x - y)
    sys.add_function(x*y - 1)

    solver = ZeroDimSolver(sys, mptype='adaptive')
    assert solver.was_randomized() is True     # ZeroDimSolver squared it up for you
    solver.solve()

    n_all = len(solver.all_solutions())        # the squared system's full path count
    assert n_all > 2
    assert len(solver.solutions()) == 2        # only the two true common solutions
    # the remaining endpoints are nonsolutions (finite, not solutions) or diverged to infinity; every endpoint
    # is exactly one of genuine solution / nonsolution / at-infinity:
    assert len(solver.solutions()) + len(solver.nonsolutions()) + len(solver.infinite_solutions()) == n_all

    return solver


def under_determined_system():
    """An under-determined system: a helpful refusal."""
    x, y = bertini.Variable('x'), bertini.Variable('y')
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([x, y]))
    sys.add_function(x**2 + y**2 - 1)          # one equation, two variables: a whole curve

    try:
        ZeroDimSolver(sys, mptype='adaptive')
        raise AssertionError("expected ZeroDimSolver to refuse an under-determined system")
    except RuntimeError as e:
        assert 'under-determined' in str(e)


def main():
    square_system()
    over_determined_system()
    under_determined_system()


if __name__ == '__main__':
    main()
