"""Finding all the solutions -- Bertini 2 tutorial.

Computes every isolated complex solution of a system with no real solutions,
and shows the family of solution accessors.

Run:  python all_solutions.py
"""

import numpy as np
import bertini as bertini
from bertini.nag_algorithm import ZeroDimSolver


def build_system():
    """A system with no real solutions: the unit circle meeting the hyperbola xy = 1."""
    x, y = bertini.Variable('x'), bertini.Variable('y')

    sys = bertini.System()
    sys.add_function(x**2 + y**2 - 1)     # the unit circle
    sys.add_function(x*y - 1)             # the hyperbola xy = 1
    sys.add_variable_group(bertini.VariableGroup([x, y]))

    return sys


def solve_and_collect(sys):
    """Solve in adaptive precision and collect the finite solutions."""
    solver = ZeroDimSolver(sys, mptype='adaptive')
    solver.solve()

    good = solver.finite_solutions()      # successful, finite endpoints -- the actual points

    assert len(good) == 4                 # all four complex solutions

    return solver, good


def check_all_complex(good):
    """Every one of them is genuinely complex, and each really is a solution."""
    for s in good:
        xv, yv = complex(s[0]), complex(s[1])
        assert abs(xv.imag) > 1e-6 or abs(yv.imag) > 1e-6   # none are real
        # and each really is a solution
        assert abs(xv*xv + yv*yv - 1) < 1e-8
        assert abs(xv*yv - 1) < 1e-8


def accessor_views(solver):
    """all_solutions() is the whole list; the categories are filtered views of it."""
    assert len(solver.all_solutions()) == 4          # one per tracked path
    assert len(solver.finite_solutions()) == 4        # the finite ones
    assert len(solver.infinite_solutions()) == 0      # the at-infinity ones (is_finite is False)


def main():
    sys = build_system()
    solver, good = solve_and_collect(sys)
    check_all_complex(good)
    accessor_views(solver)


if __name__ == '__main__':
    main()
