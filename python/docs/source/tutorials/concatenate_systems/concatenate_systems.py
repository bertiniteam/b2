"""Build a system in pieces with clone and concatenate (Bertini 2 tutorial).

Stacks two four-variable, two-function blocks into a square 4x4 system and solves it.
Run:  python concatenate_systems.py
"""

import numpy as np
import bertini


def build_pieces():
    """Build two non-square blocks (geometry + constraints) over shared variables."""
    x1, x2, x3, x4 = (bertini.Variable(n) for n in ('x1', 'x2', 'x3', 'x4'))
    group = bertini.VariableGroup([x1, x2, x3, x4])

    geometry = bertini.System()
    geometry.add_variable_group(group)
    geometry.add_function(x1*x1 - x2)        # x2 = x1²
    geometry.add_function(x2*x2 - x3)        # x3 = x2²

    constraints = bertini.System()
    constraints.add_variable_group(group)
    constraints.add_function(x3 - x4)        # x4 = x3
    constraints.add_function(x4 - x1)        # x1 = x4

    assert geometry.num_functions() == 2 and geometry.num_variables() == 4
    assert constraints.num_functions() == 2 and constraints.num_variables() == 4

    return (x1, x2, x3, x4), geometry, constraints


def stack(geometry, constraints):
    """Concatenate the two blocks into a new square system; inputs unchanged."""
    system = bertini.system.concatenate(geometry, constraints)

    assert system.num_functions() == 4 and system.num_variables() == 4
    assert geometry.num_functions() == 2        # inputs unchanged
    assert constraints.num_functions() == 2
    assert list(system.degrees()) == [2, 2, 1, 1]

    return system


def solve_and_check(system):
    """Solve the square system and verify the four solutions."""
    bertini.random.set_random_seed(1)           # reproducible start system + gamma
    zd = bertini.nag_algorithm.ZeroDimSolver(system, endgame='cauchy', mptype='adaptive', startsystem='binomial')
    zd.solve()

    solutions = [np.array([complex(c) for c in s]) for s in zd.all_solutions()]
    assert len(solutions) == 4

    # every computed point satisfies the concatenated system
    assert all(max(abs(complex(v)) for v in system.eval(s)) < 1e-7 for s in solutions)

    # the two real solutions are the zero point and the all-ones point
    real = [s for s in solutions if max(abs(c.imag) for c in s) < 1e-7]
    found = sorted(tuple(round(c.real, 4) for c in s) for s in real)
    assert found == [(0.0, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0, 1.0)]

    return solutions


def build_with_clone(variables):
    """Declare variables once, clone per block, concatenate, and solve."""
    x1, x2, x3, x4 = variables

    base = bertini.System()                     # variable setup only -- no functions yet
    base.add_variable_group(bertini.VariableGroup([x1, x2, x3, x4]))

    geometry = bertini.system.clone(base)
    geometry.add_function(x1*x1 - x2)
    geometry.add_function(x2*x2 - x3)

    constraints = bertini.system.clone(base)
    constraints.add_function(x3 - x4)
    constraints.add_function(x4 - x1)

    system = bertini.system.concatenate(geometry, constraints)
    assert system.num_functions() == 4

    bertini.random.set_random_seed(1)
    zd = bertini.nag_algorithm.ZeroDimSolver(system, endgame='cauchy', mptype='adaptive', startsystem='binomial')
    zd.solve()
    assert len(zd.all_solutions()) == 4


def main():
    variables, geometry, constraints = build_pieces()
    system = stack(geometry, constraints)
    solve_and_check(system)
    build_with_clone(variables)


if __name__ == '__main__':
    main()
