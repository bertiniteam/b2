"""Randomize an overdetermined system -- Bertini 2 tutorial.

Squares up overdetermined systems by generic randomization, then filters.
Run:  python randomize.py
"""

import numpy as np
import bertini
from bertini import linalg


def part1_one_block():
    """Part 1 -- one block of variables: three curves in the plane."""
    x, y = bertini.Variable('x'), bertini.Variable('y')

    original = bertini.System()
    original.add_variable_group(bertini.VariableGroup([x, y]))
    original.add_function(x*x + y*y - 1)      # the unit circle      (degree 2)
    original.add_function(x - y)              # the line x = y       (degree 1)
    original.add_function(2*x*x - 1)          # the lines x = ±1/√2  (degree 2)

    assert list(original.degrees()) == [2, 1, 2]   # three equations, two unknowns: overdetermined

    # randomize returns a new square system; the original is untouched
    randomized = linalg.randomize(original)

    assert randomized.num_functions() == 2          # squared up
    assert original.num_functions() == 3            # original unchanged
    assert sorted(randomized.degrees()) == [2, 2]

    # solve the square system, then filter against the original
    bertini.random.set_random_seed(1)               # reproducible generic coefficients + gamma
    zd = bertini.nag_algorithm.ZeroDimSolver(randomized, endgame='cauchy', mptype='adaptive', startsystem='binomial')
    zd.solve()
    solutions = zd.all_solutions()
    assert len(solutions) == 4                       # the two we want, plus two extraneous

    def satisfies_original(point):
        coords = np.array([complex(c) for c in point])
        return max(abs(complex(v)) for v in original.eval(coords)) < 1e-7

    true_solutions = [np.array([complex(c) for c in s]) for s in solutions if satisfies_original(s)]

    r = 1.0 / np.sqrt(2.0)
    found = sorted((round(p[0].real, 4), round(p[1].real, 4)) for p in true_solutions)
    assert found == sorted([(round(r, 4), round(r, 4)), (round(-r, 4), round(-r, 4))])


def part2_variable_groups():
    """Part 2 -- variables in separate groups: a bilinear system."""
    x, y = bertini.Variable('x'), bertini.Variable('y')

    original = bertini.System()
    original.add_variable_group(bertini.VariableGroup([x]))    # x is its own group
    original.add_variable_group(bertini.VariableGroup([y]))    # y is its own group
    original.add_function(x*y - 1)        # xy = 1
    original.add_function(x + y - 2)      # x + y = 2
    original.add_function(x - y)          # x = y

    randomized = linalg.randomize(original)
    assert randomized.num_functions() == 2

    # solve with the multihomogeneous start system, which exploits the grouping
    bertini.random.set_random_seed(3)
    zd = bertini.nag_algorithm.ZeroDimSolver(randomized, endgame='cauchy', mptype='adaptive', startsystem='mhom')
    zd.solve()

    def satisfies(point):
        coords = np.array([complex(c) for c in point])
        return max(abs(complex(v)) for v in original.eval(coords)) < 1e-7

    true_solutions = [np.array([complex(c) for c in s]) for s in zd.all_solutions() if satisfies(s)]

    assert len(true_solutions) == 1
    p = true_solutions[0]
    assert abs(p[0] - 1) < 1e-7 and abs(p[1] - 1) < 1e-7


def main():
    part1_one_block()
    part2_variable_groups()


if __name__ == '__main__':
    main()
