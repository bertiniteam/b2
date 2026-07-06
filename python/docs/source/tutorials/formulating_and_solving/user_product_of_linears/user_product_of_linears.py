"""Author your own start system: a product of linears (Bertini 2 tutorial).

Assembles the tutorial's code fragments into one runnable program.
Run:  python user_product_of_linears.py
"""

import itertools
import math

import numpy as np

import bertini
from bertini import nag_algorithm
from bertini import multiprec


def build_target():
    """A target you can check by hand: unit circle meeting a parabola."""
    x, y = bertini.Variable('x'), bertini.Variable('y')

    target = bertini.System()
    target.add_variable_group(bertini.VariableGroup([x, y]))
    target.add_function(x*x + y*y - 1)        # the unit circle
    target.add_function(y - x*x)              # the parabola y = x^2
    return target, x, y


def build_start(x, y):
    """A start system you can write down: each function a product of two linears."""
    start = bertini.System()
    start.add_variable_group(bertini.VariableGroup([x, y]))
    start.add_products_of_linears([
        [[1, 0, '-1'], [1, 0, '1']],     # s0 = (x - 1)(x + 1)
        [[0, 1, '-1'], [0, 1, '-2']],    # s1 = (y - 1)(y - 2)
    ])

    assert list(start.degrees()) == [2, 2]   # a product's degree is its number of factors
    return start


def build_start_points():
    """Start points are intersections of hyperplanes: the grid x in {1,-1} times y in {1,2}."""
    start_points = [np.array([multiprec.complex_mp(str(a)), multiprec.complex_mp(str(b))])
                    for a, b in itertools.product([1, -1], [1, 2])]
    # (1, 1), (1, 2), (-1, 1), (-1, 2)
    return start_points


def solve_homotopy(target, start, start_points):
    """Blend start into a homotopy with the gamma-trick and solve."""
    gamma = bertini.coefficient(multiprec.complex_mp('0.6', '0.8'))   # exact, off the real axis
    H = nag_algorithm.blend_homotopy(target, start, gamma=gamma)

    solver = bertini.HomotopySolver(H, start_points, target)
    solver.solve()
    solutions = solver.all_solutions()
    return solver, solutions


def check_against_known(solutions):
    """Check computed roots against the four exact answers."""
    assert len(solutions) == 4

    s5 = math.sqrt(5)
    y1, y2 = (-1 + s5) / 2, (-1 - s5) / 2
    known = [np.array([math.sqrt(y1), y1]),  np.array([-math.sqrt(y1), y1]),
             np.array([1j*math.sqrt(-y2), y2]), np.array([-1j*math.sqrt(-y2), y2])]

    computed = [np.array([complex(v) for v in p]) for p in solutions]
    for root in known:
        nearest = min(np.max(np.abs(c - root)) for c in computed)
        assert nearest < 1e-8


def classify_endpoints(solver):
    """Ask the solver's metadata which endpoints are finite, nonsingular, real."""
    assert len(solver.finite_solutions()) == 4       # all four endpoints are finite
    assert len(solver.nonsingular_solutions()) == 4  # ... and all nonsingular
    assert len(solver.real_solutions()) == 2         # two real, two with purely imaginary x


def main():
    target, x, y = build_target()
    start = build_start(x, y)
    start_points = build_start_points()
    solver, solutions = solve_homotopy(target, start, start_points)
    check_against_known(solutions)
    classify_endpoints(solver)


if __name__ == '__main__':
    main()
