"""Hauenstein distance critical points: a real point on every component (the Trott curve).

Pick a random real point p; the critical points of the squared distance to p, restricted to the
curve f=0, contain at least one real point on each connected component.  The Trott curve has four
ovals, and a generic p yields its nearest + farthest point on each -- eight real points.
"""

from fractions import Fraction

import bertini as pb
from bertini import linalg


def _trott_critical_system(px, py):
    x, y = pb.Variable('x'), pb.Variable('y')
    f = 144 * (x**4 + y**4) - 225 * (x**2 + y**2) + 350 * x**2 * y**2 + 81
    fx, fy = f.differentiate(x), f.differentiate(y)               # bertini differentiates the curve
    # gradient of f parallel to (x - p): the 2x2 determinant vanishes
    parallel = (x - linalg.coefficient(px)) * fy - (y - linalg.coefficient(py)) * fx
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_function(f)
    sys.add_function(parallel)
    return sys


def _trott(x, y):
    return 144 * (x**4 + y**4) - 225 * (x**2 + y**2) + 350 * x**2 * y**2 + 81


def test_real_point_on_every_trott_component():
    sys = _trott_critical_system(Fraction(41, 100), Fraction(23, 100))
    assert list(sys.degrees()) == [4, 4]            # 16 paths

    solver = pb.nag_algorithm.ZeroDimSolver(sys, endgame='cauchy', mptype='adaptive', startsystem='binomial')
    solver.solve()

    # keep the real solutions by the solver's own classification (is_real applies the configured
    # real_threshold) rather than a hand-picked epsilon
    sols, meta = solver.all_solutions(), solver.solution_metadata()
    reals = [(complex(s[0]).real, complex(s[1]).real)
             for s, m in zip(sols, meta) if m.is_real]

    # two critical points (nearest + farthest) on each of the four ovals
    assert len(reals) == 8

    # every real critical point actually lies on the curve
    for (a, b) in reals:
        assert abs(_trott(a, b)) < 1e-6

    # all four ovals are represented: one per (sign(x), sign(y)) quadrant region
    quadrants = {(a > 0, b > 0) for (a, b) in reals}
    assert len(quadrants) == 4
