"""Structured blocks survive a zero-dim solve.

The payoff of making the evaluation blocks homogenization-aware: a linear condition added as
a LinearFormsBlock (via add_linear / add_linear_forms) is no longer confined to
direct evaluation -- it is homogenized alongside the polynomial functions when the zero-dim
solver prepares the target, so a mixed poly + linear-forms system solves end to end.
"""

import numpy as np
import pytest

import bertini as pb


def _roots_xy(solutions):
    return sorted(
        (round(complex(s[0]).real, 4), round(complex(s[1]).real, 4)) for s in solutions
    )


def test_circle_intersect_line_via_add_linear_solves():
    # circle x^2 + y^2 - 1 (a polynomial-block function) intersect the line 2x + y - 1 = 0
    # carried as a LinearFormsBlock.  The two intersection points are (0, 1) and (0.8, -0.6).
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add(pb.VariableGroup([x, y]))
    sys.add(x * x + y * y - 1)                                   # polynomial block, degree 2
    sys.add_linear(np.array([[2, 1]]), np.array([x, y]), [-1])  # linear-forms block: 2x+y-1

    assert list(sys.degrees()) == [2, 1]

    zd = pb.ZeroDimSolver(sys, endgame='cauchy', mptype='adaptive', startsystem='binomial')
    zd.solve()
    got = _roots_xy(zd.all_solutions())

    assert len(got) == 2
    expected = sorted([(0.0, 1.0), (0.8, -0.6)])
    for (gx, gy), (ex, ey) in zip(got, expected):
        assert abs(gx - ex) < 1e-6 and abs(gy - ey) < 1e-6


def test_circle_intersect_line_via_add_linear_forms_solves():
    # same intersection, but the line is supplied as an augmented coefficient row through
    # add_linear_forms ([2, 1, -1] = 2x + 1y + (-1)).
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add(pb.VariableGroup([x, y]))
    sys.add(x * x + y * y - 1)
    sys.add_linear_forms([[2, 1, -1]])

    zd = pb.ZeroDimSolver(sys, endgame='cauchy', mptype='adaptive', startsystem='binomial')
    zd.solve()
    got = _roots_xy(zd.all_solutions())

    assert len(got) == 2
    expected = sorted([(0.0, 1.0), (0.8, -0.6)])
    for (gx, gy), (ex, ey) in zip(got, expected):
        assert abs(gx - ex) < 1e-6 and abs(gy - ey) < 1e-6
