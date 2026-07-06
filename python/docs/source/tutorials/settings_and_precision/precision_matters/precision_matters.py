"""Tutorial: When double precision is not enough (adaptive precision on cyclic-5).

Assembled from the ``.. testcode::`` blocks in index.rst.
Run:  python precision_matters.py
"""

import numpy as np
import bertini
from bertini.nag_algorithm import ZeroDimSolver


def build_cyclic_system(n=5):
    """The cyclic-n roots system, built as a total-degree-solvable System."""
    x = [bertini.Variable('x{}'.format(i)) for i in range(n)]
    w = x + x                                        # doubled, to take products that wrap around
    sys = bertini.System()
    for length in range(1, n):                       # the degree-`length` cyclic sums
        sys.add_function(np.sum([np.prod(w[s:s + length]) for s in range(n)]))
    sys.add_function(np.prod(x) - 1)                 # the product, normalized
    sys.add_variable_group(bertini.VariableGroup(x))
    return sys


def solve_adaptive(sys):
    """Solve reliably with adaptive precision, and verify against the known ground truth."""
    solver = ZeroDimSolver(sys, mptype='adaptive')
    solver.solve()
    report = solver.report()

    assert report.num_finite_solutions == 70         # every finite solution of cyclic-5
    assert report.num_failed == 0                    # no path the tracker had to give up on
    return report


def solve_double(sys):
    """The same solve in double precision -- faster, but can quietly lose a near-singular path."""
    solver_d = ZeroDimSolver(sys, mptype='double')
    solver_d.solve()
    report_d = solver_d.report()

    assert report_d.num_finite_solutions <= 70       # never too many -- but sometimes too few
    return report_d


def main():
    sys = build_cyclic_system(n=5)
    solve_adaptive(sys)
    solve_double(sys)


if __name__ == '__main__':
    main()
