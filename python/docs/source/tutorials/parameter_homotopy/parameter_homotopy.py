"""Parameter homotopy: solve once, re-solve many times (Bertini 2 tutorial).

Fixed unit circle meeting a moving horizontal line y = s/2.
Run:  python parameter_homotopy.py
"""

import bertini
from bertini import nag_algorithm

# the variables are SHARED across every member of the family: the parameter homotopy
# interpolates the members' equations, so they must be built over the same Variable objects.
x, y = bertini.Variable('x'), bertini.Variable('y')


def member(s):
    """The system { x^2 + y^2 - 1, 2y - s }: the unit circle meeting the line y = s/2."""
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([x, y]))
    sys.add_function(x*x + y*y - 1)
    sys.add_function(2*y - s)
    return sys


def solve_once():
    """Pick a generic member and solve it ab initio; its solutions are our start points."""
    generic = member(1)                                # the line y = 1/2
    first = bertini.ZeroDimSolver(generic, mptype='adaptive')
    first.solve()
    start_points = first.all_solutions()               # (+/- sqrt(3)/2, 1/2)
    return generic, start_points


def move_parameter(generic, start_points):
    """Move to another member via a parameter homotopy, without solving from scratch."""
    target = member(0)                                 # the line y = 0
    H = nag_algorithm.coefficient_parameter_homotopy(target, generic)
    moved = bertini.HomotopySolver(H, start_points, target)
    moved.solve()
    # moved.all_solutions() are now (+/- 1, 0)
    return moved


def sweep(generic, start_points):
    """Sweep many parameters, reusing the same start points, never solving from scratch."""
    for s in [0, -1, 1]:                               # lines y = 0, -1/2, 1/2
        target = member(s)
        H = nag_algorithm.coefficient_parameter_homotopy(target, generic)
        solver = bertini.HomotopySolver(H, start_points, target)
        solver.solve()
        roots = [p for p in solver.all_solutions() if len(p) == 2]
        for p in roots:
            xv, yv = complex(p[0]), complex(p[1])
            assert abs(xv*xv + yv*yv - 1) < 1e-8       # on the circle
            assert abs(2*yv - s) < 1e-8                # on the line y = s/2


def main():
    generic, start_points = solve_once()
    move_parameter(generic, start_points)
    sweep(generic, start_points)


if __name__ == '__main__':
    main()
