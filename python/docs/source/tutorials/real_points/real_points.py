"""A real point on every component of the Trott curve.

Assembles the code fragments from the ``real_points`` tutorial into one
runnable program.  We build the critical-point (squared-distance) system for
the Trott quartic, solve it, keep the real solutions -- one real point on every
component -- and plot the curve, the random point p, and the distances.

Run:  python real_points.py
"""

import os
import matplotlib
matplotlib.use('Agg')            # headless backend, must precede pyplot import
import matplotlib.pyplot as plt

import numpy as np
from fractions import Fraction
import bertini

_OUT = os.path.dirname(os.path.abspath(__file__))


def build_system():
    """Build the square critical-point system for the Trott curve."""
    x, y = bertini.Variable('x'), bertini.Variable('y')
    f = 144*(x**4 + y**4) - 225*(x**2 + y**2) + 350*x**2*y**2 + 81

    # bertini differentiates the curve symbolically -- no hand-coded gradients
    fx, fy = f.differentiate(x), f.differentiate(y)

    # a hardcoded "random" real point p, with exact rational coordinates
    PX, PY = Fraction(41, 100), Fraction(23, 100)
    parallel = (x - bertini.coefficient(PX))*fy - (y - bertini.coefficient(PY))*fx

    crit_system = bertini.System()
    crit_system.add_variable_group(bertini.VariableGroup([x, y]))
    crit_system.add_function(f)
    crit_system.add_function(parallel)

    return crit_system, PX, PY


def solve_system(crit_system):
    """Solve and keep the real critical points -- the solver decides 'real'."""
    solver = bertini.ZeroDimSolver(crit_system, mptype='adaptive')
    solver.solve()

    crit = [(complex(s[0]).real, complex(s[1]).real)
            for s in solver.real_solutions()]      # the solver decides what "real" means

    assert len(crit) == 8                  # two per oval (nearest + farthest), all four ovals hit

    return crit


def plot(crit, PX, PY):
    """Plot the curve, the point p, the critical points, and the distances."""
    px, py = float(PX), float(PY)
    gx, gy = np.meshgrid(np.linspace(-1.2, 1.2, 500), np.linspace(-1.2, 1.2, 500))
    F = 144*(gx**4 + gy**4) - 225*(gx**2 + gy**2) + 350*gx**2*gy**2 + 81

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.contour(gx, gy, F, levels=[0], colors='C0')          # the Trott curve
    for (cx, cy) in crit:
        ax.plot([px, cx], [py, cy], color='0.7', lw=0.8)    # distance to each critical point
    ax.scatter([c[0] for c in crit], [c[1] for c in crit], c='C3', s=40, label='real critical points')
    ax.scatter([px], [py], c='k', marker='*', s=200, label='random point p')
    ax.set_aspect('equal'); ax.legend(loc='upper right', fontsize=8)
    ax.set_title('A real point on every component of the Trott curve')

    fig.savefig(os.path.join(_OUT, 'real_points.png'), dpi=150)


def main():
    crit_system, PX, PY = build_system()
    crit = solve_system(crit_system)
    plot(crit, PX, PY)


if __name__ == '__main__':
    main()
