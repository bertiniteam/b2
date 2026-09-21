"""Plotting multiprecision solutions and paths.

Solves the circle meeting the parabola, then draws two figures:

    * plotting_solutions -- the real plane, and the complex x-plane of all four solutions
    * plotting_paths     -- the four homotopy paths through the complex x-plane

Run:  python plotting.py
"""

import os

import numpy as np

import matplotlib
matplotlib.use('Agg')                 # headless: no display needed
import matplotlib.pyplot as plt

import bertini
from bertini import ZeroDimSolver, SolutionPathCollector

# with ambient recording on, a repeated identical solve recalls its answers and tracks
# nothing, so the path collector would see no paths at all
bertini.recording(False)

_OUT = os.path.dirname(os.path.abspath(__file__))


def build():
    """The unit circle meeting the parabola y = x^2: two real solutions and two complex."""
    x, y = bertini.Variable('x'), bertini.Variable('y')
    system = bertini.System()
    system.add_variable_group(bertini.VariableGroup([x, y]))
    system.add_function(x * x + y * y - 1)
    system.add_function(y - x * x)
    return system


def solve_and_collect(system):
    """Solve, keeping every path.  Returns the solver and the collector."""
    bertini.random.set_random_seed(2)     # so you get exactly this picture

    solver = ZeroDimSolver(system, mptype='adaptive')
    paths = SolutionPathCollector()
    solver.add_observer(paths)
    solver.solve()
    return solver, paths


def plot_solutions(solver):
    """The real plane beside the complex x-plane; saves plotting_solutions.{svg,png}."""
    solutions = solver.all_solutions()

    # every solution is a numpy array of the multiprecision complex dtype, so the parts
    # come from bertini.real / bertini.imag -- numpy's .real/.imag would be wrong here
    xs = np.array([solution[0] for solution in solutions])
    ys = np.array([solution[1] for solution in solutions])

    x_re, x_im = bertini.real(xs), bertini.imag(xs)
    y_re = bertini.real(ys)

    real_enough = [abs(float(v)) < 1e-12 for v in bertini.imag(xs)]

    figure, (left, right) = plt.subplots(1, 2, figsize=(11, 5))

    angle = np.linspace(0, 2 * np.pi, 400)
    left.plot(np.cos(angle), np.sin(angle), color='#888888', lw=1, label='$x^2 + y^2 = 1$')
    grid = np.linspace(-1.3, 1.3, 400)
    left.plot(grid, grid ** 2, color='#4477aa', lw=1, label='$y = x^2$')
    left.scatter([v for v, keep in zip(x_re, real_enough) if keep],
                 [v for v, keep in zip(y_re, real_enough) if keep],
                 c='k', marker='o', s=70, zorder=5, label='real solutions')
    left.set_xlim(-1.4, 1.4)
    left.set_ylim(-1.2, 1.7)
    left.set_aspect(1.0)
    left.set_xlabel('x')
    left.set_ylabel('y')
    left.set_title('the real plane')
    left.legend(loc='upper left', fontsize=8)

    right.axhline(0, color='#cccccc', lw=0.8)
    right.axvline(0, color='#cccccc', lw=0.8)
    right.scatter([v for v, keep in zip(x_re, real_enough) if keep],
                  [v for v, keep in zip(x_im, real_enough) if keep],
                  c='k', marker='o', s=70, label='real')
    right.scatter([v for v, keep in zip(x_re, real_enough) if not keep],
                  [v for v, keep in zip(x_im, real_enough) if not keep],
                  c='#cc3311', marker='s', s=70, label='complex')
    right.set_aspect(1.0)
    right.set_xlabel('Re(x)')
    right.set_ylabel('Im(x)')
    right.set_title('the complex x-plane')
    right.legend(loc='upper left', fontsize=8)

    figure.tight_layout()
    figure.savefig(os.path.join(_OUT, 'plotting_solutions.svg'))
    figure.savefig(os.path.join(_OUT, 'plotting_solutions.png'), dpi=150)
    plt.close(figure)


def plot_paths(solver, paths):
    """The four homotopy paths in the complex x-plane; saves plotting_paths.{svg,png}."""
    figure, axis = plt.subplots(figsize=(6, 6))
    colours = plt.get_cmap('turbo')

    for index, path in enumerate(paths.series):
        points = path.points()                    # float64 already, one row per step
        x_affine = points[:, 1] / points[:, 0]    # dehomogenize the first coordinate
        colour = colours(index / max(len(paths.series) - 1, 1))
        axis.plot(x_affine.real, x_affine.imag, '-', color=colour, lw=1.3)
        axis.plot(x_affine.real[0], x_affine.imag[0], 'o',
                  color=colour, ms=6, mfc='white')

    xs = np.array([solution[0] for solution in solver.all_solutions()])
    axis.scatter(bertini.real(xs), bertini.imag(xs),
                 c='k', marker='*', s=160, zorder=5, label='solutions')

    axis.set_aspect(1.0)
    axis.set_box_aspect(1)
    axis.set_xlabel('Re(x)')
    axis.set_ylabel('Im(x)')
    axis.set_title('paths through the complex x-plane (open circles: start points)')
    axis.legend(loc='upper right', fontsize=8)

    figure.tight_layout()
    figure.savefig(os.path.join(_OUT, 'plotting_paths.svg'))
    figure.savefig(os.path.join(_OUT, 'plotting_paths.png'), dpi=150)
    plt.close(figure)


def main():
    system = build()
    solver, paths = solve_and_collect(system)
    plot_solutions(solver)
    plot_paths(solver, paths)


if __name__ == '__main__':
    main()
