"""A database of solutions: filtering and plotting with pandas.

Assembled from the ``.. testcode::`` blocks of the solution_database tutorial.
Solves the intersection of two plane curves (a nodal cubic and a second cubic
threaded through its node), lays the solve out as a pandas DataFrame, filters the
solutions by category, and draws two figures: the real plane and the complex planes.

Run:  python solution_database.py
"""

import os

import matplotlib
matplotlib.use('Agg')            # headless: no display needed
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

import bertini
from bertini.nag_algorithm import ZeroDimSolver

_OUT = os.path.dirname(os.path.abspath(__file__))


def two_plane_curves():
    """Build and solve the intersection of the two cubics."""
    x, y = bertini.Variable('x'), bertini.Variable('y')
    sys = bertini.System()
    sys.add_function(y*y - x**3 - x**2)     # nodal cubic   y^2 = x^3 + x^2  (node at the origin)
    sys.add_function(y - x**3 + 2*x)        # second cubic  y = x^3 - 2x     (passes through it)
    sys.add_variable_group(bertini.VariableGroup([x, y]))

    solver = ZeroDimSolver(sys, mptype='adaptive')
    solver.solve()
    return solver


def one_row_per_solution(solver):
    """Lay the solve out as a DataFrame, one row per distinct finite solution."""
    df = solver.to_dataframe()
    assert len(df) == 5                         # five distinct finite solutions
    assert {'solution', 'is_real', 'is_singular', 'multiplicity'} <= set(df.columns)
    assert len(df['solution'].iloc[0]) == 2     # each cell is the whole (x, y) point

    assert len(solver.to_dataframe(merge_multiplicities=False)) == 6     # the node counted twice

    assert df['system'].iloc[0] is df['system'].iloc[-1]   # one shared reference, not a per-row copy

    # each category is a filter on the frame
    real_simple    = df[df.is_real & ~df.is_singular]      # transverse real crossings
    complex_simple = df[~df.is_real & ~df.is_singular]      # a complex-conjugate pair
    singular       = df[df.is_singular]                     # the node, one row, multiplicity 2

    assert len(real_simple) == 2
    assert len(complex_simple) == 2
    assert len(singular) == 1 and (singular.multiplicity == 2).all()

    return df


def split_coordinates(df):
    """Pull the two coordinates out of the solution cell into their own columns."""
    df['x'] = [complex(p[0]) for p in df.solution]
    df['y'] = [complex(p[1]) for p in df.solution]
    df['x_re'], df['x_im'] = df.x.map(lambda z: z.real), df.x.map(lambda z: z.imag)
    df['y_re'], df['y_im'] = df.y.map(lambda z: z.real), df.y.map(lambda z: z.imag)

    real_simple, singular = df[df.is_real & ~df.is_singular], df[df.is_singular]
    return real_simple, singular


def plot_real_plane(real_simple, singular):
    """The two curves and the real roots, in the real (x, y) plane."""
    gx, gy = np.meshgrid(np.linspace(-2.4, 2.4, 500), np.linspace(-3.6, 3.6, 500))
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.contour(gx, gy, gy**2 - gx**3 - gx**2, levels=[0], colors='C0')   # y^2 = x^3 + x^2
    ax.contour(gx, gy, gy - gx**3 + 2*gx,     levels=[0], colors='C1')   # y = x^3 - 2x
    ax.scatter(real_simple.x_re, real_simple.y_re, c='k', marker='o', s=70, label='real, simple')
    ax.scatter(singular.x_re,    singular.y_re,    c='r', marker='*', s=300, label='real, singular')
    ax.legend(); ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_aspect('equal')
    fig.savefig(os.path.join(_OUT, 'solution_database_real_plane.svg'))
    fig.savefig(os.path.join(_OUT, 'solution_database_real_plane.png'), dpi=150)


def plot_complex_planes(df, real_simple, singular):
    """Each coordinate in its own complex plane, colored by category."""
    complex_simple = df[~df.is_real & ~df.is_singular]
    cats = [('real, simple',    real_simple,    'k',  'o',  70),
            ('complex, simple', complex_simple, 'C2', 's',  70),
            ('singular',        singular,       'r',  '*', 300)]

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, (re, im, title) in zip(axes, [('x_re', 'x_im', 'x in C'), ('y_re', 'y_im', 'y in C')]):
        ax.axhline(0, color='0.8'); ax.axvline(0, color='0.8')
        for label, sub, color, marker, size in cats:
            ax.scatter(sub[re], sub[im], c=color, marker=marker, s=size, label=label)
        ax.set_title(title); ax.set_xlabel('real part'); ax.set_ylabel('imaginary part')
    axes[0].legend()
    fig.savefig(os.path.join(_OUT, 'solution_database_complex_planes.svg'))
    fig.savefig(os.path.join(_OUT, 'solution_database_complex_planes.png'), dpi=150)


def main():
    solver = two_plane_curves()
    df = one_row_per_solution(solver)
    real_simple, singular = split_coordinates(df)
    plot_real_plane(real_simple, singular)
    plot_complex_planes(df, real_simple, singular)


if __name__ == '__main__':
    main()
