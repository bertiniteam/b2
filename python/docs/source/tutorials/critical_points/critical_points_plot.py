"""Regenerate the critical_points tutorial's 3-D curve figure (critical_points.svg/.png).

This is the standalone regenerator for the plot that the tutorial otherwise draws inline in a
``.. testcode::`` block (which only ``plt.show()``s).  It solves the same rank-deficiency system as
``critical_points.py``, then draws the full reducible curve -- the y-axis line, the two interlocking
circles, and the quartic where the cylinders meet -- with every real critical point labelled as a
smooth projection-critical point (dot) or a singular crossing (X).  The solve is seeded, so the
figure is reproducible.

The sibling ``critical_points_fibers.svg`` figure has no regenerator: it is a hand-authored
illustration with no in-repo source (see tools/refresh_doc_artifacts.py `_NO_REGENERATOR`).

Run:  python critical_points_plot.py   (or via tools/refresh_doc_artifacts.py --plots --only critical_points)
"""

import os

import matplotlib
matplotlib.use('Agg')            # headless: no display needed
import matplotlib.pyplot as plt

import numpy as np

import bertini
from bertini import linalg
from bertini.nag_algorithm import ZeroDimSolver

_OUT = os.path.dirname(os.path.abspath(__file__))


def solve_critical_points():
    """Solve the rank-deficiency system; return the finite (x, y, z) complex tuples."""
    bertini.random.set_random_seed(165)
    x, y, z = bertini.Variable('x'), bertini.Variable('y'), bertini.Variable('z')
    f = x * (x**2 + y**2 - 1)
    g = z * ((y - 1)**2 + z**2 - 1)

    pi = bertini.random_matrix(1, 3, real=True, orthonormal=False)   # a random real projection row

    J = bertini.jacobian([f, g], [x, y, z])               # 2 x 3 symbolic Jacobian of the curve
    M = np.vstack([J, linalg.as_coefficients(pi)])        # 3 x 3: J_f stacked over the projection
    v = linalg.variable_vector('v', 3)                    # the null-vector unknowns

    sys = bertini.System()
    sys.add(bertini.VariableGroup([x, y, z, *v]), f, g)   # curve equations
    linalg.add_functions(sys, M @ v)                      # M v = 0   (rank deficiency)
    sys.add_function((bertini.random_matrix(1, 3, symbolic=True) @ v)[0] - 1)   # de-zero patch h.v = 1

    solver = ZeroDimSolver(sys, endgame='cauchy', mptype='adaptive', startsystem='binomial')
    solver.solve()

    finite = []
    for s in solver.solutions():
        p = np.array(s)
        xyz = (complex(p[0]), complex(p[1]), complex(p[2]))
        if all(abs(w) < 1e6 for w in xyz):
            finite.append(xyz)
    return finite


def plot(finite):
    """Draw the full curve with its real critical points; return the figure."""
    reals = [np.array([w.real for w in p]) for p in finite if all(abs(w.imag) < 1e-7 for w in p)]
    uniq = []
    for p in reals:
        if not any(np.linalg.norm(p - q) < 1e-4 for q in uniq):
            uniq.append(p)

    def n_components(p):                                   # how many components pass through p
        X, Y, Z = p; e = 1e-5
        return sum([abs(X) < e and abs(Z) < e,                                   # line  x=z=0
                    abs(X) < e and abs((Y - 1)**2 + Z**2 - 1) < e,               # circle  x=0
                    abs(Z) < e and abs(X**2 + Y**2 - 1) < e,                     # circle  z=0
                    abs(X**2 + Y**2 - 1) < e and abs((Y - 1)**2 + Z**2 - 1) < e])# quartic

    smooth   = np.array([p for p in uniq if n_components(p) == 1])
    singular = np.array([p for p in uniq if n_components(p) >= 2])

    t = np.linspace(0, 2*np.pi, 400); tt = np.linspace(0, np.pi, 300)
    fig = plt.figure(figsize=(7.5, 6.5)); ax = fig.add_subplot(projection='3d')
    ax.plot([0, 0], [-2.5, 2.5], [0, 0], color='0.5', lw=1.2, label='line  x=z=0')
    ax.plot(np.zeros_like(t), 1 + np.cos(t), np.sin(t), 'C0', lw=1.2, label='circle  x=0')
    ax.plot(np.cos(t), np.sin(t), np.zeros_like(t), 'C1', lw=1.2, label='circle  z=0')
    for sgn in (1, -1):                                    # quartic: two z-branches
        ax.plot(np.cos(tt), np.sin(tt), sgn*np.sqrt(np.clip(np.sin(tt)*(2 - np.sin(tt)), 0, None)),
                'C2', lw=1.0, label='quartic (cylinder ∩ cylinder)' if sgn == 1 else None)
    ax.scatter(smooth[:, 0], smooth[:, 1], smooth[:, 2], c='C3', s=48, depthshade=False,
               label='projection-critical (smooth)')
    ax.scatter(singular[:, 0], singular[:, 1], singular[:, 2], marker='X', c='k', s=55,
               depthshade=False, label='singular crossings')
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.set_xlim(-2, 2); ax.set_ylim(-2.5, 2.5); ax.set_zlim(-1.6, 1.6)
    ax.legend(loc='upper left', fontsize=7.5)
    ax.set_title('All critical points: smooth (●) and singular crossings (✕)')
    return fig


def main():
    fig = plot(solve_critical_points())
    fig.savefig(os.path.join(_OUT, 'critical_points.svg'))
    fig.savefig(os.path.join(_OUT, 'critical_points.png'), dpi=150)


if __name__ == '__main__':
    main()
