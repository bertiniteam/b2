# This file is part of Bertini 2.
#
# tracking_analytic.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tracking_analytic.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with tracking_analytic.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""The tutorial's runnable source, and the regenerator for its figure.

Tracks a homotopy whose functions are not polynomials -- sine -- and draws the paths in the
complex plane.  Two frames, side by side: five paths landing on distinct roots of
``sin(x) = 1/2``, and two paths colliding at a branch point of ``sin(x) = 1``.

Nothing here is random.  gamma is the exact rational point ``(-24 + 7i)/25`` on the unit
circle, and the start points are exact, so the picture is the same every time.

Run:  python tracking_analytic.py
      (or via tools/refresh_doc_artifacts.py --plots --only tracking_analytic)
"""

import os

import matplotlib
matplotlib.use('Agg')            # headless: no display needed
import matplotlib.pyplot as plt

import numpy as np

import bertini
from bertini import SolutionPathCollector

_OUT = os.path.dirname(os.path.abspath(__file__))

PI_50 = '3.14159265358979323846264338327950288419716939937511'


def sine_homotopy(right_hand_side):
    """H = (1-t)*(sin(x) - c) + gamma*t*sin(x), and the target system it ends at.

    At t = 1 the homotopy is gamma*sin(x), whose roots are the integer multiples of pi.  At
    t = 0 it is the target, sin(x) = c.
    """
    x = bertini.Variable('x')

    target = bertini.System()
    target.add_variable_group(bertini.VariableGroup([x]))
    target.add_function(bertini.symbolics.sin(x) - bertini.coefficient(right_hand_side))

    start = bertini.System()
    start.add_variable_group(bertini.VariableGroup([x]))
    start.add_function(bertini.symbolics.sin(x))

    gamma = bertini.coefficient('-24/25') + bertini.I * bertini.coefficient('7/25')
    return bertini.system.make_homotopy(target, start, gamma=gamma), target


def zeros_of_sine(multiples):
    """Start points: the roots of sin(x) = 0 you choose to chase, as k*pi."""
    points = []
    for k in multiples:
        value = bertini.multiprec.real_mp(PI_50) * k
        points.append(np.array([bertini.multiprec.complex_mp(value)]))
    return points


def track(right_hand_side, multiples, mptype='multiple'):
    """Track the chosen start points, collecting every path.  Returns (solver, collector)."""
    homotopy, target = sine_homotopy(right_hand_side)
    solver = bertini.HomotopySolver(homotopy, zeros_of_sine(multiples), target, mptype=mptype)

    collector = SolutionPathCollector()
    solver.add_observer(collector)
    solver.solve()
    return solver, collector


def _affine(points):
    """The affine x-coordinate of a recorded path.  These systems are never homogenized, so the
    recorded points are already affine; the check keeps the helper honest if that changes."""
    return points[:, 0] if points.shape[1] == 1 else points[:, 1] / points[:, 0]


def _draw(ax, collector, solver, title):
    cmap = plt.get_cmap('viridis')
    series = collector.series
    for index, path in enumerate(series):
        x_of_t = _affine(path.points())
        color = cmap(index / max(len(series) - 1, 1))
        ax.plot(x_of_t.real, x_of_t.imag, '-', color=color, lw=1.4)
        ax.plot(x_of_t.real[0], x_of_t.imag[0], 'o', color=color, ms=7, mfc='white', zorder=4,
                label='start: a zero of sine' if index == 0 else None)

    ends = [complex(s[0]) for s in solver.all_solutions() if len(s)]
    ax.scatter([e.real for e in ends], [e.imag for e in ends],
               c='k', marker='*', s=150, zorder=5, label='endpoint')

    ax.axhline(0, color='0.85', lw=0.8, zorder=0)     # the real axis, where both ends live
    ax.set_title(title, fontsize=11)
    ax.set_xlabel('Re(x)')
    ax.set_ylabel('Im(x)')
    ax.legend(loc='upper left', fontsize=8, framealpha=0.9)


def plot():
    """Draw both frames and write the figure beside this file."""
    spread_solver, spread_paths = track('1/2', (-2, -1, 0, 1, 2))
    collide_solver, collide_paths = track(1, (0, 1))

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    _draw(axes[0], spread_paths, spread_solver,
          'sin(x) = 1/2:  five paths, five distinct roots')
    _draw(axes[1], collide_paths, collide_solver,
          'sin(x) = 1:  two paths, one double root at pi/2')

    for ax in axes:
        ax.set_box_aspect(1)

    fig.tight_layout()
    fig.savefig(os.path.join(_OUT, 'tracking_analytic.svg'))
    fig.savefig(os.path.join(_OUT, 'tracking_analytic.png'), dpi=150)
    print('wrote tracking_analytic.svg and .png to', _OUT)


def report():
    """Print what the tutorial asserts, for anyone running this file directly."""
    solver, _ = track('1/2', (-2, -1, 0, 1, 2))
    ends = sorted(complex(s[0]).real for s in solver.all_solutions() if len(s))
    print('sin(x) = 1/2 endpoints:', [round(e, 12) for e in ends])
    print('arcsin(1/2)           :', round(float(np.arcsin(0.5)), 12))

    solver, _ = track(1, (0, 1))
    print('sin(x) = 1 cycle numbers:', [m.cycle_num for m in solver.solution_metadata()])
    print('pi/2                    :', round(float(np.pi / 2), 12))


if __name__ == '__main__':
    bertini.recording(False)          # a doc figure is not a run worth archiving
    bertini.default_precision(60)
    report()
    plot()
