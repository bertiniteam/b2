# This file is part of Bertini 2.
#
# the_seam.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# the_seam.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the_seam.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""The Seam -- a line that is not in the mathematics, drawn by the paths that died on it.

For each parameter c on a grid, track

    log(x) - y = 0,      y^2 - x + c = 0

from the solution (x, y) = (1, 0) at c = 1, and keep every step of every path.  The second
variable IS the logarithm, so the tracker is continuing a logarithm rather than evaluating
one.  Some paths arrive.  Some crawl onto the negative real axis and stall there with the
step size at its floor, because that is where the numeric logarithm jumps -- a convention
inside an implementation, not a feature of the equations.

The hero frame is the x-plane, every step of every path exposed onto one buffer: arrivals in
teal, deaths in red, and the deaths pile onto the cut.  The teaching frame is the c-plane,
showing which parameters are reachable at all and how far the rest got.

Nothing here is random.  gamma is the exact rational point (-24 + 7i)/25 on the unit circle
and every grid coordinate is an exact rational, so both frames reproduce exactly.

Run:  python the_seam.py            (the full render, some minutes)
      python the_seam.py --scout    (a small one, for checking a window)
      or via tools/refresh_doc_artifacts.py --plots --only the_seam
"""

import argparse
import os
from fractions import Fraction

import matplotlib
matplotlib.use('Agg')            # headless: no display needed
import matplotlib.pyplot as plt

import numpy as np

_OUT = os.path.dirname(os.path.abspath(__file__))

# The parameter window, as exact rationals: the reachable lobe, its coastline, and the wedge
# bitten out of the lower left.
RE_RANGE = (Fraction(-14), Fraction(8))
IM_RANGE = (Fraction(-9), Fraction(9))

# The x-plane window for the exposure.  Both interesting points sit in it: the branch point of
# the logarithm at 0, and the common birthplace of every path at 1.
X_RANGE, Y_RANGE = (-2.6, 2.6), (-2.6, 2.6)

PRECISION = 40
_worker = {}


# --- the engine -----------------------------------------------------------------------------

def _init(precision):
    """Per-process setup: import bertini once, and build the constants it will reuse."""
    import bertini
    bertini.recording(False)          # tens of thousands of throwaway solves: do not archive
    bertini.default_precision(precision)
    _worker['bertini'] = bertini
    _worker['gamma'] = (bertini.coefficient('-24/25')
                        + bertini.I * bertini.coefficient('7/25'))


def _exact(re, im):
    """An exact complex coefficient node from two rationals."""
    bertini = _worker['bertini']
    return bertini.coefficient(str(re)) + bertini.I * bertini.coefficient(str(im))


def _system_for(re, im):
    bertini = _worker['bertini']
    x, y = bertini.Variable('x'), bertini.Variable('y')
    s = bertini.System()
    s.add_variable_group(bertini.VariableGroup([x, y]))
    s.add_function(bertini.symbolics.log(x) - y)      # y is the logarithm, continued not evaluated
    s.add_function(y * y - x + _exact(re, im))
    return s


def _track(re, im):
    """One parameter value.  Returns (trajectory in x, times, fate)."""
    bertini = _worker['bertini']
    target = _system_for(re, im)
    start = _system_for(Fraction(1), Fraction(0))
    H = bertini.system.make_homotopy(target, start, gamma=_worker['gamma'])

    point = np.array([bertini.multiprec.complex_mp('1', '0'),
                      bertini.multiprec.complex_mp('0', '0')])
    solver = bertini.HomotopySolver(H, [point], target, mptype='multiple')

    collector = bertini.SolutionPathCollector()
    solver.add_observer(collector)
    solver.solve()

    meta = list(solver.solution_metadata())
    fate = str(meta[0].pre_endgame_success_code).split('.')[-1] if meta else 'NeverStarted'

    if not collector.series:
        return np.zeros(0, dtype=complex), np.zeros(0, dtype=complex), fate

    path = collector.series[0]
    xs = np.array([complex(v) for v in path.points()[:, 0]])
    ts = np.array([complex(t) for t in path.times()])
    return xs, ts, fate


def _row(args):
    """One row of the parameter grid, for the pool."""
    row_index, im, reals = args
    return row_index, [(_track(re, im)) for re in reals]


def _grid(width, height):
    """Exact rational grid coordinates, so the picture does not depend on float formatting."""
    reals = [RE_RANGE[0] + (RE_RANGE[1] - RE_RANGE[0]) * Fraction(i, width - 1)
             for i in range(width)]
    imags = [IM_RANGE[1] - (IM_RANGE[1] - IM_RANGE[0]) * Fraction(j, height - 1)
             for j in range(height)]
    return reals, imags


def sweep(width, height, processes=None):
    """Track the whole grid, in parallel.  Returns rows of (xs, ts, fate)."""
    from multiprocessing import Pool

    reals, imags = _grid(width, height)
    jobs = [(j, imags[j], reals) for j in range(height)]
    rows = [None] * height

    with Pool(processes or os.cpu_count(), initializer=_init, initargs=(PRECISION,)) as pool:
        for done, (row_index, record) in enumerate(pool.imap_unordered(_row, jobs), start=1):
            rows[row_index] = record
            if done % max(height // 12, 1) == 0 or done == height:
                print(f'  rows {done}/{height}', flush=True)
    return rows


# --- the exposure ---------------------------------------------------------------------------

def _deposit(buffer, xs, weights, side, samples_per_unit):
    """Lay one trajectory onto the buffer, sampling along each segment so the trace is a
    continuous stroke rather than a row of dots, and weighting by the path time it spans --
    brightness is dwell, not step count, so refining the stepper does not change the picture."""
    if len(xs) < 2:
        return
    for index, (a, b) in enumerate(zip(xs[:-1], xs[1:])):
        steps = max(int(abs(b - a) * samples_per_unit), 2)
        walk = np.linspace(a, b, steps)
        col = (walk.real - X_RANGE[0]) / (X_RANGE[1] - X_RANGE[0]) * (side - 1)
        row = (Y_RANGE[1] - walk.imag) / (Y_RANGE[1] - Y_RANGE[0]) * (side - 1)
        keep = (col >= 0) & (col < side) & (row >= 0) & (row < side)
        if keep.any():
            np.add.at(buffer, (row[keep].astype(int), col[keep].astype(int)),
                      weights[index] / steps)


def expose(rows, side):
    """Three buffers: paths that arrived, paths that died, and where the dying ones stopped."""
    samples_per_unit = side / (X_RANGE[1] - X_RANGE[0]) * 1.6
    arrived = np.zeros((side, side))
    dying = np.zeros((side, side))
    graves = np.zeros((side, side))

    for record in rows:
        for xs, ts, fate in record:
            if len(xs) < 2:
                continue
            weights = np.abs(np.diff(np.abs(ts)))
            if fate == 'Success':
                _deposit(arrived, xs, weights, side, samples_per_unit)
            else:
                _deposit(dying, xs, weights, side, samples_per_unit)
                where = xs[-1]
                col = (where.real - X_RANGE[0]) / (X_RANGE[1] - X_RANGE[0]) * (side - 1)
                row = (Y_RANGE[1] - where.imag) / (Y_RANGE[1] - Y_RANGE[0]) * (side - 1)
                if 0 <= col < side and 0 <= row < side:
                    graves[int(row), int(col)] += 1.0
    return arrived, dying, graves


def _stretch(buffer, gain):
    """Compress the dynamic range so a faint single path and a thousand overlapping ones are
    both visible.  Normalized by the buffer's own maximum, so it is scale-free."""
    top = buffer.max()
    return buffer if top <= 0 else np.log1p(buffer * gain) / np.log1p(top * gain)


def colour(arrived, dying, graves, gain):
    live, doomed = _stretch(arrived, gain), _stretch(dying, gain)
    rest = _stretch(graves, 60)
    rgb = np.zeros(arrived.shape + (3,))
    rgb[..., 0] = doomed * 1.15 + rest * 1.40                      # deaths burn red
    rgb[..., 1] = live * 0.95 + doomed * 0.20 + rest * 0.55        # arrivals in teal
    rgb[..., 2] = live * 1.05 + doomed * 0.45 + rest * 0.30
    return np.clip(rgb, 0, 1)


# --- the frames -----------------------------------------------------------------------------

def hero(rows, side, gain, name):
    """Write the exposure.

    Saved through a 256-entry adaptive palette with Floyd-Steinberg dithering, which is
    deterministic and costs nothing visible on an image of two hues over black, but takes the
    file from about 2.2 MB to under 1 MB.  Committed figures stay small; the resolution is
    what carries the detail here, so the resolution is what is kept.
    """
    from PIL import Image

    rgb = colour(*expose(rows, side), gain=gain)
    picture = Image.fromarray((np.clip(rgb, 0, 1) * 255).astype(np.uint8), mode='RGB')
    picture = picture.quantize(colors=256, method=Image.MEDIANCUT, dither=Image.FLOYDSTEINBERG)
    picture.save(os.path.join(_OUT, name), optimize=True)
    print('  wrote', name)


def teaching(rows, width, height, name):
    """The c-plane: which parameters are reachable, and how far the rest got."""
    reached = np.zeros((height, width))
    for j, record in enumerate(rows):
        for i, (_, ts, fate) in enumerate(record):
            reached[j, i] = 1.0 if fate == 'Success' else (1.0 - abs(ts[-1]) if len(ts) else 0.0)

    fig, ax = plt.subplots(figsize=(7, 6))
    image = ax.imshow(reached, origin='upper', cmap='magma', vmin=0, vmax=1,
                      extent=[float(RE_RANGE[0]), float(RE_RANGE[1]),
                              float(IM_RANGE[0]), float(IM_RANGE[1])])
    ax.set_xlabel('Re(c)')
    ax.set_ylabel('Im(c)')
    ax.set_title('which parameters are reachable, and how far the rest got', fontsize=10)
    bar = fig.colorbar(image, ax=ax)
    bar.set_label('path time reached (1.0 = arrived)')
    fig.tight_layout()
    fig.savefig(os.path.join(_OUT, name), dpi=130)
    plt.close(fig)
    print('  wrote', name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scout', action='store_true',
                        help='a small, fast render, for checking a window')
    args = parser.parse_args()

    width, height, side, gain = (90, 74, 900, 4000) if args.scout else (240, 200, 1600, 26000)

    print(f'tracking {width * height} paths on a {width} x {height} parameter grid')
    rows = sweep(width, height)

    fates = {}
    for record in rows:
        for _, _, fate in record:
            fates[fate] = fates.get(fate, 0) + 1
    print('  fates:', fates)

    suffix = '_scout' if args.scout else ''
    hero(rows, side, gain, f'the_seam{suffix}.png')
    teaching(rows, width, height, f'the_seam_teaching{suffix}.png')


if __name__ == '__main__':
    main()
