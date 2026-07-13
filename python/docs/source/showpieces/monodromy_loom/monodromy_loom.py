"""The Monodromy Loom.

A showpiece: watch the solutions of a *parameterized* polynomial braid as its parameter walks
a closed loop around the discriminant.  For the one-parameter family

    f(x; c) = x^d - d*x - c ,

the d roots move continuously as c does; carry c once around a loop that encircles the branch
points (the c-values where two roots collide) and the roots come back **permuted** -- the
monodromy of the family, i.e. its Galois action, made visible as a literal braid.

Each strand is one root's path, plotted in 3-D as ``(Re x, Im x, theta)`` where theta is the
loop angle 0 -> 2*pi.  The whole loop is baked into a single homotopy by making the coefficient
a function of the path variable, ``c(t) = center + radius * exp(i*(theta + phi))`` with
``theta = 2*pi*(1 - t)``, so one continuous track per strand traces the entire loop and the
path variable maps straight to the loop angle.

The *texture* of each strand is the adaptive-precision tracker's own diagnostics, streamed by a
``PathDataCollector``: where the loop grazes a branch point two strands nearly collide, the step
size collapses, and the strand brightens and swells right there -- you can see the solver work
hardest exactly where the mathematics is.

Two frames are produced:

    * monodromy_loom_teaching.png -- x^3 - 3x - c, a loop around ONE branch point: two roots
      swap (a single transposition), the third rides straight.  The lesson in one picture.
    * monodromy_loom.png          -- x^5 - 5x - c, a loop around ALL FOUR branch points: a full
      5-cycle, with four pinches where the tracker sweats.  The show-off.

Regenerated through ``tools/refresh_doc_artifacts.py`` (see that tool).  This is a raster
showpiece: PNG only (a 3-D render has no meaningful SVG), so it is exempt from the tutorial
figures' png+svg rule.  It is NOT a doctest -- the docs page embeds the pre-rendered image.

Run standalone:  python monodromy_loom.py
"""

import math
import os
from fractions import Fraction

import numpy as np

import matplotlib
matplotlib.use('Agg')                 # headless: no display needed
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection

import bertini
import bertini.tracking as tracking
from bertini.multiprec import complex_mp
from bertini.tracking.observers import amp as amp_observers

_OUT = os.path.dirname(os.path.abspath(__file__))

# a vivid strand palette on near-black
_HUES = ['#3fe0ff', '#ff4fa3', '#ffd24a', '#8aff5a', '#c07bff', '#ff8a3f']
_BG = '#05060a'


# --- the engine: bake the whole parameter loop into one homotopy --------------------------------
#
# Every constant here is EXACT and nothing reaches for a float: the loop parameters are rationals
# (fractions.Fraction) and the phase is a rational multiple of pi.  Node arithmetic promotes ints,
# exposes the exact constants bertini.Pi and bertini.I, honours **, and bertini.coefficient accepts
# exact values -- so the homotopy is written verbatim, with no coercion to defeat.

def _c_of_theta(theta, center, radius):
    """The loop's moving coefficient c = center + radius * (cos theta + i sin theta), as a node.
    ``center`` / ``radius`` are exact (int or fractions.Fraction); ``theta`` is a node."""
    return (bertini.coefficient(center)
            + bertini.coefficient(radius) * (bertini.cos(theta) + bertini.I * bertini.sin(theta)))

def loom_homotopy(degree, center, radius, phi_over_pi):
    """H(x, t) = x^d - d*x - c(t),  c(t) = center + radius*exp(i*(theta + phi)),
    theta = 2pi(1-t),  phi = phi_over_pi * pi.

    At t=1 theta=phi (the start configuration); at t=0 the loop has closed.  Tracking a root of the
    start configuration from t=1 to t=0 carries it once around the loop.
    """
    x = bertini.Variable('x')
    t = bertini.Variable('t')
    theta = 2 * bertini.Pi * (1 - t) + bertini.coefficient(phi_over_pi) * bertini.Pi
    sys = bertini.System()
    sys.add_function(x**degree - degree * x - _c_of_theta(theta, center, radius))
    sys.add_path_variable(t)
    sys.add_variable_group(bertini.VariableGroup([x]))
    return sys

def start_configuration(degree, center, radius, phi_over_pi):
    """The roots of the start configuration f(x) = x^d - d*x - c(t=1), where theta = phi at t=1."""
    x = bertini.Variable('x')
    theta_start = bertini.coefficient(phi_over_pi) * bertini.Pi
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([x]))
    sys.add_function(x**degree - degree * x - _c_of_theta(theta_start, center, radius))
    solver = bertini.nag_algorithm.ZeroDimSolver(sys, mptype='adaptive')
    solver.solve()
    return [complex(s[0]) for s in solver.all_solutions()]

def track_loop(degree, center, radius, phi_over_pi, tol=1e-10):
    """Track every strand once around the loop.  Returns a list of per-strand dicts with the loop
    angle ``theta``, the complex position ``x``, and the tracker diagnostics along the path."""
    H = loom_homotopy(degree, center, radius, phi_over_pi)
    roots = start_configuration(degree, center, radius, phi_over_pi)

    tracker = bertini.AMPTracker(H)
    tracker.setup(bertini.tracking.Predictor.RK4, tol, 1e6,
                  tracking.SteppingConfig(), tracking.NewtonConfig())
    tracker.precision_setup(tracking.amp_config_from(H))

    strands = []
    for r in roots:
        collector = amp_observers.PathDataCollector()
        tracker.add_observer(collector)
        end = np.zeros(H.num_variables(), dtype=complex_mp)
        tracker.track_path(end, complex_mp(1.0, 0.0), complex_mp(0.0, 0.0),
                           np.array([complex_mp(r.real, r.imag)]))
        tracker.remove_observer(collector)

        times = collector.times()               # complex path-variable value per step
        pts = collector.points()[:, 0]           # affine x (one variable)
        dgn = collector.diagnostics()            # (n, 4): |t|, condition, precision, stepsize
        strands.append(dict(start=r, end=complex(end[0]),
                            theta=2 * math.pi * (1 - times.real), x=pts,
                            precision=dgn[:, 2], stepsize=dgn[:, 3]))
    return strands, roots


# --- rendering ----------------------------------------------------------------------------------

def _resample(strand, n):
    """Interpolate one strand onto a uniform theta grid (for smooth ribbons and pinch-finding)."""
    th = strand['theta']
    order = np.argsort(th)
    g = np.linspace(th[order].min(), th[order].max(), n)
    xr = np.interp(g, th[order], strand['x'][order].real)
    xi = np.interp(g, th[order], strand['x'][order].imag)
    logstep = np.interp(g, th[order], np.log10(np.maximum(strand['stepsize'][order], 1e-12)))
    return g, xr, xi, logstep

def _find_pinches(strands, ngrid=900, thresh=0.6):
    """The near-collisions: local minima of pairwise strand distance below ``thresh``."""
    grid = np.linspace(0, 2 * math.pi, ngrid)
    X = []
    for s in strands:
        th = s['theta']; o = np.argsort(th)
        X.append(np.interp(grid, th[o], s['x'][o].real) + 1j * np.interp(grid, th[o], s['x'][o].imag))
    X = np.array(X)
    found = []
    for i in range(len(X)):
        for j in range(i + 1, len(X)):
            d = np.abs(X[i] - X[j])
            for k in range(1, ngrid - 1):
                if d[k] < thresh and d[k] <= d[k - 1] and d[k] < d[k + 1]:
                    mid = 0.5 * (X[i, k] + X[j, k])
                    found.append((mid.real, mid.imag, grid[k]))
    found.sort(key=lambda p: p[2])
    merged = []
    for p in found:
        if not merged or abs(p[2] - merged[-1][2]) > 0.25:
            merged.append(p)
    return merged

def render(strands, out, title, subtitle):
    """Draw the braid: strands coloured by identity, brightness and width driven by tracker
    stress (step-size collapse), branch-point pinches glowing, start rings -> end stars."""
    plt.rcParams.update({'figure.facecolor': _BG, 'axes.facecolor': _BG})
    fig = plt.figure(figsize=(9, 11))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor(_BG)

    all_ls = np.concatenate([np.log10(np.maximum(s['stepsize'], 1e-12)) for s in strands])
    calm, hard = all_ls.max(), all_ls.min()          # big step = calm, small step = stressed

    def stress(ls):
        return np.clip((calm - ls) / (calm - hard + 1e-9), 0.0, 1.0)

    for i, s in enumerate(strands):
        g, xr, xi, ls = _resample(s, 700)
        st = stress(ls)
        pts = np.column_stack([xr, xi, g])
        segs = np.stack([pts[:-1], pts[1:]], axis=1)
        col = _HUES[i % len(_HUES)]
        for width, alpha, colour in ((6 + 26 * st[:-1], 0.10, col),      # outer glow
                                     (3 + 8 * st[:-1], 0.22, col),        # inner glow
                                     (0.8 + 1.6 * st[:-1], 0.92, 'white')):  # hot core
            lc = Line3DCollection(segs, colors=colour, alpha=alpha)
            lc.set_linewidth(width)
            ax.add_collection3d(lc)
        ax.scatter([xr[0]], [xi[0]], [g[0]], s=70, facecolors='none',
                   edgecolors=col, linewidths=1.6, depthshade=False)                  # start ring
        ax.scatter([xr[-1]], [xi[-1]], [g[-1]], marker='*', s=180, color=col,
                   depthshade=False, zorder=6)                                        # end star

    for (xr, xi, th) in _find_pinches(strands):
        ax.scatter([xr], [xi], [th], marker='*', s=460, color='white', depthshade=False, zorder=8)
        ax.scatter([xr], [xi], [th], marker='*', s=1500, color='#fff2c0', alpha=0.26, depthshade=False)

    ax.set_xlabel('Re(x)', color='#c9d3e0')
    ax.set_ylabel('Im(x)', color='#c9d3e0')
    ax.set_zlabel('loop angle  θ  (0 → 2π)', color='#c9d3e0')
    ax.set_zticks([0, math.pi, 2 * math.pi]); ax.set_zticklabels(['0', 'π', '2π'])
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_pane_color((0.02, 0.03, 0.06, 1.0)); a.line.set_color('#223')
    ax.tick_params(colors='#8896aa')
    ax.grid(False)
    ax.view_init(elev=16, azim=-58)
    ax.set_box_aspect((1, 1, 1.7))
    ax.set_title(title + '\n' + subtitle, color='#e8eef7', fontsize=11, pad=8)
    fig.savefig(out, dpi=150, facecolor=_BG, bbox_inches='tight')
    plt.close(fig)


# --- the two frames -----------------------------------------------------------------------------

def teaching_frame(out):
    """x^3 - 3x - c: a loop around ONE branch point (c = +2) -> a single transposition."""
    bertini.random.set_random_seed(1)     # deterministic start-root ordering -> stable render
    # x^d - d*x has branch points at c = +/-(d-1); for d=3 that is +/-2.  Centre the loop right of
    # +2 so it encircles +2, excludes -2, and its closest approach (the pinch) is at theta = pi.
    graze, radius = Fraction(1, 50), Fraction(3, 2)     # 0.02, 1.5
    center = (3 - 1) + (radius - graze)                 # 2 + (radius - graze)
    strands, _ = track_loop(3, center, radius, phi_over_pi=0)
    render(strands, out,
           'The Monodromy Loom — teaching case:  $x^3 - 3x - c$',
           'loop one branch point  →  two roots SWAP (a transposition); the third rides straight')

def showpiece_frame(out):
    """x^5 - 5x - c: a loop around ALL FOUR branch points -> a full 5-cycle, four pinches."""
    bertini.random.set_random_seed(1)     # deterministic start-root ordering -> stable render
    radius = (5 - 1) + Fraction(3, 20)     # circle |c| = R just outside the branch orbit |c| = 4
    strands, _ = track_loop(5, 0, radius, phi_over_pi=Fraction(1, 4))   # phi = pi/4 keeps pinches off the seam
    render(strands, out,
           'The Monodromy Loom — showpiece:  $x^5 - 5x - c$',
           'loop encircles four branch points  →  a full 5-cycle; four pinches where the tracker sweats')


def main():
    bertini.recording(False)          # we are here to WATCH tracking, not recall it
    teaching_frame(os.path.join(_OUT, 'monodromy_loom_teaching.png'))
    showpiece_frame(os.path.join(_OUT, 'monodromy_loom.png'))


if __name__ == '__main__':
    main()
