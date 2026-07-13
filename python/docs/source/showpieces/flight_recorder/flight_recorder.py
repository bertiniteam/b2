"""The Flight Recorder.

A showpiece: film one brutally hard homotopy path -- the descent to a highly *singular*
solution -- and read out the adaptive-precision tracker's full telemetry, step by step.

The subject is the origin (0, 0), where two rotated rose curves r = sin(m*theta) and
r = sin(n*theta) meet.  For (m, n) = (7, 5) that intersection has multiplicity 35: thirty-five
homotopy paths pile into the same point, and the ones that get there have to fight for every
digit.  We ask the solver for tight final accuracy, then watch a single path's instruments as it
goes:

    * the endgame **spiral** -- the Cauchy endgame samples a circle around the singular endpoint;
      with cycle number c the solution winds c times as t -> 0 (here c = 7), a log-radial spiral;
    * **precision** climbing 16 -> 20 -> 30 -> ... digits as adaptive precision escalates into
      mpfr to keep the accuracy the tight tolerance demands;
    * the **condition number** blowing up by many orders of magnitude as the Jacobian degenerates;
    * the **step size** sawtoothing -- grown when the going is easy, cut hard (a rejected step)
      when it is not.

All of it is real, captured by a ``PathDataCollector`` on the path's own tracker (via a
``SolutionPathCollector`` over the whole solve), then laid out as a cockpit.  The point of the
family is that it is *crankable*: raise (m, n) or tighten the tolerance and the singular point --
and the tracker's struggle -- gets arbitrarily worse.

Regenerated through ``tools/refresh_doc_artifacts.py``.  Raster (PNG) showpiece, not a doctest.

Run standalone:  python flight_recorder.py
"""

import math
import os

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection

import bertini
from bertini import ZeroDimSolver, SolutionPathCollector
from bertini.sympy_bridge import from_sympy

_OUT = os.path.dirname(os.path.abspath(__file__))
_BG = '#070a10'
_INK = '#d6e2f0'
_DIM = '#7f8da0'


# --- the crankable singular system --------------------------------------------------------------

def _rose(k):
    """The rectangular equation of the rose r = sin(k*theta): (x^2+y^2)^((k+1)/2) = Im[(x+iy)^k],
    a polynomial for odd k.  Returned as a bertini function-tree node via the sympy bridge."""
    from sympy import symbols, im, I
    xs, ys = symbols('x y', real=True)
    return from_sympy((xs**2 + ys**2)**((k + 1) // 2) - im((xs + I * ys)**k))

def system_rhodonea(m, n):
    """Two rose curves, the second rotated by a random angle so their only structured coincidence
    is the highly singular meeting at the origin.  (m, n) = (7, 5) -> multiplicity 35 at (0,0).
    Returns (system, rotation_angle) so the geometry can be drawn to match the solved system."""
    x, y = bertini.variables(list('xy'))
    f1, f2 = _rose(m), _rose(n)
    t = bertini.random_real()                       # random rotation (seed fixed by the caller)
    angle = float(t.real)
    rx = bertini.cos(t) * x + bertini.sin(t) * y
    ry = -bertini.sin(t) * x + bertini.cos(t) * y
    f2 = f2.subs({x: rx, y: ry})
    sys = bertini.System()
    sys.add_variable_group(x, y)
    sys.add([f1, f2])
    return sys, angle

def rose_curve(k, rotation=0.0, n=1400):
    """Sample the real rose r = sin(k*theta) as (X, Y), optionally rotated to match the system."""
    th = np.linspace(0, 2 * math.pi, n)
    r = np.sin(k * th)
    X, Y = r * np.cos(th), r * np.sin(th)
    c, s = math.cos(rotation), math.sin(rotation)
    return c * X - s * Y, s * X + c * Y

def _descent_spine(path):
    """Keep only the real-time continuation steps: those where |t| reaches a new minimum.  This
    drops the Cauchy endgame's circular sampling (steps at constant |t|), leaving the loop-free
    descent.  Returns (abs_t, affine) on the spine, with the singular endpoint (0,0) appended."""
    abs_t = np.abs(path.times())
    affine = path.points()[:, 1:] / path.points()[:, 0:1]
    runmin = np.minimum.accumulate(abs_t)
    keep = np.concatenate([[True], runmin[1:] < runmin[:-1]])
    at, af = abs_t[keep], affine[keep]
    return np.append(at, at[-1] * 1e-2), np.vstack([af, np.zeros((1, af.shape[1]))])


# --- fly one hard path and record everything ----------------------------------------------------

def record_hard_path(m=7, n=5, final_tolerance=1e-20, seed=2):
    """Solve system_rhodonea(m, n) with tight accuracy, collecting every path; return the richest
    singular path's telemetry (the one that escalates precision the most) plus its metadata."""
    bertini.recording(False)                        # WATCH tracking, do not recall it
    bertini.random.set_random_seed(seed)            # deterministic system + solve -> stable picture

    system, rotation = system_rhodonea(m, n)
    solver = ZeroDimSolver(system, mptype='adaptive')
    tol = solver.get_config(bertini.nag_algorithm.TolerancesConfig)
    tol.final_tolerance = final_tolerance
    solver.set_config(tol)
    cfg = solver.get_config(bertini.nag_algorithm.ZeroDimConfig)
    cfg.num_threads = 1                              # serial -> deterministic path ordering / picture
    solver.set_config(cfg)

    collector = SolutionPathCollector()
    solver.add_observer(collector)
    solver.solve()

    meta = {int(md.path_index): md for md in solver.solution_metadata()}
    cols = collector.series[0].DIAGNOSTIC_COLUMNS
    P = cols.index('precision')

    best = None
    singular_paths = []
    for path in collector.series:
        md = meta.get(path.path_index)
        if not (md and md.is_singular):
            continue
        singular_paths.append(path)
        dgn = path.diagnostics()
        key = (int(dgn[:, P].max()), len(dgn))       # most precision, then most steps
        if best is None or key > best[0]:
            best = (key, path, md)
    _, path, md = best

    dgn = path.diagnostics()
    cockpit = dict(m=m, n=n, final_tolerance=final_tolerance, path_index=int(path.path_index),
                   affine=path.points()[:, 1:] / path.points()[:, 0:1],
                   abs_t=dgn[:, 0], condition=dgn[:, 1], precision=dgn[:, 2], stepsize=dgn[:, 3],
                   cycle=int(md.cycle_num), multiplicity=int(md.multiplicity),
                   precision_digits=int(md.precision_digits), accuracy_digits=int(md.accuracy_digits))

    # system-level companion: the two roses, and every singular path's loop-free descent
    setup = dict(m=m, n=n, rotation=rotation, multiplicity=int(md.multiplicity),
                 spines=[_descent_spine(p) for p in singular_paths])
    return cockpit, setup


# --- the cockpit --------------------------------------------------------------------------------

def _style_axis(ax):
    ax.set_facecolor(_BG)
    for s in ax.spines.values():
        s.set_color('#26324a')
    ax.tick_params(colors='#8fa0bb', labelsize=8)
    ax.grid(True, color='#141b28', lw=0.7)
    ax.xaxis.label.set_color(_INK); ax.yaxis.label.set_color(_INK)

def render(rec, out):
    plt.rcParams.update({'figure.facecolor': _BG, 'axes.facecolor': _BG,
                         'font.family': 'monospace'})
    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.15, 1.0], height_ratios=[1, 1, 1],
                          hspace=0.5, wspace=0.34, left=0.06, right=0.97, top=0.86, bottom=0.09)

    n_steps = len(rec['abs_t'])
    step = np.arange(n_steps)
    eg = rec['abs_t'] < 0.1                          # the endgame portion (small |t|)
    boundary = int(np.argmax(eg)) if eg.any() else n_steps

    # ---- the endgame spiral (left, spanning all rows) ----
    axS = fig.add_subplot(gs[:, 0]); _style_axis(axS)
    xv = rec['affine'][:, 0]                          # one coordinate of the solution
    xe, te = xv[eg], rec['abs_t'][eg]
    r = np.log10(np.abs(xe)) - np.log10(np.abs(xe).min()) + 0.15    # log-radial: 0 -> singular pt
    disp = r * np.exp(1j * np.angle(xe))
    pts = np.column_stack([disp.real, disp.imag])
    segs = np.stack([pts[:-1], pts[1:]], axis=1)
    lc = LineCollection(segs, cmap='turbo',
                        norm=mcolors.LogNorm(max(te.min(), 1e-30), te.max()))
    lc.set_array(0.5 * (te[:-1] + te[1:])); lc.set_linewidth(1.7)
    axS.add_collection(lc)
    axS.scatter([0], [0], marker='*', s=320, color='white', zorder=6)
    axS.scatter([0], [0], marker='*', s=1100, color='#fff0b0', alpha=0.25, zorder=5)
    R = np.abs(disp).max() * 1.1
    axS.set_xlim(-R, R); axS.set_ylim(-R, R); axS.set_aspect(1.0)
    groups = rec['multiplicity'] // rec['cycle'] if rec['cycle'] else 0
    axS.set_title(f"Cauchy endgame spiral into the singular point  (cycle number c = {rec['cycle']})\n"
                  f"this path winds {rec['cycle']}× as t → 0; the mult-{rec['multiplicity']} point is "
                  f"{groups} cyclic groups of {rec['cycle']}  ({groups}×{rec['cycle']} = {rec['multiplicity']})",
                  color=_INK, fontsize=9.5, pad=6)
    cb = fig.colorbar(lc, ax=axS, fraction=0.035, pad=0.015)
    cb.set_label('|t| (log)', color=_INK, labelpad=-2)
    cb.ax.tick_params(colors='#8fa0bb')

    def gauge(ax, y, label, color, logy=False, drops=None):
        ax.plot(step, y, color=color, lw=1.4)
        if drops is not None:
            ax.plot(step[drops], y[drops], 'v', color='#ff5a5a', ms=4, alpha=0.8)
        if logy:
            ax.set_yscale('log')
        ax.axvline(boundary, color='#4de0c0', lw=1.0, ls=(0, (3, 3)), alpha=0.8)
        ax.set_ylabel(label); _style_axis(ax)
        ax.set_xlim(0, n_steps - 1)

    # ---- precision staircase ----
    axP = fig.add_subplot(gs[0, 1])
    gauge(axP, rec['precision'], 'precision\n(digits)', '#ffd24a')
    axP.set_ylim(min(rec['precision']) - 2, max(rec['precision']) + 6)
    axP.text(boundary, axP.get_ylim()[1], ' endgame', color='#4de0c0', fontsize=7, va='top')

    # ---- condition number ----
    axC = fig.add_subplot(gs[1, 1])
    gauge(axC, np.maximum(rec['condition'], 1.0), 'condition\nnumber', '#ff6fae', logy=True)

    # ---- step size (sawtooth), failed steps marked ----
    axZ = fig.add_subplot(gs[2, 1])
    drops = np.where(np.diff(rec['stepsize']) < 0)[0] + 1        # a cut step size = a rejected step
    gauge(axZ, np.maximum(rec['stepsize'], 1e-30), 'step size', '#5ad1ff', logy=True, drops=drops)
    axZ.set_xlabel('tracker step  (its own clock →)')
    axZ.plot([], [], 'v', color='#ff5a5a', ms=5, label='step cut')
    axZ.legend(loc='lower left', fontsize=7, facecolor=_BG, edgecolor='#26324a', labelcolor=_INK)

    # ---- title / readout bar ----
    fig.text(0.06, 0.955, "✈  FLIGHT RECORDER", color='#4de0c0', fontsize=16, fontweight='bold')
    fig.text(0.06, 0.905,
             f"one path to a multiplicity-{rec['multiplicity']} singular point  ·  "
             f"$x^{{{rec['m']}}}$-rose ∩ $x^{{{rec['n']}}}$-rose at (0,0)  ·  "
             f"final_tolerance = {rec['final_tolerance']:.0e}",
             color=_INK, fontsize=10)
    fig.text(0.97, 0.955,
             f"cycle {rec['cycle']}   ·   precision {int(rec['precision'].min())}→{int(rec['precision'].max())} digits   "
             f"·   condition ×10^{int(np.log10(rec['condition'].max()))}   ·   {n_steps} steps",
             color='#8fa0bb', fontsize=9, ha='right')

    fig.savefig(out, dpi=150, facecolor=_BG)
    plt.close(fig)
    return rec


def render_setup(setup, out):
    """The system-level companion (whole solve, not one path): the two rose curves whose crossing
    at (0,0) is the singular target, and all the homotopy paths converging on it loop-free."""
    plt.rcParams.update({'figure.facecolor': _BG, 'axes.facecolor': _BG, 'font.family': 'monospace'})
    fig, (axR, axC) = plt.subplots(1, 2, figsize=(13, 6.6))
    fig.subplots_adjust(left=0.06, right=0.95, top=0.84, bottom=0.09, wspace=0.24)

    # --- the geometry: two real rose curves meeting at the origin ---
    _style_axis(axR)
    Xa, Ya = rose_curve(setup['m'], 0.0)
    Xb, Yb = rose_curve(setup['n'], setup['rotation'])
    axR.plot(Xa, Ya, color='#ff6fae', lw=1.6, label=f"$r=\\sin({setup['m']}\\theta)$")
    axR.plot(Xb, Yb, color='#5ad1ff', lw=1.6, label=f"$r=\\sin({setup['n']}\\theta)$, rotated")
    axR.scatter([0], [0], s=200, facecolors='none', edgecolors='white', linewidths=1.3, zorder=6)  # ring, not covering
    axR.set_aspect(1.0); axR.set_xlabel('x'); axR.set_ylabel('y')
    axR.legend(loc='upper right', fontsize=9, facecolor=_BG, edgecolor='#26324a', labelcolor=_INK)
    axR.set_title(f"the problem — two rose curves meet at (0,0)\n"
                  f"a multiplicity-{setup['multiplicity']} singular intersection", color=_INK, fontsize=10, pad=6)

    # --- the solve: every singular path's loop-free descent, converging on (0,0) ---
    _style_axis(axC)
    cmap = plt.get_cmap('turbo')
    npath = len(setup['spines'])
    order = np.argsort([float(np.angle(s[1][0, 0])) for s in setup['spines']])  # hue by start angle
    for rank, i in enumerate(order):
        abs_t, affine = setup['spines'][i]
        xv = affine[:, 0]                             # one coordinate's complex plane
        pts = np.column_stack([xv.real, xv.imag])
        segs = np.stack([pts[:-1], pts[1:]], axis=1)
        rgb = cmap((rank + 0.5) / npath)              # a distinct hue per path -> 35 followable threads
        axC.add_collection(LineCollection(segs, colors=[rgb], linewidths=1.1, alpha=0.72))
        axC.plot(xv.real[0], xv.imag[0], 'o', color=rgb, ms=6, mec='white', mew=0.7, zorder=7)  # start (t≈1)
    axC.scatter([0], [0], s=170, facecolors='none', edgecolors='white', linewidths=1.4, zorder=8)  # target ring
    axC.autoscale(); axC.set_aspect(1.0)
    axC.set_xlabel('Re(x)'); axC.set_ylabel('Im(x)')
    axC.set_title(f"the solve — {npath} homotopy paths converge on it\n"
                  "real-time continuation (Cauchy endgame loops filtered out)", color=_INK, fontsize=10, pad=6)

    fig.text(0.06, 0.945, "THE SINGULAR RENDEZVOUS", color='#4de0c0', fontsize=14, fontweight='bold')
    fig.savefig(out, dpi=150, facecolor=_BG)
    plt.close(fig)


def main():
    cockpit, setup = record_hard_path(m=7, n=5, final_tolerance=1e-24, seed=2)
    render_setup(setup, os.path.join(_OUT, 'flight_recorder_setup.png'))
    render(cockpit, os.path.join(_OUT, 'flight_recorder.png'))


if __name__ == '__main__':
    main()
