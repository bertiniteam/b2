"""Watching the paths: observers and path data.

An observer is a small object you attach to a tracker (or any observable); its ``Observe``
method is called with an event every time something happens along a path.  This script builds
up from a one-line observer, to ``PathDataCollector`` recording a single path into numpy, to the
``SolutionPathCollector`` meta-observer that records every path of a whole solve -- and plots the
results.  It produces three figures:

    * observers_and_path_data  -- the six homotopy paths of z^6 - 2z^2 + 2 in the complex plane
    * cyclic3_paths            -- the cyclic-3 paths in 3-D, coloured by condition number
    * griewank_osborn_endgame  -- the Cauchy endgame loops at a singular solution, log-radial

Run:  python observers_and_path_data.py
"""

import os
from fractions import Fraction

import numpy as np

import matplotlib
matplotlib.use('Agg')                 # headless: no display needed
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection

import bertini
import bertini.tracking as tracking
from bertini import ZeroDimSolver, SolutionPathCollector
from bertini.multiprec import complex_mp
from bertini.nag_algorithm import observers as nag_observers

_OUT = os.path.dirname(os.path.abspath(__file__))


# --- Writing an observer in Python --------------------------------------------------------------

class StepPrinter(tracking.observers.amp.CustomObserver):
    def Observe(self, event):
        if isinstance(event, tracking.observers.amp.SuccessfulStep):
            trk = event.tracker()
            print("t =", complex(trk.current_time()),
                  " |z| =", trk.current_point(),
                  " cond =", float(trk.latest_condition_number()))


def custom_observer():
    """Build the y - t homotopy, attach a StepPrinter, track, and detach.  Returns the tracker."""
    y, t = bertini.Variable('y'), bertini.Variable('t')
    sys = bertini.System()
    sys.add_function(y - t)
    sys.add_path_variable(t)
    sys.add_variable_group(bertini.VariableGroup([y]))

    tracker = bertini.AMPTracker(sys)
    tracker.setup(bertini.Predictor.Euler, 1e-5, 1e5,
                  tracking.SteppingConfig(), tracking.NewtonConfig())
    tracker.precision_setup(tracking.amp_config_from(sys))

    printer = StepPrinter()
    tracker.add_observer(printer)
    end = np.zeros(sys.num_variables(), dtype=complex_mp)
    tracker.track_path(end, complex_mp(1), complex_mp(0), np.array([complex_mp(1)]))
    tracker.remove_observer(printer)

    # the ready-made CallbackObserver: attach a function without writing a class
    obs = tracking.observers.amp.CallbackObserver()
    obs.on(tracking.observers.amp.PrecisionChanged,
           lambda e: print(e.previous(), "->", e.next()))
    tracker.add_observer(obs)
    tracker.remove_observer(obs)

    return tracker, sys


# --- Collecting one path into numpy -------------------------------------------------------------

def collect_one_path(tracker, sys):
    """PathDataCollector records one path's time, points, and diagnostics into numpy arrays."""
    b = tracking.observers.amp.PathDataCollector()
    tracker.add_observer(b)
    end = np.zeros(sys.num_variables(), dtype=complex_mp)
    tracker.track_path(end, complex_mp(1), complex_mp(0), np.array([complex_mp(1)]))
    tracker.remove_observer(b)

    ts  = b.times()          # complex,  shape (n_steps,)
    zs  = b.points()         # complex,  shape (n_steps, n_vars)
    dgn = b.diagnostics()    # float,    shape (n_steps, 4): |t|, condition number, precision, stepsize

    assert len(ts) > 0 and zs.shape[1] == sys.num_variables()
    return ts, zs, dgn


# --- The meta-observer: every path of a whole solve ---------------------------------------------

def solve_degree_six():
    """Solve z^6 - 2z^2 + 2 and collect every path with a SolutionPathCollector.  Returns solver, A."""
    bertini.random.set_random_seed(2)   # so you get exactly this picture

    z = bertini.Variable('z')
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([z]))
    sys.add_function(z**6 - 2*z**2 + 2)

    solver = ZeroDimSolver(sys, mptype='adaptive')

    A = SolutionPathCollector()
    solver.add_observer(A)
    solver.solve()

    assert len(A.series) == 6          # one PathDataCollector per solution path
    return solver, A


def plot_degree_six(solver, A):
    """Draw the six paths in the complex plane; save as observers_and_path_data.{svg,png}."""
    fig, ax = plt.subplots(figsize=(6, 6))
    cmap = plt.get_cmap('turbo')
    for i, path in enumerate(A.series):
        pts  = path.points()
        zaff = pts[:, 1] / pts[:, 0]                 # dehomogenize to affine z
        color = cmap(i / max(len(A.series) - 1, 1))
        ax.plot(zaff.real, zaff.imag, '-', color=color, lw=1.3)
        ax.plot(zaff.real[0], zaff.imag[0], 'o', color=color, ms=6, mfc='white')  # start, t near 1

    sols = [complex(s[0]) for s in solver.all_solutions()]
    ax.scatter([s.real for s in sols], [s.imag for s in sols],
               c='k', marker='*', s=140, zorder=5, label='solutions')
    ax.set_aspect(1.0); ax.set_box_aspect(1)     # 1:1 data scaling, square box
    ax.set_xlabel('Re(z)'); ax.set_ylabel('Im(z)')
    ax.legend(loc='upper right', fontsize=8)
    fig.savefig(os.path.join(_OUT, 'observers_and_path_data.svg'))
    fig.savefig(os.path.join(_OUT, 'observers_and_path_data.png'), dpi=150)


# --- Build it yourself: one observer that attaches another ---------------------------------------

class PathRecorder(tracking.observers.amp.CustomObserver):
    def __init__(self, parent, path_index):
        super().__init__()
        self.parent     = parent          # the observer we report back to
        self.path_index = path_index
        self._times     = []
        self._points    = []

    def Observe(self, event):
        if isinstance(event, tracking.observers.amp.SuccessfulStep):
            trk = event.tracker()
            # copy out *now* -- the event and tracker state are valid only during this call
            self._times.append(complex(trk.current_time()))
            self._points.append([complex(z) for z in trk.current_point()])

    def report(self):
        self.parent.receive(self.path_index,
                            np.array(self._times, dtype=complex),
                            np.array(self._points, dtype=complex))


class MyPathCollector(nag_observers.CustomObserver):
    def __init__(self):
        super().__init__()
        self.paths   = {}     # path_index -> (times, points), filled in by B.report()
        self._active = {}     # path_index -> (tracker, live recorder)

    def Observe(self, event):
        if isinstance(event, nag_observers.PathStarted):
            tracker = event.tracker()      # the tracker that RUNS this path (a thread-local
            b = PathRecorder(self, event.path_index())   # clone when threaded -- not the member
            tracker.add_observer(b)                  # attach B ...   tracker)
            self._active[event.path_index()] = (tracker, b)
        elif isinstance(event, nag_observers.PathComplete):
            tracker, b = self._active.pop(event.path_index())
            tracker.remove_observer(b)               # ... and detach it when the path is done
            b.report()                               # B hands its data back to A

    def receive(self, path_index, times, points):    # B calls this
        self.paths[path_index] = (times, points)


def hand_built_collector():
    """Reconstruct SolutionPathCollector from scratch with MyPathCollector + PathRecorder."""
    bertini.random.set_random_seed(2)
    z = bertini.Variable('z')
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([z]))
    sys.add_function(z**6 - 2*z**2 + 2)

    solver = ZeroDimSolver(sys, mptype='adaptive')
    A = MyPathCollector()
    solver.add_observer(A)
    solver.solve()

    assert len(A.paths) == 6                  # same six paths SolutionPathCollector would give you
    assert all(len(times) > 0 for times, _ in A.paths.values())   # each path actually captured steps
    return A


# --- Observers, threads, and MPI: forcing a serial solve ----------------------------------------

def serial_solve():
    """Pin the solve to a single thread with num_threads = 1; SolutionPathCollector works the same."""
    z = bertini.Variable('z')
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([z]))
    sys.add_function(z**6 - 2*z**2 + 2)

    solver = ZeroDimSolver(sys, mptype='adaptive')
    cfg = solver.get_config(bertini.nag_algorithm.ZeroDimConfig)
    cfg.num_threads = 1                 # 0 = auto (all cores), 1 = serial, N = N threads
    solver.set_config(cfg)

    A = SolutionPathCollector()         # works the same in serial and threaded
    solver.add_observer(A)
    solver.solve()
    assert len(A.series) == 6
    return A


# --- A 3-D system, coloured by condition number -------------------------------------------------

def solve_cyclic3():
    """Solve the cyclic-3 system and collect every path.  Returns solver, A, COND column index."""
    bertini.random.set_random_seed(3)

    x, y, z = bertini.Variable('x'), bertini.Variable('y'), bertini.Variable('z')
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([x, y, z]))
    sys.add_function(x + y + z)
    sys.add_function(x*y + y*z + z*x)
    sys.add_function(x*y*z - 1)

    solver = ZeroDimSolver(sys, mptype='adaptive')
    A = SolutionPathCollector()
    solver.add_observer(A)
    solver.solve()

    COND = A.series[0].DIAGNOSTIC_COLUMNS.index("condition_number")
    return solver, A, COND


def plot_cyclic3(solver, A, COND):
    """Draw the cyclic-3 paths in 3-D coloured by condition number; save cyclic3_paths.{svg,png}."""
    tracks, all_cond = [], []
    for path in A.series:
        aff  = path.points()[:, 1:] / path.points()[:, 0:1]   # dehomogenize -> 3 affine coords
        real = aff.real                                        # (n, 3): Re(x), Re(y), Re(z)
        cond = path.diagnostics()[:, COND]
        tracks.append((real, cond)); all_cond.append(cond)

    norm = mcolors.LogNorm(vmin=max(min(c.min() for c in all_cond), 1.0),
                           vmax=max(c.max() for c in all_cond))

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    for real, cond in tracks:
        segs = np.stack([real[:-1], real[1:]], axis=1)         # consecutive points -> segments
        lc = Line3DCollection(segs, cmap='viridis', norm=norm)
        lc.set_array(0.5 * (cond[:-1] + cond[1:]))             # colour each segment by condition number
        lc.set_linewidth(2)
        ax.add_collection3d(lc)
        last = lc

    sols = solver.all_solutions()
    ax.scatter([complex(s[0]).real for s in sols],
               [complex(s[1]).real for s in sols],
               [complex(s[2]).real for s in sols], c='k', marker='*', s=120, depthshade=False)
    ax.set_xlabel('Re(x)'); ax.set_ylabel('Re(y)'); ax.set_zlabel('Re(z)')
    ax.set_box_aspect((1, 1, 1))                               # 1:1:1 data aspect ratio
    fig.colorbar(last, ax=ax, shrink=0.6, pad=0.1, label='condition number (log)')
    fig.savefig(os.path.join(_OUT, 'cyclic3_paths.svg'))
    fig.savefig(os.path.join(_OUT, 'cyclic3_paths.png'), dpi=150)


# --- Watching the Cauchy endgame at a singular solution -----------------------------------------

def solve_griewank_osborn():
    """Solve Griewank-Osborn and pick out the singular paths.  Returns the singular path list."""
    bertini.random.set_random_seed(1)

    x, y = bertini.Variable('x'), bertini.Variable('y')
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([x, y]))
    sys.add_function(bertini.coefficient(Fraction(29, 16)) * x**3 - 2*x*y)  # exact rational coeff
    sys.add_function(y - x**2)

    solver = ZeroDimSolver(sys, mptype='adaptive')
    A = SolutionPathCollector()
    solver.add_observer(A)
    solver.solve()

    # let the solver classify which paths ended at a singular solution
    singular_idx = {int(m.path_index) for m in solver.solution_metadata() if m.is_singular}
    singular = [p for p in A.series if p.path_index in singular_idx]
    return singular


def plot_griewank_osborn(singular):
    """Draw the Cauchy endgame loops on a log-radial scale; save griewank_osborn_endgame.{svg,png}."""
    path = max(singular, key=len)                 # one representative singular path
    aff  = path.points()[:, 1:] / path.points()[:, 0:1]
    xv   = aff[:, 0]
    t    = np.abs(path.times())
    eg   = t < 0.1                                 # the endgame portion (small |t|)
    xv, t = xv[eg], t[eg]

    r    = np.log10(np.abs(xv)) - np.log10(np.abs(xv).min()) + 0.1   # log-radial coordinate
    disp = r * np.exp(1j * np.angle(xv))

    fig, ax = plt.subplots(figsize=(6, 6))
    pts  = np.column_stack([disp.real, disp.imag])
    segs = np.stack([pts[:-1], pts[1:]], axis=1)
    lc = LineCollection(segs, cmap='plasma', norm=mcolors.LogNorm(t.min(), t.max()))
    lc.set_array(0.5 * (t[:-1] + t[1:]))
    ax.add_collection(lc)
    ax.scatter([0], [0], c='k', marker='*', s=160, zorder=5)  # the singular solution, at the centre
    R = -(np.log10(np.abs(xv).min())) + 0.3
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_aspect(1.0)
    fig.colorbar(lc, ax=ax, label='|t|  (log)')
    fig.savefig(os.path.join(_OUT, 'griewank_osborn_endgame.svg'))
    fig.savefig(os.path.join(_OUT, 'griewank_osborn_endgame.png'), dpi=150)


def main():
    # observer basics + one-path collection on the y - t homotopy
    tracker, sys = custom_observer()
    collect_one_path(tracker, sys)

    # meta-observer: every path of the degree-six solve, then plot (figure 1)
    solver, A = solve_degree_six()
    plot_degree_six(solver, A)

    # hand-built collector reconstructing SolutionPathCollector, and the serial-solve variant
    hand_built_collector()
    serial_solve()

    # cyclic-3 in 3-D coloured by condition number (figure 2)
    c3_solver, c3_A, COND = solve_cyclic3()
    plot_cyclic3(c3_solver, c3_A, COND)

    # Cauchy endgame at a singular Griewank-Osborn solution (figure 3)
    singular = solve_griewank_osborn()
    plot_griewank_osborn(singular)


if __name__ == '__main__':
    main()
