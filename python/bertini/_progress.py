# This file is part of Bertini 2.
#
# python/bertini/_progress.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/_progress.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this file.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team

"""Say how far a solve has got, while it is still running.

A solve announces each path as it starts and finishes, so watching those events gives a
count without the solver having to know anything about display.
:class:`bertini.ProgressReport` is an observer that turns them into a bar:

    paths  41/64 [00:12<00:19]  38 ok, 2 diverged, 1 failed

It is on by default, and arranges to be invisible in every case where a bar would be an
intrusion rather than a help:

* **Not a terminal, nothing drawn.**  tqdm's own ``disable=None`` turns it off when the
  stream is not a TTY, so a piped script, a CI run, the test suite and the docs build all
  see nothing.  It writes to stderr, so even forced on it cannot contaminate captured
  output.
* **Fast solves stay silent.**  Nothing is drawn for the first
  :data:`DELAY_SECONDS`, so the many small solves inside a loop, a tutorial or a test
  finish before a bar would ever appear.
* **One rank draws.**  Under MPI only the manager rank reports; the workers would
  otherwise each draw their own.

**Redrawing is throttled, not the events.**  At a million paths, printing on every one is
the problem b1 solved with its ``printpathprogress`` setting; here the redraw interval
defaults to whatever gives a couple of hundred updates across the whole solve, at any
scale, and can be set to an exact number of paths.  Note what this does *not* save: the
two python callbacks per path still happen, because the solver emits the events either
way.  That is immaterial at thousands of paths and would not be at hundreds of millions,
where the answer is an observer that stays in C++ (issue filed).

**A recalled path is not a tracked path.**  A solve that finds its answers in the records
runs no paths at all and emits no path events, so the bar would sit at zero and then
vanish.  Recalled paths are reported on their own line instead, and the bar sizes itself
to the paths actually being tracked.
"""

import sys as _sys
import time as _time
from contextlib import contextmanager as _contextmanager

from bertini._pybertini.nag_algorithms import observers as _observers


DELAY_SECONDS = 5.0
"""How long a solve must run before anything is drawn."""

TARGET_UPDATES = 200
"""About how many times to redraw over a whole solve, when the interval is not given."""


def _tqdm():
    """The tqdm class, or None when tqdm is not installed."""
    try:
        from tqdm.auto import tqdm
    except ImportError:
        return None
    return tqdm


def _this_rank_reports():
    """Only the manager rank draws; a worker rank would draw over it."""
    try:
        import bertini.parallel as parallel
        return bool(parallel.is_manager())
    except Exception:
        return True


def _outcome_name(code):
    """Which tally a finished path belongs to."""
    name = getattr(code, 'name', str(code))
    if name == 'Success':
        return 'ok'
    if name in ('SecurityMaxNormReached', 'GoingToInfinity'):
        return 'diverged'
    if name in ('NeverStarted', 'ExternallyTerminated', 'WallClockLimitReached'):
        return 'stopped'
    return 'failed'


class ProgressReport(_observers.CustomObserver):
    """Watch a solve and report how far it has got.

    Attach it to a solver and it reports for one solve; :func:`bertini.solve` and
    ``ZeroDimSolver.solve`` do that for you unless told not to.

    Parameters
    ----------
    solver : the solver to watch
        Read for its path counts, which the events do not carry.
    every : int, optional
        Redraw once per this many finished paths.  The default divides the path count to
        give about :data:`TARGET_UPDATES` redraws over the solve.
    delay : float, optional
        Draw nothing until the solve has run this long.  Defaults to
        :data:`DELAY_SECONDS`.
    stream : file, optional
        Where to draw.  Defaults to stderr, so it never mixes with a program's output.
    force : bool, optional
        Draw even when the stream is not a terminal.  Off by default, which is what keeps
        a bar out of logs and captured output.
    """

    def __init__(self, solver, every=None, delay=None, stream=None, force=False):
        super().__init__()
        self.solver = solver
        self.every = every
        self.delay = DELAY_SECONDS if delay is None else delay
        self.stream = stream or _sys.stderr
        self.force = force

        self.total = 0
        self.recalled = 0
        self.finished = 0
        self.tallies = {}
        self.began = None

        self._bar = None
        self._since_redraw = 0
        self._sized = False

    # ----- the observer seam ---------------------------------------------------------

    def Observe(self, event):
        """Handle one solve event.  Called by the solver, possibly from a worker thread."""
        try:
            self._observe(event)
        except Exception:
            pass          # a progress display must never be the reason a solve fails

    def _observe(self, event):
        kind = type(event).__name__

        if kind == 'AlgorithmStarted':
            self._start()
        elif kind == 'PathComplete':
            self._path_done(event)
        elif kind == 'AlgorithmComplete':
            self._finish()

    # ----- the phases ----------------------------------------------------------------

    def _start(self):
        self.began = _time.time()
        self.total = int(self._ask('num_paths', 0))
        self.tallies = {}
        self.finished = 0
        self._sized = False

    def _size_to_the_work(self):
        """Once recall has been resolved, the bar counts only the paths being tracked."""
        self._sized = True
        self.recalled = int(self._ask('num_paths_recalled', 0))

        tracking = max(self.total - self.recalled, 0)
        if self.every is None:
            self.every = max(1, tracking // TARGET_UPDATES)

        tqdm = _tqdm()
        if tqdm is None or not _this_rank_reports() or tracking == 0:
            return

        self._bar = tqdm(total=tracking, desc='paths', unit='path',
                         disable=False if self.force else None,   # None: off when not a terminal
                         delay=self.delay, file=self.stream,
                         leave=True, smoothing=0.1)

    def _path_done(self, event):
        if not self._sized:
            self._size_to_the_work()

        self.finished += 1
        outcome = _outcome_name(self._outcome_of(event))
        self.tallies[outcome] = self.tallies.get(outcome, 0) + 1

        self._since_redraw += 1
        if self._bar is not None and self._since_redraw >= self.every:
            self._bar.update(self._since_redraw)
            self._bar.set_postfix_str(self.tally_text(), refresh=False)
            self._since_redraw = 0

    def _finish(self):
        if not self._sized:
            self._size_to_the_work()

        if self._bar is not None:
            if self._since_redraw:
                self._bar.update(self._since_redraw)
                self._since_redraw = 0
            self._bar.set_postfix_str(self.tally_text(), refresh=False)
            self._bar.close()
            self._bar = None

        # the recalled paths never reached the bar, having never been tracked
        if self.recalled and _this_rank_reports():
            self._say('recalled %d of %d paths from the records'
                      % (self.recalled, self.total))

    def _say(self, message):
        """One line, on the same stream and under the same conditions as the bar."""
        try:
            if self.force or self.stream.isatty():
                self.stream.write(message + '\n')
                self.stream.flush()
        except Exception:
            pass

    # ----- the parts a caller might want ---------------------------------------------

    def tally_text(self):
        """`38 ok, 2 diverged, 1 failed`, in a fixed order."""
        order = ('ok', 'diverged', 'failed', 'stopped')
        return ', '.join('%d %s' % (self.tallies[name], name)
                         for name in order if self.tallies.get(name))

    def close(self):
        """Take down the bar, if one is up.  Safe to call more than once."""
        if self._bar is not None:
            self._bar.close()
            self._bar = None

    def _ask(self, name, default):
        try:
            member = getattr(self.solver, name)
        except AttributeError:
            return default
        try:
            return member() if callable(member) else member
        except Exception:
            return default

    @staticmethod
    def _outcome_of(event):
        try:
            return event.outcome()
        except Exception:
            return None


@_contextmanager
def watching(solver, show_progress=True, **kwargs):
    """Report on the solve that runs inside this block.

    Yields the :class:`ProgressReport`, or None when asked for no progress.  Holds a
    reference to it for as long as it is attached: the solver keeps observers without
    owning them, so letting one be collected mid-solve would be a crash.
    """
    if not show_progress:
        yield None
        return

    report = ProgressReport(solver, **kwargs)
    solver.add_observer(report)
    try:
        yield report
    finally:
        try:
            solver.remove_observer(report)
        finally:
            report.close()
