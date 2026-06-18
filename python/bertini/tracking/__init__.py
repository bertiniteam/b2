# This file is part of Bertini 2.
# 
# python/bertini/tracking/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/bertini/tracking/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/bertini/tracking/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Copyright(C) Bertini2 Development Team
# 
#  See <http://www.gnu.org/licenses/> for a copy of the license, 
#  as well as COPYING.  Bertini2 is provided with permitted 
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
# 
#  silviana amethyst
#  UWEC
#  Spring 2018
# 





"""
Tracking-specific things -- trackers, configs
"""
from bertini._pybertini import tracking as _pybtracking
from bertini._pybertini.tracking import *

# trackers gain configure()/config_names() built on their type-list config interface;
# configs gain update()/to_dict()/from_dict()/repr/eq.
from ..config import _enhance_owners, _enhance_all
_enhance_owners(_pybtracking)
_enhance_all(_pybtracking)

__all__ = dir(_pybtracking)


# What an observer may return from Observe() to control its own subscription.
# Returning None (the usual case) means KeepObserving.
from bertini._pybertini.detail import ObserveResult

_obs_amp = _pybtracking.observers.amp
_obs_dbl = _pybtracking.observers.double
_obs_mul = _pybtracking.observers.multiple

AMPTracker.observers = _obs_amp
DoublePrecisionTracker.observers = _obs_dbl
MultiplePrecisionTracker.observers = _obs_mul


def _make_callback_observer(AbstractClass):
    class CallbackObserver(AbstractClass):
        """Observer that routes events to registered Python callables.

        Usage::

            obs = CallbackObserver()
            obs.on(observers.amp.TrackingStarted, lambda e: print("tracking started"))
            obs.on(observers.amp.PrecisionChanged,
                   lambda e: print(e.previous(), "->", e.next()))
            tracker.add_observer(obs)
            tracker.track_path(...)
            tracker.remove_observer(obs)
        """
        def __init__(self):
            super().__init__()
            self._callbacks = {}

        def on(self, event_type, callback):
            """Register *callback* to be called when an event of *event_type* arrives.

            Returns self for chaining.
            """
            self._callbacks.setdefault(event_type, []).append(callback)
            return self

        def Observe(self, event):
            for event_type, callbacks in self._callbacks.items():
                if isinstance(event, event_type):
                    for cb in callbacks:
                        cb(event)

    return CallbackObserver


_obs_amp.CallbackObserver = _make_callback_observer(_obs_amp.CustomObserver)
_obs_dbl.CallbackObserver = _make_callback_observer(_obs_dbl.CustomObserver)
_obs_mul.CallbackObserver = _make_callback_observer(_obs_mul.CustomObserver)


def _make_path_observers(obs_mod):
    """Build the path-collecting observers for one precision's observer module.

    Returns ``(PathDataCollector, PathCollectionObserver)``.  ``obs_mod`` is e.g.
    ``bertini.tracking.observers.amp`` and supplies ``.CustomObserver`` plus the event
    classes ``SuccessfulStep`` / ``TrackingStarted`` / ``TrackingEnded``.
    """

    class PathDataCollector(obs_mod.CustomObserver):
        """Collects one tracked path into a time series, for plotting.

        Attach to a tracker; on every successful step it records the time, the
        space point, and a few step diagnostics.  Because an adaptive-precision
        tracker hands back arbitrary-precision (mpfr) numbers -- which do not pack
        into a single numpy array -- every value is cast to a plain python
        ``complex``/``float`` as it is collected (double precision is plenty for a
        picture).  The collected data is offered as several typed numpy arrays
        (or, optionally, a pandas DataFrame)::

            b = PathDataCollector()
            tracker.add_observer(b)
            tracker.track_path(...)
            tracker.remove_observer(b)
            t   = b.times()          # complex, shape (n_steps,)
            z   = b.points()         # complex, shape (n_steps, n_vars)
            dgn = b.diagnostics()    # float,   shape (n_steps, 4): |t|, cond, prec, stepsize
        """

        #: column order of the float array returned by :meth:`diagnostics`.
        DIAGNOSTIC_COLUMNS = ("abs_t", "condition_number", "precision", "stepsize")

        def __init__(self):
            super().__init__()
            #: the time this track started at, if known (set by PathCollectionObserver
            #: from the TrackingStarted event).  Lets you tell a main homotopy path
            #: (starts at the global start time) from an endgame sub-track (starts
            #: near the endgame boundary), since a solver reuses one tracker for both.
            self.start_time = None
            self._t = []        # list[complex]
            self._points = []   # list[list[complex]]
            self._diag = []     # list[[abs_t, cond, prec, stepsize]]

        def Observe(self, event):
            if not isinstance(event, obs_mod.SuccessfulStep):
                return
            trk = event.tracker()
            tval = complex(trk.current_time())
            self._t.append(tval)
            self._points.append([complex(z) for z in trk.current_point()])
            self._diag.append([
                abs(tval),
                float(trk.latest_condition_number()),
                float(trk.current_precision()),
                float(trk.current_stepsize()),
            ])

        def __len__(self):
            return len(self._t)

        def times(self):
            import numpy as np
            return np.array(self._t, dtype=complex)

        def points(self):
            import numpy as np
            return np.array(self._points, dtype=complex)

        def diagnostics(self):
            import numpy as np
            return np.array(self._diag, dtype=float).reshape(-1, len(self.DIAGNOSTIC_COLUMNS))

        def as_dataframe(self):
            """Return the whole path as a pandas DataFrame with named columns.

            Requires pandas (``pip install pandas``).  Columns: ``t`` (complex),
            one ``z{i}`` (complex) per variable, then the diagnostic columns.
            """
            try:
                import pandas as pd
            except ImportError as e:  # pragma: no cover - exercised only without pandas
                raise ImportError(
                    "PathDataCollector.as_dataframe() needs pandas; "
                    "install it, or use times()/points()/diagnostics() instead."
                ) from e
            import numpy as np
            data = {"t": self.times()}
            pts = self.points()
            for i in range(pts.shape[1] if pts.size else 0):
                data["z{}".format(i)] = pts[:, i]
            diag = self.diagnostics()
            for j, name in enumerate(self.DIAGNOSTIC_COLUMNS):
                data[name] = diag[:, j] if diag.size else np.empty(0)
            return pd.DataFrame(data)

    class PathCollectionObserver(obs_mod.CustomObserver):
        """A meta-observer: collects *every* path of a multi-path run.

        Attach one of these to a tracker (e.g. ``zd.get_tracker()``) before a
        solve.  Each time the tracker starts a path it spins up a fresh
        :class:`PathDataCollector`, attaches it, and -- when the path ends --
        harvests it into :attr:`series` and detaches it.  Attaching/detaching
        happens from inside ``Observe`` and relies on the observable deferring
        those mutations until the current notification finishes.

        After the run, :attr:`series` is a list of finished ``PathDataCollector``
        objects, one per path, in the order the tracker ran them::

            a = PathCollectionObserver()
            zd.get_tracker().add_observer(a)
            zd.solve()
            for path in a.series:
                t, z = path.times(), path.points()
                ...
        """

        def __init__(self):
            super().__init__()
            self.series = []      # list[PathDataCollector], one finished path each
            self._active = None
            self._tracker = None

        def Observe(self, event):
            if isinstance(event, obs_mod.TrackingStarted):
                trk = event.tracker()
                collector = PathDataCollector()
                collector.start_time = complex(trk.current_time())
                self._active = collector
                self._tracker = trk
                trk.add_observer(collector)        # deferred; starts on the next event
            elif isinstance(event, obs_mod.TrackingEnded):
                if self._active is not None:
                    self._tracker.remove_observer(self._active)   # deferred
                    self.series.append(self._active)
                    self._active = None
                    self._tracker = None

        def __len__(self):
            return len(self.series)

    return PathDataCollector, PathCollectionObserver


for _m in (_obs_amp, _obs_dbl, _obs_mul):
    _m.PathDataCollector, _m.PathCollectionObserver = _make_path_observers(_m)