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
from ..config import enhance_owners, enhance_all
enhance_owners(_pybtracking)
enhance_all(_pybtracking)

__all__ = dir(_pybtracking)


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


_obs_amp.CallbackObserver = _make_callback_observer(_obs_amp.Abstract)
_obs_dbl.CallbackObserver = _make_callback_observer(_obs_dbl.Abstract)
_obs_mul.CallbackObserver = _make_callback_observer(_obs_mul.Abstract)