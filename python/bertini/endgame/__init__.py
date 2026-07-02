# This file is part of Bertini 2.
# 
# python/bertini/endgame/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/bertini/endgame/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/bertini/endgame/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
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
Endgame-specific things -- endgames, configs
***********************************************

Endgames allow the computation of singular endpoints.  

Flavors
=========

There are two basic flavors of endgame implemented:

1. Power Series, commonly written PS or PSEG
2. Cauchy

Both estimate the cycle number and use it to compute a root at a time which is never tracked to.  PSEG uses Hermite interpolation and extrapolation, and Cauchy uses loops around the target time coupled with the `Cauchy integral formula <https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula>`_.  Both continue until two successive approximations of the root match to a given tolerance (:py:attr:`bertini.endgame.EndgameConfig.final_tolerance`).

The implementations of the endgames go with a particular tracker, hence there are six provided endgame types.  Choose the one that goes with your selected tracker type.  Adaptive Multiple Precision is a good choice.


Implementation reference
=========================

AMP Endgames
-------------

* :class:`~bertini.endgame.AMPCauchyEndgame`
* :class:`~bertini.endgame.AMPPowerSeriesEndgame`

Fixed Double Precision Endgames
---------------------------------

* :class:`~bertini.endgame.FixedDoubleCauchyEndgame`
* :class:`~bertini.endgame.FixedDoublePowerSeriesEndgame`

Fixed Multiple Precision  Endgames
-------------------------------------

* :class:`~bertini.endgame.FixedMultipleCauchyEndgame`
* :class:`~bertini.endgame.FixedMultiplePowerSeriesEndgame`

"""

import bertini._pybertini.endgame as _pybe

from bertini._pybertini.endgame import observers
from bertini._pybertini.endgame import (
    CauchyConfig,
    EndgameConfig,
    PowerSeriesConfig,
    SecurityConfig,
)

# configs gain update()/to_dict()/from_dict()/repr/eq
from ..config import _enhance_all
_enhance_all(_pybe)


class _EndgameBase:
    """Thin wrapper that stores the endgame boundary time at construction.

    The C++ Run() re-precisions the stored boundary time to match whatever
    precision the incoming point is at, so callers never need to manage that.
    """

    def __init__(self, cpp_class, tracker, boundary_time):
        object.__setattr__(self, '_eg', cpp_class(tracker))
        self._eg.set_boundary_time(boundary_time)

    def run(self, point):
        return self._eg.run(point)

    def __getattr__(self, name):
        if name == '_eg':
            raise AttributeError('_eg not initialized')
        return getattr(self._eg, name)


class AMPCauchyEndgame(_EndgameBase):
    def __init__(self, tracker, boundary_time):
        super().__init__(_pybe.AMPCauchyEndgame, tracker, boundary_time)

class AMPPowerSeriesEndgame(_EndgameBase):
    def __init__(self, tracker, boundary_time):
        super().__init__(_pybe.AMPPowerSeriesEndgame, tracker, boundary_time)

class FixedDoubleCauchyEndgame(_EndgameBase):
    def __init__(self, tracker, boundary_time):
        super().__init__(_pybe.FixedDoubleCauchyEndgame, tracker, boundary_time)

class FixedDoublePowerSeriesEndgame(_EndgameBase):
    def __init__(self, tracker, boundary_time):
        super().__init__(_pybe.FixedDoublePowerSeriesEndgame, tracker, boundary_time)

class FixedMultipleCauchyEndgame(_EndgameBase):
    def __init__(self, tracker, boundary_time):
        super().__init__(_pybe.FixedMultipleCauchyEndgame, tracker, boundary_time)

class FixedMultiplePowerSeriesEndgame(_EndgameBase):
    def __init__(self, tracker, boundary_time):
        super().__init__(_pybe.FixedMultiplePowerSeriesEndgame, tracker, boundary_time)


__all__ = [
    'AMPCauchyEndgame',
    'AMPPowerSeriesEndgame',
    'CauchyConfig',
    'EndgameConfig',
    'FixedDoubleCauchyEndgame',
    'FixedDoublePowerSeriesEndgame',
    'FixedMultipleCauchyEndgame',
    'FixedMultiplePowerSeriesEndgame',
    'PowerSeriesConfig',
    'SecurityConfig',
    'observers',
]
