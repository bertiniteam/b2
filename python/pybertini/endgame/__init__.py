"""
Endgame-specific things -- endgames, configs
***********************************************

Endgames allow the computation of singular endpoints.  

Flavors
=========

There are two basic flavors of endgame implemented:

1. Power Series, commonly written PS or PSEG
2. Cauchy

Both estimate the cycle number and use it to compute a root at a time which is never tracked to.  PSEG uses Hermite interpolation and extrapolation, and Cauchy uses loops around the target time coupled with the `Cauchy integral formula <https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula>`_.  Both continue until two successive approximations of the root match to a given tolerance (:py:attr:`pybertini.endgame.config.Endgame.final_tolerance`).

The implementations of the endgames go with a particular tracker, hence there are six provided endgame types.  Choose the one that goes with your selected tracker type.  Adaptive Multiple Precision is a good choice.


Implementation reference
=========================

AMP Endgames
-------------

* :class:`~pybertini.endgame.AMPCauchyEG`
* :class:`~pybertini.endgame.AMPPSEG`

Fixed Double Precision Endgames
---------------------------------

* :class:`~pybertini.endgame.FixedDoublePSEG`
* :class:`~pybertini.endgame.FixedDoublePSEG`

Fixed Multiple Precision  Endgames
-------------------------------------

* :class:`~pybertini.endgame.FixedMultiplePSEG`
* :class:`~pybertini.endgame.FixedMultiplePSEG`

"""

import _pybertini.endgame

from _pybertini.endgame import *

__all__ = ['AMPCauchyEG',
 'AMPPSEG',
 'FixedDoubleCauchyEG',
 'FixedDoublePSEG',
 'FixedMultipleCauchyEG',
 'FixedMultiplePSEG',
 '__doc__',
 '__loader__',
 '__name__',
 '__package__',
 '__spec__',
 'config',
 'observers']
