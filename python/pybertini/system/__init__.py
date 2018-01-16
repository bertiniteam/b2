


"""
Provides utilities for working with systems of functions -- polynomials are intended, although you can work with functions involving things like trig functions, arbitrary powers, etc.

Making a new `System` is the starting point you want, probably:

::

	sys = pybertini.system.System()

-----------

There are also things available in the `start_system` submodule.

🏞

"""

import _pybertini

from _pybertini import System # brings the type System
from _pybertini import start_system

__all__ = ['System','start_system']
