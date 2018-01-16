


"""
Provides utilities for working with systems of functions -- polynomials are intended, although you can work with functions involving things like trig functions, arbitrary powers, etc.

Making a new `System` is the starting point you want, probably:

::

	sys = pybertini.system.System()

-----------

There are also things available in the `start_system` submodule.

🏞

"""

import _pybertini.system

from _pybertini.system import * # brings the type System
from _pybertini.system import start_system

__all__ = dir(_pybertini.system)
__all__.extend(['start_system'])
