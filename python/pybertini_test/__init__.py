"""
PyBertini -- Python bindings for Bertini2.

This code is licensed under the GNU Public License, Version 3.
"""


### this __init__.py is strongly inspired by that for GalSim 
### https://github.com/GalSim-developers/GalSim

from ._version import __version__, __version_info__


version = __version__


# put stuff in the pybertini namespace

from .system import System
