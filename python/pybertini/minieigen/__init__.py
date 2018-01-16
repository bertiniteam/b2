"""
The `MiniEigen <https://github.com/eudoxos/minieigen>`_ submodule, for PyBertini

Bertini2 uses `Eigen <https://eigen.tuxfamily.org/>`_ for linear algebra, enabling expression templates on linear algebraic objects and operations for arbitrary types -- most importantly the multiprecision complex numbers needed for working near singularities in algebraic geometry.

This module exposes some functionality.  Not near all of Eigen is exposed.  We could use help on this.  If you know some Boost.Python and some Eigen, please consider helping expose more of Eigen through MiniEigen.  Bertini2 is not the host of MiniEigen -- user Eudoxus on GitHub provided MiniEigen.

If you encounter any problems with this functionality, please ask for help via a `GitHub issue <https://github.com/bertiniteam/b2/issues>`_.

"""


import _pybertini
import _pybertini.minieigen

from _pybertini.minieigen import *

__all__ = dir(_pybertini.minieigen)

