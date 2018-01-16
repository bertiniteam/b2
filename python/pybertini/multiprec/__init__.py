"""
Multiprecision types


"""

import _pybertini.mpfr

from _pybertini.mpfr import *

import _pybertini.minieigen

Vector = _pybertini.minieigen.VectorXmp
Matrix = _pybertini.minieigen.MatrixXmp


__all__ = dir(_pybertini.mpfr)
__all__.extend(['Vector','Matrix'])
