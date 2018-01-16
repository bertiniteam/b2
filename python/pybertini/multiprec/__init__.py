"""
Multiprecision types


"""

import _pybertini.multiprec

from _pybertini.multiprec import *

import _pybertini.minieigen

Vector = _pybertini.minieigen.VectorXmp
Matrix = _pybertini.minieigen.MatrixXmp


__all__ = dir(_pybertini.multiprec)
__all__.extend(['Vector','Matrix'])
