"""
Mutliprecision types


"""

import _pybertini.mpfr

from _pybertini.mpfr import *

Vector = _pybertini.VectorXmp
Matrix = _pybertini.MatrixXmp


__all__ = dir(_pybertini.mpfr)
__all__.extend(['Vector','Matrix'])
