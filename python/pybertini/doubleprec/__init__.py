"""
Double precision types, beyond those built in.


"""

import _pybertini.minieigen

Vector = _pybertini.minieigen.VectorXd
Matrix = _pybertini.minieigen.MatrixXd

__all__ = ['Vector','Matrix']
