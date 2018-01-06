import _pybertini

from _pybertini import mpfr as multiprec

multiprec.Vector = _pybertini.VectorXmp
multiprec.Matrix = _pybertini.MatrixXmp
