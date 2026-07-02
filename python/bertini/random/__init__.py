import numpy as _np

from bertini._pybertini import random as _pybrand
from bertini._pybertini.random import *


def random_matrix(rows, cols, real=False, units=False, orthonormal=True, symbolic=False):
    """A ``rows`` x ``cols`` random matrix, at the current default precision.

    The building block for a random linear form / projection / slice / patch -- they are all just
    linear functions, and "what you use them for is up to you."  Returns a numpy array of
    arbitrary-precision values (``bertini.multiprec.Complex``); pass ``symbolic=True`` to instead
    get coefficient *nodes* ready to drop straight into function-tree expressions.

    Parameters
    ----------
    rows, cols : int
        The shape of the matrix.
    real : bool, default False
        When True the entries are real (zero imaginary part); otherwise complex.
    units : bool, default False
        When True (and not ``orthonormal``) each entry has modulus 1 (a ``real`` unit is +/-1);
        otherwise each entry is drawn with bounded modulus (box-uniform in [-1,1] per component).
    orthonormal : bool, default True
        When True the rows are conjugate-orthonormal (QR-factored from a square matrix of units,
        then truncated -- perfectly conditioned).  ``real=True`` gives a real orthogonal matrix.
        ``orthonormal`` takes precedence over ``units`` (orthonormal rows are already normalized).
    symbolic : bool, default False
        When True, return a numpy object array of coefficient *nodes*
        (via :func:`bertini.coefficients`); otherwise numeric ``multiprec.complex_mp``.

    Notes
    -----
    Reproducible via :func:`bertini.random.set_random_seed`.  A projection is just a linear
    functional with zero constant term, so its gradient row is exactly ``random_matrix(1, n)``.
    """
    if orthonormal:
        M = _np.array(_pybrand.conjugate_orthonormal_matrix(rows, cols, real))
    else:
        if units:
            draw = _pybrand.real_unit if real else _pybrand.complex_unit
        else:
            draw = _pybrand.real_bounded_modulus if real else _pybrand.complex_bounded_modulus
        M = _np.array([[draw() for _ in range(cols)] for _ in range(rows)], dtype=object)

    # eigenpy collapses a single-row matrix to 1-D; keep the (rows, cols) shape predictable.
    M = M.reshape(rows, cols)

    if symbolic:
        from bertini._coefficients import coefficients
        return coefficients(M)
    return M


__all__ = dir(_pybrand) + ['random_matrix']
