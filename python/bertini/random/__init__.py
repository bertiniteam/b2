import numpy as _np

from bertini._pybertini import random as _pybrand
from bertini._pybertini.random import *


def random_real(symbolic=False):
    """A random real number of bounded modulus, at the current default precision.

    Box-uniform in ``[-1, 1]`` and pulled away from 0 -- the Bertini genericity draw.
    Reproducible via :func:`set_random_seed`.

    Parameters
    ----------
    symbolic : bool, default False
        When False (default) return a numeric ``multiprec.real_mp`` value.  When True return a
        constant function-tree *node* (a real :class:`~bertini.symbolics.Complex` leaf, imaginary
        part 0) ready to drop straight into an expression -- the spelling to reach for when you
        want a random coefficient rather than a random value.  ``bertini.symbolics.random_real``
        is the always-symbolic shortcut.

    Notes
    -----
    The node carries exactly the digits drawn, i.e. the current default precision -- set
    ``bertini.default_precision(1000)`` first for a 1000-digit constant.  (A random float can be no
    more precise than its draw; for a *precision-independent* exact random constant use
    ``bertini.symbolics.Rational.rand_real`` instead.)
    """
    v = _pybrand.random_real()
    if symbolic:
        from bertini._coefficients import coefficient
        return coefficient(v)
    return v


def random_complex(symbolic=False):
    """A random complex number of bounded modulus, at the current default precision.

    Modulus pulled toward 1 (away from 0 and infinity) -- the Bertini genericity draw.
    Reproducible via :func:`set_random_seed`.

    Parameters
    ----------
    symbolic : bool, default False
        When False (default) return a numeric ``multiprec.complex_mp`` value.  When True return a
        constant function-tree *node* (a :class:`~bertini.symbolics.Complex` leaf) ready to drop
        straight into an expression.  ``bertini.symbolics.random_complex`` is the always-symbolic
        shortcut.

    Notes
    -----
    The node carries exactly the digits drawn, i.e. the current default precision -- set
    ``bertini.default_precision(1000)`` first for a 1000-digit constant.  For a
    *precision-independent* exact random constant use ``bertini.symbolics.Rational.rand`` instead.
    """
    v = _pybrand.random_complex()
    if symbolic:
        from bertini._coefficients import coefficient
        return coefficient(v)
    return v


def random_vector(size, real=False, symbolic=False):
    """A random length-``size`` vector of bounded-modulus numbers, at the current default precision.

    The natural random projection / linear-functional coefficient vector (generic and
    seed-reproducible, unlike the quantized orthonormal :func:`random_matrix`).

    Parameters
    ----------
    size : int
        The length of the vector.
    real : bool, default False
        When True the entries are real (``real_mp``); otherwise complex (``complex_mp``).
    symbolic : bool, default False
        When True return a numpy object array of constant coefficient *nodes* (via
        :func:`bertini.coefficients`) rather than numeric multiprecision values -- each entry a
        :class:`~bertini.symbolics.Complex` leaf.

    Notes
    -----
    As for :func:`random_real`, a symbolic entry carries the current default precision's worth of
    digits; set ``bertini.default_precision`` first if you need more.
    """
    v = _pybrand.random_vector(size, real)
    if symbolic:
        from bertini._coefficients import coefficients
        return coefficients(v)
    return v


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

    For a **generic real** direction (e.g. a real projection), prefer :func:`bertini.random_vector`
    with ``real=True``: a real *orthonormal* matrix is QR-factored from a matrix of real units (+/-1),
    so its entries are **quantized** (a real ``random_matrix(3, 1)`` has entries +/-1/sqrt(3)) and look
    seed-independent -- fine for conditioning, but not a generic direction.  ``random_vector`` draws
    continuous bounded-modulus reals instead (generic and seed-reproducible).
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
