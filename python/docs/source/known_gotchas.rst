⚠️ Known gotchas
*****************************

A few sharp edges fall out of how Bertini 2's multiprecision numbers are exposed to
NumPy.  They are collected here with the idiom that works and the idiom that bites.

.. _gotcha-numpy-reductions:

NumPy reductions over multiprecision arrays
============================================

Bertini 2 exposes :class:`~bertini.real_mp` (variable-precision real) and
:class:`~bertini.complex_mp` (variable-precision complex) as **custom NumPy
dtypes** (via eigenpy).  Element-wise math, indexing, ``@`` / :func:`numpy.dot`,
:func:`numpy.linalg.norm`, :func:`numpy.cumsum` and friends all work on arrays of these
dtypes.

The sharp edge is the **identity-seeded reductions** -- :func:`numpy.sum`,
:func:`numpy.prod`, :func:`numpy.mean`, and a bare ``ufunc.reduce`` with no ``initial=``:

.. code-block:: python

    >>> import numpy as np
    >>> from bertini.multiprec import complex_mp
    >>> v = np.array([complex_mp(3), complex_mp(4)], dtype=complex_mp)
    >>> np.sum(v)                                    # ✗ DON'T -- may crash
    Traceback (most recent call last):
        ...
    SystemError: <method 'reduce' of 'numpy.ufunc' objects> returned NULL without setting an exception

Why this happens
----------------

To reduce an array, NumPy needs the reduction's **identity element** -- ``0`` for a sum,
``1`` for a product -- materialized *in the array's dtype*.  For these custom dtypes,
NumPy (2.x) tries to build that identity from a Python ``int`` and fails *inside its own
machinery*, before ever calling Bertini's dtype code, returning ``NULL`` without setting a
Python exception (hence the bare ``SystemError``).  It is a limitation of the NumPy ↔
eigenpy boundary for legacy user-defined dtypes, **not** something Bertini can patch in its
bindings, and whether it triggers depends on the exact NumPy / eigenpy / Eigen build -- so
it may "work" on one machine and crash on another.  Treat it as always-unsupported.

What to do instead
------------------

Give the reduction a value to start from, or use an operation that seeds itself from the
data (``dot`` / ``norm`` / a plain Python loop).  All of these are stable across builds:

.. doctest::

    >>> import numpy as np
    >>> from bertini.multiprec import real_mp, complex_mp
    >>> v = np.array([complex_mp(3), complex_mp(4)], dtype=complex_mp)
    >>> w = np.array([real_mp(1), real_mp(2), real_mp(3)], dtype=real_mp)

    >>> # ✓ sum: hand ufunc.reduce an explicit, correctly-typed identity
    >>> complex(np.add.reduce(v, initial=complex_mp(0)))
    (7+0j)
    >>> float(np.add.reduce(w, initial=real_mp(0)))
    6.0

    >>> # ✓ sum: or just use Python's built-in sum()
    >>> float(sum(w))
    6.0

    >>> # ✓ Euclidean norm works directly (it routes through the dtype's dot slot)
    >>> float(np.linalg.norm(v))
    5.0

    >>> # ✓ sum of squares: np.dot does not conjugate, so dot(v, v) == sum(v**2)
    >>> complex(np.dot(v, v))
    (25+0j)

    >>> # ✓ mean: reduce with an identity, then divide by the count
    >>> float(np.add.reduce(w, initial=real_mp(0)) / w.size)
    2.0

In short: anywhere you would reach for ``np.sum(a)`` / ``np.prod(a)`` / ``np.mean(a)`` on a
``real_mp`` or ``complex_mp`` array, reach for ``np.add.reduce(a, initial=...)`` (or
``np.multiply.reduce(a, initial=...)``), ``np.dot`` / ``np.linalg.norm``, or a plain Python
``sum`` instead.

The bertini helpers do this for you
-----------------------------------

The same boundary makes ``arr.real`` / ``np.abs(arr)`` / ``np.round(arr)`` unreliable on these
dtypes.  So bertini ships elementwise helpers that operate over a scalar, a list, or a numpy array and
return **mp-native** results (never a float64 collapse):

.. code-block:: python

    import bertini

    bertini.real(pt)      # real parts, as real_mp        (replaces arr.real)
    bertini.imag(pt)      # imaginary parts, as real_mp
    bertini.abs(pt)       # magnitudes, as real_mp         (replaces np.abs)
    bertini.conj(pt)      # complex conjugates
    bertini.round(pt, 8)  # rounded, staying mp            (replaces np.round)
    bertini.sum(pt)       # sum, staying mp                (sidesteps the reduction gotcha above)
    bertini.norm(pt)      # Euclidean norm, as real_mp
    bertini.is_real(pt)   # is every coordinate real (|imag| < tol)?  -> bool

They live at the top level (``bertini.abs``, ...); the builtin-shadowing names (``abs``, ``round``,
``sum``) are deliberately kept out of ``from bertini import *``, so a star-import never clobbers the
Python builtins.  For the tolerance point comparison behind de-duplication, see
:func:`bertini.is_distinct_up_to`.
