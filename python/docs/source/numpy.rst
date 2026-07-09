Multiprecision numbers and NumPy
********************************

Bertini 2 exposes :class:`~bertini.real_mp` (variable-precision real) and
:class:`~bertini.complex_mp` (variable-precision complex) as **native NumPy dtypes**
(via eigenpy).  A plain numpy array *is* the multiprecision vector -- no wrapper types,
no conversion step -- and the numpy you already know works on it:

* element-wise arithmetic, comparisons, ``@`` / :func:`numpy.dot` /
  :func:`numpy.linalg.norm`,
* the ufunc family: :func:`numpy.abs`, :func:`numpy.conj`, ``exp`` / ``log`` /
  the trigonometric and hyperbolic functions and their inverses, ``power``, ``sign``,
  ``sqrt``, the rounding family (``floor`` / ``ceil`` / ``trunc`` / ``rint``, with
  ``rint`` / ``round`` on complex rounding each component), ``minimum`` / ``maximum``,
  ``isnan`` / ``isinf`` / ``isfinite``, and friends,
* sorting and order statistics on the real type: :func:`numpy.sort`,
  :func:`numpy.argsort`, :func:`numpy.searchsorted`, :func:`numpy.argmax` /
  :func:`numpy.argmin`, :func:`numpy.median`,
* reductions: :func:`numpy.sum`, :func:`numpy.prod`, :func:`numpy.mean`,
  :func:`numpy.cumsum`, ``ufunc.reduce``.

The contract behind all of it: every ufunc loop calls the same multiprecision function
that the scalar :mod:`bertini.multiprec` function of the same name binds, at the
operands' precision.  ``np.exp(a)[i]`` is exactly ``mp.exp(a[i])`` -- numpy is a fast
way to say the same thing, never a detour through float64.

.. _numpy-float64-boundary:

Your digits are protected: the float64 boundary
================================================

A Python ``float`` is a 53-bit binary number: ``0.1`` as a float is *not* the decimal
``0.1`` you typed.  So Bertini draws a deliberate line -- a float64 **value** never
silently enters a multiprecision computation.  Values come from strings
(``real_mp('0.1')``) or integers (exact), and high-precision results only drop to
double when you explicitly ask (``float(x)``, ``complex(z)``, ``arr.astype(float)``,
``arr.astype(complex)``).

**Comparisons are a different story** -- checking a residual against a tolerance never
pollutes anything, so it just works, floats and all:

.. doctest::

    >>> import numpy as np
    >>> from bertini.multiprec import real_mp, complex_mp
    >>> a = np.array([complex_mp(3), complex_mp(4)])
    >>> b = np.array([complex_mp(3), complex_mp(4)])

    >>> # the closeness check: compare against a double tolerance, get bools back
    >>> bool(np.all(np.abs(a - b) < 1e-10))
    True

    >>> # the comparison is exact -- 1e-22 is not "rounded away" against 1e-30
    >>> bool(np.all(np.array([real_mp('1e-22')]) < 1e-30))
    False

This mirrors the C++ solvers, whose tolerances are doubles too.

Where the boundary shows as an error, that is it doing its job:

* ``a + 0.1`` raises (``ufunc ... not supported``) -- a float value would have entered
  the computation.  Say ``a + real_mp('0.1')``.
* ``a == 0.1`` raises -- exact equality against a float literal is the trap the
  boundary exists for.  Compare against ``real_mp('0.1')``, or use a tolerance.
* ``np.isclose`` / ``np.allclose`` raise ``DTypePromotionError`` -- they *compute*
  ``atol + rtol*np.abs(b)`` with float64 internally.  Use the one-liner above.
* ``np.round(a, decimals)`` with nonzero ``decimals`` fails (it scales by a float
  power of ten internally); plain ``np.round(a)`` / :func:`numpy.rint` work, and
  :func:`bertini.round` does decimal-digit rounding at full precision.

.. _numpy-complex-components:

Component access on complex arrays
===================================

The one place numpy itself cannot be taught about a user-defined dtype: the ndarray
attributes ``.real`` and ``.imag``.  numpy hardwires them to its own built-in complex
types -- on any other dtype ``.real`` returns the array itself and ``.imag`` returns
zeros, with no hook for the bindings to fix or even detect it.

Bertini defends every spelling it can reach:

* **Solution points are simply correct.**  The arrays returned by ``bertini.solve``
  override ``.real`` / ``.imag`` at the subclass level, so ``sol.real`` / ``sol.imag``
  give the true components at full precision.
* **The numpy functions fail loudly.**  Importing bertini wraps :func:`numpy.real` /
  :func:`numpy.imag` / :func:`numpy.angle`: on a plain mp-complex array they raise a
  ``TypeError`` naming the right tool, instead of silently returning wrong values.
  All other inputs pass straight through to numpy.

.. warning::

   The raw ``.real`` / ``.imag`` **attributes on a plain ndarray** are the one
   spelling nothing can reach: on a ``complex_mp`` array you built yourself (not a
   Solution), they return wrong values silently.  Complex *scalars* are fine --
   ``w[0].imag`` is exact.

The component accessors -- same results, full precision, array-capable:

.. doctest::

    >>> import numpy as np
    >>> import bertini.multiprec as mp
    >>> from bertini.multiprec import complex_mp
    >>> w = np.array([complex_mp(1, 2), complex_mp(3, 4)])

    >>> [float(x) for x in mp.real(w)]
    [1.0, 3.0]
    >>> [float(x) for x in mp.imag(w)]
    [2.0, 4.0]

    >>> # the np.angle equivalent
    >>> [round(float(x), 4) for x in mp.arg(w)]
    [1.1071, 0.9273]

Convenience helpers
====================

For everyday work with solution points, the top-level helpers accept a scalar, a list,
or a numpy array, and return multiprecision-native results:

.. code-block:: python

    import bertini

    bertini.real(pt)      # real parts, as real_mp
    bertini.imag(pt)      # imaginary parts, as real_mp
    bertini.abs(pt)       # magnitudes, as real_mp
    bertini.conj(pt)      # complex conjugates
    bertini.round(pt, 8)  # rounded to 8 DECIMAL digits, staying mp
    bertini.sum(pt)       # sum, staying mp
    bertini.norm(pt)      # Euclidean norm, as real_mp
    bertini.is_real(pt)   # is every coordinate real (|imag| < tol)?  -> bool

On multiprecision-dtype arrays these ride the native numpy loops; on lists and mixed
input they work element-wise.  The builtin-shadowing names (``abs``, ``round``,
``sum``) are deliberately kept out of ``from bertini import *``, so a star-import
never clobbers the Python builtins.  For the tolerance point comparison behind
de-duplication, see :func:`bertini.is_distinct_up_to`.

And for the whole vocabulary in one go -- the elementary functions *and* these
helpers, working on symbolic expressions, multiprecision numbers, and numpy
containers alike, dispatched per argument::

    from bertini.operators import *

    f = sin(x) + Pi*y            # symbolic (x, y Variables)
    v = sin(real_mp('0.5'))      # numeric, full precision
    m = abs(solutions[0])        # numpy containers, mp dtypes included
    t = arg(complex_mp(1, 1))    # components: arg/real/imag/conj

This star-import *does* shadow ``abs``/``round``/``sum`` in your namespace -- that is
its point (they fall back to builtin behavior on plain python input), and it is
opt-in.

Boundaries by design
=====================

* **Complex numbers are unordered.**  ``complex_mp`` arrays do not support
  :func:`numpy.sort`, :func:`numpy.argmax`, ``minimum`` / ``maximum``, or ``<``.
  (numpy's ``complex128`` sorts lexicographically for historical reasons; Bertini
  does not reproduce that.)  Sort a derived real quantity: ``np.argsort(np.abs(w))``.
* **No casts to integer types.**  ``arr.astype(int)`` fails; go through double
  explicitly if you truly want it.
