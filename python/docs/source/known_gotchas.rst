⚠️ Known gotchas
*****************************

A few sharp edges fall out of how Bertini 2's multiprecision numbers are exposed to
NumPy.  They are collected here with the idiom that works and the idiom that bites.

Bertini 2 exposes :class:`~bertini.real_mp` (variable-precision real) and
:class:`~bertini.complex_mp` (variable-precision complex) as **custom NumPy dtypes**
(via eigenpy), and nearly everything works on arrays of them:

* element-wise arithmetic, comparisons, ``@`` / :func:`numpy.dot` /
  :func:`numpy.linalg.norm`,
* the whole ufunc family: :func:`numpy.abs`, :func:`numpy.conj`, ``exp`` / ``log`` /
  the trigonometric and hyperbolic functions and their inverses, ``power``, ``sign``,
  ``sqrt``, the rounding family (``floor`` / ``ceil`` / ``trunc`` / ``rint``),
  ``minimum`` / ``maximum``, ``isnan`` / ``isinf`` / ``isfinite``, and friends,
* sorting and order statistics on the real type: :func:`numpy.sort`,
  :func:`numpy.argsort`, :func:`numpy.searchsorted`, :func:`numpy.argmax` /
  :func:`numpy.argmin`, :func:`numpy.median`,
* the identity-seeded reductions: :func:`numpy.sum`, :func:`numpy.prod`,
  :func:`numpy.mean`, :func:`numpy.cumsum`, bare ``ufunc.reduce``.

Every ufunc loop calls the same multiprecision function that the scalar
:mod:`bertini.multiprec` function of the same name binds, at the operands' precision --
``np.exp(a)[i]`` is exactly ``mp.exp(a[i])``.

The sharp edges that remain are below.  They are *by design* (the float64 boundary) or
*by NumPy limitation* (component access on a user dtype), not bugs to be waited out.

.. _gotcha-float64-boundary:

Mixing float64 into multiprecision arrays
==========================================

A Python ``float`` (or NumPy ``float64``) does **not** silently promote into a
multiprecision *value*.  This is deliberate: ``0.1`` as a float64 is not the number you
typed -- it carries binary noise past the 16th digit -- so anywhere a double could
sneak into a high-precision computation, Bertini makes you say what you mean
(construct from a *string*: ``real_mp('0.1')``).

The line is drawn at whether the float's *value* flows into the computation:

* **Tolerance comparisons against a float are allowed.**  ``np.abs(a - b) < 1e-10``
  works: an ordering comparison yields a ``bool``, so no float ever enters a
  multiprecision value -- and the comparison is *exact* (the mp value is compared
  against the double, not rounded down to double first).  This matches the C++
  solvers, whose tolerances are doubles.
* **Mixed arithmetic is blocked** (``a + 0.1`` raises ``ufunc ... not supported``),
  and so is **mixed equality** (``a == 0.1`` -- exact equality against a float
  literal is precisely the trap the boundary exists for).
* ``np.isclose(a, b)`` / ``np.allclose(a, b)`` raise ``DTypePromotionError``: they
  *compute* ``atol + rtol*np.abs(b)`` with float64 values internally, which is
  arithmetic, not comparison.
* ``np.round(a, decimals)`` with nonzero ``decimals`` fails for the same reason (it
  scales by a float power of ten internally).  Plain ``np.round(a)`` /
  :func:`numpy.rint` work.

The closeness idiom is therefore just:

.. doctest::

    >>> import numpy as np
    >>> from bertini.multiprec import real_mp, complex_mp
    >>> a = np.array([complex_mp(3), complex_mp(4)])
    >>> b = np.array([complex_mp(3), complex_mp(4)])

    >>> # ✓ the replacement for np.allclose(a, b): compare against a double tolerance
    >>> bool(np.all(np.abs(a - b) < 1e-10))
    True

    >>> # ✓ the comparison is exact -- 1e-22 is not lost against a 1e-30 tolerance
    >>> bool(np.all(np.array([real_mp('1e-22')]) < 1e-30))
    False

    >>> # ✗ np.allclose itself cannot work: it computes with float64 tolerances
    >>> np.allclose(a, b)
    Traceback (most recent call last):
        ...
    numpy.exceptions.DTypePromotionError: The DType <class 'numpy.dtype[complex_mp]'> could not be promoted by <class 'numpy.dtype[float64]'>. This means that no common DType exists for the given inputs. For example they cannot be stored in a single array unless the dtype is `object`. The full list of DTypes is: (<class 'numpy.dtype[complex_mp]'>, <class 'numpy.dtype[float64]'>)

Integers are fine in both directions of intent -- they are exact, so they convert and
promote freely (``a * 3``, ``np.power(a, 2)``, ``real_mp(7)``).  And converting mp
*down* to double is available when you ask for it explicitly (``float(x)``,
``complex(z)``, ``arr.astype(float)``, ``arr.astype(complex)``) -- you are consciously
truncating.

.. _gotcha-complex-components:

Component access on complex arrays: never ``.real`` / ``.imag``
================================================================

NumPy does not know a user-defined dtype is complex-like, so the ndarray attributes
``.real`` and ``.imag`` (and :func:`numpy.real`, :func:`numpy.imag`,
:func:`numpy.angle`, which route through them) return **silently wrong values** on
``complex_mp`` arrays: ``.real`` returns the complex values themselves and ``.imag``
returns zeros.  This is a NumPy limitation for legacy user dtypes; Bertini cannot hook
those attributes.

Use the :mod:`bertini.multiprec` component functions, which accept arrays:

.. doctest::

    >>> import numpy as np
    >>> import bertini.multiprec as mp
    >>> from bertini.multiprec import complex_mp
    >>> w = np.array([complex_mp(1, 2), complex_mp(3, 4)])

    >>> # ✓ arrays of the real/imaginary parts, as real_mp
    >>> [float(x) for x in mp.real(w)]
    [1.0, 3.0]
    >>> [float(x) for x in mp.imag(w)]
    [2.0, 4.0]

    >>> # ✓ the np.angle replacement
    >>> [round(float(x), 4) for x in mp.arg(w)]
    [1.1071, 0.9273]

    >>> # ✗ the ndarray attribute is wrong (numpy returns zeros -- not Bertini's doing)
    >>> complex(w.imag[0])
    0j

On complex *scalars* the properties are correct (``w[0].imag`` is exact); it is only
the *array* attributes that lie.

.. _gotcha-complex-ordering:

Complex is unordered
=====================

There is no ``<`` on complex numbers, so ``complex_mp`` arrays do not support
:func:`numpy.sort`, :func:`numpy.argmax`, :func:`numpy.minimum` / ``maximum``, or the
ordering comparisons.  (NumPy's built-in ``complex128`` sorts lexicographically for
historical reasons; Bertini deliberately does not reproduce that.)  Sort a derived real
quantity instead -- e.g. ``np.argsort(np.abs(w))``.

Casts that do not exist
========================

There is deliberately **no** cast from the multiprecision types to any integer type.
``arr.astype(int)`` on an mp array will fail; go through double first if you truly
want it, accepting the truncation.

.. _gotcha-numpy-reductions:

A historical note on reductions
================================

Earlier versions of this page declared :func:`numpy.sum` / :func:`numpy.prod` /
:func:`numpy.mean` unsupported on mp arrays: on some older NumPy builds the
identity-seeded reduce failed inside NumPy with ``SystemError: ... returned NULL
without setting an exception``.  With current NumPy (verified on 2.3 and 2.4 series;
regression-tested in CI on all platforms) they work:

.. doctest::

    >>> import numpy as np
    >>> from bertini.multiprec import real_mp, complex_mp
    >>> v = np.array([complex_mp(3), complex_mp(4)])
    >>> w = np.array([real_mp(1), real_mp(2), real_mp(3)])

    >>> complex(np.sum(v))
    (7+0j)
    >>> float(np.mean(w))
    2.0

If you are pinned to an old NumPy and see that ``SystemError``, the old idioms all
still work and remain the portable fallback:

.. doctest::

    >>> # explicit, correctly-typed identity
    >>> complex(np.add.reduce(v, initial=complex_mp(0)))
    (7+0j)
    >>> # python's own sum
    >>> float(sum(w))
    6.0
    >>> # norm / sum of squares route through the dtype's dot slot
    >>> float(np.linalg.norm(v))
    5.0
    >>> complex(np.dot(v, v))
    (25+0j)
