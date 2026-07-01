🪢 Author your own start system: a product of linears
*********************************************************

.. testsetup:: *

   import bertini

Every other tutorial lets bertini *generate* the start system for you -- a total-degree start
system, or the multihomogeneous one -- or (in :ref:`the parameter-homotopy tutorial
<tutorials>`) reuses the solutions of an earlier solve.  This one is different: **you** author
the start system, by hand, and drive the homotopy yourself.

The vehicle is a **product of linear forms**.  A linear form is :math:`c\cdot[x;1]` (the trailing
``1`` carries the constant term), and a start function is a product of them:

.. math::

   s_i(x) \;=\; \prod_r \bigl( c_{i,r} \cdot [x;1] \bigr).

Each factor :math:`c_{i,r}\cdot[x;1] = 0` is a **hyperplane**, so a start solution is just an
intersection of one hyperplane per function -- something you can write down and check by eye.  No
opaque generated coefficients: you see exactly why the start points are what they are.  The
machinery underneath is the same first-class evaluation block bertini's own multihomogeneous
start system uses, exposed through :mod:`bertini.linalg`.

A target you can check by hand
==============================

Take a unit circle meeting a parabola -- two quadratics in two variables, so Bézout says four
solutions:

.. testcode::

    import numpy as np
    import bertini
    from bertini import linalg, nag_algorithm
    from bertini import multiprec

    x, y = bertini.Variable('x'), bertini.Variable('y')

    target = bertini.System()
    target.add_variable_group(bertini.VariableGroup([x, y]))
    target.add_function(x*x + y*y - 1)        # the unit circle
    target.add_function(y - x*x)              # the parabola y = x^2

Substituting :math:`y = x^2` into :math:`x^2 + y^2 - 1` gives :math:`y^2 + y - 1 = 0`, so
:math:`y = \tfrac{-1 \pm \sqrt 5}{2}`.  The larger root is positive and yields a **real** pair
:math:`(\pm\sqrt{y}, y)`; the smaller is negative and yields a **purely imaginary** :math:`x`
pair.  Keep those four exact answers in your pocket -- we will check against them at the end.

A start system you can write down
=================================

For a target of degrees :math:`(2, 2)` we need a start system of the same degrees with solutions
we already know.  Make each start function a product of **two** linear forms, chosen so the
factors are coordinate-aligned:

.. testcode::

    start = bertini.System()
    start.add_variable_group(bertini.VariableGroup([x, y]))
    linalg.add_products_of_linears(start, [
        [[1, 0, '-1'], [1, 0, '1']],     # s0 = (x - 1)(x + 1)
        [[0, 1, '-1'], [0, 1, '-2']],    # s1 = (y - 1)(y - 2)
    ])

    assert list(start.degrees()) == [2, 2]   # a product's degree is its number of factors

Each entry of the list is one function's coefficient matrix: **one row per linear factor**, the
trailing column being that factor's constant term.  So ``[[1, 0, '-1'], [1, 0, '1']]`` is
:math:`(1\,x + 0\,y - 1)(1\,x + 0\,y + 1) = (x-1)(x+1)`.

.. note::

   **Coefficients must be exact.**  :func:`~bertini.linalg.add_products_of_linears` refuses
   Python floats: a 64-bit literal carries only ~16 digits and would cap the precision of every
   downstream computation.  Pass ints, :class:`fractions.Fraction`, exact strings (``'-1'``,
   ``'3/4'``), or :mod:`bertini.multiprec` values.

   This is a genuine *product* (degree = number of factors), the first-class C++
   ``ProductsOfLinearsBlock``.  It is **not** the same as
   :func:`~bertini.linalg.add_linear_forms`, which adds a *stack* of degree-1 linear forms.

Start points are intersections of hyperplanes
==============================================

Because :math:`s_0` vanishes when :math:`x = \pm 1` and :math:`s_1` vanishes when
:math:`y \in \{1, 2\}`, a start solution picks one factor (one hyperplane) from each function and
solves the resulting linear system.  Here that is just the grid :math:`x \in \{1, -1\}` times
:math:`y \in \{1, 2\}` -- four points, written down by hand:

.. testcode::

    import itertools
    start_points = [np.array([multiprec.Complex(str(a)), multiprec.Complex(str(b))])
                    for a, b in itertools.product([1, -1], [1, 2])]
    # (1, 1), (1, 2), (-1, 1), (-1, 2)

There are :math:`2 \times 2 = 4` of them -- the Bézout number of the start system, as it must be.

Blend into a homotopy and solve
===============================

Now couple the start system to the target with the gamma-trick straight-line homotopy
:math:`H = (1-t)\,\text{target} + \gamma\,t\,\text{start}`.  Because the start system carries a
structured evaluation block (the product of linears), it cannot be fused by ordinary node
arithmetic; :func:`~bertini.nag_algorithm.blend_homotopy` combines the two whole systems with a
blend block instead:

.. testcode::

    gamma = linalg.coefficient(multiprec.Complex('0.6', '0.8'))   # exact, off the real axis
    H = nag_algorithm.blend_homotopy(target, start, gamma=gamma)

    solver = nag_algorithm.HomotopySolver(H, start_points, target)
    solver.solve()
    solutions = solver.all_solutions()

At :math:`t=1` the homotopy is :math:`\gamma\,\text{start}`, so our four points are its roots; at
:math:`t=0` it is the target.  :func:`~bertini.nag_algorithm.HomotopySolver` runs the full zero-dim
pipeline -- pre-endgame tracking, the midpath check, the endgame, post-processing -- on the
homotopy and start points you supplied, rather than ones it generated.

Choosing :math:`\gamma` off the real axis makes the straight-line path miss the (measure-zero)
singular locus where solutions collide.  We pin it to the exact value :math:`0.6 + 0.8i` so the
example is reproducible; a random complex :math:`\gamma` (omit the argument) works just as well.

Check against the known answers
===============================

Guard against an empty result, then measure the **distance** from each known root to the nearest
computed solution, with an infinity-norm so the check is faithful to scale rather than a
scale-naive residual:

.. testcode::

    assert len(solutions) == 4

    import math
    s5 = math.sqrt(5)
    y1, y2 = (-1 + s5) / 2, (-1 - s5) / 2
    known = [np.array([math.sqrt(y1), y1]),  np.array([-math.sqrt(y1), y1]),
             np.array([1j*math.sqrt(-y2), y2]), np.array([-1j*math.sqrt(-y2), y2])]

    computed = [np.array([complex(v) for v in p]) for p in solutions]
    for root in known:
        nearest = min(np.max(np.abs(c - root)) for c in computed)
        assert nearest < 1e-8

Classify the endpoints
======================

Rather than hand-rolling cutoffs, ask the solver's metadata which endpoints are finite, real, or
singular.  Two of our four solutions are real (the :math:`(\pm\sqrt{y_1}, y_1)` pair); the other
two have purely imaginary :math:`x`:

.. testcode::

    assert len(solver.finite_solutions()) == 4       # all four endpoints are finite
    assert len(solver.nonsingular_solutions()) == 4  # ... and all nonsingular
    assert len(solver.real_solutions()) == 2         # two real, two with purely imaginary x

That is the whole arc: an exact, hand-authored product-of-linears start system, four start points
you wrote down yourself, blended into a homotopy and tracked through the same solver bertini uses
for its generated start systems -- landing on the target's four roots, neatly split into real and
complex.

Complete example
================

The whole tutorial as one runnable script -- assemble nothing, just run it:

.. literalinclude:: user_product_of_linears.py
   :language: python
   :caption: user_product_of_linears.py
