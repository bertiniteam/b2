🌐 Finding *all* the solutions
*********************************************

.. testsetup:: *

   import bertini

Most numerical solvers find *a* solution near where you start them.  Bertini can compute **all
isolated complex solutions** of a polynomial system -- including the complex ones a
real-valued solver can never see -- and it knows in advance an upper bound on how many to expect.

A system with no real solutions
===============================

Take the unit circle and the hyperbola :math:`xy = 1`

.. testcode::

    import numpy as np
    import bertini as bertini
    from bertini.nag_algorithm import ZeroDim

    x, y = bertini.Variable('x'), bertini.Variable('y')

    sys = bertini.System()
    sys.add_function(x**2 + y**2 - 1)     # the unit circle
    sys.add_function(x*y - 1)             # the hyperbola xy = 1
    sys.add_variable_group(bertini.VariableGroup([x, y]))

The circle and that hyperbola do **not** meet anywhere in the real plane, so a real
Newton solver returns nothing useful.  But over the complex numbers they meet in exactly
:math:`2 \times 2 = 4` points (the product of the degrees -- the total-degree Bézout
number).

Solve, and collect the finite solutions.  ``ZeroDim`` is a small factory over the bound
solver classes: ``ZeroDim(sys)`` is the Cauchy endgame in adaptive precision with a total-degree
start system, and you pick the rest with strings -- ``endgame=`` (``'cauchy'`` / ``'powerseries'``),
``mptype=`` (``'double'`` / ``'multiple'`` / ``'adaptive'``), and ``startsystem=`` (``'totaldegree'``
/ ``'mhom'``).  Here we ask for adaptive precision explicitly so ill-conditioned paths still succeed:

.. testcode::

    solver = ZeroDim(sys, mptype='adaptive')
    solver.solve()

    good = solver.finite_solutions()      # successful, finite endpoints -- the actual points

    assert len(good) == 4                 # all four complex solutions

Every one of them is genuinely complex:

.. testcode::

    for s in good:
        xv, yv = complex(s[0]), complex(s[1])
        assert abs(xv.imag) > 1e-6 or abs(yv.imag) > 1e-6   # none are real
        # and each really is a solution
        assert abs(xv*xv + yv*yv - 1) < 1e-8
        assert abs(xv*yv - 1) < 1e-8

There are no real solutions to filter out here -- the whole solution set lives off the
real plane, which is exactly the point: homotopy continuation found a structure that real
methods cannot.

Why the count is guaranteed
===========================

Bertini did not stumble onto four solutions by luck.  It built a *start system* with
exactly the total-degree Bézout number of start points (four), and tracked one homotopy
path from each.  Paths that diverge to infinity are reported as failures rather than
silently dropped, so the successful endpoints are the complete set of finite isolated
solutions.

Every solution, and the at-infinity ones
========================================

``finite_solutions()`` above is one of a family of accessors.  ``all_solutions()`` is the whole
list -- **one entry per tracked path**, finite and at-infinity alike -- and the categories are
filtered views of it.  Because this system has four finite solutions and the total-degree start
tracked exactly four paths, none diverged, so the finite set *is* the whole set here:

.. testcode::

    assert len(solver.all_solutions()) == 4          # one per tracked path
    assert len(solver.finite_solutions()) == 4        # the finite ones
    assert len(solver.infinite_solutions()) == 0      # the at-infinity ones (is_finite is False)

``infinite_solutions()`` is the complement of ``finite_solutions()`` within ``all_solutions()`` --
the endpoints the endgame resolved as diverging to infinity.  The other filtered views are
``real_solutions()``, ``nonsingular_solutions()``, and ``singular_solutions()``; every accessor
takes ``user_coords=False`` to hand back the solver's internal homogenized coordinates instead of
your variables'.

The database of solutions
=========================

For bookkeeping across a whole solve -- or many solves -- ``to_dataframe()`` returns the solve as a
:mod:`pandas` DataFrame: **one row per solution**, with a column per coordinate (``x0``, ``x1``,
...) followed by every field of the per-solution metadata (``is_finite``, ``is_real``,
``is_singular``, ``multiplicity``, ``condition_number``, ``endgame_success``, ...).  Each category
is then a one-line filter on the frame, and you have the metadata right alongside the points::

    df = solver.to_dataframe()                  # finite solutions only, by default
    df[df.is_real & ~df.is_singular]            # the nonsingular real ones, with their metadata

    everything = solver.to_dataframe(omit_infinite=False)   # every tracked path
    everything[~everything.is_finite]                       # the endpoints at infinity

``omit_infinite`` is ``True`` by default, so the frame holds just the genuine finite solutions.
:mod:`pandas` is an optional dependency: ``to_dataframe()`` raises a helpful :class:`ImportError`
if it is not installed, and the point accessors above are the no-pandas path.

A note on precision
===================

We used the **adaptive-precision** tracker.  On a well-conditioned path it works in fast
double precision; when a path approaches a singularity and double precision would lose
the solution, it transparently raises the working precision (and lowers it again
afterwards).  That is what lets the same call reliably solve both gentle and nasty
systems.  You can inspect or set the baseline precision with
:func:`bertini.default_precision`; solutions can be read back at that precision through the
multiprecision interface rather than rounded to ``float``.

When a system has structure -- several groups of variables that each appear with low
degree -- the total-degree count above is wasteful.  The
:doc:`eigenvalues_by_homotopy` tutorial shows the multihomogeneous start system, which
tracks one path per *actual* solution instead.
