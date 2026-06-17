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

::

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

Solve, and collect the successful endpoints.  ``ZeroDim`` is a small factory over the bound
solver classes: ``ZeroDim(sys)`` is the Cauchy endgame in multiple precision with a total-degree
start system, and you pick the rest with strings -- ``endgame=`` (``'cauchy'`` / ``'powerseries'``),
``mptype=`` (``'double'`` / ``'multiple'`` / ``'adaptive'``), and ``startsystem=`` (``'totaldegree'``
/ ``'mhom'``).  Here we ask for adaptive precision so ill-conditioned paths still succeed::

    solver = ZeroDim(sys, mptype='adaptive')
    solver.solve()

    OK = int(bertini.tracking.SuccessCode.Success)
    sols = solver.solutions()
    meta = solver.solution_metadata()
    good = [sols[i] for i in range(len(sols))
            if int(meta[i].endgame_success) == OK]

    assert len(good) == 4                 # all four complex solutions

Every one of them is genuinely complex::

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
