🔁 Parameter homotopy: solve once, re-solve many times
**********************************************************

.. testsetup:: *

   import bertini
   _docs_ambient = bertini.records_dir()   # restored in testcleanup

.. testcleanup:: *

   import shutil
   shutil.rmtree("circle_sweep_records", ignore_errors=True)
   bertini.records_dir(_docs_ambient)

The most expensive part of homotopy continuation is the *first* solve -- the one that has to
find every solution from scratch.  But a great many problems are really a **family** of systems
that differ only in their coefficients: the same equations, evaluated at different parameter
values.  The eigenvalues of :math:`A(s)` as :math:`s` varies; the intersections of a fixed curve
with a moving line; a design swept across a range of settings.

For such a family you should solve **once**, at a generic parameter, and then **reuse** those
solutions as the start points of a *parameter homotopy* that slides the coefficients to whatever
value you actually care about -- as many times as you like, each move tracking just the solutions
you already have.  This is bertini's "evaluate as little as possible" in action, and it is
pleasantly parallel across the parameter values.

The tools are :func:`bertini.nag_algorithm.coefficient_parameter_homotopy` (builds the homotopy)
and :func:`bertini.HomotopySolver` (runs the full zero-dim pipeline -- pre-endgame
tracking, the midpath check, the endgame, post-processing -- from a homotopy you constructed and
a list of start points you already have).

A family of systems
===================

Take a fixed unit circle intersected with a horizontal line whose height is the parameter:

.. testcode::

    import bertini
    from bertini import nag_algorithm

    bertini.records_dir("circle_sweep_records")   # one line: record everything below
    bertini.random.set_random_seed(42)            # same seed => reruns resume, not recompute

    # the variables are SHARED across every member of the family: the parameter homotopy
    # interpolates the members' equations, so they must be built over the same Variable objects.
    x, y = bertini.Variable('x'), bertini.Variable('y')

    def member(s):
        """The system { x^2 + y^2 - 1, 2y - s }: the unit circle meeting the line y = s/2."""
        sys = bertini.System()
        sys.add_variable_group(bertini.VariableGroup([x, y]))
        sys.add_function(x*x + y*y - 1)
        sys.add_function(2*y - s)
        return sys

Solve once, ab initio
=====================

Pick a generic member and solve it the usual way (a total-degree start system).  Its two
solutions are the start points we will reuse forever after:

.. testcode::

    generic = member(1)                                # the line y = 1/2
    first = bertini.ZeroDimSolver(generic, mptype='adaptive')
    first.solve()
    start_points = first.all_solutions()                   # (+/- sqrt(3)/2, 1/2)

Move the parameter -- without solving again
===========================================

To reach another member, build the parameter homotopy from that member back to the generic one
and track the start points through it:

.. testcode::

    target = member(0)                                 # the line y = 0
    H = nag_algorithm.coefficient_parameter_homotopy(target, generic)
    moved = bertini.HomotopySolver(H, start_points, target)
    moved.solve()
    # moved.all_solutions() are now (+/- 1, 0)

``coefficient_parameter_homotopy(target, generic)`` is just :math:`(1-t)\,\text{target} +
t\,\text{generic}` with ``t`` as the path variable: at :math:`t=1` it is ``generic`` (so its
solutions are our start points) and at :math:`t=0` it is ``target``.

Now the payoff -- sweep as many parameters as you want, reusing the *same* start points, never
solving from scratch again:

.. testcode::

    for s in [0, -1, 1]:                               # lines y = 0, -1/2, 1/2
        target = member(s)
        H = nag_algorithm.coefficient_parameter_homotopy(target, generic)
        solver = bertini.HomotopySolver(H, start_points, target)
        solver.solve()
        roots = [p for p in solver.all_solutions() if len(p) == 2]
        for p in roots:
            xv, yv = complex(p[0]), complex(p[1])
            assert abs(xv*xv + yv*yv - 1) < 1e-8       # on the circle
            assert abs(2*yv - s) < 1e-8                # on the line y = s/2

Each iteration tracks only the two solutions we already have, to the new line -- not a fresh
total-degree solve.  Across a large sweep that is the difference between tracking a handful of
paths per parameter and tracking the full Bézout count every time.

The sweep, on the record
========================

Those two lines back at the top were not decoration.  Naming the ambient directory with
``records_dir`` turns recording on for **every** solver in the process -- the bare
:class:`~bertini.ZeroDimSolver` and :class:`~bertini.HomotopySolver` used here included,
not just :func:`bertini.solve`.  Every solve above wrote durable records of what it
computed, and every solve *consults* the records before computing.  The pinned seed is
what makes that pay: with the same seed a rerun of the script rebuilds the *same*
homotopies (randomness is seed-rooted), so every already-answered solve recalls from
the records instead of tracking again.  Kill a thousand-member sweep at member 700 and
rerun -- the first 700 come back instantly and the sweep continues where it died.
Without a pinned seed each run draws fresh randomness, and there is nothing to resume
*from*.

Read the sweep back with the navigation tools (see :doc:`../your_records/index` for the
full story):

.. testcode::

    runs = bertini.runs("circle_sweep_records")
    print(len(runs) >= 5, all(runs['num_paths'] == 2))

.. testoutput::

    True True

One ab-initio solve, then nothing but two-path parameter moves: exactly the yoga this
tutorial preaches, now auditable after the fact.

Choosing the generic member
===========================

A parameter homotopy works because the *singular* parameter values -- where solutions collide or
run off to infinity -- form a measure-zero set, so a straight-line path between two generic
values misses them.  In the example above the real path stays safely inside :math:`|s| < 2`
(the circle and line stay transverse), so a real generic value is fine.  In general, choose the
generic member's coefficients to be **generic complex numbers**; then the path avoids the bad
set with probability one.  And keep coefficients *exact* -- :mod:`bertini.linalg` refuses python
floats so low-precision literals cannot silently cap your precision.

Complete example
================

The whole tutorial as one runnable script -- assemble nothing, just run it:

.. literalinclude:: parameter_homotopy.py
   :language: python
   :caption: parameter_homotopy.py
