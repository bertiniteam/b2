🔁 Parameter homotopy: solve once, re-solve many times
**********************************************************

.. testsetup:: *

   import bertini

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
and :func:`bertini.nag_algorithm.HomotopySolver` (runs the full zero-dim pipeline -- pre-endgame
tracking, the midpath check, the endgame, post-processing -- from a homotopy you constructed and
a list of start points you already have).

A family of systems
===================

Take a fixed unit circle intersected with a horizontal line whose height is the parameter:

.. testcode::

    import bertini
    from bertini import nag_algorithm

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
    first = nag_algorithm.ZeroDimSolver(generic, mptype='adaptive')
    first.solve()
    start_points = first.all_solutions()                   # (+/- sqrt(3)/2, 1/2)

Move the parameter -- without solving again
===========================================

To reach another member, build the parameter homotopy from that member back to the generic one
and track the start points through it:

.. testcode::

    target = member(0)                                 # the line y = 0
    H = nag_algorithm.coefficient_parameter_homotopy(target, generic)
    moved = nag_algorithm.HomotopySolver(H, start_points, target)
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
        solver = nag_algorithm.HomotopySolver(H, start_points, target)
        solver.solve()
        roots = [p for p in solver.all_solutions() if len(p) == 2]
        for p in roots:
            xv, yv = complex(p[0]), complex(p[1])
            assert abs(xv*xv + yv*yv - 1) < 1e-8       # on the circle
            assert abs(2*yv - s) < 1e-8                # on the line y = s/2

Each iteration tracks only the two solutions we already have, to the new line -- not a fresh
total-degree solve.  Across a large sweep that is the difference between tracking a handful of
paths per parameter and tracking the full Bézout count every time.

Choosing the generic member
===========================

A parameter homotopy works because the *singular* parameter values -- where solutions collide or
run off to infinity -- form a measure-zero set, so a straight-line path between two generic
values misses them.  In the example above the real path stays safely inside :math:`|s| < 2`
(the circle and line stay transverse), so a real generic value is fine.  In general, choose the
generic member's coefficients to be **generic complex numbers**; then the path avoids the bad
set with probability one.  And keep coefficients *exact* -- :mod:`bertini.linalg` refuses python
floats so low-precision literals cannot silently cap your precision.
