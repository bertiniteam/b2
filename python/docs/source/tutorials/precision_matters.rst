🎚️ When double precision is not enough
*****************************************

.. testsetup:: *

   import bertini

Tracking a path is a numerical tightrope walk.  Most paths are easy, but a few can pass close to a
singularity, where the system's Jacobian is nearly rank-deficient and a tiny numerical error gets
amplified enormously.  In ``double`` precision -- 53 bits, about 16 digits -- the tracker can run
out of room on such a path: it shrinks the step size to keep the error in check, hits the floor,
and **gives up**.  When that happens you do not get a wrong answer, you get a *missing* one: a
genuine solution silently absent from the results.

The cure is **adaptive precision**.  Instead of failing, the tracker raises the working precision
on exactly the hard stretches of exactly the hard paths -- 100, 200 digits if it needs them -- and
follows the path to the end.  This tutorial shows the difference on a system where it matters.

A system with near-singular paths
=================================

The cyclic-:math:`n` roots system is a classic benchmark.  For :math:`n = 5` it has a number of
finite solutions that is *known exactly* -- **70** -- so we have a ground truth to check against
(G. Björck and R. Fröberg, *J. Symbolic Comput.* 12(3), 1991).  Its total-degree homotopy tracks
:math:`5! = 120` paths; 70 converge to finite roots and the other 50 diverge to infinity.  Several
of the 70 sit in tight, near-coincident clusters (the system has a dihedral symmetry), and the
paths leading to them are ill-conditioned -- the kind double precision struggles with.

.. testcode::

   import numpy as np
   import bertini
   from bertini.nag_algorithm import ZeroDim

   n = 5
   x = [bertini.Variable('x{}'.format(i)) for i in range(n)]
   w = x + x                                        # doubled, to take products that wrap around
   sys = bertini.System()
   for length in range(1, n):                       # the degree-`length` cyclic sums
       sys.add_function(np.sum([np.prod(w[s:s + length]) for s in range(n)]))
   sys.add_function(np.prod(x) - 1)                 # the product, normalized
   sys.add_variable_group(bertini.VariableGroup(x))

Solve it reliably, with adaptive precision
==========================================

Ask for adaptive precision with ``mptype='adaptive'`` and solve.  We check the answer the same way
a careful user should: not by trusting the raw endpoint list, but by asking the solver how it
**classified** each path, using the *named* success codes:

.. testcode::

   SC = bertini.tracking.SuccessCode

   solver = ZeroDim(sys, mptype='adaptive')
   solver.solve()
   md = solver.solution_metadata()

   # A genuine solution is a path whose endgame SUCCEEDED and whose endpoint is FINITE.  Count
   # DISTINCT points by summing 1/multiplicity, so a true multiple root is counted once.
   finite = [m for m in md if m.endgame_success == SC.Success and m.is_finite]
   distinct = round(sum(1.0 / m.multiplicity for m in finite))
   assert distinct == 70                            # every finite solution of cyclic-5

Just as important: **no path was lost**.  Every one of the 120 paths reached a definite outcome --
either it succeeded, or it cleanly diverged (``GoingToInfinity`` is a *result*, not a failure).  A
path that ended any other way is one the tracker could not follow, and any root it was heading for
is missing:

.. testcode::

   lost = [m for m in md
           if m.endgame_success not in (SC.Success, SC.GoingToInfinity)]
   assert lost == []                                # nothing fell off the tightrope

How a count can lie
===================

The same solve in ``double`` precision is faster, and *most* of the time it also finds all 70.  But
about one run in ten, one of those near-singular paths fails -- the tracker hits its minimum step
size and abandons it -- and the distinct count quietly comes back **69**:

.. testcode::

   solver_d = ZeroDim(sys, mptype='double')
   solver_d.solve()
   md_d = solver_d.solution_metadata()

   distinct_d = round(sum(1.0 / m.multiplicity
                          for m in md_d if m.endgame_success == SC.Success and m.is_finite))
   assert distinct_d <= 70                          # never too many -- but sometimes too few

   lost_d = [m for m in md_d
             if m.endgame_success not in (SC.Success, SC.GoingToInfinity)]

The crucial point is that the solver *knows*.  When a path is dropped, it is right there in the
metadata with a name on it:

.. code-block:: python

   for m in lost_d:
       print(m.solution_index, m.endgame_success)

   # on a run that drops a root, this prints something like:
   #   99  SuccessCode.MinStepSizeReached

``SuccessCode.MinStepSizeReached`` is the tracker telling you, in plain language, that it could not
follow that path far enough in double precision.  A program that only looks at ``len(solutions())``
or the distinct count never sees it -- it just reports 69 and moves on, silently wrong.

The lesson
==========

Two habits keep a solve honest:

#. **Reach for adaptive precision when paths are hard.**  ``mptype='adaptive'`` costs more
   arithmetic than ``'double'``, but it spends that cost *only* where the geometry demands it, and
   it turns "sometimes 69" into "always 70".  Pinning a random seed only makes a flaky run
   *reproducible*; it is never the fix.
#. **Check the success codes, by name, not just the count.**  ``[m for m in
   solver.solution_metadata() if m.endgame_success not in (SC.Success, SC.GoingToInfinity)]`` is the
   list of paths the solver could not resolve.  If it is empty, you found everything; if it is not,
   the names tell you what went wrong.
