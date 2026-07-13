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
   from bertini.nag_algorithm import ZeroDimSolver

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

Ask for adaptive precision with ``mptype='adaptive'``, solve, and check the answer the way a careful
user should -- not by trusting the raw endpoint count, but by asking the solver for its **report**:

.. testcode::

   solver = ZeroDimSolver(sys, mptype='adaptive')
   solver.solve()
   report = solver.report()

   assert report.num_finite_solutions == 70         # every finite solution of cyclic-5
   assert report.num_failed == 0                    # no path the tracker had to give up on

The report is a one-call summary of how every path ended up.  ``print(report)`` shows it::

   zero-dim solve -- 120 paths tracked
     finite solutions    70   (distinct)
     diverged            50
     FAILED               0
     ----
     ...
     all paths resolved? yes

``num_finite_solutions`` counts *distinct* points (it sums ``1/multiplicity``, so a genuine multiple
root is counted once).  ``num_failed`` is the **precision** headline: it counts paths the tracker had
to *give up* on -- the ``MinStepSizeReached`` failures that adaptive precision exists to prevent.
Zero of them here.  (The report also exposes ``all_paths_resolved``, a stricter all-in-one flag that
*additionally* requires the midpath check to have cleared every path crossing -- a separate,
predictor-driven concern covered in :doc:`/tutorials/settings_and_precision/crossed_paths/index`, not a precision one.)

How a count can lie
===================

The same solve in ``double`` precision is faster, and *most* of the time it also finds all 70.  But
about one run in ten, one of those near-singular paths fails -- the tracker hits its minimum step
size and abandons it -- and the finite count quietly comes back **69**:

.. testcode::

   solver_d = ZeroDimSolver(sys, mptype='double')
   solver_d.solve()
   report_d = solver_d.report()

   assert report_d.num_finite_solutions <= 70       # never too many -- but sometimes too few

The crucial point is that the solver *knows*; a short count is never silent.  When a path is lost,
``report_d.all_paths_resolved`` is ``False`` and the reason is named in ``report_d.failures_by_reason``.
On a run that drops a root, ``print(report_d)`` shows it plainly:

.. code-block:: text

   zero-dim solve -- 120 paths tracked
     finite solutions    69   (distinct)
     diverged            50
     FAILED               1    1 x MinStepSizeReached
     ----
     ...
     all paths resolved? NO

``MinStepSizeReached`` is the tracker telling you, in plain language, that it could not follow that
path far enough in double precision.  A program that trusts only the count never sees it -- it
reports 69 and moves on, silently wrong.  ``report.all_paths_resolved`` does.

The lesson
==========

Two habits keep a solve honest:

#. **Reach for adaptive precision when paths are hard.**  ``mptype='adaptive'`` costs more
   arithmetic than ``'double'``, but it spends that cost *only* where the geometry demands it, and it
   turns "sometimes 69" into "always 70".  Pinning a random seed only makes a flaky run
   *reproducible*; it is never the fix.
#. **Check** ``report().num_failed``\ **, not just the count.**  The solve report classifies every
   path; if any *failed to track*, ``num_failed`` is nonzero and ``failures_by_reason`` names the
   reason.  A count alone can hide a root the tracker silently lost.

Complete example
================

The whole tutorial as one runnable script -- assemble nothing, just run it:

.. literalinclude:: precision_matters.py
   :language: python
   :caption: precision_matters.py
