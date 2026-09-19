Auxiliary coordinates
*********************

Some systems carry coordinates that exist for the construction rather than for the answer.  A
critical-point system in null-vector form is ``f; v^T M; patch``, where ``v`` needs its own patch
or the system would admit the zero null vector.  That patch fixes a normalization, so ``v``'s
direction is the meaningful object and its magnitude is whatever the patch made it.

Bertini judges a point by its largest coordinate.  It does this three times: the tracker truncates
a path whose point exceeds ``path_truncation_threshold``, the endgame abandons one that exceeds
``Security.max_norm``, and the solver classifies an endpoint as at infinity when it exceeds
``endpoint_finite_threshold``.  If the largest coordinate belongs to a block like ``v``, all three
are reading a number nobody chose.  No threshold is the right one, because the quantity is not
governed.

The fix is to say which coordinates the question is about.

Declaring them
==============

A system is told which of its coordinates are auxiliary, by variable group or by index:

.. testcode::

   import numpy as np
   import bertini as pb

   x, v = pb.Variable('x'), pb.Variable('v')

   sys = pb.System()
   sys.add_variable_group(pb.VariableGroup([x]))
   sys.add_variable_group(pb.VariableGroup([v]))
   sys.add_function(x*x - 1)
   sys.add_function(v - 1)

   point = np.array([complex(1e3), complex(1e9)], dtype=complex)
   print(sys.is_finite(point, 1e5))

   sys.set_auxiliary_variable_groups([1])
   print(sys.is_finite(point, 1e5))

.. testoutput::

   False
   True

The point is the same in both calls.  In the first, ``v = 1e9`` decides the verdict.  In the
second, the system has been told that ``v`` is not what finiteness is about, so the verdict rests
on ``x = 1e3``.

``set_auxiliary_coordinates`` takes indices instead of groups, for a system whose grouping does
not separate the two kinds of coordinate.  Indices are into a point in your own coordinates, in
variable-group order, so they still mean the same thing after the solver homogenizes and patches.
The two lists are unioned; either can be emptied by passing ``[]``.

What it changes in a solve
==========================

``{x^2 - 1, v - 100}`` has the solutions ``(1, 100)`` and ``(-1, 100)``.  With the finiteness
threshold at 10, ``v`` puts both of them past it:

.. testcode::

   def solve_it(auxiliary):
       s = pb.System()
       s.add_variable_group(pb.VariableGroup([x, v]))
       s.add_function(x*x - 1)
       s.add_function(v - 100)
       if auxiliary:
           s.set_auxiliary_coordinates([1])

       solver = pb.nag_algorithm.ZeroDimSolver(s)
       solver.update(endpoint_finite_threshold=10)
       solver.solve()
       md = solver.solution_metadata()
       return sum(1 for m in md if m.is_finite), sum(1 for m in md if m.is_real)

   print(solve_it(auxiliary=False))
   print(solve_it(auxiliary=True))

.. testoutput::

   (0, 0)
   (2, 2)

Both runs track the same two paths to the same two points, and both report ``Success``.  What
changes is the verdict: without the declaration the solver calls two ordinary roots infinite, and
calls them complex as well, because realness is measured over the same coordinates.

The half you cannot work around
===============================

Classification happens after the solve, so a caller who disagrees with it can always compute
their own.  Truncation happens during tracking.  A path the tracker abandoned partway, because a
coordinate that was never the question grew large, is simply not there afterwards, and no
post-processing recovers it.  That is the reason this belongs to the system rather than to a
setting applied at the end: the same declaration reaches the tracker and the endgame, which is
where it matters most.

What it does not mean
=====================

Bertini attaches no further meaning to the word.  An auxiliary coordinate is tracked like any
other, is stored in the solution, is written to records, and is returned to you.  It is excluded
from two questions and nothing else.

Marking every coordinate auxiliary is refused.  A point with nothing left to judge would be
finite and real unconditionally, which is not a verdict worth rendering.

It is part of the system's identity
===================================

Two systems that differ only in which coordinates are auxiliary have different content digests:

.. testcode::

   plain = pb.System()
   plain.add_variable_group(pb.VariableGroup([x, v]))
   plain.add_function(x*x - 1)
   plain.add_function(v - 100)

   declared = pb.System.from_canonical(plain.canonical_encoding())
   declared.set_auxiliary_coordinates([1])

   print(plain.content_digest() == declared.content_digest())

.. testoutput::

   False

This is deliberate.  The tracker truncates on the coordinates that are *not* auxiliary, so the two
systems do not compute the same thing, and a records directory must not answer one with the other.
Adding or removing a declaration makes a new ask, which is recomputed rather than recalled.

.. seealso::

   :doc:`/detailed/configuration` for the thresholds named here, and
   :doc:`../crossed_paths/index` for the other per-path verdict a solve can render.
