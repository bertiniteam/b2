⏱️ Stopping a solve, and giving paths a budget
**********************************************

Path costs vary widely.  In a sweep of tens of thousands of paths a few may cost a thousand
times the median, and those few set the wall-clock time of the whole sweep.  This page covers
stopping a solve, giving paths a wall-clock budget, and what a path records when it is
abandoned.

The mechanism is the same throughout: a path can be stopped between steps.  At the top of
every predictor-corrector step the tracker checks whether a stop has been requested and
whether its wall-clock deadline has passed, and if so returns there with its position intact.
The overrun is at most the step in flight.

Ctrl-C works
============

Press Ctrl-C during a solve and it stops at its next step.  You get a ``KeyboardInterrupt``,
as you would from any Python call, and the solver is left exactly as the stop found it:

.. code-block:: python

    solver = bertini.ZeroDimSolver(sys)
    try:
        solver.solve()
    except KeyboardInterrupt:
        pass
    solver.was_stopped_early()          # True
    solver.all_solutions()              # what finished, one entry per path
    solver.solution_metadata()          # which paths succeeded, were stopped, or never began
    solver.solve()                      # resumes: finished paths recalled, the rest tracked

Nothing is destroyed by the exception.  The paths that finished are readable, the paths that
were in flight are marked ``ExternallyTerminated``, and the paths not yet begun stay
``NeverStarted``.  Because every finished path is on record (see
:doc:`../../record_keeping/automatic_record_keeping/index`), solving again picks up where you
left off.

The same stop can be requested from code.  :func:`bertini.request_stop` sets the flag that
Ctrl-C sets, and the solve returns normally instead of raising.  A signal cannot be delivered
from inside a tutorial, so here the request comes from an observer that fires once two paths
have completed:

.. testcode::

   import bertini as pb

   x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')

   # a system with enough paths (64) that a stop lands while paths remain
   big = pb.System()
   big.add_variable_group(pb.VariableGroup([x, y, z]))
   big.add_function(x**4 + y**4 + z**4 + x*y*z - 1)
   big.add_function(x**4 - 2*y**4 + z**3 + x*y - 3)
   big.add_function(x**3*z + y**4 - z**4 + x + 2*y - 5)

   class StopAfter(pb.nag_algorithm.observers.CustomObserver):
       def __init__(self, n):
           super().__init__()
           self.n, self.completed = n, 0
       def Observe(self, event):
           if isinstance(event, pb.nag_algorithm.observers.PathComplete):
               self.completed += 1
               if self.completed == self.n:
                   pb.request_stop()

   solver = pb.ZeroDimSolver(big)
   trigger = StopAfter(2)
   solver.add_observer(trigger)
   solver.solve()

   codes = [m.pre_endgame_success_code for m in solver.solution_metadata()]
   print('stopped early:', solver.was_stopped_early())
   print('never started:', solver.num_paths_never_started() == codes.count(pb.SuccessCode.NeverStarted))
   print('every path accounted for:',
         codes.count(pb.SuccessCode.Success)
         + codes.count(pb.SuccessCode.ExternallyTerminated)
         + codes.count(pb.SuccessCode.NeverStarted) == len(codes))

.. testoutput::

   stopped early: True
   never started: True
   every path accounted for: True

How many paths were in flight when the stop landed, versus not yet begun, depends on how many
threads were running, which is why the tutorial prints the invariants and not the counts.
Being stopped is not sticky.  Take the observer off and solve again, and the same solver runs
to completion:

.. testcode::

   solver.remove_observer(trigger)
   solver.solve()
   print('stopped early:', solver.was_stopped_early(),
         '  never started:', solver.num_paths_never_started())

.. testoutput::

   stopped early: False   never started: 0

What an abandoned path leaves behind
====================================

A path that did not succeed -- whatever the reason -- says where it got to.  Its metadata
carries ``final_time_used``, the time it reached; ``last_point``, the point it was at, in the
solver's internal coordinates; ``num_successful_steps`` and ``num_failed_steps``, the
predictor-corrector steps over the whole path; and ``max_precision_used``.  A successful
path carries the same fields except ``last_point``, which is empty: its endpoint is the
solution itself.

To see the stamp on a path that got partway, give the tracker a budget of five steps.  The
system from here on is two quadrics with four well-separated roots, so four paths:

.. testcode::

   s = pb.System()
   s.add_variable_group(pb.VariableGroup([x, y]))
   s.add_function(x**2 - 1)
   s.add_function(y**2 - 1)

   solver = pb.ZeroDimSolver(s)
   solver.get_tracker().get_stepping().update(max_num_steps=5)
   solver.solve()

   m = solver.solution_metadata()[0]
   print(m.pre_endgame_success_code)
   print('steps taken:', m.num_successful_steps + m.num_failed_steps)
   print('stopped short of the start at t=1:', 0 < abs(m.final_time_used) < 1)
   print('coordinates in the last point:', len(m.last_point))

.. testoutput::

   MaxNumStepsTaken
   steps taken: 5
   stopped short of the start at t=1: True
   coordinates in the last point: 3

The last point has three coordinates, not two: the solver tracks in homogeneous coordinates,
so the point is the two variables plus the homogenizing coordinate.  A solution's coordinates
are dehomogenized for you; a last point is left as the tracker had it.

The stamp is recorded with the path.  The time reached, the steps taken and the precision
describe how far a path got and what it cost in terms that do not depend on the machine.

A budget per path
=================

``max_path_wall_clock_duration``, in seconds, gives every path a wall-clock budget.  A path
that has not finished when its budget runs out is abandoned between steps with
``WallClockLimitReached``, stamped as above, and recorded with the budget that stopped it.
The budget spans the whole path, pre-endgame tracking and endgame together.

A tutorial has to be the same on every machine, so the budgets below are either a nanosecond,
which is gone before the first step anywhere, or an hour, which is never reached.  In real
use you would pick something like ten times your median path time.

.. testcode::

   solver = pb.ZeroDimSolver(s)
   solver.update(max_path_wall_clock_duration=1e-9)
   solver.solve()

   print({m.pre_endgame_success_code.name for m in solver.solution_metadata()})
   print('budget recorded with each path:',
         all(m.wall_clock_limit_seconds == 1e-9 for m in solver.solution_metadata()))
   print('stopped early:', solver.was_stopped_early())

.. testoutput::

   {'WallClockLimitReached'}
   budget recorded with each path: True
   stopped early: False

Running out of a per-path budget does not count as being stopped early: the solve was not
interrupted, each path used up its own allowance, and every path was asked.  With a larger
budget the same paths finish:

.. testcode::

   solver.update(max_path_wall_clock_duration=3600)
   solver.solve()
   print({m.endgame_success_code.name for m in solver.solution_metadata()})

.. testoutput::

   {'Success'}

A bare tracker can be limited on its own, with no solver anywhere:
``tracker.set_max_wall_clock_duration(seconds)`` arms a deadline that holds across
``track_path`` calls until ``tracker.clear_max_wall_clock_time()``.  It is a deadline rather
than a duration because an endgame issues hundreds of tracking calls for one path, and a
budget that restarted with each call would never bite.

A budget for the whole solve
============================

``max_solve_wall_clock_duration`` budgets the whole ``solve()`` call instead.  Once it runs out
no further path is started and any path in flight is abandoned, and the solve reads exactly as
if you had pressed Ctrl-C:

.. testcode::

   solver = pb.ZeroDimSolver(s)
   solver.update(max_solve_wall_clock_duration=1e-9)
   solver.solve()
   print('stopped early:', solver.was_stopped_early(),
         '  never started:', solver.num_paths_never_started())

.. testoutput::

   stopped early: True   never started: 4

The per-path budget is the finer tool: it guarantees progress path by path, and a path it
abandons is a fact about that path.  The whole-solve budget is for thinking in terms of the
call, and a path it abandons is treated as interrupted, because the reason had nothing to do
with the path.

Budgets, records, and what counts as answered
=============================================

Neither budget is part of a run's identity.  A budget says how long to wait for an answer,
not what the answer is, so a sixty-second run and a sixty-one-second run are the same ask and
share their records.  That raises a question the records system has to answer: when a path
was abandoned last time, may its record stand in for tracking it this time?

The rule is that a record is reusable when the run that would reuse it asks no more of that
path than the run that produced it.  A path that completed -- succeeded, diverged, or failed
under settings that are part of the ask -- always qualifies.  An abandoned path was ended for
a reason kept out of the ask, so whether to reuse it is your call, made through
``RecordsConfig.recall``:

* ``RecallPolicy.Completed`` (the default) reuses completed paths and re-tracks abandoned
  ones, with one refinement: a path a per-path budget cut off is reused while the current
  budget is no larger than the one that abandoned it, and re-tracked once the budget is
  raised.  An interrupted path is always re-tracked.
* ``RecallPolicy.Everything`` reuses every recorded outcome as recorded, abandonments
  included.
* ``RecallPolicy.Nothing`` tracks every path fresh, records or not: for observers,
  benchmarking, and re-verification.  ``recall=False`` means the same thing.

The default, with the seed fixed so that every solve below is the same ask:

.. testcode::

   import tempfile, os
   where = os.path.join(tempfile.mkdtemp(), 'records')

   def solve_with(path_budget):
       pb.random.set_random_seed(42)
       solver = pb.ZeroDimSolver(s)
       solver.record_to(where)
       solver.update(max_path_wall_clock_duration=path_budget)
       solver.solve()
       return solver

   first = solve_with(1e-9)                  # four abandonments, recorded with their budget
   print('recalled', first.num_paths_recalled())

   same = solve_with(1e-9)                   # same patience asked: the abandonments stand in
   print('recalled', same.num_paths_recalled())

   more = solve_with(3600)                   # more patience: re-tracked, and they finish
   print('recalled', more.num_paths_recalled(),
         ' finished', sum(m.endgame_success_code == pb.SuccessCode.Success
                          for m in more.solution_metadata()))

   again = solve_with(3600)                  # now they are completed paths: reused
   print('recalled', again.num_paths_recalled())

.. testoutput::

   recalled 0
   recalled 4
   recalled 0  finished 4
   recalled 4

Each of those runs recorded its paths, abandonments included.  The record of an abandoned
path holds its code, where it got to, and the budget it ran under; the policy decides how a
later run uses it.

A wall clock is machine dependent, so the same budget means the same amount of work only on
the same machine.  The comparison above is exact on one machine and approximate across
machines.  The steps, precision and time reached in the stamp do not depend on the machine.
