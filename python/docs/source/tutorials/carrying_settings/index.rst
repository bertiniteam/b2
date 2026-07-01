⚙️ Carrying settings across related solves
*******************************************

.. testsetup:: *

   import numpy as np
   import bertini

A real computation is rarely one solve. A numerical irreducible decomposition, a parameter sweep, a
multi-stage pipeline -- each runs *many* related solves, and you usually want the **same** tracking
settings on every one of them. Re-typing tolerances and step controls for each solver is both tedious
and a place for them to silently drift apart.

Bertini gives every **config owner** -- every tracker and every solver -- the same small interface for
this. (Owning configs is a capability, not a class: a tracker is not a solver, but both carry configs,
so both get these methods.)

Set fields by name
==================

Every numeric config field accepts a **string** -- the exact, noise-free way to write a number
(``0.05`` the Python float is not :math:`1/20`; ``"0.05"`` is). And you can set fields right on the
solver without naming which config struct they live in: :meth:`update` routes each field to whichever
config owns it.

.. testcode::

   x, y = bertini.Variable('x'), bertini.Variable('y')
   system = bertini.System()
   system.add_function(x*x + y*y - 1)
   system.add_function(x + y)
   system.add_variable_group(bertini.VariableGroup([x, y]))

   solver = bertini.nag_algorithm.ZeroDimSolver(system)
   solver.update(final_tolerance="1e-11",                 # -> TolerancesConfig
                 max_num_crossed_path_resolve_attempts=3) # -> ZeroDimConfig

   from bertini.nag_algorithm import TolerancesConfig, ZeroDimConfig
   assert solver.get_config(TolerancesConfig).final_tolerance == 1e-11
   assert solver.get_config(ZeroDimConfig).max_num_crossed_path_resolve_attempts == 3

A misspelled field, or one that lives on a different owner, raises immediately rather than silently
doing nothing -- ``max_step_size`` is a *tracker* setting, so the solver rejects it:

.. testcode::

   try:
       solver.update(max_step_size="0.05")     # that field is on the tracker, not the solver
       raise AssertionError("should have raised")
   except AttributeError:
       pass

   solver.get_tracker().update(max_step_size="0.05")   # set it where it lives
   from bertini.tracking import SteppingConfig
   assert solver.get_tracker().get_config(SteppingConfig).max_step_size == bertini.multiprec.Float("0.05")

Carry a whole bundle
====================

:meth:`get_settings` hands you the owner's entire configuration as a plain dict ``{name: config}`` --
independent, picklable copies. :meth:`set_settings` applies one back. So you configure **once** and
stamp every subsequent solver with the same settings.

.. testcode::

   reference = bertini.nag_algorithm.ZeroDimSolver(system)
   reference.update(final_tolerance="1e-11", newton_before_endgame="1e-6")
   settings = reference.get_settings()
   assert set(settings) == set(reference.config_names())     # one entry per config

   # ... later, for each related solve ...
   next_solver = bertini.nag_algorithm.ZeroDimSolver(system)
   next_solver.set_settings(settings)
   assert next_solver.get_config(TolerancesConfig).final_tolerance == 1e-11

The bundle is an ordinary picklable value, so it can be stored or sent to another process -- the same
settings then drive a worker's solves as drive the manager's.

.. testcode::

   import pickle
   carried = pickle.loads(pickle.dumps(settings))
   worker_solver = bertini.nag_algorithm.ZeroDimSolver(system)
   worker_solver.set_settings(carried)
   assert worker_solver.get_config(TolerancesConfig).final_tolerance == 1e-11

The settings are precision-agnostic
===================================

A bundle built from one precision model applies **unchanged** to another. This is what lets a
multi-stage workflow mix precision models -- a quick double-precision pass, then an adaptive
re-solve -- while carrying one set of tracking tolerances through all of them.

.. testcode::

   tuned = bertini.nag_algorithm.ZeroDimSolver(system, mptype='double')
   tuned.update(final_tolerance="1e-10")
   bundle = tuned.get_settings()

   for mptype in ('multiple', 'adaptive'):
       solver = bertini.nag_algorithm.ZeroDimSolver(system, mptype=mptype)
       solver.set_settings(bundle)                       # drops on cleanly, any precision model
       assert solver.get_config(TolerancesConfig).final_tolerance == 1e-10

By default :meth:`set_settings` applies only the configs the target actually has and skips the rest,
so a bundle also moves between *different kinds* of owner (a tracker has no solver-level tolerances).
Pass ``strict=True`` to instead require every config in the bundle to be applicable.

.. testcode::

   from bertini.tracking import AMPTracker
   tracker = AMPTracker(system)
   tracker.set_settings(settings)             # silently skips the solver-only configs
   try:
       tracker.set_settings(settings, strict=True)
       raise AssertionError("should have raised")
   except KeyError:
       pass

Put together, the pattern for any multi-solve workflow is: **tune one owner, ``get_settings()``, and
``set_settings()`` on each of the rest.**

Complete example
================

The whole tutorial as one runnable script -- assemble nothing, just run it:

.. literalinclude:: carrying_settings.py
   :language: python
   :caption: carrying_settings.py
