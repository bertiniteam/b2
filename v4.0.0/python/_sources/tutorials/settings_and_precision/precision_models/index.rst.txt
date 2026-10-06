🔬 Precision models: double, multiple, adaptive
*************************************************

.. testsetup:: *

   import numpy as np
   import bertini

Bertini 2 can track in three precision models, and ``ZeroDimSolver`` selects between them with ``mptype``:

* ``'double'`` -- hardware double precision (≈16 digits). Fastest; fine when the paths are
  well-conditioned.
* ``'multiple'`` -- a fixed multiprecision, the same number of digits throughout.
* ``'adaptive'`` -- adaptive multiprecision (AMP, the default): the precision rises and falls along
  each path as the conditioning demands. The most robust, and what makes a near-singular path solvable.

This tutorial is about how the three models differ *in the interface* -- the types you get back and
how you set the precision. For *why and when* you reach for adaptive precision (the conditioning
story), see :doc:`/tutorials/settings_and_precision/precision_matters/index`.

.. testcode::

   x, y = bertini.Variable('x'), bertini.Variable('y')
   system = bertini.System()
   system.add_function(x*x + y*y - 1)
   system.add_function(x + y)
   system.add_variable_group(bertini.VariableGroup([x, y]))

   names = {mptype: type(bertini.ZeroDimSolver(system, mptype=mptype)).__name__
            for mptype in ('double', 'multiple', 'adaptive')}
   assert names['double']   == 'ZeroDimSolverPowerSeriesDoublePrecision'
   assert names['multiple'] == 'ZeroDimSolverPowerSeriesFixedMultiplePrecision'
   assert names['adaptive'] == 'ZeroDimSolverPowerSeriesAdaptivePrecision'

Reading solutions: convert with ``complex()``
=============================================

The model changes the **type** of the numbers you get back. A double solve returns NumPy
``complex128``; a multiprecision solve returns arrays of :class:`bertini.complex_mp` (NumPy
arrays with ``dtype`` ``complex_mp``). A very portable habit is to convert
each coordinate with :func:`complex`:

.. testcode::

   bertini.random.set_random_seed(2)

   dbl = bertini.ZeroDimSolver(system, mptype='double'); dbl.solve()
   assert dbl.all_solutions()[0].dtype == np.complex128

   amp = bertini.ZeroDimSolver(system, mptype='adaptive'); amp.solve()
   assert str(amp.all_solutions()[0].dtype) == 'complex_mp'     # bertini.complex_mp

   # the same code reads either one:
   def to_python(solution):
       return np.array([complex(c) for c in solution])

   for solver in (dbl, amp):
       pts = sorted(tuple(np.round(to_python(s).real, 4)) for s in solver.all_solutions())
       assert pts == [(-0.7071, 0.7071), (0.7071, -0.7071)]

But it potentially loses precision.  You choose.  I don't know what the right thing to do is.

Setting the precision
=====================

A **fixed multiple** solve works at one precision everywhere: the tracker and every point it
produces carry the same number of digits, taken from the default precision when the solver is
constructed. There is nothing to set on the system. A ``System`` carries no precision of its own --
it evaluates at the precision of whatever point it is handed -- so the same ``system`` object serves
a 16-digit solve, a 40-digit solve and an adaptive solve unchanged:

.. testcode::

   bertini.default_precision(40)                  # 40 digits for this solve
   m = bertini.ZeroDimSolver(system, mptype='multiple')
   m.solve()
   assert len(m.all_solutions()) == 2
   assert m.all_solutions()[0][0].precision == 40  # the points carry the digits

   bertini.default_precision(30)                  # restore a modest default

An **adaptive** solve manages precision itself; its knobs live in the AMP config -- most usefully
``maximum_precision``, the ceiling above which a path is declared to have failed:

.. testcode::

   from bertini.tracking import AMPConfig
   amp = bertini.ZeroDimSolver(system, mptype='adaptive')
   assert amp.get_tracker().get_config(AMPConfig).maximum_precision == 300
   amp.get_tracker().update(maximum_precision=200)     # tighten the ceiling
   assert amp.get_tracker().get_config(AMPConfig).maximum_precision == 200

The config surface is the same shape across models
==================================================

Apart from the precision-specific tracker knobs (the ``amp`` config exists only for an adaptive
tracker), the configs are the same across precision models -- ``tolerances``, ``zero_dim``, stepping,
and so on are not precision-stamped. That is exactly why a settings bundle carries between models (see
:doc:`/tutorials/settings_and_precision/carrying_settings/index`):

.. testcode::

   shared = {'tolerances', 'zero_dim', 'post_processing', 'auto_retrack'}
   for mptype in ('double', 'multiple', 'adaptive'):
       names = set(bertini.ZeroDimSolver(system, mptype=mptype).config_names())
       assert shared <= names

A rule of thumb: reach for ``'adaptive'`` by default; drop to ``'double'`` for speed when the
problem is well-conditioned; use ``'multiple'`` when you want a fixed, known precision throughout.

Complete example
================

.. literalinclude:: precision_models.py
   :language: python
   :caption: precision_models.py
