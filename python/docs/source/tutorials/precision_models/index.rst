🔬 Precision models: double, multiple, adaptive
*************************************************

.. testsetup:: *

   import numpy as np
   import bertini

Bertini can track in three precision models, and ``ZeroDimSolver`` selects between them with ``mptype``:

* ``'double'`` -- hardware double precision (≈16 digits). Fastest; fine when the paths are
  well-conditioned.
* ``'multiple'`` -- a fixed multiprecision, the same number of digits throughout.
* ``'adaptive'`` -- adaptive multiprecision (AMP, the default): the precision rises and falls along
  each path as the conditioning demands. The most robust, and what makes a near-singular path solvable.

This tutorial is about how the three models differ *in the interface* -- the types you get back and
how you set the precision. For *why and when* you reach for adaptive precision (the conditioning
story), see :doc:`/tutorials/precision_matters/index`.

.. testcode::

   x, y = bertini.Variable('x'), bertini.Variable('y')
   system = bertini.System()
   system.add_function(x*x + y*y - 1)
   system.add_function(x + y)
   system.add_variable_group(bertini.VariableGroup([x, y]))

   names = {mptype: type(bertini.nag_algorithm.ZeroDimSolver(system, mptype=mptype)).__name__
            for mptype in ('double', 'multiple', 'adaptive')}
   assert names['double']   == 'ZeroDimSolverCauchyDoublePrecision'
   assert names['multiple'] == 'ZeroDimSolverCauchyFixedMultiplePrecision'
   assert names['adaptive'] == 'ZeroDimSolverCauchyAdaptivePrecision'

   # The solver type encodes the endgame and the precision model, but NOT the start system -- that is
   # a construction choice now (the algorithm holds the start system polymorphically), so it no longer
   # appears as a suffix on the class name.

Reading solutions: convert with ``complex()``
=============================================

The model changes the **type** of the numbers you get back. A double solve returns NumPy
``complex128``; a multiprecision solve returns arrays of :class:`bertini.multiprec.Complex` (NumPy
arrays with ``dtype`` ``Complex``). The portable habit -- works in every model -- is to convert each
coordinate with :func:`complex`:

.. testcode::

   bertini.random.set_random_seed(2)

   dbl = bertini.nag_algorithm.ZeroDimSolver(system, mptype='double'); dbl.solve()
   assert dbl.all_solutions()[0].dtype == np.complex128

   amp = bertini.nag_algorithm.ZeroDimSolver(system, mptype='adaptive'); amp.solve()
   assert str(amp.all_solutions()[0].dtype) == 'Complex'        # bertini.multiprec.Complex

   # the same code reads either one:
   def to_python(solution):
       return np.array([complex(c) for c in solution])

   for solver in (dbl, amp):
       pts = sorted(tuple(np.round(to_python(s).real, 4)) for s in solver.all_solutions())
       assert pts == [(-0.7071, 0.7071), (0.7071, -0.7071)]

Setting the precision
=====================

A **fixed multiple** solve works at one precision *everywhere*: the system, the tracker, and the
points must all agree. So set the default precision and lift the system to it before constructing the
solver:

.. testcode::

   bertini.default_precision(40)                  # 40 digits for this solve
   system.precision(40)                           # the system must match
   m = bertini.nag_algorithm.ZeroDimSolver(system, mptype='multiple')
   m.solve()
   assert len(m.all_solutions()) == 2

   bertini.default_precision(30)                  # restore a modest default
   system.precision(30)

(Mismatch is reported, not silently tolerated: constructing a multiple-precision solve whose system
sits at a different precision raises rather than tracking at the wrong precision.)

An **adaptive** solve manages precision itself; its knobs live in the AMP config -- most usefully
``maximum_precision``, the ceiling above which a path is declared to have failed:

.. testcode::

   from bertini.tracking import AMPConfig
   amp = bertini.nag_algorithm.ZeroDimSolver(system, mptype='adaptive')
   assert amp.get_tracker().get_config(AMPConfig).maximum_precision == 300
   amp.get_tracker().update(maximum_precision=200)     # tighten the ceiling
   assert amp.get_tracker().get_config(AMPConfig).maximum_precision == 200

The config surface is the same shape across models
==================================================

Apart from the precision-specific tracker knobs (the ``amp`` config exists only for an adaptive
tracker), the configs are the same across precision models -- ``tolerances``, ``zero_dim``, stepping,
and so on are not precision-stamped. That is exactly why a settings bundle carries between models (see
:doc:`/tutorials/carrying_settings/index`):

.. testcode::

   shared = {'tolerances', 'zero_dim', 'post_processing', 'auto_retrack'}
   for mptype in ('double', 'multiple', 'adaptive'):
       names = set(bertini.nag_algorithm.ZeroDimSolver(system, mptype=mptype).config_names())
       assert shared <= names

So the rule of thumb: reach for ``'adaptive'`` by default; drop to ``'double'`` for speed when the
problem is well-conditioned; use ``'multiple'`` when you want a fixed, known precision throughout. And
read every solution through :func:`complex` so your code does not care which model produced it.
