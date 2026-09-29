🧱 Build a system in pieces with clone and concatenate
*********************************************************

.. testsetup:: *

   import numpy as np
   import bertini

Sometimes a polynomial system arrives in **pieces**. One block of equations describes some
geometry; another block adds constraints. It is natural to build each block as its own
:class:`~bertini.system.System` and then stack them into one square system to solve.

:func:`~bertini.system.concatenate` does exactly that: it appends the second system's functions
onto a copy of the first, producing a new system. The two blocks share the same unknowns, so the
result is one system over those unknowns with all the functions together. This tutorial stacks two
**four-variable, two-function** blocks into a square :math:`4 \times 4` system and solves it -- and
shows how :func:`~bertini.system.clone` lets you declare the variables once and reuse them for every
block.

The pieces
==========

Four unknowns :math:`x_1, x_2, x_3, x_4`. The first block is some **geometry** -- two equations
relating the unknowns -- and the second block adds two **constraints**. Each block is a perfectly
good :class:`~bertini.system.System` on its own; neither is square (two equations, four unknowns),
so neither can be solved alone.

.. testcode::

   x1, x2, x3, x4 = (bertini.Variable(n) for n in ('x1', 'x2', 'x3', 'x4'))
   group = bertini.VariableGroup([x1, x2, x3, x4])

   geometry = bertini.System()
   geometry.add_variable_group(group)
   geometry.add_function(x1*x1 - x2)        # x2 = x1²
   geometry.add_function(x2*x2 - x3)        # x3 = x2²

   constraints = bertini.System()
   constraints.add_variable_group(group)
   constraints.add_function(x3 - x4)        # x4 = x3
   constraints.add_function(x4 - x1)        # x1 = x4

   assert geometry.num_functions() == 2 and geometry.num_variables() == 4
   assert constraints.num_functions() == 2 and constraints.num_variables() == 4

Both blocks are written over the **same** variables. That is the only requirement for concatenation:
the two systems must share their variable ordering. Reusing the same ``group`` (and the same
``Variable`` objects) guarantees it -- and even rebuilding the variables by name would, since
variables in Bertini are canonical by name.

Stack them and solve
=====================

:func:`~bertini.system.concatenate` returns a **new** square system; the two blocks are left
untouched.

.. testcode::

   system = bertini.system.concatenate(geometry, constraints)

   assert system.num_functions() == 4 and system.num_variables() == 4
   assert geometry.num_functions() == 2        # inputs unchanged
   assert constraints.num_functions() == 2
   assert system.degrees() == [2, 2, 1, 1]

The combined system is square, so a total-degree start system can solve it. Its path count is the
product of the degrees, :math:`2 \times 2 \times 1 \times 1 = 4`.

.. testcode::

   bertini.random.set_random_seed(1)           # reproducible start system + gamma
   zd = bertini.ZeroDimSolver(system, endgame='cauchy', mptype='adaptive', startsystem='binomial')
   zd.solve()

   solutions = [np.array([complex(c) for c in s]) for s in zd.all_solutions()]
   assert len(solutions) == 4

   # every computed point satisfies the concatenated system
   assert all(max(abs(complex(v)) for v in system.eval(s)) < 1e-7 for s in solutions)

By hand: the four equations force :math:`x_4 = x_1`, :math:`x_3 = x_4 = x_1`, :math:`x_2 = x_1^2`,
and :math:`x_2^2 = x_3 = x_1`, so :math:`x_1^4 = x_1`. That gives :math:`x_1 \in \{0, 1, \omega,
\omega^2\}` (the cube roots of unity, and zero) -- four solutions, of which two are real.

.. testcode::

   real = [s for s in solutions if max(abs(c.imag) for c in s) < 1e-7]
   found = sorted(tuple(round(c.real, 4) for c in s) for s in real)
   assert found == [(0.0, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0, 1.0)]

Declare the variables once, with clone
=======================================

Above, each block re-declared the variable group. When you are assembling several blocks, it is
tidier to set up the variables **once** and :func:`~bertini.system.clone` that setup for each block.

A clone shares the original's variables (and gets its own evaluation memory), so blocks cloned from
a common base are guaranteed to share variable ordering -- exactly what concatenation needs. Adding
functions to one clone does not touch the others.

.. testcode::

   base = bertini.System()                     # variable setup only -- no functions yet
   base.add_variable_group(bertini.VariableGroup([x1, x2, x3, x4]))

   geometry = bertini.system.clone(base)
   geometry.add_function(x1*x1 - x2)
   geometry.add_function(x2*x2 - x3)

   constraints = bertini.system.clone(base)
   constraints.add_function(x3 - x4)
   constraints.add_function(x4 - x1)

   system = bertini.system.concatenate(geometry, constraints)
   assert system.num_functions() == 4

   bertini.random.set_random_seed(1)
   zd = bertini.ZeroDimSolver(system, endgame='cauchy', mptype='adaptive', startsystem='binomial')
   zd.solve()
   assert len(zd.all_solutions()) == 4

Same four solutions, assembled from cloned blocks. The pattern -- **set up variables once, clone per
block, concatenate, solve** -- scales to as many blocks as your problem comes in.

.. note::

   What makes the blocks line up is that they are over the same variables, and a variable is its
   name: the clone's :math:`x_1` *is* the base's :math:`x_1`. The same holds for a system copied with
   ``copy.deepcopy`` or sent through :mod:`pickle` (to another process, say), so those combine with
   the original too. ``clone`` is the direct way, since it has no archive to go through.

Two systems are combined only when they are over the same variables, in the same groups and the same
order. A point is a list of numbers, and the order of the variables is what says which number is
which, so systems over :math:`(x, y)` and :math:`(y, x)` are refused, not quietly matched up by
position:

.. testcode::

   u, v = bertini.Variable('u'), bertini.Variable('v')

   one_way = bertini.System()
   one_way.add_variable_group(bertini.VariableGroup([u, v]))
   one_way.add_function(u - 2)

   other_way = bertini.System()
   other_way.add_variable_group(bertini.VariableGroup([v, u]))
   other_way.add_function(v - 3)

   try:
       bertini.system.concatenate(one_way, other_way)
   except RuntimeError as refusal:
       print('same variables in a different order' in str(refusal))

.. testoutput::

   True

A system that declares no variables at all takes the variables of the one it is combined with.

Complete example
================

.. literalinclude:: concatenate_systems.py
   :language: python
   :caption: concatenate_systems.py
