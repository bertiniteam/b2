🛷 Move a slice, hold the system fixed 
****************************************

.. testsetup:: *

   import numpy as np
   import bertini
   from bertini import linalg, nag_algorithm

Regeneration and witness-set work share a shape: most equations are **fixed** -- the polynomial
system, plus "below" linear slices that cut the dimension -- and only a small part **moves** along
the path variable. The fixed equations should be evaluated **once** per step: never duplicated,
never scaled by the path coefficient, never differentiated in :math:`t`.

:func:`~bertini.nag_algorithm.moving_homotopy` builds exactly that. You hand it the **fixed** system
and the two endpoints of the **moving** rows, and it returns

.. math::

   H \;=\; \bigl[\; \text{fixed's blocks} \;;\; (1-t)\,\text{end} + \gamma\,t\,\text{start} \;\bigr],

keeping the fixed equations as their own evaluation blocks and moving only the rest. Pair it with
:func:`~bertini.nag_algorithm.HomotopySolver` and the start points you already know.

Move one slice
==============

A fixed unit circle, sliced by a line that slides from the x-axis (:math:`y=0`) to the diagonal
(:math:`y-x=0`). At :math:`t=1` the slice is the x-axis, so the start points are the circle's
intersections with it, :math:`(\pm 1, 0)`; at :math:`t=0` the slice is the diagonal:

.. testcode::

    x, y = bertini.Variable('x'), bertini.Variable('y')

    fixed = bertini.System()
    fixed.add_variable_group(bertini.VariableGroup([x, y]))
    fixed.add_function(x*x + y*y - 1)                 # the fixed unit circle

    start_moving = bertini.System()                   # the slice at t=1: the x-axis y = 0
    start_moving.add_variable_group(bertini.VariableGroup([x, y]))
    start_moving.add_function(y)

    end_moving = bertini.System()                     # the slice at t=0: the diagonal y - x = 0
    end_moving.add_variable_group(bertini.VariableGroup([x, y]))
    end_moving.add_function(y - x)

.. testcode::

    gamma = linalg.coefficient(bertini.multiprec.Complex('0.6', '0.8'))   # off the real axis
    H = nag_algorithm.moving_homotopy(fixed, start_moving, end_moving, gamma=gamma)

    target = bertini.system.concatenate(fixed, end_moving)   # the t=0 system: circle + diagonal
    start_points = [np.array([bertini.multiprec.Complex('1'),  bertini.multiprec.Complex('0')]),
                    np.array([bertini.multiprec.Complex('-1'), bertini.multiprec.Complex('0')])]

    solver = nag_algorithm.HomotopySolver(H, start_points, target)
    solver.solve()
    roots = sorted((round(complex(s[0]).real, 4), round(complex(s[1]).real, 4))
                   for s in solver.all_solutions())

    r = round(1 / np.sqrt(2), 4)
    assert roots == sorted([(r, r), (-r, -r)])        # circle ∩ diagonal = (±1/√2, ±1/√2)

The circle never moved; only the slice did. The two start points slid along the circle to the two
diagonal intersections.

A static slice and a moving slice
=====================================

Now in three variables, the genuinely regeneration-flavored case: a fixed unit **sphere** (a
surface) cut to dimension zero by **two** slices -- a **static** one ``z = 0`` and a **moving** one
that slides ``y = 0`` :math:`\to` ``y - x = 0``. Two blocks are fixed (the sphere and the static
slice); only the moving slice carries :math:`t`:

.. testcode::

    x, y, z = bertini.Variable('x'), bertini.Variable('y'), bertini.Variable('z')

    fixed = bertini.System()
    fixed.add_variable_group(bertini.VariableGroup([x, y, z]))
    fixed.add_function(x*x + y*y + z*z - 1)                                # the sphere
    linalg.add_linear(fixed, np.array([[0, 0, 1]]), np.array([x, y, z]))   # static slice z = 0

    start_moving = bertini.System(); start_moving.add_variable_group(bertini.VariableGroup([x, y, z]))
    start_moving.add_function(y)                                           # moving slice at t=1
    end_moving = bertini.System(); end_moving.add_variable_group(bertini.VariableGroup([x, y, z]))
    end_moving.add_function(y - x)                                         # moving slice at t=0

.. testcode::

    H = nag_algorithm.moving_homotopy(fixed, start_moving, end_moving, gamma=gamma)
    target = bertini.system.concatenate(fixed, end_moving)
    start_points = [np.array([bertini.multiprec.Complex(str(a)), bertini.multiprec.Complex('0'),
                              bertini.multiprec.Complex('0')]) for a in (1, -1)]

    solver = nag_algorithm.HomotopySolver(H, start_points, target)
    solver.solve()
    roots = sorted((round(complex(s[0]).real, 4), round(complex(s[1]).real, 4), round(complex(s[2]).real, 4))
                   for s in solver.all_solutions())
    r = round(1 / np.sqrt(2), 4)
    assert roots == sorted([(r, r, 0.0), (-r, -r, 0.0)])

**The fixed system is left out of the motion.** Ask the homotopy for :math:`dH/dt`: the rows of the
two fixed blocks (the sphere and the static slice) are exactly zero -- they are evaluated once and
never differentiated as the slice moves -- while only the moving row is nonzero:

.. testcode::

    dHdt = H.eval_time_derivative(np.array([bertini.multiprec.Complex('0.3'),
                                            bertini.multiprec.Complex('0.4'),
                                            bertini.multiprec.Complex('0.5')]),
                                  bertini.multiprec.Complex('0.5'))
    assert abs(complex(dHdt[0])) == 0.0      # sphere row: out of dH/dt
    assert abs(complex(dHdt[1])) == 0.0      # static slice row: out of dH/dt
    assert abs(complex(dHdt[2])) > 0.0       # only the moving slice carries t

Deform a product of linears into a polynomial
=============================================

The regeneration "add a degree" step itself is the same construction with a different moving pair:
deform a **product of linears** into the actual polynomial, holding a slice fixed. Here the product
:math:`(x-1)(x+1)` deforms into the circle, with the slice ``y = 1/2`` static. At :math:`t=1` the
moving row is :math:`\gamma\,(x-1)(x+1)`, whose roots on the slice are :math:`(\pm 1, 1/2)`:

.. testcode::

    x, y = bertini.Variable('x'), bertini.Variable('y')

    fixed = bertini.System(); fixed.add_variable_group(bertini.VariableGroup([x, y]))
    linalg.add_linear(fixed, np.array([[0, 1]]), np.array([x, y]), ['-1/2'])    # static slice y = 1/2

    start_moving = bertini.System(); start_moving.add_variable_group(bertini.VariableGroup([x, y]))
    linalg.add_products_of_linears(start_moving, [[[1, 0, -1], [1, 0, 1]]])      # (x-1)(x+1), a structured block

    end_moving = bertini.System(); end_moving.add_variable_group(bertini.VariableGroup([x, y]))
    end_moving.add_function(x*x + y*y - 1)                                       # the polynomial

.. testcode::

    H = nag_algorithm.moving_homotopy(fixed, start_moving, end_moving, gamma=gamma)
    target = bertini.system.concatenate(fixed, end_moving)
    start_points = [np.array([bertini.multiprec.Complex(str(a)), bertini.multiprec.Complex('0.5')])
                    for a in (1, -1)]

    solver = nag_algorithm.HomotopySolver(H, start_points, target)
    solver.solve()
    roots = sorted((round(complex(s[0]).real, 4), round(complex(s[1]).real, 4))
                   for s in solver.all_solutions())
    s3 = round(np.sqrt(3) / 2, 4)
    assert roots == sorted([(s3, 0.5), (-s3, 0.5)])     # circle ∩ {y = 1/2}

The product-of-linears is a first-class structured block (it is never expanded into nodes), the
static slice stays put, and only the one regenerating row carries :math:`t`. That is the unit of
work a regeneration cascade repeats -- and at every step the fixed equations are evaluated once.
