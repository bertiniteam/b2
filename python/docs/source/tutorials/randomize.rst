🎲 Randomize an overdetermined system
*************************************************

.. testsetup:: *

   import numpy as np
   import bertini
   from bertini import linalg

Most tutorials solve a **square** system -- as many equations as unknowns. But polynomial models
are often **overdetermined**: more equations than unknowns. The common zeros are still (at most)
isolated points, but you cannot hand such a system to a total-degree start system, which needs a
square target.

The fix is **randomization**: replace the :math:`N` functions with :math:`n` generic combinations
(:math:`n` = the number of unknowns). Every common zero of the original system is a zero of each
combination, so the genuine solutions survive -- but the square randomized system has *extra*
solutions too. The workflow is always: **randomize, solve the square system, then keep only the
points that also satisfy the original system.**

This tutorial does that twice, to show two different things: first an ordinary system in one block
of variables, then a system whose variables split into groups -- where the grouping cuts the path
count.

Part 1 -- one block of variables
=================================

Three curves in the plane: the unit circle, the line :math:`x = y`, and the vertical lines
:math:`2x^2 = 1`. Along :math:`x = y` the circle gives :math:`2x^2 = 1`, which is the third
equation, so all three meet at :math:`(\pm\tfrac{1}{\sqrt2}, \pm\tfrac{1}{\sqrt2})` -- two real
points you can read off by hand.

.. testcode::

   x, y = bertini.Variable('x'), bertini.Variable('y')

   original = bertini.System()
   original.add_variable_group(bertini.VariableGroup([x, y]))
   original.add_function(x*x + y*y - 1)      # the unit circle      (degree 2)
   original.add_function(x - y)              # the line x = y       (degree 1)
   original.add_function(2*x*x - 1)          # the lines x = ±1/√2  (degree 2)

   assert list(original.degrees()) == [2, 1, 2]   # three equations, two unknowns: overdetermined

:func:`~bertini.linalg.randomize` returns a **new** square system; the original is left untouched.
For a single variable group it sorts the functions by descending degree and forms :math:`R = [\,I
\mid C\,]` with a random block :math:`C`, so the randomized degrees are the :math:`n` **largest**
original degrees and the total-degree path count is their product -- here :math:`2 \times 2 = 4`,
the minimum a randomization can achieve.

.. testcode::

   randomized = linalg.randomize(original)

   assert randomized.num_functions() == 2          # squared up
   assert original.num_functions() == 3            # original unchanged
   assert sorted(randomized.degrees()) == [2, 2]

.. note::

   The descending sort is what keeps the count minimal: the degree-1 line is folded into the
   degree-2 rows rather than becoming a degree-1 target that a folded degree-2 function would
   inflate. The payoff grows with the degree spread -- for degrees :math:`(3, 2, 2)` randomization
   tracks :math:`3 \times 2 = 6` paths, where a degree-blind squaring (every row pushed to the
   maximum degree 3) would track :math:`3^2 = 9`. You may also pass your own exact matrix,
   ``linalg.randomize(original, R)``, in which case the functions are kept in their given order.

Solve the square system, then filter: evaluate the **original** system at each computed point and
keep the ones with a tiny residual.

.. testcode::

   bertini.random.set_random_seed(1)               # reproducible generic coefficients + gamma
   zd = bertini.nag_algorithm.ZeroDimSolver(randomized, endgame='cauchy', mptype='adaptive', startsystem='rootsofunity')
   zd.solve()
   solutions = zd.all_solutions()
   assert len(solutions) == 4                       # the two we want, plus two extraneous

   def satisfies_original(point):
       coords = np.array([complex(c) for c in point])
       return max(abs(complex(v)) for v in original.eval(coords)) < 1e-7

   true_solutions = [np.array([complex(c) for c in s]) for s in solutions if satisfies_original(s)]

   r = 1.0 / np.sqrt(2.0)
   found = sorted((round(p[0].real, 4), round(p[1].real, 4)) for p in true_solutions)
   assert found == sorted([(round(r, 4), round(r, 4)), (round(-r, 4), round(-r, 4))])

Two of the four points survive the filter: :math:`(\tfrac{1}{\sqrt2}, \tfrac{1}{\sqrt2})` and
:math:`(-\tfrac{1}{\sqrt2}, -\tfrac{1}{\sqrt2})` -- exactly where all three curves meet. The other
two solve the random combination but not the original equations, so the filter drops them.

Part 2 -- variables in separate groups
=======================================

Now something different: when the unknowns split into separate **variable groups**, randomization
works group by group and the *multihomogeneous* start system tracks far fewer paths.

The vehicle is a **bilinear** system -- each function is degree 1 in :math:`x` and degree 1 in
:math:`y` -- with :math:`x` and :math:`y` in their own groups. Three equations whose only common
solution is :math:`(1, 1)`: from :math:`x = y` and :math:`x + y = 2` we get :math:`x = y = 1`, and
indeed :math:`xy = 1`.

.. testcode::

   x, y = bertini.Variable('x'), bertini.Variable('y')

   original = bertini.System()
   original.add_variable_group(bertini.VariableGroup([x]))    # x is its own group
   original.add_variable_group(bertini.VariableGroup([y]))    # y is its own group
   original.add_function(x*y - 1)        # xy = 1
   original.add_function(x + y - 2)      # x + y = 2
   original.add_function(x - y)          # x = y

   randomized = linalg.randomize(original)
   assert randomized.num_functions() == 2

Solve with the **multihomogeneous** start system, which exploits the grouping. Two bilinear
equations on :math:`\mathbb{P}^1 \times \mathbb{P}^1` meet in the multihomogeneous Bézout number
of paths -- here :math:`2` (the coefficient of :math:`z_x z_y` in :math:`(z_x + z_y)^2`) -- whereas
a total-degree start would track :math:`2^2 = 4`. Same answer, half the work.

.. testcode::

   bertini.random.set_random_seed(3)
   zd = bertini.nag_algorithm.ZeroDimSolver(randomized, endgame='cauchy', mptype='adaptive', startsystem='mhom')
   zd.solve()

   def satisfies(point):
       coords = np.array([complex(c) for c in point])
       return max(abs(complex(v)) for v in original.eval(coords)) < 1e-7

   true_solutions = [np.array([complex(c) for c in s]) for s in zd.all_solutions() if satisfies(s)]

   assert len(true_solutions) == 1
   p = true_solutions[0]
   assert abs(p[0] - 1) < 1e-7 and abs(p[1] - 1) < 1e-7

One solution survives: :math:`(1, 1)`, exactly as read off by hand. The difference from Part 1 is
the **grouping** -- declaring ``x`` and ``y`` as separate variable groups let the multihomogeneous
start system see the bilinear structure and track fewer paths to the same answer.
