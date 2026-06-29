🎯 Solving a system with ZeroDimSolver
***************************************************

``ZeroDimSolver`` is the algorithm for the most common task: given a polynomial system, find all of
its **isolated complex solutions**.  You hand it a system; it forms a start system and a homotopy,
tracks the paths, and classifies the endpoints.  (If you instead already have a homotopy and a list
of start points -- the parameter-homotopy workflow -- reach for :func:`~bertini.nag_algorithm.HomotopySolver`.)

This tutorial shows ``ZeroDimSolver`` on the three shapes of input you will meet: a well-posed
**square** system, an **over-determined** system (more equations than unknowns), and an
**under-determined** one (fewer).  The point is that ``ZeroDimSolver`` does the right thing -- and
tells you clearly -- in each case.

.. testsetup:: *

   import bertini

A pleasant square system
========================

A *square* system has exactly one equation per variable, and generically finitely many isolated
solutions.  Take the unit circle meeting the line :math:`x = y`:

.. testcode::

    import bertini
    from bertini.nag_algorithm import ZeroDimSolver

    x, y = bertini.Variable('x'), bertini.Variable('y')

    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([x, y]))
    sys.add_function(x**2 + y**2 - 1)     # the unit circle
    sys.add_function(x - y)               # the line x = y

    solver = ZeroDimSolver(sys, mptype='adaptive')
    solver.solve()

There are two intersection points, :math:`\pm(1/\sqrt2,\, 1/\sqrt2)`:

.. testcode::

    assert solver.was_randomized() is False    # already square -- nothing to do
    assert len(solver.finite_solutions()) == 2

An over-determined system: extra equations, extraneous solutions filtered out
=============================================================================

An *over-determined* system has **more equations than variables**.  It can still have isolated
solutions -- the points where *all* the equations vanish at once -- but a start system needs a
square target to track.  ``ZeroDimSolver`` handles this for you: it **randomizes** the system down
to square (replacing the :math:`N` equations by :math:`n` generic combinations of them), solves the
square system, and then **filters out the extraneous solutions** that the squaring introduces by
re-checking each one against your original equations.

Three equations in two variables, whose only common solutions are :math:`(1, 1)` and
:math:`(-1, -1)`:

.. testcode::

    x, y = bertini.Variable('x'), bertini.Variable('y')
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([x, y]))
    sys.add_function(x**2 + y**2 - 2)
    sys.add_function(x - y)
    sys.add_function(x*y - 1)

    solver = ZeroDimSolver(sys, mptype='adaptive')
    assert solver.was_randomized() is True     # ZeroDimSolver squared it up for you
    solver.solve()

The squared system has more endpoints than there are genuine solutions (the extras satisfy the
random combinations but not the original equations).  ``all_solutions()`` shows them all;
``solutions()`` returns only the genuine roots, and the extraneous ones -- flagged
``is_nonsolution`` in the metadata -- are surfaced separately by ``nonsolutions()``:

.. testcode::

    n_all = len(solver.all_solutions())        # the squared system's full path count
    assert n_all > 2
    assert len(solver.solutions()) == 2        # only the two true common solutions
    assert len(solver.nonsolutions()) == n_all - 2   # the extraneous, junk points

``solutions()`` is the general getter: by default the finite genuine solutions, with keywords to
filter by category -- ``solutions(real=False)`` for the complex ones, ``solutions(singular=False)``
for the nonsingular ones, ``solutions(infinite=True)`` to also include the at-infinity endpoints,
``solutions(nonsolution=True)`` to opt the junk back in.  You can inspect the exact combinations the
squaring used with ``solver.randomization_matrix()``.

An under-determined system: a helpful refusal
=============================================

An *under-determined* system has **fewer equations than variables**, so its solution set is
**positive-dimensional** -- a curve, a surface, ... -- not a finite set of points.  ``ZeroDimSolver``
computes isolated solutions only, so it declines, with an error that says why rather than quietly
returning garbage:

.. testcode::

    x, y = bertini.Variable('x'), bertini.Variable('y')
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([x, y]))
    sys.add_function(x**2 + y**2 - 1)          # one equation, two variables: a whole curve

    try:
        ZeroDimSolver(sys, mptype='adaptive')
        raise AssertionError("expected ZeroDimSolver to refuse an under-determined system")
    except RuntimeError as e:
        assert 'under-determined' in str(e)

The same refusal protects you from a subtler case: a system that is square *by equation count* but
whose solution set is still positive-dimensional (for instance a repeated equation).  ``ZeroDimSolver``
checks the Jacobian rank at a generic point and raises if the solution set cannot be
zero-dimensional.
