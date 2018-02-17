🎮 Using an endgame to compute singular endpoints 
*********************************************************



Background
==============

Polynomial systems often have singular solutions.  In numerical algebraic geometry, we want to compute all solutions, even the challenging singular ones.  The normal method of homotopy continuation with straight-line tracking fails to compute such roots, because tracking to a place where the Jacobian is non-invertible using methods that require inverting the Jacobian is doomed to fail [#]_.  

So, if we can't track to a singular solution, but we still want to track to compute them, what are we to do?  We track around them, or near them, but not actually to them.  These methods are collectively called *endgames*, a term coined to evoke a sense of chess :cite:`morgan1990computing` :cite:`morgan1992computing` :cite:`morgan1992power`.  Thanks, Andrew Sommese, Charles Wampler, and Alexander Morgan, for everything you have given our community.




Application
=============

There are two implemented endgames in Bertini:

#. Power series -- uses `Hermite interpolation <https://en.wikipedia.org/wiki/Hermite_interpolation>`_ across a sequence of geometrically-spaced points (in time) to extrapolate to a target time.
#. Cauchy -- uses `Cauchy's integral formula <https://en.wikipedia.org/wiki/Cauchy's_integral_formula>`_

Each is provided in the three precision modes, double, fixed multiple, and adaptive.  Since we are using the adaptive tracker in this tutorial, we will of course use the adaptive endgame.  I really like the Cauchy endgame, so we're in the land of the ``pybertini.endgame.AMPCauchyEG``.

To make an endgame, we need to feed it the tracker that is used to run.  There are also config structs to play with, that control the way things are computed.

::

	eg = pybertini.endgame.AMPCauchyEG(tr)

Since the endgame hasn't been run yet things are empty and default::

	assert(eg.cycle_number()==0)
	assert(eg.final_approximation()==pybertini.VectorXmp())

The endgames are used by invoking ``run``, feeding it the point we are tracking on, the time we are at, and the time we want to track to.












Further reading
=================

The following three papers laid the foundation for endgames and computation of singular endpoints:

* Computing singular solutions to nonlinear analytic systems :cite:`morgan1990computing`
* Computing singular solutions to polynomial systems :cite:`morgan1992computing` 
* A power series method for computing singular solutions to nonlinear analytic systems :cite:`morgan1992power`.

Footnotes
---------

.. [#]  No, we don't actually invert the Jacobian in practice while solving the Davidenko differential equation, but numerical issues exist no matter which method you use to solve the system.



