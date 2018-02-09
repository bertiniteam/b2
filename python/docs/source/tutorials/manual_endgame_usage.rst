Using an endgame to compute singular endpoints
***************************************************



Background
==============




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



