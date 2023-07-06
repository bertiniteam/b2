🎮 Using an endgame to compute singular endpoints 
*********************************************************



Background
==============

Polynomial systems often have singular solutions.  In numerical algebraic geometry, we want to compute all solutions, even the challenging singular ones.  The normal method of homotopy continuation with straight-line tracking fails to compute such roots, because tracking to a place where the Jacobian is non-invertible using methods that require inverting the Jacobian is doomed to fail [#]_.  

So, if we can't track to a singular solution, but we still want to track to compute them, what are we to do?  We track around them, or near them, but not actually to them.  These methods are collectively called *endgames*, a term coined to evoke a sense of chess :cite:`morgan1990computing` :cite:`morgan1992computing` :cite:`morgan1992power`.  Thanks, Andrew Sommese, Charles Wampler, and Alexander Morgan, for everything you have given our community.

Endgames represent a way to finish a tracking of a path, when the endpoint is possibly singular.  Rather than track all the way to the endtime, you instead run an endgame that uses mathematical theory to compute the root.

Endgames in PyBertini
==========================

An endgame is a computational tool that one does in the final stage of a path track to a possibly singular root.  There are two implemented endgames in Bertini:

#. Power series (PSEG) -- uses `Hermite interpolation <https://en.wikipedia.org/wiki/Hermite_interpolation>`_ across a sequence of geometrically-spaced points (in time) to extrapolate to a target time :cite:`morgan1992power`. 
#. Cauchy (CauchyEG)-- uses `Cauchy's integral formula <https://en.wikipedia.org/wiki/Cauchy's_integral_formula>`_ in a sequence of circles about the root you are computing.  

Both try to compute the cycle number :math:`c` for the root.  In PSEG, :math:`c` is used as the degree of a Hermite interpolant used to extrapolate to 0.  In CauchyEG,  it is used for the number of cycles to walk before doing a trapezoid-rule integral.

Each is provided in the three precision modes, double, fixed multiple, and adaptive.  Since we are using the :class:`~pybertini.tracking.AMPTracker` in this tutorial, we will of course use the adaptive endgame.  I really like the Cauchy endgame, so we're in the land of the :class:`~pybertini.endgame.AMPCauchyEG`.


Example
----------


Form a system
~~~~~~~~~~~~~~~~

The Griewank-Osborne system has one multiplicity-three singular solution at the origin :cite:`griewank1983analysis`.  It comes pre-built for us as part of Bertini2's C++ core, and is accessible by peeking into the `precon` module.  

.. todo::

	expose the precon namespace.  it's a 1-hour task, and danielle 😈 should do it.

Let's build it from scratch, for the practice.

:: 

	import pybertini

	gw = pybertini.System()

	x = pybertini.Variable("x")
	y = pybertini.Variable("y")

	vg = pybertini.VariableGroup()
	vg.append(x)
	vg.append(y)
	gw.add_variable_group(vg)

	gw.add_function(pybertini.multiprec.Rational(29,16)*x**3 - 2*x*y)
	gw.add_function(y - x**2)


Form a start system and homotopy 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next, we make the total degree start system for `gw`, and couple it using the gamma trick :cite:`morgan1987homotopy` and a path variable.

::

	t = pybertini.Variable('t')
	td = pybertini.system.start_system.TotalDegree(gw)
	gamma = pybertini.function_tree.symbol.Rational.rand()
	hom = (1-t)*gw + t*gamma*td
	hom.add_path_variable(t)



🛤 Track to the endgame boundary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make a tracker.  I use adaptive precision a lot, so we'll roll with that.  There are also double and fixed-multiple versions.  See the other tutorials or the detailed documentation.

::

	tr = pybertini.tracking.AMPTracker(hom)

	start_time = pybertini.multiprec.Complex("1")
	eg_boundary = pybertini.multiprec.Complex("0.1")

	midpath_points = [None]*td.num_start_points()
	for ii in range(td.num_start_points()):
		midpath_points[ii] = pybertini.multiprec.Vector()
		code = tr.track_path(result=midpath_points[ii], start_time=start_time, end_time=eg_boundary, start_point=td.start_point_mp(ii))
		if code != pybertini.tracking.SuccessCode.Success:
			print('uh oh, tracking a path before the endgame boundary failed, successcode ' + code)




🎮 Use the endgame
~~~~~~~~~~~~~~~~~~~~


To make an endgame, we need to feed it the tracker that is used to run.  There are also config structs to play with, that control the way things are computed.

::

	eg = pybertini.endgame.AMPCauchyEG(tr)

	# make an observer to be able to see what's going on inside
	ob = pybertini.endgame.observers.amp_cauchy.GoryDetailLogger()

	eg.add_observer(ob)

Since the endgame hasn't been run yet things are empty and default::

	assert(eg.cycle_number()==0)
	assert(eg.final_approximation()==np.array([]))

The endgames are used by invoking ``run``, feeding it the point we are tracking on, the time we are at, and the time we want to track to. ::

	final_points = []


	target_time = pybertini.multiprec.Complex(0)
	codes = []
	for ii in range(td.num_start_points()):
		eg_boundary.precision( midpath_points[ii][0].precision())
		target_time.precision( midpath_points[ii][0].precision())
		print('before {} {} {}'.format(eg_boundary.precision(), target_time.precision(), midpath_points[ii][0].precision()))
		codes.append(eg.run(start_time=eg_boundary, target_time=target_time, start_point=midpath_points[ii]))
		print('path {} -- code {}'.format(ii,codes[-1]))
		print(eg.final_approximation())
		# final_points.append(copy.deep_copy(eg.final_approximation()))
		print('after {} {} {}'.format(eg_boundary.precision(), target_time.precision(), midpath_points[ii][0].precision()))

.. todo::

	the endgame returns its `final_approximation` by reference, so capturing its value into a list makes many references to this internal variable, not copies of the point.  so, one should take deepcopy's of the vector, but they are not currently pickleable due to the complex multiprecision class.  an issue has been filed (#148) and this issue will be solved shortly (danielle, 20180227)

Conclusion
============


Using a singular endgame, we can compute singular endpoints of homotopy paths.  What an age to live in!  🌌





📚 Further reading
========================

The following three papers (cited above) laid the foundation for endgames and computation of singular endpoints:

* Computing singular solutions to nonlinear analytic systems :cite:`morgan1990computing`
* Computing singular solutions to polynomial systems :cite:`morgan1992computing` 
* A power series method for computing singular solutions to nonlinear analytic systems :cite:`morgan1992power`.

👣 Footnotes
-------------

.. [#]  No, we don't actually invert the Jacobian in practice while solving the Davidenko differential equation, but numerical issues exist no matter which method you use to solve the system.



