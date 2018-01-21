Tracking to nonsingular endpoints
**********************************************

PyBertini works by setting up systems, setting up algorithms to use those systems, and doing something with the output.

Forming a system
=================


Let's make a couple of variables::

	x = pybertini.function_tree.symbol.Variable("x") #yes, you can make a variable not match its name...
	y = pybertini.function_tree.symbol.Variable("y")

Now, make a few symbolic expressions out of them::

	f = x**2 + y**2 -1
	g = x+y

There's no need to "set them equal to 0" -- expressions used as functions in a system in Bertini are taken to be equal to zero.  If you have an equality that's not zero, move one side to the other.

Let's make an empty system, then build into it::

	sys = pybertini.System()
	sys.add_function(f)
	sys.add_function(g)

``sys`` doesn't know its variables yet, so let's group them into an affine variable group [#]_, and stuff it into ``sys``::

	grp = pybertini.VariableGroup()
	grp.append(x)
	grp.append(y)
	sys.add_variable_group(grp)

Let's check that the degrees of our functions are correct::

	d = sys.degrees()
	assert(d[0]==2)
	assert(d[1]==1)

What happens if we add a non-polynomial function to our system?

::

	sys.add_function(x**-1)
	sys.add_function( pybertini.function_tree.sin(x) )
	d = sys.degrees()
	assert(d[2]==-1) # unsurprising, but actually a coincidence
	assert(d[3]==-1) # also -1.  anything non-polynomial is a negative number.  sin has no degree


Forming a homotopy
==================

A homotopy in Numerical Algebraic Geometry glues together a start system and a target system.  Above, we formed a target system, ``sys``.  Now, let's make a start system ``td``, and couple it to ``sys``.

The most basic, easiest to form and solve, start system is the Total Degree (TD) start system.  It is implemented as a first-class object in Bertini and PyBertini.  It takes in a polynomial system as its argument, and self-forms.::

	del sys #we mal-formed our system above by adding too many functions, and non-polynomial functions to it.
	# so, we start over
	sys = pybertini.System()
	sys.add_variable_group(grp)
	sys.add_function(f)
	sys.add_function(g)

	td = pybertini.TotalDegree(sys)

Wonderful, now we have an easy-to-solve system, the structure of which mirrors that of our target system.  Every start system comes with a method for generating its start points, by integer index.::

	# generates the 1th (0-based offsets in python) start point
	# at double precision
	td.start_point_d(1)

	# generate the 1th point at current multiple precision
	sp = td.start_point_mp(1)
	assert(pybertini.default_precision() == sp[1].precision())

Finally, we couple ``sys`` and ``td``::

	t = pybertini.function_tree.symbol.Variable("t")
	homotopy = (1-t)*sys + t*td
	homotopy.add_path_variable(t)

Now, we have the minimum theoretical ingredients for solving a polynomial system using Numerical Algebraic Geometry: a homotopy, a target system, and a start system.

Tracking a single path
======================

There are three basic trackers available in PyBertini:

#. Fixed double precision: ``pybertini.tracking.DoublePrecisionTracker``
#. Fixed multiple precision: ``pybertini.tracking.MultiplePrecisionTracker``
#. Adaptive precision: ``pybertini.tracking.AMPTracker``

Each brings its own advantages and disadvantages.  And, each has its ambient numeric type.

Let's use the adaptive one, since adaptivity is generally a good trait to have.  ``AMPTracker`` uses variable-precision vectors and matrices in its ambient work -- that is, you feed it multiprecisions, and get back multiprecisions.  Internally, it will use double precision when it can, and higher when it has to.

We associate a system with a tracker when we make it.  You cannot make a tracker without telling the tracker which system it will be tracking...

::

	tr = pybertini.tracking.AMPTracker(homotopy)
	tr.set_tolerance(1e-5) # track the path to 5 digits or so

	# adjust some stepping settings
	stepping = pybertini.tracking.config.SteppingConfig()
	stepping.max_step_size = pybertini.multiprec.rational(1,13)

	#then, set the config into the tracker.


Once we feel comfortable with the configs (of which there are many, see the book or elsewhere in this site, perhaps), we can track a path.

::

	result = pybertini.VectorXmp()
	tr.track_path(result, pybertini.multiprec.complex(1), pybertini.multiprec.complex(0), td.start_point_mp(0))

Let's generate a log of what was computed along the way, first making an observer, and then attaching it to the tracker.

::

	#make observer

	#attach

Re-running it, you should find the logfile ``bertini#.log``.

Using an endgame to compute singular endpoints
===============================================

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



A complete tracking of paths
=============================

::
	
	import pybertini

	x = pybertini.function_tree.symbol.Variable("x") #yes, you can make a variable not match its name...
	y = pybertini.function_tree.symbol.Variable("y")
	f = x**2 + y**2 -1
	g = x+y

	sys = pybertini.System()
	sys.add_function(f)
	sys.add_function(g)

	grp = pybertini.VariableGroup()
	grp.append(x)
	grp.append(y)
	sys.add_variable_group(grp)

	td = pybertini.start_system.TotalDegree(sys)

	t = pybertini.function_tree.symbol.Variable("t")
	homotopy = (1-t)*sys + t*td
	homotopy.add_path_variable(t)

	tr = pybertini.tracking.AMPTracker(homotopy)

	g = pybertini.tracking.observers.amp.GoryDetailLogger()

	tr.add_observer(g)
	tr.tracking_tolerance(1e-5) # track the path to 5 digits or so
	tr.infinite_truncation_tolerance(1e5)
	# tr.predictor(pybertini.tracking.Predictor.RK4)
	stepping = pybertini.tracking.config.SteppingConfig()
	# stepping.max_step_size = pybertini.multiprec.rational(1,13)

	results = []

	for ii in range(td.num_start_points()):
		results.append(pybertini.multiprec.Vector())
		tr.track_path(result=results[-1], start_time=pybertini.multiprec.complex(1), end_time=pybertini.multiprec.complex(0), start_point=td.start_point_mp(ii))

	tr.remove_observer(g)

Footnotes
---------

.. [#]  Affinely-grouped variables live together in the same complex space, :math:`\mathbb{C}^N`.  The alternative is projectively-grouped variables, which live in a copy of :math:`\mathbb{P}^N`.
