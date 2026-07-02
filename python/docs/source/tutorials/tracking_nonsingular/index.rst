🛤 Tracking to nonsingular endpoints 
**********************************************

.. testsetup:: *

   import bertini

Bertini works by setting up systems, setting up algorithms to use those systems, and doing something with the output.

Forming a system
=================

First, gain access to bertini:

.. testcode::

    import bertini
    import numpy as np

Let's make a couple of :class:`~bertini.symbolics.Variable`'s:

.. testcode::

	x = bertini.symbolics.Variable("x") #yes, you can make a variable not match its name...
	y = bertini.symbolics.Variable("y")

Now, make a few symbolic expressions out of them:

.. testcode::

	f = x**2 + y**2 -1  # ** is exponentiation in Python.
	g = x+y

There's no need to "set them equal to 0" -- expressions used as functions in a system in Bertini are taken to be equal to zero.  If you have an equality that's not zero, move one side to the other.

Let's make an empty :class:`~bertini.system.System`, then build into it:

.. testcode::

	sys = bertini.System()
	sys.add_function(f)
	sys.add_function(g)

``sys`` doesn't know its variables yet, so let's group them into an affine :class:`~bertini.VariableGroup` [#]_, and stuff it into ``sys``:

.. testcode::

	grp = bertini.VariableGroup()
	grp.append(x)
	grp.append(y)
	sys.add_variable_group(grp)

Let's check that the degrees of our functions are correct:

.. testcode::

	d = sys.degrees()
	assert(d[0]==2)  # f is degree 2 (highest power in any term is 2)
	assert(d[1]==1)  # g is degree 1 (highest power in any term is 2)


Aside -- a brief exploration into non-algebraic things
---------------------------------------------------------


What happens if we add a non-polynomial function to our system?

.. testcode::

	sys.add_function(x**-1)  # happily accepts a non-polynomial function.
	sys.add_function(bertini.symbolics.sin(x) )
	d = sys.degrees()
	assert(d[2]==-1) # unsurprising, but actually a coincidence
	assert(d[3]==-1) # also -1.  anything non-polynomial is a negative number.  
	# sin has no well-defined degree

	# bertini uses negative degree to indicate non-polynomial


correcting our system -- a return to algebraicness
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

We can indeed do homotopy continuation with a non-algebraic systems.  What we cannot do is form a start system that we can guarantee will track to all solutions of the target system.  (because of things like :math:`\sin(x)` having infinitely many solutions, etc)

.. testcode::

	del sys #we mal-formed our system above by adding too many functions, and non-polynomial functions to it.
	# so, we start over
	sys = bertini.System()
	sys.add_variable_group(grp)
	sys.add_function(f)
	sys.add_function(g)



Forming a start system
=========================

To solve our algebraic system ``sys``, we need a corresponding start system -- one with related structure, but that is actually solvable without too much trouble.  Bertini2 has several implemented options.  The most basic (easiest to form and solve) start system is the Total Degree (TD) start system.  It is implemented as a first-class object in Bertini.  It takes in a polynomial system as its argument, and self-forms.


Above, we formed a target system, ``sys``.  Now, let's make a start system ``td``.  Later, we will couple it to ``sys``.
It's trivial to make a total degree start system (:class:`~bertini.system.start_system.TotalDegreeLinearProduct`):

.. testcode::

	td = bertini.system.start_system.TotalDegreeLinearProduct(sys)

Note that you have to pass in the target system into the constructor of the total degree, or you get an error.


Wonderful, now we have an easy-to-solve system ``td``, the structure of which mirrors that of our target system.  Every start system comes with a method ``start_point_*`` for generating its start points, by integer index.

.. testcode::

	# generate the 1th (0-based offsets in python) start point
	sp_d = td.start_point_d(1)# at double precision

	sp_mp = td.start_point_mp(1) # generate the 1th point at current default multiple precision
	assert(bertini.default_precision() == sp_mp[1].precision)


Forming a homotopy
==================


We turn next to the act of path tracking.  This is the core computational method of numerical algebraic geometry, and it requires a continuous deformation between systems, called a "homotopy".  

A homotopy in Numerical Algebraic Geometry glues together a start system and a target system, such that we can later "continue" from one into the other.   Observe:


We couple ``sys`` and ``td``:

.. testcode::

	t = bertini.Variable("t")     # make a path variable
	homotopy = (1-t)*sys + t*td     # glue
	homotopy.add_path_variable(t)   # indicate the path var

Now, we have the minimum theoretical ingredients for solving a polynomial system using Numerical Algebraic Geometry: 

#. a homotopy ``homotopy``, 
#. a target system ``sys``, 
#. and a start system ``td``.

as well as a few other incidentals which will be implicitly used, such as a path variable ``t``.


Tracking a single path
======================

There are three basic trackers available in Bertini 2:


#. Fixed double precision: :class:`~bertini.DoublePrecisionTracker`
#. Fixed multiple precision: :class:`~bertini.MultiplePrecisionTracker`
#. Adaptive precision: :class:`~bertini.AMPTracker`

Each brings its own advantages and disadvantages.  And, each has its ambient numeric type.

Let's use the adaptive one, since adaptivity is generally a good trait to have.  ``AMPTracker`` uses variable-precision vectors and matrices in its ambient work -- that is, you feed it multiprecisions, and get back multiprecisions.  Internally, it will use double precision when it can, and higher when it has to.

We associate a system with a tracker when we make it.  You cannot make a tracker without telling the tracker which system it will be tracking...

.. testcode::

	tr = bertini.AMPTracker(homotopy)
	tr.tracking_tolerance(1e-5) # track the path to 5 digits or so

	# adjust some stepping settings
	stepping = bertini.tracking.SteppingConfig()
	stepping.max_step_size = bertini.multiprec.real_mp(1)/bertini.multiprec.real_mp(13)

	#then, set the config into the tracker.
	tr.set_stepping(stepping)


Once we feel comfortable with the configs (of which there are many, see the book or elsewhere in this site, perhaps), we can track a path.

.. testcode::

	result = np.zeros((2,), dtype=bertini.multiprec.complex_mp)
	tr.track_path(result,bertini.multiprec.complex_mp(1),bertini.multiprec.complex_mp(0), td.start_point_mp(0))

Logging to inspect the path that was tracked
---------------------------------------------


Let's generate a log of what was computed along the way, first making an :mod:`observer <bertini.tracking.observers>`, and then attaching it to the tracker.

.. testcode::

	#make observer
	g = bertini.tracking.observers.amp.GoryDetailLogger()
	
	#attach
	tr.add_observer(g)

Re-running it, you should find a ton of stuff printed to the screen.

.. testcode::

	tr.track_path(result,bertini.multiprec.complex_mp(1),bertini.multiprec.complex_mp(0), td.start_point_mp(0))

If you are going to keep tracking, but want to turn off the logging, remove the observer.

.. testcode::

	tr.remove_observer(g)


A complete tracking of paths
=============================


Now that we've tracked a single path, you might want to loop over all start points.  Awesome!  The next blob takes all the above, and puts it into a single blob.  Enjoy!


.. testcode:: tracking_nonsingular_main
	
	import bertini
	import numpy as np

	x = bertini.symbolics.Variable("x") #yes, you can make a variable not match its name...
	y = bertini.symbolics.Variable("y")
	f = x**2 + y**2 -1
	g = x+y

	sys = bertini.System()
	sys.add_function(f)
	sys.add_function(g)

	grp = bertini.VariableGroup()
	grp.append(x)
	grp.append(y)
	sys.add_variable_group(grp)

	td = bertini.system.start_system.TotalDegreeLinearProduct(sys)

	t = bertini.Variable("t")
	homotopy = (1-t)*sys + t*td
	homotopy.add_path_variable(t)

	tr = bertini.AMPTracker(homotopy)

	#commented out for screen-saving.
	#g = bertini.tracking.observers.amp.GoryDetailLogger()
	#tr.add_observer(g)  
	# one could alsobertini.logging.init() and set a file name, 
	# so it gets piped there instead of wherever Boost.Log goes by default.

	tr.tracking_tolerance(1e-5) # track the path to 5 digits or so
	tr.infinite_truncation_tolerance(1e5)
	tr.predictor(bertini.Predictor.RK4)
	stepping = bertini.tracking.SteppingConfig()
	stepping.max_step_size = bertini.multiprec.real_mp(1)/bertini.multiprec.real_mp(13)

	# set the config into the tracker
	tr.set_stepping(stepping)

	results = [] # make an empty list into which to put the results
	expected_code = bertini.SuccessCode.Success
	codes = []
	for ii in range(td.num_start_points()):
		results.append(np.zeros((2,),dtype=bertini.multiprec.complex_mp))
		codes.append(tr.track_path(result=results[-1], start_time=bertini.multiprec.complex_mp(1), end_time=bertini.multiprec.complex_mp(0), start_point=td.start_point_mp(ii)))

	#tr.remove_observer(g)

	# the tracked endpoints are now in the list ``results``
	print(codes == [expected_code]*2)

.. testoutput:: tracking_nonsingular_main

	True





Footnotes
---------

.. [#]  Affinely-grouped variables live together in the same complex space, :math:`\mathbb{C}^N`.  The alternative is projectively-grouped variables, which live in a copy of :math:`\mathbb{P}^N`.

Complete example
================

The whole tutorial as one runnable script -- assemble nothing, just run it:

.. literalinclude:: tracking_nonsingular.py
   :language: python
   :caption: tracking_nonsingular.py
