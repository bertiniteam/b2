Tutorial - Evaluation of cyclic-:math:`n` polynomials
*******************************************************

Bertini is software for algebraic geometry.  This means we work with systems of polynomials, a critical component of which is system and function evaluation.

Bertini2 allows us to set up many kinds of functions, and thus systems, by exploting operator overloading.

Make some symbols
==================

Let's start by making some variables, programmatically [1]_.  

::

	import pybertini
	import numpy

	num_vars = 10
	x = [None] * num_vars
	for ii in range(num_vars):
	    x[ii] = pybertini.Variable('x' + str(ii))

Huzzah, we have `num_vars` variables!  This was hard to do in Bertini 1's classic style input files.  Now we can do it directly! 🎯

Write a function to produce the cyclic :math:`n` polynomials :cite:`cyclic_n`.

::

	def cyclic(vars):
	    n = len(vars)
	    f = [None] * len(vars)
	    y = []
	    for ii in range(2):
	        for x in vars:
	            y.append(x)
	        
	    for ii in range(n):
	        f[ii] = numpy.sum( [numpy.prod(y[jj:jj+ii+1]) for jj in range(n)] ) 
	    
	    # the last one is minus one
	    f[-1] = f[-1]-1
	    return f

Now we will make a System, and put the cyclic polynomials into it.

::

	sys = pybertini.System()

	for f in cyclic(x):
	    sys.add_function(f)
	    
	print(sys) # long screen output, i know


.. [1] This is one of the reasons we wrote Bertini2's symbolic C++ core and exposed it to Python.
.. [2] https://www.sciencedirect.com/science/article/pii/S0377042702006982
