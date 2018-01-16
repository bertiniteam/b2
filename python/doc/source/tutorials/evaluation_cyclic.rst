Evaluation of cyclic-:math:`n` polynomials
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

We also need to associate the variables with the system.  Unassociated variables are left unknown, and retain their value until elsewhere set.

::
	
	vg = pybertini.VariableGroup()
	for var in x:
		vg.append(var)
	sys.add_variable_group(vg)

Let's simplify this.  It will modify elements of the constructed function tree, even those held externally -- Bertini uses shared pointers under the hood, so pay attention to where you re-use parts of your functions, because later modification of them without deep cloning will cause ... modification elsewhere, too.  

::

	pybertini.system.simplify(sys)

Now, let's evaluate it at the origin -- all zero's (0 is the default value for multiprecision complex numbers in Bertini2).  The returned value should be all zero's except the last entry, which should be -1.

::

	s = pybertini.multiprec.Vector() # todo allow int in constructor
	s.resize(num_vars)
	sys.eval(s)

Yay, all zeros, except the last one is -1.  Huzzah.

Let's change the values of our vector, and re-evaluate.

::

	for ii in range(num_vars):
		s[ii] = pybertini.multiprec.complex(ii)
	sys.eval(s)


There is much more one can do, too!  Please write the authors, particularly Dani Brake, for more.

.. [1] This is one of the reasons we wrote Bertini2's symbolic C++ core and exposed it to Python.