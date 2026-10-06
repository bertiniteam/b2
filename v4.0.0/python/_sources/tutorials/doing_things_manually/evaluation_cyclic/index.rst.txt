♻️ Evaluation of cyclic-:math:`n` polynomials
*******************************************************

Bertini is software for algebraic geometry.  This means we work with systems of polynomials, a critical component of which is system and function evaluation.

Bertini2 allows us to set up many kinds of functions, and thus systems, by exploting operator overloading.

Make some symbols
==================

Let's start by making some variables, programmatically [1]_.

.. testcode::

    import bertini
    import numpy

    num_vars = 10
    x = [None] * num_vars # preallocate the list
    for ii in range(num_vars):
        x[ii] = bertini.Variable('x' + str(ii))

Huzzah, we have ``num_vars`` variables!  This was hard to do in Bertini 1's classic style input files.  Now we can do it directly! 🎯

Write a function to produce the cyclic :math:`n` polynomials :cite:`cyclic_n`.

.. testcode::

    def cyclic(vars):
        n = len(vars)
        f = [None] * len(vars)
        y = []
        for ii in range(2):
            for x in vars:
                y.append(x)

        for ii in range(n-1):
            f[ii] = numpy.sum( [numpy.prod(y[jj:jj+ii+1]) for jj in range(n)] )

        # the last one is minus one
        f[-1] = numpy.prod(vars)-1
        return f

Now we will make a ``System``, and put the cyclic polynomials into it.

.. testcode::

    sys = bertini.System()

    for f in cyclic(x):
        sys.add_function(f)

    print(sys) # long screen output, i know

.. testoutput::
   :options: +ELLIPSIS

   ...
   10 functions:
     f_0 = x0+x1+x2+x3+x4+x5+x6+x7+x8+x9
   ...
     f_9 = x0*x1*x2*x3*x4*x5*x6*x7*x8*x9-1
   ...

We also need to associate the variables with the system.  Unassociated variables are left unknown, and retain their value until elsewhere set.

.. testcode::

    vg = bertini.VariableGroup()
    for var in x:
        vg.append(var)
    sys.add_variable_group(vg)



Now, let's evaluate it at the origin -- all zero's.  The returned value should be all zero's except the last entry, which should be -1.

.. testcode::

    s = numpy.zeros((10,), dtype=bertini.multiprec.complex_mp)
    result = sys.eval(s)
    assert complex(result[-1]) == -1                       # last cyclic function is (prod x) - 1
    assert all(complex(v) == 0 for v in result[:-1])       # the rest vanish at the origin

Yay, all zeros, except the last one is :math:`-1`.  Huzzah.

Let's change the values of our vector, and re-evaluate.

.. testcode::

    for ii in range(num_vars):
        s[ii] = bertini.multiprec.complex_mp(ii)
    result = sys.eval(s)


There is much more one can do, too!  Check out the other tutorials.  This is one of the older ones.

.. [1] This is one of the reasons we wrote Bertini2's symbolic C++ core and exposed it to Python.

Complete example
================

The whole tutorial as one runnable script:

.. literalinclude:: evaluation_cyclic.py
   :language: python
   :caption: evaluation_cyclic.py
