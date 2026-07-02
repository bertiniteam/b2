🎮 Using an endgame to compute singular endpoints 
*********************************************************


.. note::

    This tutorial is useful to those who want very detailed control of running pieces of various algorithms in Bertini 2.  If you just want to be able to solve a system, this tutorial is probably not for you.

Background
==============

Polynomial systems often have singular solutions.  In numerical algebraic geometry, we want to compute all solutions, even the challenging singular ones.  The normal method of homotopy continuation with straight-line tracking fails to compute such roots, because tracking to a place where the Jacobian is non-invertible using methods that require inverting the Jacobian is doomed to fail [#]_.  

So, if we can't track to a singular solution, but we still want to track to compute them, what are we to do?  We track around them, or near them, but not actually to them.  These methods are collectively called *endgames*, a term coined to evoke a sense of chess :cite:`morgan1990computing` :cite:`morgan1992computing` :cite:`morgan1992power`.  Thanks, Andrew Sommese, Charles Wampler, and Alexander Morgan, for everything you have given our community.

Endgames represent a way to finish a tracking of a path, when the endpoint is possibly singular.  Rather than track all the way to the endtime, you instead run an endgame that uses mathematical theory to compute the root.

Endgames in Bertini 2
==========================

An endgame is a computational tool that one does in the final stage of a path track to a possibly singular root.  There are two implemented endgames in Bertini:

#. Power series -- uses `Hermite interpolation <https://en.wikipedia.org/wiki/Hermite_interpolation>`_ across a sequence of geometrically-spaced points (in time) to extrapolate to a target time :cite:`morgan1992power`.
#. Cauchy -- uses `Cauchy's integral formula <https://en.wikipedia.org/wiki/Cauchy's_integral_formula>`_ in a sequence of circles about the root you are computing.

Both try to compute the cycle number :math:`c` for the root.  In the power series endgame, :math:`c` is used as the degree of a Hermite interpolant used to extrapolate to 0.  In the Cauchy endgame,  it is used for the number of cycles to walk before returning to the same point, computing a trapezoid-rule integral along the way.

Each is provided in the three precision modes, double, fixed multiple, and adaptive.  Since we are using the :class:`~bertini.AMPTracker` in this tutorial, we will of course use the adaptive endgame.  I really like the Cauchy endgame, so we're in the land of the :class:`~bertini.endgame.AMPCauchyEndgame`.


Example
----------


Form a system
~~~~~~~~~~~~~~~~

The Griewank-Osborne system has one multiplicity-three singular solution at the origin :cite:`griewank1983analysis`.  Let's build it from scratch, for the practice.

.. testcode::

    import bertini

    gw = bertini.System()

    x = bertini.Variable("x")
    y = bertini.Variable("y")

    vg = bertini.VariableGroup()
    vg.append(x)
    vg.append(y)
    gw.add_variable_group(vg)

    gw.add_function(bertini.multiprec.rational_mp(29,16)*x**3 - 2*x*y)
    gw.add_function(y - x**2)


Form a start system and homotopy 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next, we make the total degree start system for `gw`, and couple it using the gamma trick :cite:`morgan1987homotopy` and a path variable.

.. testcode::

    t = bertini.Variable('t')
    td = bertini.system.start_system.TotalDegreeLinearProduct(gw)
    gamma = bertini.symbolics.Rational.rand()
    hom = (1-t)*gw + t*gamma*td
    hom.add_path_variable(t)



🛤 Track to the endgame boundary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make a tracker.  I use adaptive precision a lot, so we'll roll with that.  There are also double and fixed-multiple versions.  See the other tutorials or the detailed documentation.

.. testcode::

    tr = bertini.AMPTracker(hom)

    start_time = bertini.multiprec.complex_mp("1")
    eg_boundary = bertini.multiprec.complex_mp("0.1")

    midpath_points = [None]*td.num_start_points()
    for ii in range(td.num_start_points()):
        midpath_points[ii] = bertini.multiprec.Vector(gw.num_variables())   # result must be pre-sized
        code = tr.track_path(result=midpath_points[ii], start_time=start_time, end_time=eg_boundary, start_point=td.start_point_mp(ii))
        assert code == bertini.SuccessCode.Success                 # every path reaches the boundary




🎮 Use the endgame
~~~~~~~~~~~~~~~~~~~~


To make an endgame, we feed it the tracker it should drive and the **endgame-boundary time**
(where the endgame takes over from straight-line tracking -- the :math:`t=0.1` we tracked to
above).  There are also config structs to play with, that control the way things are computed.

.. testcode::

    eg = bertini.endgame.AMPCauchyEndgame(tr, eg_boundary)

    # make an observer to be able to see what's going on inside
    ob = bertini.endgame.observers.amp_cauchy.GoryDetailLogger()

    eg.add_observer(ob)

The ``GoryDetailLogger`` writes a running, step-by-step account of the endgame to Bertini's log as
it runs -- handy when you want to watch the loops being walked.  Since the endgame hasn't been run
yet, things are empty and default:

.. testcode::

    assert eg.cycle_number() == 0
    assert len(eg.final_approximation()) == 0    # nothing computed yet

The endgame is used by invoking ``run``, feeding it just the boundary point to refine: the
endgame-boundary time and the target time (:math:`t=0`) were both fixed when we constructed the
endgame, so all ``run`` needs is where to start.

.. testcode::

    final_points = []
    codes = []
    for ii in range(td.num_start_points()):
        codes.append(eg.run(midpath_points[ii]))   # refine from the boundary down to t = 0
        final_points.append(eg.final_approximation())

The Griewank-Osborne system's only finite solution is the triple point at the origin, so exactly
three of the six homotopy paths converge there (the other three run off to infinity):

.. testcode::

    origin_hits = sum(1 for fa in final_points
                      if len(fa) and max(abs(complex(v)) for v in fa) < 1e-6)
    print('paths landing on the triple point:', origin_hits)

.. testoutput::

    paths landing on the triple point: 3


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




Complete example
================

The whole tutorial as one runnable script -- assemble nothing, just run it:

.. literalinclude:: manual_endgame_usage.py
   :language: python
   :caption: manual_endgame_usage.py
