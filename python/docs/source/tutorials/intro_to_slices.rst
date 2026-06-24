✂️ Intro to slices
*********************************************

.. testsetup:: *

   import numpy as np
   import bertini
   from bertini import linalg
   from bertini import multiprec as mp
   from bertini.nag_algorithm import Slice, WitnessSetMultiplePrecision

A **witness set** is how numerical algebraic geometry represents a positive-dimensional component
of a variety -- a curve, a surface, and so on. It is a triple:

* a **system** (the equations whose solution set contains the component),
* a generic linear **slice** (enough hyperplanes to cut the component down to points), and
* the **points** where the component meets that slice.

The two headline invariants fall straight out: the **dimension** of the component is the number of
hyperplanes in the slice, and the **degree** is the number of witness points.

The system and the points are ordinary objects you already know. The **slice** is the new one, and
it carries a fair amount of semantics -- so most of this tutorial is about slices.

The slice: a stack of linear forms
==================================

A :class:`~bertini.nag_algorithm.Slice` is a stack of linear forms :math:`M\,[x ; 1]` -- one row per
form, the trailing column of each row holding that form's constant term. You build one from exact
coefficients with :func:`~bertini.linalg.slice_from_coefficients`:

.. testcode::

   x, y = bertini.Variable('x'), bertini.Variable('y')

   # two linear forms on (x, y):  2x + 3y + 1   and   x - y + 4
   s = linalg.slice_from_coefficients([[2, 3, 1], [1, -1, 4]], [x, y])

   assert s.dimension() == 2        # two forms -> cuts a 2-dimensional component
   assert s.num_variables() == 2

or generate a generic one (the usual case in practice -- a witness set's slice is random):

.. testcode::

   x0, x1, x2 = bertini.Variable('x0'), bertini.Variable('x1'), bertini.Variable('x2')
   vg = bertini.VariableGroup([x0, x1, x2])

   generic = Slice.random_complex(vg, 2)    # two random hyperplanes in 3-space
   assert generic.dimension() == 2
   assert generic.num_variables() == 3

.. note::

   Coefficients must be **exact** (ints, ``fractions.Fraction``, exact strings like ``'3/4'``, or
   ``bertini.multiprec`` values). A Python ``float`` is refused -- its ~16 digits would cap the
   precision of the arbitrary-precision tree. See :func:`~bertini.linalg.coefficient`.

A slice is a sequence of forms -- and the shape rules matter
============================================================

This is the part with "lots of semantics running around," so it is worth being precise. A slice
behaves like a **Python sequence of its linear forms**, with the usual list semantics:

* an **integer index** selects an *element* -- the i-th form's coefficient **vector**;
* a **slice** ``[i:j]`` (or a list of indices) selects a *sub-collection* -- a new (sub-)``Slice``;
* the whole coefficient matrix is :meth:`~bertini.nag_algorithm.Slice.coefficients`, **always 2-D**.

.. testcode::

   # an ELEMENT: slice[i] is the i-th form's coefficient VECTOR (1-D), "a line is a vector"
   v = np.asarray(s[0])
   assert v.shape == (3,)
   assert [complex(c).real for c in v] == [2.0, 3.0, 1.0]      # 2x + 3y + 1

   # iterating yields the form vectors
   assert len(list(s)) == 2

   # a SUB-COLLECTION: slice[i:j] is a new Slice
   first = s[:1]
   assert first.dimension() == 1

   # the whole matrix: always 2-D, (num_forms, num_variables + 1)
   assert s.coefficients().shape == (2, 3)

The "always 2-D" is the important guarantee. The underlying binding library (eigenpy today) returns
a *one-row* matrix to Python as a **1-D** array -- a data-dependent shape that silently breaks code
written for the many-row case. This bites constantly in practice: a **curve** (dimension 1) has a
**single-form** slice, so its coefficient matrix is exactly the one that would collapse.

bertini owns the shape so you do not have to think about it: ``coefficients()`` is **always** 2-D,
even for one form, and the 1-D "vector" view is something you *ask for* by indexing an element
(``s[i]``) -- never something the form count hands you by surprise.

.. testcode::

   one_form = linalg.slice_from_coefficients([[2, 3, 1]], [x, y])
   assert one_form.coefficients().shape == (1, 3)     # NOT (3,) -- it never collapses
   assert np.asarray(one_form[0]).shape == (3,)       # the vector view is explicit

(The "why" -- and the list of accessors that still carry the raw hazard -- is recorded in
``docs/adr/0033``.)

A slice rides along with a system
=================================

A slice does **not** know about homogenization or patches; the *system* owns those. A slice is just
linear forms, meant to be carried alongside a system. Two ways to put them together:

.. testcode::

   # as_system(): a standalone System of just the slice's forms
   only_slice = s.as_system()
   assert only_slice.num_functions() == 2

   # add_to(): append the slice's forms to an existing system (over the same variables)
   sphere = bertini.System()
   sphere.add_variable_group(vg)
   sphere.add_function(x0*x0 + x1*x1 + x2*x2 - 1)      # the unit sphere: a surface in 3-space
   generic.add_to(sphere)
   assert sphere.num_functions() == 3                  # 1 sphere equation + 2 slice forms

.. note::

   ``add_to`` is homogenization-aware. If the system has already been homogenized, the slice's
   constant term is folded onto the homogenizing variable (exactly as the system did to its own
   functions), so the appended forms stay consistent. You never homogenize the slice yourself.

The witness set
===============

Now assemble the triple. You can build a witness set all at once, or incrementally -- add points as
you find them:

.. testcode::

   def pt(*entries):
       return np.array([mp.Complex(str(e)) for e in entries], dtype=mp.Complex)

   sys = bertini.System()
   sys.add_variable_group(vg)
   sys.add_function(x0*x0 + x1*x1 + x2*x2 - 1)         # the sphere again (a 2-dim component)
   slice2 = Slice.random_complex(vg, 2)

   # all at once: points + slice + system
   w = WitnessSetMultiplePrecision([pt(1, 0, 0), pt(0, 1, 0)], slice2, sys)
   assert w.degree() == 2          # two witness points
   assert w.dimension() == 2       # a surface
   assert w.is_consistent()        # 3 variables - 1 equation == slice dimension 2

   # or incrementally
   w2 = WitnessSetMultiplePrecision()
   w2.set_system(sys)
   w2.set_slice(slice2)
   w2.add_point(pt(1, 0, 0))
   w2.add_point(pt(0, 1, 0))
   assert w2.degree() == 2

A witness set prints a readable summary:

.. testcode::

   print(repr(w))

.. testoutput::

   WitnessSet -- dimension 2, degree 2 (consistent)
     slice:  2 linear forms on 3 variables
     system: 1 function on 3 variables

Carrying it around, and emitting to Bertini 1
=============================================

A witness set (and a slice) serializes, so you can pickle one, send it to another process, or
checkpoint it -- the deserialized system is restored ready to evaluate:

.. testcode::

   import pickle
   w_again = pickle.loads(pickle.dumps(w))
   assert w_again.degree() == 2
   assert w_again.is_consistent()

You can also emit the **witness (square) system** -- the system with the slice's forms appended --
as a Bertini 1 classic input file, for cross-validation against the original solver:

.. testcode::

   square = w.witness_system()          # system + slice, the square system to track
   assert square.num_functions() == 3
   text = w.to_classic_input()
   assert 'CONFIG' in text and 'INPUT' in text

Building start systems from slices (regeneration)
=================================================

Because a slice's rows share their layout with a products-of-linears block, a list of slices drops
straight into a product-of-linears start system -- the bridge regeneration is built on. Each slice
becomes one function that is the product of its hyperplanes:

.. testcode::

   start = bertini.System()
   start.add_variable_group(bertini.VariableGroup([x, y]))
   sa = linalg.slice_from_coefficients([[1, 0, -1], [1, 0, 1]], [x, y])    # (x - 1)(x + 1)
   sb = linalg.slice_from_coefficients([[0, 1, -1], [0, 1, -2]], [x, y])   # (y - 1)(y - 2)
   linalg.add_slices_as_products(start, [sa, sb])
   assert list(start.degrees()) == [2, 2]
