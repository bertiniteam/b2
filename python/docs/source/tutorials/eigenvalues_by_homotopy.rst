🧮 Eigenvalues by homotopy continuation
*********************************************

.. testsetup:: *

   import bertini

The eigenvalue problem :math:`A x = \lambda x` is the cleanest place to see *why*
numerical algebraic geometry cares about **variable-group structure**.  Written as a
polynomial system it is

.. math::

   (A - \lambda I)\, x = 0,

a system that is degree 1 in the eigenvector coordinates :math:`x` and degree 1 in the
eigenvalue :math:`\lambda`.  Treating :math:`x` and :math:`\lambda` as *separate* groups
of variables, every equation has bidegree :math:`(1,1)`, and the **multihomogeneous
Bézout number** is :math:`n` -- exactly the number of eigenvalues.  The naive total-degree
count would be :math:`2^n`.  For an :math:`8\times 8` matrix that is 8 paths to track
instead of 256.

This is the payoff of Bertini's multihomogeneous start system: it tracks one path per
*actual* solution, not one per total-degree phantom.

Setting up the system
=====================

We need ``numpy`` for the matrix and the cross-check, ``bertini`` for the solve, and the
:mod:`bertini.linalg` layer so we can write the equations as actual linear algebra::

    import numpy as np
    import bertini as bertini
    from bertini import linalg
    from bertini.nag_algorithm import ZeroDimCauchyAdaptivePrecisionMHomogeneous

Pick a small symmetric matrix (real, distinct eigenvalues make the check easy to read)::

    A = np.array([[2, 1, 0],
                  [1, 3, 1],
                  [0, 1, 4]])
    n = A.shape[0]

An eigenvector is only defined up to scale, so :math:`(A-\lambda I)x = 0` alone has a whole
line of solutions for each eigenvalue (and the trivial :math:`x=0`).  We pin the scale with
one generic linear normalization :math:`c\cdot x = 1`, which selects a single representative
on each eigenline and excludes :math:`x=0`::

    c = np.array([5, 8, 3])   # any generic integer vector works

Make a **vector of variables** for the eigenvector, and a scalar for the eigenvalue::

    x = linalg.variable_vector('x', n)        # array([x0, x1, x2], dtype=object)
    lam = bertini.Variable('lam')

Now the equations *are* linear algebra.  :math:`(A-\lambda I)x` is just ``A @ x - lam*x``,
a length-:math:`n` array of expressions, and the normalization is ``c @ x - 1``::

    sys = bertini.System()
    linalg.add_functions(sys, A @ x - lam * x)     # the rows of (A - lam I) x
    sys.add_function(c @ x - 1)                     # fixes the eigenvector scale

The eigenvector coordinates and the eigenvalue go into **separate variable groups** -- that
is what makes the problem multihomogeneous::

    sys.add_variable_group(bertini.VariableGroup(list(x)))   # the eigenvector group
    sys.add_variable_group(bertini.VariableGroup([lam]))     # the eigenvalue group

.. note::

   The matrix ``A`` and the vector ``c`` here are **integers**, which enter the function
   tree exactly.  ``bertini.linalg`` deliberately refuses plain Python ``float`` entries --
   a 64-bit float would cap the precision of the whole arbitrary-precision computation at
   ~16 digits.  For non-integer coefficients pass exact values
   (``linalg.as_coefficients([['5/2', '1'], ...])`` accepts fractions, exact decimal
   strings, and bertini multiprecision values).

Solving
=======

Use the adaptive-precision, Cauchy-endgame zero-dimensional solver with the
**multihomogeneous** start system::

    solver = ZeroDimCauchyAdaptivePrecisionMHomogeneous(sys)
    solver.solve()

The solver tracks :math:`n` paths -- the multihomogeneous Bézout number -- and each one
ends at an eigenpair.  Collect the successfully-tracked solutions::

    OK = int(bertini.tracking.SuccessCode.Success)
    sols = solver.solutions()                 # user coordinates: [x0, x1, x2, lam]
    meta = solver.solution_metadata()
    good = [sols[i] for i in range(len(sols))
            if int(meta[i].endgame_success) == OK]

    assert len(good) == n                     # one solution per eigenvalue

Reading off the eigenvalues
===========================

The eigenvalue is the last coordinate of each solution; the rest are the (normalized)
eigenvector.  Compare against ``numpy`` -- they agree::

    recovered = sorted(complex(s[n]).real for s in good)
    expected  = sorted(np.linalg.eigvals(A).real)

    for got, want in zip(recovered, expected):
        assert abs(got - want) < 1e-6

Why this matters, and where it is going
=======================================

Two things are worth noticing.  First, the eigenvalues came out of *path tracking*, not a
dedicated eigensolver -- the same machinery solves any polynomial system, and here the
multihomogeneous structure made it efficient.  Second, we never built :math:`(A-\lambda
I)x` entry by entry: with :math:`x` a vector of variables, the condition reads ``A @ x -
lam*x`` -- the linear algebra you would write on paper.  That is :mod:`bertini.linalg`,
whose motivating example is exactly this problem.  Coefficients on variables stay exact, so
the precision of the solve is never silently capped by a stray floating-point literal.

The same idea extends to **matrices of variables** (``linalg.variable_matrix``), which is
where problems like rank conditions and regeneration are headed.  Today the vector
equations are expanded into scalar polynomials; a future evaluation block will carry the
matrix structure all the way into the numerics for larger problems.
