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

We need ``numpy`` for the matrix and the cross-check, and ``bertini`` for the solve::

    import numpy as np
    import bertini as bertini
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

    c = np.array([5, 8, 3])   # any generic vector works

Now build the variables and the equations.  The eigenvector coordinates and the eigenvalue
go into **separate variable groups** -- that is what makes this multihomogeneous::

    xs = [bertini.Variable(f'x{i}') for i in range(n)]
    lam = bertini.Variable('lam')

    sys = bertini.System()

    # one equation per row of (A - lam*I) x
    for i in range(n):
        row = -lam * xs[i]
        for j in range(n):
            if A[i, j]:
                row = row + int(A[i, j]) * xs[j]
        sys.add_function(row)

    # the scale-fixing normalization c . x - 1
    norm = -1
    for j in range(n):
        norm = norm + int(c[j]) * xs[j]
    sys.add_function(norm)

    sys.add_variable_group(bertini.VariableGroup(xs))      # the eigenvector group
    sys.add_variable_group(bertini.VariableGroup([lam]))   # the eigenvalue group

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
multihomogeneous structure made it efficient.  Second, building :math:`(A-\lambda I)x` by
hand, entry by entry, is verbose.  The natural way to write this is

.. math::

   A x - \lambda x = 0,

with :math:`x` a *vector of variables* and :math:`A` a matrix.  A succinct
linear-algebra layer over systems -- expressing conditions on vectors and matrices of
variables directly -- is in progress, and the eigenvalue problem is its motivating example.
