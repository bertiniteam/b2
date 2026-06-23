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
:mod:`bertini.linalg` layer so we can write the equations as actual linear algebra:

.. testcode::

    import numpy as np
    import bertini as bertini
    from bertini import linalg
    from bertini.nag_algorithm import ZeroDim

Pick a small symmetric matrix (real, distinct eigenvalues make the check easy to read):

.. testcode::

    A = np.array([[2, 1, 0],
                  [1, 3, 1],
                  [0, 1, 4]])
    n = A.shape[0]

With the linear-algebra layer the equations *are* linear algebra: make a **vector of
variables** for the eigenvector and a scalar for the eigenvalue, and :math:`(A-\lambda I)x`
is just ``A @ x - lam*x``:

.. testcode::

    x = linalg.variable_vector('x', n)        # array([x0, x1, x2], dtype=object)
    lam = bertini.Variable('lam')

.. note::

   The matrix ``A`` here is **integer**, so it enters the function tree exactly.
   ``bertini.linalg`` deliberately refuses plain Python ``float`` entries -- a 64-bit float
   would cap the precision of the whole arbitrary-precision computation at ~16 digits.  For
   non-integer coefficients pass exact values
   (``linalg.as_coefficients([['5/2', '1'], ...])`` accepts fractions, exact decimal
   strings, and bertini multiprecision values).

Two ways to make it zero-dimensional
====================================

An eigenvector is only defined up to scale, so :math:`(A-\lambda I)x = 0` alone has a whole
*line* of solutions for each eigenvalue (and the trivial :math:`x=0`).  There are two clean
ways to cut that down to isolated points, and Bertini supports both.

**Affine, with a normalization.**  Keep :math:`x` an ordinary (affine) variable group and
add one generic linear equation :math:`c\cdot x = 1` to pin the scale:

.. testcode::

    def build_affine():
        x = linalg.variable_vector('x', n)
        lam = bertini.Variable('lam')
        c = np.array([5, 8, 3])                          # any generic integer vector
        sys = bertini.System()
        linalg.add_functions(sys, A @ x - lam * x)       # the rows of (A - lam I) x
        sys.add_function(c @ x - 1)                      # fix the eigenvector scale
        sys.add_variable_group(bertini.VariableGroup(list(x)))
        sys.add_variable_group(bertini.VariableGroup([lam]))
        return sys

**Projective.**  Declare :math:`x` a *projective* (homogeneous) variable group: the
eigenvector then lives in :math:`\mathbb{P}^{n-1}` natively, where scale is already
quotiented out, so **no normalization equation is needed**:

.. testcode::

    def build_projective():
        x = linalg.variable_vector('x', n)
        lam = bertini.Variable('lam')
        sys = bertini.System()
        linalg.add_functions(sys, A @ x - lam * x)       # nothing else!
        sys.add_hom_variable_group(bertini.VariableGroup(list(x)))   # x in P^{n-1}
        sys.add_variable_group(bertini.VariableGroup([lam]))
        return sys

Both are multihomogeneous with the same Bézout number :math:`n`, so both track exactly
:math:`n` paths -- the projective version just expresses the geometry directly instead of
bolting on a normalization.

Solving, and reading off the eigenvalues
========================================

The solve is identical for either system.  A small helper runs it and pulls out the
eigenvalues (``lam`` is the last user coordinate of each solution; the rest are the
eigenvector):

.. testcode::

    def eigenvalues_of(system):
        solver = ZeroDim(system, mptype='adaptive', startsystem='mhom')
        solver.solve()
        good = solver.finite_solutions()                 # the n eigenpairs (lam is the last coord)
        assert len(good) == n                            # one path per eigenvalue
        return sorted(complex(s[len(s) - 1]).real for s in good)

Both formulations recover ``numpy``'s eigenvalues:

.. testcode::

    expected = sorted(np.linalg.eigvals(A).real)
    for build in (build_affine, build_projective):
        got = eigenvalues_of(build())
        for g, e in zip(got, expected):
            assert abs(g - e) < 1e-6

Comparing the two
=================

Both track the same number of paths, but the projective system is **leaner** -- one fewer
equation, and the eigenvector group needs no homogenizing coordinate -- so it is usually
the faster of the two.  Time them and see:

.. testcode::

    import time

    for name, build in (("affine+normalization", build_affine),
                        ("projective          ", build_projective)):
        t0 = time.perf_counter()
        got = eigenvalues_of(build())
        dt = time.perf_counter() - t0
        print(f"{name}: {got}  in {dt:.3f}s")

.. testoutput::
   :options: +ELLIPSIS

   affine+normalization: [...]  in ...s
   projective          : [...]  in ...s

The projective system is the more direct statement of the problem and typically solves
faster; the affine one is handy when you would rather not think projectively.  Either way
the solver returns the same eigenvalues -- pick whichever reads better for your problem.

Why this matters, and where it is going
=======================================

Two things are worth noticing.  First, the eigenvalues came out of *path tracking*, not a
dedicated eigensolver -- the same machinery solves any polynomial system, and here the
multihomogeneous structure made it efficient (``n`` paths, not :math:`2^n`).  Second, we
never built :math:`(A-\lambda I)x` entry by entry: with :math:`x` a vector of variables, the
condition reads ``A @ x - lam*x`` -- the linear algebra you would write on paper.  That is
:mod:`bertini.linalg`, whose motivating example is exactly this problem.  Coefficients on
variables stay exact, so the precision of the solve is never silently capped by a stray
floating-point literal.

The same idea extends to **matrices of variables** (``linalg.variable_matrix``), which is
where problems like rank conditions and regeneration are headed.  Today the vector
equations are expanded into scalar polynomials; a future evaluation block will carry the
matrix structure all the way into the numerics for larger problems.
