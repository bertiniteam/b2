"""Eigenvalues by numerical algebraic geometry, multihomogeneously.

The eigenproblem A x = lam x is the canonical example where multihomogeneous structure
beats total degree.  Writing it as a polynomial system

    (A - lam I) x = 0          (n equations, each bidegree (1,1) in {x} and {lam})
    c . x - 1     = 0          (a generic linear normalization fixing the eigenvector scale)

over two variable groups {x_0..x_{n-1}} and {lam}, the m-homogeneous Bezout number is n
(exactly the eigenvalue count) -- versus the total-degree bound 2^n.

This exercises the full block-composed path under the adaptive-precision tracker: the
products-of-linears MHom start, the blend-block homotopy, and (crucially) precision
staying in lockstep through the Cauchy endgame on the harder MHom paths.  We solve and
cross-check against numpy.linalg.eigvals -- the test is self-validating.
"""

import numpy as np
import pytest

import bertini as pb
from bertini import linalg as bla
from bertini.nag_algorithm import ZeroDimCauchyAdaptivePrecisionMHomogeneous

OK = int(pb.tracking.SuccessCode.Success)


def _eigen_system(A, c):
    """Build (A - lam I)x = 0, c.x - 1 = 0 with {x} and {lam} as variable groups.

    Written with the linear-algebra layer: x is a vector of variables, so the n eigen-
    equations are the single vector expression ``A @ x - lam*x``.
    """
    n = A.shape[0]
    x = bla.variable_vector('x', n)
    lam = pb.Variable('lam')
    sys = pb.System()
    bla.add_functions(sys, A @ x - lam * x)          # the rows of (A - lam I) x
    sys.add_function(c @ x - 1)                       # c . x - 1  (fixes eigenvector scale)
    sys.add_variable_group(pb.VariableGroup(list(x)))    # eigenvector group
    sys.add_variable_group(pb.VariableGroup([lam]))      # eigenvalue group
    return sys, lam


def _recovered_eigenvalues(solver, lam_index):
    sols = solver.solutions()
    md = solver.solution_metadata()
    return sorted(
        complex(sols[i][lam_index]).real
        for i in range(len(sols))
        if int(md[i].endgame_success) == OK and len(sols[i]) > 0
    )


def test_symmetric_3x3_eigenvalues_match_numpy():
    A = np.array([[2, 1, 0], [1, 3, 1], [0, 1, 4]])  # distinct real eigenvalues
    c = np.array([5, 8, 3])                            # a generic normalization
    sys, _ = _eigen_system(A, c)

    solver = ZeroDimCauchyAdaptivePrecisionMHomogeneous(sys)
    solver.solve()

    # user (dehomogenized) coordinates are [x0, x1, x2, lam]; lam is last.
    got = _recovered_eigenvalues(solver, lam_index=3)
    expected = sorted(np.linalg.eigvals(A).real)

    assert len(got) == A.shape[0]                      # MHom Bezout n == #eigenvalues
    for g, e in zip(got, expected):
        assert abs(g - e) < 1e-6
