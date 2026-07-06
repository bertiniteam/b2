"""Eigenvalues by numerical algebraic geometry, multihomogeneously.

The eigenproblem A x = lam x is the canonical example where multihomogeneous structure
beats total degree.  There are two equally valid ways to make it a zero-dimensional system,
and both give the m-homogeneous Bezout number n (exactly the eigenvalue count, versus the
total-degree bound 2^n):

  * affine + normalization: x an affine group, plus a generic linear form c.x - 1 = 0 to
    fix the eigenvector scale (otherwise each eigenline is a curve of solutions).
  * projective: x a projective (homogeneous) variable group -- the eigenvector lives in
    P^{n-1} natively, so no normalization equation is needed.

Both exercise the full block-composed path under the adaptive-precision tracker (the
products-of-linears MHom start, the blend-block homotopy, precision lockstep through the
Cauchy endgame), and the projective case additionally exercises projective-group MHom
start-point generation.  We cross-check against numpy.linalg.eigvals -- self-validating.
"""

import numpy as np
import pytest

import bertini as pb
from bertini import ZeroDimSolver

OK = int(pb.SuccessCode.Success)


def _eigen_system_affine(A, c):
    """(A - lam I)x = 0, c.x - 1 = 0 with x and lam as AFFINE variable groups."""
    n = A.shape[0]
    x = np.array(pb.variables('x', n), dtype=object)
    lam = pb.Variable('lam')
    sys = pb.System()
    sys.add_functions(A @ x - lam * x)          # the rows of (A - lam I) x
    sys.add_function(c @ x - 1)                       # fix the eigenvector scale
    sys.add_variable_group(pb.VariableGroup(list(x)))
    sys.add_variable_group(pb.VariableGroup([lam]))
    return sys


def _eigen_system_projective(A):
    """(A - lam I)x = 0 with x a PROJECTIVE group and lam affine -- no normalization."""
    n = A.shape[0]
    x = np.array(pb.variables('x', n), dtype=object)
    lam = pb.Variable('lam')
    sys = pb.System()
    sys.add_functions(A @ x - lam * x)
    sys.add_hom_variable_group(pb.VariableGroup(list(x)))   # eigenvector in P^{n-1}
    sys.add_variable_group(pb.VariableGroup([lam]))
    return sys


def _recovered_eigenvalues(solver):
    # lam is the last user (dehomogenized) coordinate in both formulations.
    sols = solver.all_solutions()
    md = solver.solution_metadata()
    return sorted(
        complex(sols[i][len(sols[i]) - 1]).real
        for i in range(len(sols))
        if int(md[i].endgame_success_code) == OK and len(sols[i]) > 0
    )


@pytest.mark.parametrize("make_system", [
    lambda A: _eigen_system_affine(A, np.array([5, 8, 3])),
    _eigen_system_projective,
])
def test_symmetric_3x3_eigenvalues_match_numpy(make_system):
    A = np.array([[2, 1, 0], [1, 3, 1], [0, 1, 4]])  # distinct real eigenvalues

    solver = ZeroDimSolver(make_system(A), endgame='cauchy', mptype='adaptive', startsystem='mhom')
    solver.solve()

    got = _recovered_eigenvalues(solver)
    expected = sorted(np.linalg.eigvals(A).real)

    assert len(got) == A.shape[0]                      # MHom Bezout n == #eigenvalues
    for g, e in zip(got, expected):
        assert abs(g - e) < 1e-6
