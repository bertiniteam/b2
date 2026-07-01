"""Eigenvalues by homotopy continuation -- Bertini 2 tutorial.

Solve the eigenvalue problem (A - lam I) x = 0 as a multihomogeneous polynomial
system, recovering the eigenvalues by path tracking.

Run:  python eigenvalues_by_homotopy.py
"""

import time

import numpy as np
import bertini as bertini
from bertini import linalg
from bertini.nag_algorithm import ZeroDimSolver


def setup_matrix():
    # small symmetric matrix: real, distinct eigenvalues make the check easy to read
    A = np.array([[2, 1, 0],
                  [1, 3, 1],
                  [0, 1, 4]])
    n = A.shape[0]
    return A, n


def build_affine(A, n):
    x = linalg.variable_vector('x', n)
    lam = bertini.Variable('lam')
    c = np.array([5, 8, 3])                          # any generic integer vector
    sys = bertini.System()
    linalg.add_functions(sys, A @ x - lam * x)       # the rows of (A - lam I) x
    sys.add_function(c @ x - 1)                      # fix the eigenvector scale
    sys.add_variable_group(bertini.VariableGroup(list(x)))
    sys.add_variable_group(bertini.VariableGroup([lam]))
    return sys


def build_projective(A, n):
    x = linalg.variable_vector('x', n)
    lam = bertini.Variable('lam')
    sys = bertini.System()
    linalg.add_functions(sys, A @ x - lam * x)       # nothing else!
    sys.add_hom_variable_group(bertini.VariableGroup(list(x)))   # x in P^{n-1}
    sys.add_variable_group(bertini.VariableGroup([lam]))
    return sys


def eigenvalues_of(system, n):
    solver = ZeroDimSolver(system, mptype='adaptive', startsystem='mhom')
    solver.solve()
    good = solver.finite_solutions()                 # the n eigenpairs (lam is the last coord)
    assert len(good) == n                            # one path per eigenvalue
    return sorted(complex(s[len(s) - 1]).real for s in good)


def check_both(A, n):
    # both formulations recover numpy's eigenvalues
    expected = sorted(np.linalg.eigvals(A).real)
    for build in (build_affine, build_projective):
        got = eigenvalues_of(build(A, n), n)
        for g, e in zip(got, expected):
            assert abs(g - e) < 1e-6


def compare_timing(A, n):
    for name, build in (("affine+normalization", build_affine),
                        ("projective          ", build_projective)):
        t0 = time.perf_counter()
        got = eigenvalues_of(build(A, n), n)
        dt = time.perf_counter() - t0
        print(f"{name}: {got}  in {dt:.3f}s")


def main():
    A, n = setup_matrix()
    check_both(A, n)
    compare_timing(A, n)


if __name__ == '__main__':
    main()
