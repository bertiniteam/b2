#!/usr/bin/env python
"""Find ALL eigenvalues of a matrix as a distributed multihomogeneous homotopy solve.

The eigenproblem ``(A - lam I) x = 0`` is the textbook case where multihomogeneous structure
crushes total degree: write the eigenvector ``x`` in projective space ``P^{n-1}`` and ``lam`` as
an affine unknown, and the multihomogeneous Bezout number is exactly ``n`` -- one path per
eigenvalue -- versus the total-degree bound ``2^n``.  So an ``n x n`` matrix gives ``n`` paths,
and those ``n`` independent paths spread across MPI workers (make ``n`` large enough that there
is real work to hand out).

Run it::

    python solve_eigenvalues.py --size 24                            # serial baseline
    mpirun -n 5 python solve_eigenvalues.py --size 24                # 1 manager + 4 workers
    OMP_NUM_THREADS=4 mpirun -n 3 --bind-to none python solve_eigenvalues.py --size 24

Only rank 0 prints; it checks the recovered eigenvalues against ``numpy.linalg.eigvals`` --
correctness, not just a path count.
"""

import argparse
import os
import time

import numpy as np
from mpi4py import MPI

import bertini as pb
from bertini import linalg
from bertini.nag_algorithm import ZeroDimSolver

OK = int(pb.tracking.SuccessCode.Success)


def random_symmetric_integer_matrix(n, seed, spread=9):
    """A reproducible n x n symmetric integer matrix (real, generically distinct spectrum).

    The same ``seed`` builds the same matrix on every rank -- essential, since each rank
    constructs its own copy of the system to solve.
    """
    rng = np.random.default_rng(seed)
    M = rng.integers(-spread, spread + 1, size=(n, n))
    return M + M.T                       # symmetric => real eigenvalues


def eigen_system(A):
    """(A - lam I) x = 0 with x a PROJECTIVE group and lam affine.

    Projective ``x`` lives in ``P^{n-1}`` natively, so the eigenvector scale is already fixed --
    no normalization equation, and the multihomogeneous Bezout number is exactly ``n``.
    """
    n = A.shape[0]
    x = linalg.variable_vector('x', n)
    lam = pb.Variable('lam')
    sys = pb.System()
    linalg.add_functions(sys, A @ x - lam * x)                  # the rows of (A - lam I) x
    sys.add_hom_variable_group(pb.VariableGroup(list(x)))    # eigenvector in P^{n-1}
    sys.add_variable_group(pb.VariableGroup([lam]))
    return sys


def recovered_eigenvalues(solver):
    """lam is the last dehomogenized coordinate; collect it from each successful, finite path.

    We trust the solver's own classification: a real eigenvalue is the lam coordinate of a path
    whose endgame succeeded and that the library calls finite (is_finite) -- no hand-rolled cutoff.
    """
    sols = solver.all_solutions()
    md = solver.solution_metadata()
    return [
        complex(sols[i][len(sols[i]) - 1]).real
        for i in range(len(sols))
        if int(md[i].endgame_success_code) == OK and md[i].is_finite and len(sols[i]) > 0
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', type=int, default=24, help='matrix size n (default 24)')
    parser.add_argument('--seed', type=int, default=20260614,
                        help='RNG seed for the matrix AND the homotopy (default fixed, so runs are '
                             'reproducible and serial/distributed are directly comparable)')
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    pb.random.set_random_seed(args.seed)        # reproducible homotopy (same on every rank)
    A = random_symmetric_integer_matrix(args.size, args.seed)
    solver = ZeroDimSolver(eigen_system(A), endgame='cauchy', mptype='adaptive', startsystem='mhom')

    start = time.time()
    if comm.Get_size() > 1:
        solver.solve(communicator=comm)         # manager (rank 0) hands paths to the workers
    else:
        solver.solve()                          # serial: a lone manager would have no workers
    elapsed = time.time() - start

    if not pb.parallel.is_manager():
        return

    got = recovered_eigenvalues(solver)
    expected = sorted(np.linalg.eigvals(A).real)

    # Distinct recovered eigenvalues: a path can occasionally duplicate one, so dedupe (the
    # eigenvalues are real and well separated) before counting -- the meaningful quantity is the
    # number of distinct eigenvalues, which is deterministic.
    distinct = []
    for g in sorted(got):
        if not distinct or abs(g - distinct[-1]) >= 1e-6:
            distinct.append(g)

    print('eigen-{}:  ranks={}  threads/rank={}  paths={}  recovered={}  wall={:.1f}s'.format(
        args.size, comm.Get_size(), os.environ.get('OMP_NUM_THREADS', '1'),
        args.size, len(distinct), elapsed))

    # Correctness: every numpy eigenvalue is matched by a homotopy-recovered one (closest-match,
    # robust to ordering and tiny imaginary parts), and we recover all n distinct eigenvalues.
    assert len(distinct) == args.size, \
        'expected {} distinct eigenvalues, recovered {}'.format(args.size, len(distinct))
    max_err = max(min(abs(g - e) for g in distinct) for e in expected)
    assert max_err < 1e-6, \
        'recovered eigenvalues disagree with numpy by {:.2e}'.format(max_err)
    print('  OK: all {} eigenvalues match numpy.linalg.eigvals (max error {:.2e})'.format(
        args.size, max_err))


if __name__ == '__main__':
    main()
