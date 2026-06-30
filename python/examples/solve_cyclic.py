#!/usr/bin/env python
"""Solve the cyclic-n system, distributed across MPI ranks (with optional threads per rank).

The cyclic-n polynomials are a classic benchmark for polynomial system solving.  Total-degree
homotopy tracks one path per Bezout count (the product of the degrees), and the work is
pleasantly parallel: each path is independent.  Bertini hands the paths out from rank 0 (the
manager) to the worker ranks, and -- if OMP_NUM_THREADS > 1 -- each worker runs several tracking
threads.

Run it::

    python solve_cyclic.py --n 6                              # serial baseline
    mpirun -n 5 python solve_cyclic.py --n 6                  # 1 manager + 4 workers
    OMP_NUM_THREADS=4 mpirun -n 3 --bind-to none python solve_cyclic.py --n 6   # 2 workers x 4 threads

Only rank 0 prints; it checks the distinct finite-solution count -- taken from the solver's own
solution metadata -- against the known cyclic-n value.
"""

import argparse
import math
import os
import time

import numpy as np
from mpi4py import MPI

import bertini as pb

# distinct finite (nonsingular) solution counts of the cyclic-n system
KNOWN_FINITE = {4: 16, 5: 70, 6: 156, 7: 924, 8: 1152}


def cyclic_system(n):
    """The cyclic-n system: the n-1 cyclic sums of products, plus (prod - 1)."""
    x = [pb.Variable('x{}'.format(i)) for i in range(n)]
    w = x + x  # doubled, to take products that wrap around
    sys = pb.System()
    for length in range(1, n):                       # degree-`length` cyclic sum
        sys.add_function(np.sum([np.prod(w[start:start + length]) for start in range(n)]))
    sys.add_function(np.prod(x) - 1)                 # the product, normalized
    sys.add_variable_group(pb.VariableGroup(x))
    return sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=6, help='which cyclic-n (default 6)')
    parser.add_argument('--seed', type=int, default=None,
                        help='optional RNG seed for a reproducible run (default: random each run); '
                             'the same seed on every rank fixes the homotopy, so serial and '
                             'distributed runs are directly comparable')
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    if args.seed is not None:
        pb.random.set_random_seed(args.seed)
    system = cyclic_system(args.n)
    solver = pb.nag_algorithm.ZeroDimSolver(system, endgame='cauchy', mptype='adaptive', startsystem='binomial')

    # Adaptive precision is what makes this reliable.  cyclic-n has a few near-singular paths that
    # sit right at the edge of double-precision tracking tolerance; in double precision one of them
    # occasionally fails the endgame (MinStepSizeReached) and a genuine root is silently lost -- so
    # the distinct count comes back 1 or 2 short of the known value.  Adaptive precision raises the
    # working precision on exactly those paths, so every run finds all the finite solutions.

    start = time.time()
    if comm.Get_size() > 1:
        solver.solve(communicator=comm)           # manager (rank 0) hands paths to the workers
    else:
        solver.solve()                            # serial: a lone manager would have no workers
    elapsed = time.time() - start

    if not pb.parallel.is_manager():
        return

    # Correctness, not just bookkeeping.  The solver's own report classifies every path; we trust it
    # rather than rolling our own cutoff.  We print it whenever anything is flagged, so problems are
    # surfaced -- but we gate correctness on the ground truth, distinguishing two very different flags.
    report = solver.report()

    print('cyclic-{}:  ranks={}  threads/rank={}  paths tracked={}  finite solutions={}  wall={:.1f}s'.format(
        args.n, comm.Get_size(), os.environ.get('OMP_NUM_THREADS', '1'),
        report.num_paths_tracked, report.num_finite_solutions, elapsed))

    if not report.all_paths_resolved:        # a failed path OR an unresolved crossing -- show what and why
        print(report)

    assert report.num_paths_tracked == math.factorial(args.n), \
        'expected {} paths, got {}'.format(math.factorial(args.n), report.num_paths_tracked)
    if args.n in KNOWN_FINITE:
        # A FAILED path (the tracker gave up) means a genuinely lost root -- a hard error.  An
        # UNRESOLVED CROSSING is a conservative "may be wrong" warning the midpath check couldn't
        # clear; the count can still be exactly right (and here it is), so it is reported above, not
        # asserted.  Gate on no failed paths + the known finite count.
        assert report.num_failed == 0, 'some paths failed to track -- see the report above'
        assert report.num_finite_solutions == KNOWN_FINITE[args.n], \
            'expected {} finite solutions, got {}'.format(KNOWN_FINITE[args.n], report.num_finite_solutions)
        print('  OK: {} paths, {} finite solutions -- matches the known cyclic-{} count'.format(
            report.num_paths_tracked, report.num_finite_solutions, args.n))


if __name__ == '__main__':
    main()
