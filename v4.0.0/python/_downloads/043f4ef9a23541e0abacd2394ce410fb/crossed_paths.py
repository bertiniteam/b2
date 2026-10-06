#!/usr/bin/env python
"""Provoke a path crossing on purpose, and watch Bertini 2 detect and re-track it.

Two *distinct* homotopy paths reaching the **same** point at the endgame boundary is a
probability-0 event -- the homotopy's gamma is random, so generic paths never meet.  When it does
happen it is a **path crossing**: a symptom that the tracker under-resolved the paths (most often a
too-low-order predictor or a too-loose tolerance), not a benign duplicate.  The zero-dim solver
checks for this at the endgame boundary and, by default, re-tracks the offending paths with
tightened settings.

To see the machinery work we deliberately break the tracking: the **Euler** predictor (order 1)
with a **loose** pre-endgame tolerance, on the cyclic-5 system.  The fixed seed only makes the bad
behaviour reproducible -- it is NOT the fix.  The fix is the solver's re-tracking (and, in normal
use, the better default predictor / a tighter tolerance).

Run it::

    python crossed_paths.py            # uses a seed known to cross on the author's machine
    python crossed_paths.py --seed 7   # try another seed

Which seed crosses depends on platform numerics; if the chosen one does not cross, the script says
so and suggests trying another.
"""

import argparse
import numpy as np

import bertini as pb

CYCLIC_N = 5
KNOWN_FINITE = 70            # distinct finite solutions of cyclic-5
OK = int(pb.SuccessCode.Success)


def cyclic_system(n):
    """The cyclic-n system: the n-1 cyclic sums of products, plus (prod - 1)."""
    x = [pb.Variable('x{}'.format(i)) for i in range(n)]
    w = x + x
    sys = pb.System()
    for length in range(1, n):
        sys.add_function(np.sum([np.prod(w[start:start + length]) for start in range(n)]))
    sys.add_function(np.prod(x) - 1)
    sys.add_variable_group(pb.VariableGroup(x))
    return sys


def solve(seed, resolve_attempts, tol=1e-4):
    """Solve cyclic-5 with the Euler predictor + a loose tolerance at a fixed seed.

    ``resolve_attempts`` is how many times the solver may re-track crossed paths at the endgame
    boundary; 0 means "detect and report only".  Returns (report, distinct_finite_count).
    """
    pb.random.set_random_seed(seed)
    solver = pb.ZeroDimSolver(cyclic_system(CYCLIC_N), endgame='cauchy', mptype='double', startsystem='binomial')

    # The deliberate culprit: a first-order predictor.
    solver.get_tracker().predictor(pb.Predictor.Euler)

    tols = solver.get_config(pb.nag_algorithm.TolerancesConfig)
    tols.newton_before_endgame = tol
    tols.newton_during_endgame = tol / 10
    solver.set_config(tols)

    zd = solver.get_config(pb.nag_algorithm.ZeroDimConfig)
    zd.max_num_crossed_path_resolve_attempts = resolve_attempts
    solver.set_config(zd)

    solver.solve()

    report = solver.endgame_boundary_metadata()
    finite = [m for m in solver.solution_metadata()
              if int(m.endgame_success_code) == OK and m.is_finite]
    distinct = round(sum(1.0 / m.multiplicity for m in finite))
    return report, distinct


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--seed', type=int, default=3,
                        help='RNG seed (default 3; chosen to cross on the author\'s machine)')
    args = parser.parse_args()

    # First, turn OFF re-tracking, so we can see the damage a crossing does on its own.
    report, distinct = solve(args.seed, resolve_attempts=0)
    print('Euler + loose tolerance, re-tracking DISABLED (max_num_crossed_path_resolve_attempts=0):')
    print('  crossings detected : {}'.format(report.num_crossings_detected))
    print('  crossed path indices: {}'.format(list(report.crossed_path_indices)))
    print('  midpath check passed: {}'.format(report.passed))
    print('  distinct finite solutions: {} / {}'.format(distinct, KNOWN_FINITE))

    if report.num_crossings_detected == 0:
        print('\nNo crossing was provoked at seed {} on this machine.  Try another --seed.'.format(args.seed))
        return

    print('\n  --> the crossed paths merged into one, so a true solution was lost'
          if distinct < KNOWN_FINITE else
          '\n  --> a crossing was detected (the endpoints happened to still be recoverable)')

    # Now the default behaviour: re-track the crossed paths with tightened settings.
    report, distinct = solve(args.seed, resolve_attempts=2)
    print('\nSame solve, re-tracking ENABLED (the default):')
    print('  crossings detected : {}'.format(report.num_crossings_detected))
    print('  re-track attempts  : {}'.format(report.num_resolve_attempts))
    print('  midpath check passed: {}'.format(report.passed))
    print('  distinct finite solutions: {} / {}'.format(distinct, KNOWN_FINITE))

    assert distinct == KNOWN_FINITE, \
        'expected the re-track to recover all {} solutions'.format(KNOWN_FINITE)
    print('\n  --> the crossing was resolved and the lost solution recovered.')
    print('      The real lesson: use a higher-order predictor (the default RKF45) or a tighter')
    print('      tolerance, and the crossing never happens in the first place.')


if __name__ == '__main__':
    main()
