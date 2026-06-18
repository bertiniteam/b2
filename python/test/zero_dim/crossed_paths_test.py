"""Path-crossing detection and the endgame-boundary re-track machinery, at the Python layer.

Two *distinct* paths landing on the same point at the endgame boundary is a probability-0 event:
it is a path crossing (under-resolved tracking), which the zero-dim solver is supposed to detect
and re-track.  Mirrors the C++ tests in ``core/test/nag_algorithms/zero_dim.cpp``; here we exercise
the whole pipeline and the ``endgame_boundary_metadata`` report.

To *provoke* a crossing we deliberately use the low-order Euler predictor with a loose pre-endgame
tolerance on the cyclic-5 system -- the regime documented in ADR-0017.  The fixed seed only makes
the bad behaviour reproducible; the *fix* is the solver's re-tracking (and, in real use, the better
default predictor / tighter tolerance), never the seed (see ADR-0017).  Because exactly which seed
crosses depends on platform numerics, the test searches a small range for one rather than hard-coding
it.
"""

import numpy as np
import pytest

import bertini as pb
from bertini.nag_algorithm import ZeroDimCauchyDoublePrecisionTotalDegree

OK = int(pb.tracking.SuccessCode.Success)

CYCLIC_N = 5
KNOWN_FINITE = 70           # distinct finite solutions of cyclic-5
LOOSE_TOL = 1e-4            # loose enough (with Euler) to occasionally cross a pair of paths
SEEDS_TO_TRY = range(1, 13)


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


def _solve(seed, resolve_attempts):
    """Solve cyclic-5 with Euler + loose tolerance at a fixed seed.

    Returns (report, distinct_finite_count).
    """
    pb.random.set_random_seed(seed)
    solver = ZeroDimCauchyDoublePrecisionTotalDegree(cyclic_system(CYCLIC_N))

    # Euler (order 1) is the root cause of the crossing we want to provoke.
    solver.get_tracker().predictor(pb.tracking.Predictor.Euler)

    tol = solver.get_config(pb.nag_algorithm.TolerancesConfig)
    tol.newton_before_endgame = LOOSE_TOL
    tol.newton_during_endgame = LOOSE_TOL / 10
    solver.set_config(tol)

    zd = solver.get_config(pb.nag_algorithm.ZeroDimConfigDoublePrec)
    zd.max_num_crossed_path_resolve_attempts = resolve_attempts
    solver.set_config(zd)

    solver.solve()

    report = solver.endgame_boundary_metadata()
    finite = [m for m in solver.solution_metadata()
              if int(m.endgame_success) == OK and m.is_finite]
    distinct = round(sum(1.0 / m.multiplicity for m in finite))
    return report, distinct


def _find_crossing_seed():
    """Find a seed whose report-only solve detects an (unresolved) crossing; skip if none."""
    for seed in SEEDS_TO_TRY:
        report, distinct = _solve(seed, resolve_attempts=0)
        if report.num_crossings_detected > 0:
            return seed, report, distinct
    pytest.skip("no path crossing provoked in the searched seed range on this platform; "
                "the detection/resolve logic is still covered by the deterministic C++ tests")


@pytest.mark.timeout(900)  # several cyclic-5 Euler solves; loose Euler tracking is slow, and the
                           # manylinux test container is markedly slower per-core than 360s allowed
                           # (it timed out there while passing on macOS) -- this is a safety net
                           # against a true hang, not a performance assertion.
def test_crossing_is_detected_then_resolved():
    """Report-only loses a solution to the crossing; re-tracking recovers it."""
    seed, report_only, distinct_unresolved = _find_crossing_seed()

    # With re-tracking disabled the crossing is detected but left unresolved, and the two merged
    # paths cost us a solution: the distinct finite count comes up short.
    assert not report_only.passed
    assert report_only.num_resolve_attempts == 0
    assert len(report_only.crossed_path_indices) == report_only.num_crossings_detected
    assert distinct_unresolved < KNOWN_FINITE

    # Re-enable the default re-tracking: the same crossing is detected, then resolved, and the lost
    # solution is recovered -- the count is now exactly right.
    resolved, distinct_resolved = _solve(seed, resolve_attempts=2)
    assert resolved.num_crossings_detected > 0
    assert resolved.num_resolve_attempts >= 1
    assert resolved.passed
    assert distinct_resolved == KNOWN_FINITE


def test_clean_solve_reports_no_crossings():
    """With the good default predictor and a sane tolerance, no crossing is detected."""
    pb.random.set_random_seed(1)
    solver = ZeroDimCauchyDoublePrecisionTotalDegree(cyclic_system(CYCLIC_N))
    tol = solver.get_config(pb.nag_algorithm.TolerancesConfig)
    tol.newton_before_endgame = 1e-6
    tol.newton_during_endgame = 1e-7
    solver.set_config(tol)
    solver.solve()

    report = solver.endgame_boundary_metadata()
    assert report.passed
    assert report.num_crossings_detected == 0
    assert report.num_resolve_attempts == 0
    assert len(solver.endgame_boundary_solutions()) > 0
