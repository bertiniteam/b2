"""Cross-software regression: bertini2's cyclic-5 finite solutions match an independent oracle.

The oracle is data/cyclic5_finite_solutions.txt -- the 70 finite solutions computed by Bertini 1.7.0
(a separate solver) with slightly-tighter-than-default tolerances.  The finite solutions are points
of the variety; they are invariant under the choice of random homotopy, so two solvers' results are
the same SET even though their paths (different random gammas) are not.  Comparing the set is what
makes this a meaningful cross-software check: it catches a silently lost root (the count comes back
short) or a wrong endpoint -- the exact failure ADR-0017 / the cyclic-5 diagnostic warned about.

To regenerate the oracle, see the provenance header in the data file.
"""

from pathlib import Path

import numpy as np
import pytest

import bertini as pb
from bertini.nag_algorithm import TolerancesConfig


N = 5
KNOWN_FINITE = 70
# The two solvers each track to their final tolerance, so the same root differs between them by about
# that much.  Match within 10 * (our final_tolerance) under the infinity norm -- tight (~1e-10 at the
# default final_tolerance of 1e-11), so a wrong or mis-tracked endpoint fails, not just a missing one.
FINAL_TOL_MATCH_FACTOR = 10
DATA = Path(__file__).parent / "data" / "cyclic5_finite_solutions.txt"


def cyclic_system(n):
    """The cyclic-n system (same construction as examples/solve_cyclic.cyclic_system)."""
    x = [pb.Variable('x{}'.format(i)) for i in range(n)]
    w = x + x
    sys = pb.System()
    for length in range(1, n):
        sys.add_function(np.sum([np.prod(w[start:start + length]) for start in range(n)]))
    sys.add_function(np.prod(x) - 1)
    sys.add_variable_group(pb.VariableGroup(x))
    return sys


def load_trusted():
    rows = []
    for line in DATA.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        f = [float(v) for v in line.split()]
        rows.append([f[2 * k] + 1j * f[2 * k + 1] for k in range(N)])
    return np.array(rows)


def solve_finite():
    pb.random.set_random_seed(1)   # deterministic homotopy for a stable test
    solver = pb.nag_algorithm.ZeroDimSolver(cyclic_system(N), endgame='cauchy', mptype='adaptive', startsystem='binomial')
    final_tol = float(solver.get_config(TolerancesConfig).final_tolerance)
    solver.solve()
    ours = np.array([[complex(c) for c in s] for s in solver.finite_solutions()])
    return ours, FINAL_TOL_MATCH_FACTOR * final_tol


def test_trusted_set_is_intact():
    trusted = load_trusted()
    assert trusted.shape == (KNOWN_FINITE, N)


# The adaptive cyclic-5 solve is the same heavy solve as the solve_cyclic example: on a slow Windows
# CI box it can exceed the global 180s pytest-timeout, so override it here (mirrors the example test).
@pytest.mark.timeout(600)
def test_b2_finite_solutions_match_the_oracle():
    trusted = load_trusted()
    ours, match_tol = solve_finite()

    assert len(ours) == KNOWN_FINITE, \
        "expected {} finite solutions, got {}".format(KNOWN_FINITE, len(ours))

    # Bijection: each oracle point is matched by exactly one of ours, within match_tol (= 10 *
    # final_tolerance, infinity norm).  Greedy nearest matching is valid because the 70 roots are
    # well separated (>> match_tol apart).
    unused = list(range(len(ours)))
    for t in trusted:
        dists = [max(abs(ours[j] - t)) for j in unused]   # infinity norm over the 5 coordinates
        k = int(np.argmin(dists))
        assert dists[k] < match_tol, \
            "no computed solution within {:.1e} of an oracle solution (closest {:.2e})".format(match_tol, dists[k])
        unused.pop(k)
    assert not unused, "computed solutions left over after matching the oracle -- duplicates or extras"
