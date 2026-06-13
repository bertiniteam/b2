"""Multihomogeneous zero-dim solving from Python.

Exercises the products-of-linears MHom start system + blend-block homotopy through
the bound ZeroDim solver -- the Python side of the block-composed System work.

The fixed-double tracker is conditioning-fragile for MHom paths (a given random gamma
drives both paths to MinStepSize maybe two times in three), and RandomMp -- which feeds
the gamma and the MHom start coefficients -- is not yet reseedable, so the solve is not
seed-deterministic.  These are therefore *capability* tests: retry over gamma draws and
assert that MHom can solve, and that when it does the solutions are genuine roots.
"""

import pytest

import bertini as pb
from bertini.nag_algorithm import (
    ZeroDimCauchyDoublePrecisionMHomogeneous,
    ZeroDimPowerSeriesDoublePrecisionMHomogeneous,
    ZeroDimCauchyAdaptivePrecisionMHomogeneous,
)

# Compare success codes by integer value: the enhance machinery's config __eq__ can be
# reached through nested comparisons and trips over enum members, so avoid enum-vs-enum ==.
OK = int(pb.tracking.SuccessCode.Success)


def _two_group_system():
    """x*y - 1 = 0, x + y = 0 over variable groups {x}, {y}.

    y = -x gives -x^2 - 1 = 0, so x = +/- i -> the two solutions (i,-i), (-i,i).
    Bidegrees (1,1),(1,1) -> m-homogeneous Bezout number 2, below the total-degree
    Bezout number 4: MHom tracks 2 paths, not 4.
    """
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_function(x * y - 1)
    sys.add_function(x + y)
    sys.add_variable_group(pb.VariableGroup([x]))
    sys.add_variable_group(pb.VariableGroup([y]))
    return sys


def _successful_roots(solver):
    """User-coordinate solutions of paths that reached the endgame successfully."""
    sols = solver.solutions()
    md = solver.solution_metadata()
    return [sols[i] for i in range(len(sols))
            if int(md[i].endgame_success) == OK and len(sols[i]) == 2]


def _solve_with_retry(solver_cls, attempts=40):
    for _ in range(attempts):
        solver = solver_cls(_two_group_system())
        solver.solve()
        good = _successful_roots(solver)
        if len(good) == 2:
            return good
    return None


@pytest.mark.parametrize("solver_cls", [
    ZeroDimCauchyDoublePrecisionMHomogeneous,
    ZeroDimPowerSeriesDoublePrecisionMHomogeneous,
    ZeroDimCauchyAdaptivePrecisionMHomogeneous,  # adaptive: handles the harder MHom paths
])
def test_mhom_solves_two_variable_group_system(solver_cls):
    good = _solve_with_retry(solver_cls)
    assert good is not None, "MHom did not solve in the retry budget"

    # both paths converged to genuine roots, and to the two *distinct* solutions
    for s in good:
        xv, yv = complex(s[0]), complex(s[1])
        assert abs(xv * yv - 1) < 1e-7
        assert abs(xv + yv) < 1e-7
    assert abs(complex(good[0][0]) - complex(good[1][0])) > 1e-3
