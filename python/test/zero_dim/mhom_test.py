"""Multihomogeneous zero-dim solving from Python.

Exercises the products-of-linears MHom start system + blend-block homotopy through the bound
ZeroDim solver -- the Python side of the block-composed System work.

These pin a fixed seed (``set_random_seed`` -- RandomMp is reseedable now), so the homotopy is the
same every run: the tests pass or fail deterministically, with no gamma-retry loop.  (An earlier
40x retry-until-success loop, and then a smaller one, ballooned Windows CI.)
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


def _solve(solver_cls, seed=1):
    """Solve the two-group system with a FIXED seed, so the homotopy is reproducible."""
    pb.random.set_random_seed(seed)          # gamma + MHom start coefficients are now seedable
    solver = solver_cls(_two_group_system())
    solver.solve()
    return solver


def _successful_roots(solver):
    """User-coordinate solutions of paths that reached the endgame successfully."""
    sols = solver.all_solutions()
    md = solver.solution_metadata()
    return [sols[i] for i in range(len(sols))
            if int(md[i].endgame_success) == OK and len(sols[i]) == 2]


def test_mhom_solves_adaptive_precision():
    """AMP is the robust MHom path; with a fixed seed it solves this system deterministically."""
    good = _successful_roots(_solve(ZeroDimCauchyAdaptivePrecisionMHomogeneous))
    assert len(good) == 2, "adaptive-precision MHom did not solve"

    # both paths converged to genuine roots, and to the two *distinct* solutions
    for s in good:
        xv, yv = complex(s[0]), complex(s[1])
        assert abs(xv * yv - 1) < 1e-7
        assert abs(xv + yv) < 1e-7
    assert abs(complex(good[0][0]) - complex(good[1][0])) > 1e-3


@pytest.mark.parametrize("solver_cls", [
    ZeroDimCauchyDoublePrecisionMHomogeneous,
    ZeroDimPowerSeriesDoublePrecisionMHomogeneous,
])
def test_fixed_double_mhom_solves(solver_cls):
    """Fixed-double MHom solves this system at the fixed seed -- both paths reach the two genuine
    roots.  This used to only smoke-test that the binding ran, tolerating a tracking RuntimeError;
    that tolerance was masking the uninitialized patch-Jacobian bug (garbage Jacobian -> diverging
    corrector -> the tracker giving up).  With the patch Jacobian fixed, fixed double tracks
    cleanly -- so we assert the roots, no guard."""
    good = _successful_roots(_solve(solver_cls))
    assert len(good) == 2, "fixed-double MHom did not solve"

    for s in good:
        xv, yv = complex(s[0]), complex(s[1])
        assert abs(xv * yv - 1) < 1e-7
        assert abs(xv + yv) < 1e-7
    assert abs(complex(good[0][0]) - complex(good[1][0])) > 1e-3
