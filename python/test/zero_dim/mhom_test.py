"""Multihomogeneous zero-dim solving from Python.

Exercises the products-of-linears MHom start system + blend-block homotopy through
the bound ZeroDim solver -- the Python side of the block-composed System work.

The fixed-double tracker is conditioning-fragile for MHom paths (a given random gamma drives
both paths to MinStepSize maybe two times in three), and RandomMp -- which feeds the gamma and
the MHom start coefficients -- is not yet reseedable.  So *adaptive precision* (AMP) is the path
we require to solve; the fixed-double bindings only get a cheap "it runs" smoke.  We deliberately
do NOT retry over many gamma draws here: a 40x retry loop ran ~40 full solves per variant and
made the Windows CI run balloon.
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


def _solve(solver_cls):
    solver = solver_cls(_two_group_system())
    solver.solve()
    return solver


def _successful_roots(solver):
    """User-coordinate solutions of paths that reached the endgame successfully."""
    sols = solver.solutions()
    md = solver.solution_metadata()
    return [sols[i] for i in range(len(sols))
            if int(md[i].endgame_success) == OK and len(sols[i]) == 2]


def test_mhom_solves_adaptive_precision():
    """AMP is the robust MHom path -- require it to solve.  A tiny budget (3) absorbs the rare
    unlucky gamma without the retry storms that made fixed-double MHom dominate CI."""
    good = []
    for _ in range(3):
        good = _successful_roots(_solve(ZeroDimCauchyAdaptivePrecisionMHomogeneous))
        if len(good) == 2:
            break
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
def test_fixed_double_mhom_runs(solver_cls):
    """Fixed-double MHom is conditioning-fragile, so we do NOT assert it finds the roots -- that's
    AMP's job above.  This is a cheap binding/run smoke: a single solve (no gamma retries).  An
    unlucky gamma can make the fixed-double tracker fail (it raises) -- which is itself fine for a
    smoke test (the binding ran), so we accept either a clean attempt of the two MHom paths or a
    tracking failure; we only require it not to crash the interpreter."""
    solver = solver_cls(_two_group_system())
    try:
        solver.solve()
    except RuntimeError:
        return  # fixed-double tracking gave up on this gamma; the binding still ran fine
    assert len(solver.solutions()) == 2  # the m-homogeneous Bezout number: two paths attempted
