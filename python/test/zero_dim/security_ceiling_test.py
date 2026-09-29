"""The security ceiling is one rule, applied alike by both endgames.

At security level 0, two consecutive endpoint approximations whose dehomogenized infinity norm is
above ``max_norm`` truncate the path.  The Cauchy endgame used to test acceptance before security,
so it returned Success for finite roots above the ceiling, which the power series endgame
truncated.  The two endgames differ in how they estimate the root and in nothing else.

``{x^2 - 100, y - 2x}`` has the finite roots ``(+-10, +-20)``; the ceiling is lowered to 5 to put
them above it, because the double precision tracker is not reliable on a badly scaled system (it
runs out of steps before the endgame on ``x^2 = 1e6`` always, and on ``x^2 = 1e4`` in some
draws).  The last test is the case as a user meets it: roots at ``(+-1e4, +-2e4)`` and the
default ceiling, under adaptive precision.
"""

import pytest

import bertini as pb
from bertini import ZeroDimSolver

ENDGAMES = ['powerseries', 'cauchy']
PRECISION_MODELS = ['adaptive', 'double']


def _solve(endgame, mptype, scale=100, **settings):
    x, y = pb.variables(['x', 'y'])
    system = pb.System()
    system.add_variable_group([x, y])
    system.add_functions([x ** 2 - pb.coefficient(scale), y - pb.coefficient(2) * x])
    solver = ZeroDimSolver(system, endgame=endgame, mptype=mptype)
    if settings:
        solver.update(**settings)
    solver.solve(show_progress=False)
    codes = [str(m.endgame_success_code) for m in solver.solution_metadata()]
    return codes, solver.solutions()


@pytest.mark.parametrize('mptype', PRECISION_MODELS)
@pytest.mark.parametrize('endgame', ENDGAMES)
def test_roots_above_the_ceiling_are_truncated(endgame, mptype):
    codes, finite = _solve(endgame, mptype, max_norm=5)
    assert codes == ['SecurityMaxNormReached'] * 2
    assert len(finite) == 0


@pytest.mark.parametrize('mptype', PRECISION_MODELS)
@pytest.mark.parametrize('endgame', ENDGAMES)
def test_roots_under_the_ceiling_are_kept(endgame, mptype):
    codes, finite = _solve(endgame, mptype)                 # default max_norm, 1e4
    assert codes == ['Success'] * 2
    assert len(finite) == 2
    assert sorted(round(abs(complex(p[0]))) for p in finite) == [10, 10]


@pytest.mark.parametrize('mptype', PRECISION_MODELS)
@pytest.mark.parametrize('endgame', ENDGAMES)
def test_security_level_one_does_not_truncate(endgame, mptype):
    codes, finite = _solve(endgame, mptype, max_norm=5, level=1)
    assert codes == ['Success'] * 2
    assert len(finite) == 2


@pytest.mark.parametrize('endgame', ENDGAMES)
def test_roots_above_the_default_ceiling_are_truncated(endgame):
    codes, finite = _solve(endgame, 'adaptive', scale=100000000)   # roots at (+-1e4, +-2e4)
    assert codes == ['SecurityMaxNormReached'] * 2
    assert len(finite) == 0
    codes, finite = _solve(endgame, 'adaptive', scale=100000000, max_norm=1e6)
    assert codes == ['Success'] * 2
    assert len(finite) == 2
