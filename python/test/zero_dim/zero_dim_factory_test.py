"""The friendly ZeroDim(...) factory selects the right bound solver class by string.

Instead of typing ZeroDimCauchyFixedMultiplePrecisionTotalDegree, you say ZeroDim(system) (the
defaults) or ZeroDim(system, endgame=..., mptype=..., startsystem=...).
"""

import pytest

import bertini as pb
from bertini.nag_algorithm import ZeroDim
from bertini._pybertini import nag_algorithms as _n


def _system():
    # circle x^2 + y^2 - 1 meeting the line x = y -> two solutions (+/- 1/sqrt 2, +/- 1/sqrt 2).
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_function(x * x + y * y - 1)
    sys.add_function(x - y)
    return sys


def test_defaults_are_cauchy_multiple_totaldegree():
    solver = ZeroDim(_system())
    assert isinstance(solver, _n.ZeroDimCauchyFixedMultiplePrecisionTotalDegree)


@pytest.mark.parametrize("kwargs, expected", [
    (dict(),                                                   'ZeroDimCauchyFixedMultiplePrecisionTotalDegree'),
    (dict(mptype='double'),                                   'ZeroDimCauchyDoublePrecisionTotalDegree'),
    (dict(mptype='dbl'),                                      'ZeroDimCauchyDoublePrecisionTotalDegree'),
    (dict(mptype='multiple'),                                 'ZeroDimCauchyFixedMultiplePrecisionTotalDegree'),
    (dict(mptype='amp'),                                      'ZeroDimCauchyAdaptivePrecisionTotalDegree'),
    (dict(mptype='adaptive'),                                 'ZeroDimCauchyAdaptivePrecisionTotalDegree'),
    (dict(endgame='powerseries'),                             'ZeroDimPowerSeriesFixedMultiplePrecisionTotalDegree'),
    (dict(endgame='power_series'),                            'ZeroDimPowerSeriesFixedMultiplePrecisionTotalDegree'),
    (dict(startsystem='mhom'),                                'ZeroDimCauchyFixedMultiplePrecisionMHomogeneous'),
    (dict(startsystem='td'),                                  'ZeroDimCauchyFixedMultiplePrecisionTotalDegree'),
    # the docstring's headline example, and a fully-specified power-series/double/total-degree:
    (dict(endgame='cauchy', mptype='amp', startsystem='mhom'),'ZeroDimCauchyAdaptivePrecisionMHomogeneous'),
    (dict(endgame='power_series', mptype='dbl', startsystem='td'),
                                                              'ZeroDimPowerSeriesDoublePrecisionTotalDegree'),
])
def test_factory_selects_expected_class(kwargs, expected):
    solver = ZeroDim(_system(), **kwargs)
    assert type(solver).__name__ == expected
    assert isinstance(solver, getattr(_n, expected))


@pytest.mark.parametrize("bad", [
    dict(endgame='nope'),
    dict(mptype='quad'),
    dict(startsystem='bogus'),
    dict(startsystem='user'),   # user-homotopy is not built by ZeroDim -> use user_homotopy()
])
def test_factory_rejects_unknown_selectors(bad):
    with pytest.raises(ValueError):
        ZeroDim(_system(), **bad)


def test_factory_returns_a_real_solver():
    # the returned object is a genuine solver, not just a class lookup: solve end to end.
    solver = ZeroDim(_system(), mptype='amp')
    solver.solve()
    assert len(solver.solutions()) == 2


def test_precision_is_an_alias_for_mptype():
    assert isinstance(ZeroDim(_system(), precision='amp'),
                      _n.ZeroDimCauchyAdaptivePrecisionTotalDegree)
    # precision overrides mptype when both are given
    assert isinstance(ZeroDim(_system(), mptype='double', precision='adaptive'),
                      _n.ZeroDimCauchyAdaptivePrecisionTotalDegree)


def test_user_startsystem_points_at_user_homotopy():
    with pytest.raises(ValueError, match='user_homotopy'):
        ZeroDim(_system(), startsystem='user')
