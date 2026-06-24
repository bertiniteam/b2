"""The friendly ZeroDim(...) factory selects the right bound solver class by string.

Instead of typing ZeroDimCauchyAdaptivePrecisionTotalDegree, you say ZeroDim(system) (the
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


def test_defaults_are_cauchy_adaptive_totaldegree():
    solver = ZeroDim(_system())
    assert isinstance(solver, _n.ZeroDimCauchyAdaptivePrecisionTotalDegree)


@pytest.mark.parametrize("kwargs, expected", [
    (dict(),                                                   'ZeroDimCauchyAdaptivePrecisionTotalDegree'),
    (dict(mptype='double'),                                   'ZeroDimCauchyDoublePrecisionTotalDegree'),
    (dict(mptype='dbl'),                                      'ZeroDimCauchyDoublePrecisionTotalDegree'),
    (dict(mptype='multiple'),                                 'ZeroDimCauchyFixedMultiplePrecisionTotalDegree'),
    (dict(mptype='amp'),                                      'ZeroDimCauchyAdaptivePrecisionTotalDegree'),
    (dict(mptype='adaptive'),                                 'ZeroDimCauchyAdaptivePrecisionTotalDegree'),
    (dict(endgame='powerseries'),                             'ZeroDimPowerSeriesAdaptivePrecisionTotalDegree'),
    (dict(endgame='power_series'),                            'ZeroDimPowerSeriesAdaptivePrecisionTotalDegree'),
    (dict(startsystem='mhom'),                                'ZeroDimCauchyAdaptivePrecisionMHomogeneous'),
    (dict(startsystem='td'),                                  'ZeroDimCauchyAdaptivePrecisionTotalDegree'),
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
    assert len(solver.all_solutions()) == 2


def test_precision_is_an_alias_for_mptype():
    assert isinstance(ZeroDim(_system(), precision='amp'),
                      _n.ZeroDimCauchyAdaptivePrecisionTotalDegree)
    # precision overrides mptype when both are given
    assert isinstance(ZeroDim(_system(), mptype='double', precision='adaptive'),
                      _n.ZeroDimCauchyAdaptivePrecisionTotalDegree)


def test_user_startsystem_points_at_user_homotopy():
    with pytest.raises(ValueError, match='user_homotopy'):
        ZeroDim(_system(), startsystem='user')


def _two_affine_group_system():
    # x*y - 1, x + y over groups {x}, {y}: a multihomogeneous system (total degree would throw).
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x]))
    sys.add_variable_group(pb.VariableGroup([y]))
    sys.add_function(x * y - 1)
    sys.add_function(x + y)
    return sys


def _projective_system():
    x0, x1 = pb.Variable('x0'), pb.Variable('x1')
    sys = pb.System()
    sys.add_hom_variable_group(pb.VariableGroup([x0, x1]))
    sys.add_function(x0 * x0 - x1 * x1)
    return sys


def test_startsystem_inferred_from_variable_groups():
    # default startsystem='infer' mirrors the C++ blackbox InferStartType: a single affine group is
    # total degree; two-or-more groups, or any projective group, is multihomogeneous.  This matters
    # because total degree THROWS ("more than one affine variable group") on a multi-group system.
    assert isinstance(ZeroDim(_system()),
                      _n.ZeroDimCauchyAdaptivePrecisionTotalDegree)
    assert isinstance(ZeroDim(_two_affine_group_system()),
                      _n.ZeroDimCauchyAdaptivePrecisionMHomogeneous)
    assert isinstance(ZeroDim(_projective_system()),
                      _n.ZeroDimCauchyAdaptivePrecisionMHomogeneous)


def test_fixed_multiple_precision_solves():
    # Regression: a fixed-multiple zero-dim solve used to throw at the start of tracking --
    # "start point ... has differing precision from default (20!=16)" -- because the start points
    # (default precision) and the config-driven thread precision (double) disagreed.  The ambient
    # precision is now uniform (sourced from the multiprecision default), so a plain multiple-
    # precision solve completes.  See ZeroDimConfig::initial_ambient_precision.
    solver = ZeroDim(_system(), mptype='multiple')
    solver.solve()
    assert len(solver.finite_solutions()) == 2
