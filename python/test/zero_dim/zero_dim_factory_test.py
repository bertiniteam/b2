"""The friendly ZeroDimSolver(...) factory selects the right bound solver class by string.

Instead of typing ZeroDimSolverCauchyAdaptivePrecision, you say ZeroDimSolver(system) (the defaults) or
ZeroDimSolver(system, endgame=..., mptype=..., startsystem=...).

Note: after ZeroDimSolver was de-templated off the start-system type, the bound solver class encodes only
the endgame and precision (e.g. ``ZeroDimSolverCauchyAdaptivePrecision``); the start system is chosen at
construction and held polymorphically, so it is NO LONGER part of the class name.  Start-system
selection is therefore verified behaviorally here, not by the type.
"""

import pytest

import bertini as pb
from bertini import ZeroDimSolver
from bertini._pybertini import nag_algorithms as _n


def _system():
    # circle x^2 + y^2 - 1 meeting the line x = y -> two solutions (+/- 1/sqrt 2, +/- 1/sqrt 2).
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_function(x * x + y * y - 1)
    sys.add_function(x - y)
    return sys


def test_defaults_are_cauchy_adaptive():
    # default endgame + precision; start system (total degree) is not in the type name anymore.
    solver = ZeroDimSolver(_system())
    assert isinstance(solver, _n.ZeroDimSolverCauchyAdaptivePrecision)


@pytest.mark.parametrize("kwargs, expected", [
    (dict(),                                                   'ZeroDimSolverCauchyAdaptivePrecision'),
    (dict(mptype='double'),                                   'ZeroDimSolverCauchyDoublePrecision'),
    (dict(mptype='dbl'),                                      'ZeroDimSolverCauchyDoublePrecision'),
    (dict(mptype='multiple'),                                 'ZeroDimSolverCauchyFixedMultiplePrecision'),
    (dict(mptype='amp'),                                      'ZeroDimSolverCauchyAdaptivePrecision'),
    (dict(mptype='adaptive'),                                 'ZeroDimSolverCauchyAdaptivePrecision'),
    (dict(endgame='powerseries'),                             'ZeroDimSolverPowerSeriesAdaptivePrecision'),
    (dict(endgame='power_series'),                            'ZeroDimSolverPowerSeriesAdaptivePrecision'),
    # the start system no longer changes the type -- only endgame + precision do:
    (dict(startsystem='mhom'),                                'ZeroDimSolverCauchyAdaptivePrecision'),
    (dict(startsystem='binomial'),                            'ZeroDimSolverCauchyAdaptivePrecision'),
    (dict(endgame='cauchy', mptype='amp', startsystem='mhom'),'ZeroDimSolverCauchyAdaptivePrecision'),
    (dict(endgame='power_series', mptype='dbl', startsystem='linearproduct'),
                                                              'ZeroDimSolverPowerSeriesDoublePrecision'),
])
def test_factory_selects_expected_class(kwargs, expected):
    solver = ZeroDimSolver(_system(), **kwargs)
    assert type(solver).__name__ == expected
    assert isinstance(solver, getattr(_n, expected))


@pytest.mark.parametrize("bad", [
    dict(endgame='nope'),
    dict(mptype='quad'),
    dict(startsystem='bogus'),
    dict(startsystem='user'),   # user-homotopy is not built by ZeroDimSolver -> use user_homotopy()
])
def test_factory_rejects_unknown_selectors(bad):
    with pytest.raises(ValueError):
        ZeroDimSolver(_system(), **bad)


def test_factory_returns_a_real_solver():
    # the returned object is a genuine solver, not just a class lookup: solve end to end.
    solver = ZeroDimSolver(_system(), mptype='amp')
    solver.solve()
    assert len(solver.all_solutions()) == 2


def test_precision_is_an_integer_number_of_digits():
    import bertini as b
    # precision= is now a DIGIT COUNT, applied via default_precision -- not the model selector
    b.default_precision(16)
    solver = ZeroDimSolver(_system(), mptype='multiple', precision=80)
    assert b.default_precision() == 80                     # the digit count took effect
    assert isinstance(solver, _n.ZeroDimSolverCauchyFixedMultiplePrecision)   # mptype picked the class


def test_precision_as_a_string_is_the_deprecated_model_alias():
    import warnings
    # a STRING precision is the old model-selector alias: honored, but warns
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        solver = ZeroDimSolver(_system(), precision='amp')
    assert isinstance(solver, _n.ZeroDimSolverCauchyAdaptivePrecision)
    assert any(issubclass(x.category, DeprecationWarning) for x in w)


def test_user_startsystem_points_at_homotopy_solver():
    with pytest.raises(ValueError, match='HomotopySolver'):
        ZeroDimSolver(_system(), startsystem='user')


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
    # total degree; two-or-more groups, or any projective group, is multihomogeneous.  Since the start
    # system is no longer in the solver type, we verify the inference BEHAVIORALLY: a multi-group /
    # projective system would THROW ("more than one affine variable group") if total degree were
    # (wrongly) inferred, so successfully constructing AND solving it proves MHom was chosen.
    single = ZeroDimSolver(_system())
    single.solve()
    assert len(single.finite_solutions()) == 2

    multi = ZeroDimSolver(_two_affine_group_system())     # total degree would throw here
    multi.solve()
    assert len(multi.finite_solutions()) == 2       # x*y=1, x+y=0 -> (i,-i) and (-i,i)

    proj = ZeroDimSolver(_projective_system())            # projective group -> MHom
    proj.solve()
    assert len(proj.all_solutions()) >= 1


def test_fixed_multiple_precision_solves():
    # Regression: a fixed-multiple zero-dim solve used to throw at the start of tracking --
    # "start point ... has differing precision from default (20!=16)" -- because the start points
    # (default precision) and the config-driven thread precision (double) disagreed.  The ambient
    # precision is now uniform (sourced from the multiprecision default), so a plain multiple-
    # precision solve completes.  See ZeroDimConfig::initial_ambient_precision.
    solver = ZeroDimSolver(_system(), mptype='multiple')
    solver.solve()
    assert len(solver.finite_solutions()) == 2
