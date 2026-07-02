"""Solution-metadata classification and the PostProcessing config knobs that drive it.

Mirrors the C++ tests in ``core/test/nag_algorithms/zero_dim.cpp`` at the Python layer: the
per-solution flags (``is_finite`` / ``is_real`` / ``is_singular`` / ``multiplicity``) follow
Bertini 1's rules, and the ``PostProcessingConfig`` settings actually drive them -- which is the
regression that motivated the work (the flags used to never be computed, and the thresholds were
never applied).
"""

import pytest

import bertini as pb
from bertini import ZeroDimSolver

OK = int(pb.SuccessCode.Success)


def _one_var_solver(build_function):
    """A solver for the single-variable-group system { build_function(x) }."""
    x = pb.Variable('x')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x]))
    sys.add_function(build_function(x))
    return ZeroDimSolver(sys, endgame='cauchy', mptype='double', startsystem='binomial')


def _tally(solver):
    succ = [m for m in solver.solution_metadata() if int(m.endgame_success_code) == OK]
    return {
        'success': len(succ),
        'finite': sum(1 for m in succ if m.is_finite),
        'real': sum(1 for m in succ if m.is_real),
        'singular': sum(1 for m in succ if m.is_singular),
    }


def test_finite_real_nonsingular():
    solver = _one_var_solver(lambda x: x * x - 1)        # roots +/- 1
    solver.solve()
    assert _tally(solver) == {'success': 2, 'finite': 2, 'real': 2, 'singular': 0}
    assert all(m.multiplicity == 1 for m in solver.solution_metadata()
               if int(m.endgame_success_code) == OK)


def test_finite_complex_not_real():
    solver = _one_var_solver(lambda x: x * x + 1)        # roots +/- i
    solver.solve()
    c = _tally(solver)
    assert c['success'] == 2 and c['finite'] == 2 and c['real'] == 0 and c['singular'] == 0


def test_singular_double_root():
    solver = _one_var_solver(lambda x: x * x)            # double root at 0
    solver.solve()
    succ = [m for m in solver.solution_metadata() if int(m.endgame_success_code) == OK]
    assert len(succ) >= 1
    assert all(m.is_finite for m in succ)
    assert all(m.is_singular for m in succ)             # multiple and/or ill-conditioned


def test_endpoint_finite_threshold_is_applied():
    solver = _one_var_solver(lambda x: x * x - 1)        # roots +/- 1, infinity norm 1
    pp = solver.get_config(pb.nag_algorithm.PostProcessingConfig)
    pp.endpoint_finite_threshold = 0.5                   # 1 > 0.5 -> at infinity
    solver.set_config(pp)
    solver.solve()
    c = _tally(solver)
    assert c['success'] == 2 and c['finite'] == 0        # reclassified infinite by the lowered cutoff


def test_condition_number_threshold_is_applied():
    solver = _one_var_solver(lambda x: x * x - 1)        # simple, well-conditioned roots
    pp = solver.get_config(pb.nag_algorithm.PostProcessingConfig)
    pp.condition_number_threshold = 1e-3                 # any condition number exceeds this
    solver.set_config(pp)
    solver.solve()
    c = _tally(solver)
    assert c['success'] == 2 and c['singular'] == 2      # reclassified singular by the lowered threshold


def test_postprocessing_config_roundtrip():
    solver = _one_var_solver(lambda x: x * x - 1)
    pp = solver.get_config(pb.nag_algorithm.PostProcessingConfig)
    pp.endpoint_finite_threshold = 12345.0
    pp.same_point_tolerance_multiplier = 7.0
    pp.condition_number_threshold = 99.0
    pp.real_threshold = 1e-3
    solver.set_config(pp)
    back = solver.get_config(pb.nag_algorithm.PostProcessingConfig)
    assert back.endpoint_finite_threshold == 12345.0
    assert back.same_point_tolerance_multiplier == 7.0
    assert back.condition_number_threshold == 99.0
    assert back.real_threshold == 1e-3


def test_postprocessing_config_defaults_match_bertini1():
    pp = pb.nag_algorithm.PostProcessingConfig()
    assert pp.endpoint_finite_threshold == pytest.approx(1e5)
    assert pp.same_point_tolerance_multiplier == pytest.approx(10.0)
    assert pp.condition_number_threshold == pytest.approx(1e8)
    assert pp.real_threshold == pytest.approx(1e-8)
