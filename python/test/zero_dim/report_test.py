"""solver.report(): the end-of-solve diagnostic summary (SolveReport).

Correctness of the bucketing/failure logic is pinned in C++ (zero_dim.cpp,
solve_report_buckets_metadata_and_flags_failures).  Here we verify the binding: a clean solve's
fields, the human-readable __str__, and that failures_by_reason is keyed by the named SuccessCode.
"""

import bertini as pb
from bertini.nag_algorithm import ZeroDim


def _two_circles_solver():
    # x^2 - 1 = 0 and y^2 - 1 = 0: total degree 4, all four roots finite, none at infinity.
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_function(x**2 - 1)
    sys.add_function(y**2 - 1)
    sys.add_variable_group(pb.VariableGroup([x, y]))
    solver = ZeroDim(sys, mptype='double')
    solver.solve()
    return solver


def test_report_on_a_clean_solve():
    r = _two_circles_solver().report()

    assert r.num_paths_tracked == 4
    assert r.num_finite_solutions == 4
    assert r.num_diverged == 0
    assert r.num_failed == 0
    assert r.all_paths_resolved
    assert r.failures_by_reason == {}
    assert r.midpath.passed

    # the three buckets partition every path
    assert r.num_finite_endpoints + r.num_diverged + r.num_failed == r.num_paths_tracked


def test_report_str_is_human_readable():
    text = str(_two_circles_solver().report())
    assert text
    assert 'finite solutions' in text
    assert 'FAILED' in text
    assert 'all paths resolved' in text


def test_failures_keyed_by_named_success_code():
    # a clean solve has no failures, but the dict is keyed by the NAMED SuccessCode enum, never ints
    fails = _two_circles_solver().report().failures_by_reason
    assert isinstance(fails, dict)
    for code in fails:
        assert isinstance(code, pb.tracking.SuccessCode)
