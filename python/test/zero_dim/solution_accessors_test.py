"""The filtered-solution convenience accessors: finite / real / nonsingular / singular_solutions.

The filtering logic is pinned in C++ (zero_dim.cpp, filtered_solution_accessors).  Here we verify
the binding: the counts, that each kind is a subset of the finite set, and that they agree with
solver.report().
"""

import bertini as pb
from bertini.nag_algorithm import ZeroDimSolver


def _one_var_solver(build):
    x = pb.Variable('x')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x]))
    sys.add_function(build(x))
    s = ZeroDimSolver(sys, mptype='double')
    s.solve()
    return s


def test_finite_real_nonsingular_roots():
    s = _one_var_solver(lambda x: x * x - 1)          # roots +/-1: finite, real, nonsingular
    assert len(s.finite_solutions()) == 2
    assert len(s.real_solutions()) == 2
    assert len(s.nonsingular_solutions()) == 2
    assert len(s.singular_solutions()) == 0


def test_complex_roots_are_finite_but_not_real():
    s = _one_var_solver(lambda x: x * x + 1)          # roots +/-i: finite, not real
    assert len(s.finite_solutions()) == 2
    assert len(s.real_solutions()) == 0
    assert len(s.nonsingular_solutions()) == 2
    assert len(s.singular_solutions()) == 0


def test_singular_double_root():
    s = _one_var_solver(lambda x: x * x)              # double root at 0: singular, real
    finite = len(s.finite_solutions())
    assert finite >= 1
    assert len(s.singular_solutions()) == finite
    assert len(s.nonsingular_solutions()) == 0
    assert len(s.real_solutions()) == finite          # 0 is real


def test_subsets_partition_and_agree_with_report():
    s = _one_var_solver(lambda x: x * x - 1)
    r = s.report()
    # singular and nonsingular partition the finite set
    assert len(s.singular_solutions()) + len(s.nonsingular_solutions()) == len(s.finite_solutions())
    # the accessors agree with the report's counts
    assert len(s.finite_solutions()) == r.num_finite_endpoints
    assert len(s.real_solutions()) == r.num_real
    assert len(s.singular_solutions()) == r.num_singular


def test_user_vs_internal_coordinates():
    s = _one_var_solver(lambda x: x * x - 1)
    # same count in either representation (the coordinates themselves generally differ)
    assert len(s.finite_solutions(user_coords=False)) == len(s.finite_solutions())
