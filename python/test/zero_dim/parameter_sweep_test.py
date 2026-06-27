"""parameter_sweep: one ab-initio solve, track to many coefficient values.

Serial/threaded coverage here; the MPI (comm=) path is exercised by
examples/parallel_parameter_homotopy.py under mpirun (see the MPI CI smoke job).
"""

import pytest

import bertini as pb
from bertini.nag_algorithm import parameter_sweep

C = pb.multiprec.Complex

# Shared variables across the family -- the homotopy interpolates coefficients, so every member
# must be built over the same Variable objects.
x, y = pb.Variable('x'), pb.Variable('y')


def _make_system(params):
    a, b = params
    s = pb.System()
    s.add_variable_group(pb.VariableGroup([x, y]))
    s.add_function(x * x - a)        # x^2 = a
    s.add_function(y * y - b)        # y^2 = b
    return s


_GENERIC = (C("0.371", "0.59"), C("-0.83", "0.21"))     # generic complex


def test_returns_one_solver_per_target():
    targets = [(C("4", "0"), C("9", "0")), (C("2", "0"), C("5", "0"))]
    solvers = parameter_sweep(_make_system, _GENERIC, targets)
    assert len(solvers) == len(targets)
    for solver in solvers:
        assert len(solver.all_solutions()) == 4        # x^2=a, y^2=b -> 4 solutions


def test_collect_reduces_each_solve():
    targets = [(C("4", "0"), C("9", "0")), (C("2", "0"), C("5", "0"))]
    counts = parameter_sweep(_make_system, _GENERIC, targets,
                             collect=lambda s: len(s.all_solutions()))
    assert counts == [4, 4]                             # collect-form returns values, not solvers


def test_sweep_matches_a_direct_solve():
    """A swept target must give the same solution set as solving that member from scratch."""
    from bertini.nag_algorithm import ZeroDim

    target = (C("4", "0"), C("9", "0"))
    swept = parameter_sweep(_make_system, _GENERIC, [target])[0]

    direct = ZeroDim(_make_system(target), mptype='adaptive')
    direct.solve()

    def key(solver):
        ks = set()
        for v in solver.all_solutions():
            ks.add(tuple(sorted((round(complex(c).real, 5), round(complex(c).imag, 5)) for c in v)))
        return ks

    assert key(swept) == key(direct)


def test_comm_without_collect_is_an_error():
    """Passing a communicator but no collect must fail fast (solvers can't cross MPI ranks)."""
    class _FakeComm:                          # never actually used -- the guard fires first
        def Get_rank(self): return 0
        def Get_size(self): return 1
    with pytest.raises(ValueError):
        parameter_sweep(_make_system, _GENERIC, [(C("4", "0"), C("9", "0"))], comm=_FakeComm())
