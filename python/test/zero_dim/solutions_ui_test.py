"""Solution-surface UI: dedup by default (#299), metadata-by-point (#302), tolerance point
comparison (#304), and group projection of solutions (Cluster G)."""

import numpy as np
import pytest

import bertini as pb
from bertini import ZeroDimSolver


def _solve(functions, groups):
    sys = pb.System()
    for g in groups:
        sys.add_variable_group(g)
    sys.add_functions(functions)
    solver = ZeroDimSolver(sys)
    solver.solve()
    return sys, solver


@pytest.fixture
def double_root():
    # {x^2, y^2} -> the single solution (0,0) of multiplicity 4
    x, y = pb.variables(['x', 'y'])
    sys, solver = _solve([x * x, y * y], [[x, y]])
    return sys, solver


@pytest.fixture
def four_simple_roots():
    # x^2 - 1, y^2 - 1 -> four distinct simple roots (+-1, +-1)
    x, y = pb.variables(['x', 'y'])
    sys, solver = _solve([x * x - 1, y * y - 1], [[x, y]])
    return sys, solver


# --- #299: dedup multiplicities by default ---------------------------------------------------

def test_solutions_dedup_multiplicities_by_default(double_root):
    _, solver = double_root
    assert len(solver.finite_solutions()) == 1                          # merged (the default)
    assert len(solver.finite_solutions(merge_multiplicities=False)) == 4
    assert len(solver.solutions()) == 1
    assert len(solver.solutions(merge_multiplicities=False)) == 4


def test_simple_roots_unaffected_by_merge(four_simple_roots):
    _, solver = four_simple_roots
    assert len(solver.finite_solutions()) == 4
    assert len(solver.finite_solutions(merge_multiplicities=False)) == 4


# --- #302: metadata_for(point) ---------------------------------------------------------------

def test_metadata_for_point_returns_representative(double_root):
    _, solver = double_root
    pt = solver.finite_solutions()[0]                                   # (0,0)
    md = solver.metadata_for(pt, tol=1e-5)
    assert md.multiplicity == 4                                         # a single record, learns m
    assert md.multiplicity_representative


def test_metadata_for_point_coincident_is_a_list(double_root):
    _, solver = double_root
    pt = solver.finite_solutions()[0]
    all_md = solver.metadata_for(pt, tol=1e-5, coincident=True)
    assert isinstance(all_md, list)
    assert len(all_md) == 4                                             # all coincident copies


def test_metadata_for_missing_point_raises(double_root):
    _, solver = double_root
    far = solver.finite_solutions()[0].copy()
    for i in range(len(far)):
        far[i] = far[i] + 1000
    with pytest.raises(RuntimeError):
        solver.metadata_for(far, tol=1e-5)


# --- #304: is_distinct_up_to on solution points ----------------------------------------------

def test_is_distinct_up_to_tells_solutions_apart(four_simple_roots):
    _, solver = four_simple_roots
    sols = solver.finite_solutions()
    assert pb.is_distinct_up_to(sols[0], sols[1], 1e-6)                 # different roots
    assert not pb.is_distinct_up_to(sols[0], sols[0], 1e-6)            # same point


# --- Cluster G: project solutions onto one variable group ------------------------------------

def test_group_projection_of_solutions():
    # two affine groups [x] and [y]; x^2-1, y^2-1 -> (+-1, +-1)
    x, y = pb.variables(['x', 'y'])
    sys, solver = _solve([x * x - 1, y * y - 1], [[x], [y]])
    full = solver.finite_solutions()
    assert len(full) == 4 and all(len(p) == 2 for p in full)

    # by FIFO index
    g0 = solver.finite_solutions(group=0)
    assert all(len(p) == 1 for p in g0)                                # just the x coordinate

    # by VariableGroup object
    vg1 = sys.variable_groups()[1]
    g1 = solver.finite_solutions(group=vg1)
    assert all(len(p) == 1 for p in g1)

    # the x-coordinates (deduped over the group) are {+1, -1}
    xs = sorted({round(float(p[0].real)) for p in g0})
    assert xs == [-1, 1]


def test_group_projection_in_to_dataframe():
    pd = pytest.importorskip('pandas')
    x, y = pb.variables(['x', 'y'])
    sys, solver = _solve([x * x - 1, y * y - 1], [[x], [y]])
    df = solver.to_dataframe(group=0)
    assert len(df) == 4
    assert all(len(sol) == 1 for sol in df.solution)                   # solution cell holds only group 0
