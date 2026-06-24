"""Solution-access ergonomics on ZeroDim: all_solutions / infinite_solutions / to_dataframe.

These are binding-layer / Python-sugar interface tests (the numerical correctness of the solve
itself is gated by the C++ suite).  The system {x*y - 1, x - 1} has Bezout number 2 but exactly
ONE affine solution (x=1, y=1); the second total-degree path necessarily diverges, so it gives a
reliable mix of one finite and one at-infinity endpoint to exercise the category accessors.
"""

import sys

import pytest

import bertini as pb
from bertini.nag_algorithm import ZeroDim


@pytest.fixture
def two_circles_solver():
    """x^2 + y^2 - 1 = 0 and x + y = 0: two finite solutions, none at infinity."""
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_function(x**2 + y**2 - 1)
    sys.add_function(x + y)
    sys.add_variable_group(pb.VariableGroup([x, y]))
    solver = ZeroDim(sys)
    solver.solve()
    return solver


@pytest.fixture
def one_finite_one_infinite_solver():
    """x*y - 1 = 0 and x - 1 = 0: Bezout 2, one finite solution (1,1), one path to infinity."""
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_function(x * y - 1)
    sys.add_function(x - 1)
    sys.add_variable_group(pb.VariableGroup([x, y]))
    solver = ZeroDim(sys)
    solver.solve()
    return solver


# --- the rename: solutions() -> all_solutions() is a hard break (no deprecated alias) ---

def test_all_solutions_is_the_accessor(two_circles_solver):
    assert len(two_circles_solver.all_solutions()) == 2


def test_old_solutions_name_is_gone(two_circles_solver):
    """The rename is a hard break: the old name must not silently still work."""
    assert not hasattr(two_circles_solver, 'solutions')


# --- infinite_solutions: the at-infinity complement of finite_solutions ---

def test_infinite_solutions_complement(one_finite_one_infinite_solver):
    s = one_finite_one_infinite_solver
    assert len(s.all_solutions()) == 2
    assert len(s.finite_solutions()) == 1
    assert len(s.infinite_solutions()) == 1
    # no failed paths here, so finite + infinite partitions the whole list
    assert len(s.finite_solutions()) + len(s.infinite_solutions()) == len(s.all_solutions())


def test_no_infinite_solutions_when_all_finite(two_circles_solver):
    assert len(two_circles_solver.infinite_solutions()) == 0
    assert len(two_circles_solver.finite_solutions()) == 2


# --- to_dataframe: the "database of solutions" ---

def test_to_dataframe_finite_only_by_default(one_finite_one_infinite_solver):
    pd = pytest.importorskip("pandas")
    df = one_finite_one_infinite_solver.to_dataframe()
    # omit_infinite defaults True: just the one genuine finite solution
    assert len(df) == 1
    assert bool(df['is_finite'].all())
    # coordinates first, then metadata columns
    assert 'x0' in df.columns and 'x1' in df.columns
    for col in ('is_finite', 'is_real', 'is_singular', 'multiplicity', 'endgame_success'):
        assert col in df.columns


def test_to_dataframe_all_paths(one_finite_one_infinite_solver):
    pd = pytest.importorskip("pandas")
    df = one_finite_one_infinite_solver.to_dataframe(omit_infinite=False)
    assert len(df) == 2                       # every tracked path
    assert int(df['is_finite'].sum()) == 1    # one finite, one at infinity


def test_to_dataframe_filters_match_accessors(one_finite_one_infinite_solver):
    pd = pytest.importorskip("pandas")
    s = one_finite_one_infinite_solver
    df = s.to_dataframe(omit_infinite=False)
    assert len(df[df.is_finite]) == len(s.finite_solutions())
    assert len(df[~df.is_finite]) == len(s.infinite_solutions())


def test_to_dataframe_coords_match_points(two_circles_solver):
    pd = pytest.importorskip("pandas")
    df = two_circles_solver.to_dataframe()
    pts = two_circles_solver.finite_solutions()
    assert len(df) == len(pts)
    # the x0 column reproduces each finite solution's first coordinate
    key = lambda z: (z.real, z.imag)
    df_x0 = sorted((complex(v) for v in df['x0']), key=key)
    pt_x0 = sorted((complex(p[0]) for p in pts), key=key)
    for a, b in zip(df_x0, pt_x0):
        assert a == pytest.approx(b)


def test_to_dataframe_without_pandas_raises(two_circles_solver, monkeypatch):
    """When pandas is absent, to_dataframe raises a helpful ImportError (not the point accessors)."""
    monkeypatch.setitem(sys.modules, 'pandas', None)   # makes `import pandas` raise ImportError
    with pytest.raises(ImportError, match="pandas"):
        two_circles_solver.to_dataframe()
