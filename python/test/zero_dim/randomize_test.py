"""Randomizing an overdetermined system, solving the square result, and filtering.

An overdetermined system (more functions than variables) cannot be solved directly -- a
total-degree start system needs a square system.  bertini.linalg.randomize (System.Randomize)
squares it up with a generic combination R*F whose isolated solutions still contain the
original's; we solve the square system, then keep only the computed solutions that also satisfy
the original (overdetermined) system, discarding the extraneous ones the randomization introduces.

Covers: the single-affine-group path (with the descending-degree sort + [I|C] and its
homogenizing-variable power deficits exercised through the real solver), the multihomogeneous
path, no-mutation of the original, the matrix getter, the user-supplied matrix, and the
exact-coefficient discipline.
"""

import numpy as np
import pytest

import bertini as pb
from bertini import linalg


def _to_complex(solution, num_vars):
    """A solution point (in user coordinates) as a length-num_vars list of Python complex."""
    return [complex(solution[k]) for k in range(num_vars)]


def _residuals(system, point):
    """Max-abs residual of the (overdetermined) system at an affine point in user coordinates."""
    vals = system.eval(np.array(point, dtype=complex))   # double-precision eval
    return max(abs(complex(vals[i])) for i in range(system.num_functions()))


def _true_solutions(original, solutions, num_vars=2, tol=1e-7):
    """Computed solutions of the randomized system that also satisfy the original system."""
    out = []
    for s in solutions:
        p = _to_complex(s, num_vars)
        if _residuals(original, p) < tol:
            out.append(p)
    return out


def _round_pairs(points, ndig=4):
    return sorted((round(p[0].real, ndig), round(p[1].real, ndig)) for p in points)


def test_unequal_degree_single_group_solves_and_filters():
    # Overdetermined, UNEQUAL degrees (2, 1, 2): the circle, the line x=y, and 2x^2-1.
    # Common zeros are exactly (+/- 1/sqrt2, +/- 1/sqrt2) along x=y -- two real points.
    # Randomization keeps target degrees (2, 2) [the two largest], so the degree-1 function is
    # folded in with a homogenizing-variable power: this drives the h-power path through the solver.
    x, y = pb.Variable('x'), pb.Variable('y')
    original = pb.System()
    original.add_variable_group(pb.VariableGroup([x, y]))
    original.add_function(x * x + y * y - 1)
    original.add_function(x - y)
    original.add_function(2 * x * x - 1)

    pb.random.set_random_seed(1)
    randomized = linalg.randomize(original)
    assert randomized.num_functions() == 2
    assert sorted(randomized.degrees()) == [2, 2]          # product 4 = the path count
    assert original.num_functions() == 3                    # original untouched

    zd = pb.nag_algorithm.ZeroDimSolver(randomized, endgame='cauchy', mptype='adaptive', startsystem='rootsofunity')
    zd.solve()
    sols = zd.all_solutions()
    assert len(sols) == 4                                   # Bezout 2*2: two true + two extraneous

    true = _true_solutions(original, sols)
    r = 1.0 / np.sqrt(2.0)
    assert _round_pairs(true) == sorted([(round(r, 4), round(r, 4)),
                                         (round(-r, 4), round(-r, 4))])


def test_equal_degree_conics_solves_and_filters():
    # Three conics through exactly (1, 0) and (0, 1): circle, xy, and x^2+y^2-x-y.
    x, y = pb.Variable('x'), pb.Variable('y')
    original = pb.System()
    original.add_variable_group(pb.VariableGroup([x, y]))
    original.add_function(x * x + y * y - 1)
    original.add_function(x * y)
    original.add_function(x * x + y * y - x - y)

    pb.random.set_random_seed(2)
    randomized = linalg.randomize(original)
    assert randomized.num_functions() == 2

    zd = pb.nag_algorithm.ZeroDimSolver(randomized, endgame='cauchy', mptype='adaptive', startsystem='rootsofunity')
    zd.solve()
    sols = zd.all_solutions()
    assert len(sols) == 4

    true = _true_solutions(original, sols)
    assert _round_pairs(true) == sorted([(1.0, 0.0), (0.0, 1.0)])


def test_multihomogeneous_bilinear_solves_and_filters():
    # Two affine variable groups {x}, {y}; three bilinear (multidegree (1,1)) functions sharing
    # the single common solution (1, 1).  Randomize to two; the multihomogeneous start tracks
    # the bilinear Bezout (2 paths) rather than the total degree (4).
    x, y = pb.Variable('x'), pb.Variable('y')
    original = pb.System()
    original.add_variable_group(pb.VariableGroup([x]))
    original.add_variable_group(pb.VariableGroup([y]))
    original.add_function(x * y - 1)
    original.add_function(x + y - 2)
    original.add_function(x - y)

    pb.random.set_random_seed(3)
    randomized = linalg.randomize(original)
    assert randomized.num_functions() == 2

    zd = pb.nag_algorithm.ZeroDimSolver(randomized, endgame='cauchy', mptype='adaptive', startsystem='mhom')
    zd.solve()
    sols = zd.all_solutions()
    assert len(sols) >= 1

    true = _true_solutions(original, sols)
    assert _round_pairs(true) == [(1.0, 1.0)]


def test_user_supplied_matrix_combination():
    # g0 = f0 ; g1 = 3 f0 + 5 f2, with f0=x^2+y^2-1 (deg 2) and f2=2x^2-1 (deg 2) in author order.
    x, y = pb.Variable('x'), pb.Variable('y')
    original = pb.System()
    original.add_variable_group(pb.VariableGroup([x, y]))
    original.add_function(x * x + y * y - 1)
    original.add_function(x - y)
    original.add_function(2 * x * x - 1)

    randomized = linalg.randomize(original, [[1, 0, 0], [3, 0, 5]])
    assert randomized.num_functions() == 2

    # at (x, y) = (2, 3): f0 = 12, f2 = 7 -> g0 = 12, g1 = 3*12 + 5*7 = 71.
    g = randomized.eval(np.array([pb.multiprec.Complex('2'), pb.multiprec.Complex('3')]))
    g = [complex(v) for v in g]
    assert abs(g[0] - 12) < 1e-9
    assert abs(g[1] - 71) < 1e-9


def test_randomization_matrix_round_trips():
    x, y = pb.Variable('x'), pb.Variable('y')
    original = pb.System()
    original.add_variable_group(pb.VariableGroup([x, y]))
    original.add_function(x * x + y * y - 1)
    original.add_function(x - y)
    original.add_function(2 * x * x - 1)

    R = linalg.randomize(original, [[1, 0, 0], [3, 0, 5]]).randomization_matrix()
    assert R.shape == (2, 3)
    assert abs(complex(R[0, 0]) - 1) < 1e-12
    assert abs(complex(R[1, 0]) - 3) < 1e-12
    assert abs(complex(R[1, 2]) - 5) < 1e-12


def test_float_coefficient_refused():
    x, y = pb.Variable('x'), pb.Variable('y')
    original = pb.System()
    original.add_variable_group(pb.VariableGroup([x, y]))
    original.add_function(x * x + y * y - 1)
    original.add_function(x - y)
    original.add_function(2 * x * x - 1)
    with pytest.raises(TypeError):
        linalg.randomize(original, [[1.5, 0, 0], [0, 1, 0]])   # Python float -> refused
