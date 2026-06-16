"""A user-authored product-of-linears start system, solved end to end.

Unlike the generated (total-degree / multihomogeneous) start systems, here the *user* authors
the start system as an explicit product of linear forms with exact coefficients.  Each factor
c.[x;1] = 0 is a hyperplane, so the start solutions are exact intersections of one hyperplane
per function -- writable by hand.  We blend that start system into a homotopy
(nag_algorithm.blend_homotopy) and track its start points to the target's roots
(nag_algorithm.user_homotopy), exercising the first-class C++ ProductsOfLinearsBlock through the
whole zero-dim pipeline.
"""

import itertools
import math

import numpy as np

import bertini as pb
from bertini import linalg
from bertini import multiprec as mp
from bertini import nag_algorithm


def _target():
    # T = { x^2 + y^2 - 1, y - x^2 }: a unit circle meeting a parabola, Bezout number 4.
    x, y = pb.Variable('x'), pb.Variable('y')
    T = pb.System()
    T.add_variable_group(pb.VariableGroup([x, y]))
    T.add_function(x * x + y * y - 1)
    T.add_function(y - x * x)
    return T


def _start():
    # S = { (x - 1)(x + 1), (y - 1)(y - 2) }: one product of two linear forms per function,
    # written as exact augmented coefficient rows ([coeff_x, coeff_y, constant]).
    x, y = pb.Variable('x'), pb.Variable('y')
    S = pb.System()
    S.add_variable_group(pb.VariableGroup([x, y]))
    linalg.add_products_of_linears(S, [
        [[1, 0, '-1'], [1, 0, '1']],   # (x - 1)(x + 1)
        [[0, 1, '-1'], [0, 1, '-2']],  # (y - 1)(y - 2)
    ])
    return S


def _start_points():
    # one hyperplane per slot: x in {1, -1}, y in {1, 2} -- the four intersections, by hand.
    return [np.array([mp.Complex(str(a)), mp.Complex(str(b))])
            for a, b in itertools.product([1, -1], [1, 2])]


def _known_solutions():
    # Substituting y = x^2 into x^2 + y^2 - 1 gives y^2 + y - 1 = 0, so y = (-1 +/- sqrt 5)/2.
    # y1 > 0 -> a real x pair (x = +/- sqrt y1); y2 < 0 -> a purely imaginary x pair.
    s5 = math.sqrt(5)
    y1, y2 = (-1 + s5) / 2, (-1 - s5) / 2
    xr, xc = math.sqrt(y1), math.sqrt(-y2)
    return [np.array([xr, y1]), np.array([-xr, y1]),
            np.array([1j * xc, y2]), np.array([-1j * xc, y2])]


# exact gamma off the real axis (|gamma| = 1) -> a reproducible straight-line path that misses
# the measure-zero singular locus.
_GAMMA = mp.Complex('0.6', '0.8')


def test_user_authored_product_of_linears_solves_to_known_roots():
    T, S = _target(), _start()
    H = nag_algorithm.blend_homotopy(T, S, gamma=linalg.coefficient(_GAMMA))
    solver = nag_algorithm.user_homotopy(H, _start_points(), T)
    solver.solve()
    sols = solver.solutions()

    assert len(sols) == 4                                   # guard against an empty/short result

    # distance to the known solutions with an infinity-norm (scale-faithful), not a scale-naive
    # residual: every known root must be matched by some computed solution.
    computed = [np.array([complex(v) for v in p]) for p in sols]
    for known in _known_solutions():
        nearest = min(np.max(np.abs(c - known)) for c in computed)
        assert nearest < 1e-8


def test_metadata_splits_real_and_complex():
    T, S = _target(), _start()
    H = nag_algorithm.blend_homotopy(T, S, gamma=linalg.coefficient(_GAMMA))
    solver = nag_algorithm.user_homotopy(H, _start_points(), T)
    solver.solve()

    md = solver.solution_metadata()
    assert len(md) == 4
    assert all(m.is_finite for m in md)
    assert not any(m.is_singular for m in md)
    assert sum(1 for m in md if m.is_real) == 2             # the (+/- sqrt y1, y1) real pair
    assert sum(1 for m in md if not m.is_real) == 2         # the purely-imaginary-x pair


def test_coefficient_parameter_homotopy_does_not_drop_a_structured_start():
    # The footgun (ADR-0020): System node arithmetic only combines the polynomial block, so a
    # products-of-linears start would be silently dropped.  coefficient_parameter_homotopy must
    # instead blend, so that H at t=1 IS the start system and vanishes at the start points.  We
    # check this by direct evaluation rather than by tracking: the no-gamma-trick real path is
    # conditioning-fragile for this hand-picked example (which is exactly why blend_homotopy's
    # off-axis gamma exists), so a track here would be flaky -- but the homotopy is still correct.
    T, S = _target(), _start()
    H = nag_algorithm.coefficient_parameter_homotopy(T, S)
    assert H.have_path_variable()
    assert H.num_functions() == 2

    # H|t=1 = 1 * S, so it must vanish at every start point (it would not if S were dropped).
    H.set_path_variable(complex(1))
    for a, b in itertools.product([1, -1], [1, 2]):
        H.set_variables(np.array([complex(a), complex(b)]))
        assert max(abs(complex(v)) for v in H.eval()) < 1e-10

    # H|t=0 = T, so it must vanish at the target's known roots.
    H.set_path_variable(complex(0))
    for root in _known_solutions():
        H.set_variables(np.array([complex(root[0]), complex(root[1])]))
        assert max(abs(complex(v)) for v in H.eval()) < 1e-9
