# This file is part of Bertini 2.
#
# python/test/random/seeding_test.py is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/random/seeding_test.py is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this file.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""Determinism contract for the global RNG seed.

set_random_seed() controls EVERY random type, including the multiprecision draws
(complex_in_minus_one_to_one / complex_unit / real_as_complex, which feed gamma
and patch coefficients).  These used to slip through an unseedable private engine;
they now share the one per-thread engine, so the same seed reproduces a whole run.
"""

import numpy as np
import pytest

import bertini as pb
from bertini.nag_algorithm import ZeroDimCauchyAdaptivePrecisionTotalDegree


def _draw_sequence(n=6):
    """A sequence of multiprecision random draws, as full-precision strings."""
    return [repr(pb.random.complex_in_minus_one_to_one()) for _ in range(n)]


def test_get_returns_what_was_set():
    pb.random.set_random_seed(424242)
    assert pb.random.get_random_seed() == 424242


def test_same_seed_reproduces_multiprecision_draws():
    pb.random.set_random_seed(1234)
    first = _draw_sequence()

    pb.random.set_random_seed(1234)
    second = _draw_sequence()

    assert first == second


def test_different_seed_changes_multiprecision_draws():
    pb.random.set_random_seed(1234)
    first = _draw_sequence()

    pb.random.set_random_seed(5678)
    other = _draw_sequence()

    assert first != other


def test_unit_draws_also_honor_the_seed():
    pb.random.set_random_seed(99)
    first = [repr(pb.random.complex_unit()) for _ in range(4)]

    pb.random.set_random_seed(99)
    second = [repr(pb.random.complex_unit()) for _ in range(4)]

    assert first == second


# --- end-to-end: a whole solve is reproducible under a fixed seed ---
# The solver's target system carries a random patch; solutions(user_coords=False)
# exposes those patched (internal) coordinates.  Same seed -> identical patch ->
# identical internal coordinates.  Either seed -> the same user solutions, since
# the patch is just a representation.

INV_SQRT2 = 1 / np.sqrt(2)
KNOWN = (
    np.array([complex(INV_SQRT2), complex(-INV_SQRT2)]),
    np.array([complex(-INV_SQRT2), complex(INV_SQRT2)]),
)


def _solve_circle_line(seed):
    """Set the seed, then build and solve x^2+y^2-1, x+y from scratch."""
    pb.random.set_random_seed(seed)
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_function(x**2 + y**2 - 1)
    sys.add_function(x + y)
    sys.add_variable_group(pb.VariableGroup([x, y]))
    solver = ZeroDimCauchyAdaptivePrecisionTotalDegree(sys)
    solver.solve()
    return solver


def _as_complex(pt):
    return np.array([complex(pt[i]) for i in range(len(pt))])


def _sorted_internal(solver):
    """Internal (patched) coords, sorted to a canonical order for comparison."""
    internal = [_as_complex(s) for s in solver.all_solutions(user_coords=False)]
    return sorted(internal, key=lambda v: (round(v[1].real, 9), round(v[1].imag, 9)))


def test_same_seed_reproduces_patched_solution_coordinates():
    a = _sorted_internal(_solve_circle_line(31415))
    b = _sorted_internal(_solve_circle_line(31415))

    assert len(a) == len(b) == 2
    for va, vb in zip(a, b):
        # identical seed -> identical random patch -> identical internal coords
        assert np.linalg.norm(va - vb) < 1e-10


def test_user_solutions_are_seed_independent():
    """The patch differs by seed, but the actual solutions do not."""
    for seed in (31415, 27182):
        solver = _solve_circle_line(seed)
        sols = solver.all_solutions()
        assert len(sols) == 2
        for s in sols:
            d = min(np.linalg.norm(_as_complex(s) - k) for k in KNOWN)
            assert d < 1e-8
