# This file is part of Bertini 2.
#
# python/test/tracking/analytic_homotopy_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/tracking/analytic_homotopy_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/tracking/analytic_homotopy_test.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""The interface for tracking a homotopy whose functions are not polynomials.

There is no start system for such a target -- sine has infinitely many roots, and no degree
bounds them -- but a homotopy needs no start system, only start points.  Tracking works, and
so does the endgame.  Adaptive precision does not, because its error bounds are derived from
a degree the system does not have, and it now says so instead of failing obscurely.

Correctness lives in the C++ suite; what is checked here is the seam.
"""

import numpy as np
import pytest

import bertini
import bertini.tracking  # noqa: F401  (registers the tracking submodule)


GAMMA_RE, GAMMA_IM = '-24/25', '7/25'      # an exact point on the unit circle: 7-24-25
ZEROS_OF_SINE = ['-3.14159265358979323846264338327950288',
                 '0',
                 '3.14159265358979323846264338327950288']


def sine_homotopy(target_value='1/2'):
    """H = (1-t)*(sin(x) - value) + gamma*t*sin(x), and its target system."""
    x = bertini.Variable('x')

    target = bertini.System()
    target.add_variable_group(bertini.VariableGroup([x]))
    target.add_function(bertini.symbolics.sin(x) - bertini.coefficient(target_value))

    start = bertini.System()
    start.add_variable_group(bertini.VariableGroup([x]))
    start.add_function(bertini.symbolics.sin(x))

    gamma = (bertini.coefficient(GAMMA_RE) + bertini.I * bertini.coefficient(GAMMA_IM))
    return bertini.system.make_homotopy(target, start, gamma=gamma), target


def start_points():
    return [np.array([bertini.multiprec.complex_mp(v)]) for v in ZEROS_OF_SINE]


def test_an_analytic_homotopy_reports_no_degree():
    H, _ = sine_homotopy()
    assert not H.is_polynomial()
    assert list(H.degrees()) == [-1]


@pytest.mark.parametrize('mptype', ['double', 'multiple'])
def test_fixed_precision_tracks_an_analytic_homotopy(mptype):
    """Every path from a zero of sine lands on a root of sin(x) = 1/2, which is
    arcsin(1/2) + 2*pi*k or pi - arcsin(1/2) + 2*pi*k."""
    H, target = sine_homotopy()
    solver = bertini.HomotopySolver(H, start_points(), target, mptype=mptype)
    solver.solve()

    ends = [complex(s[0]) for s in solver.all_solutions() if len(s)]
    assert len(ends) == len(ZEROS_OF_SINE)

    exact = [np.pi / 6 + 2 * np.pi * k for k in (-1, 0, 1)] + \
            [5 * np.pi / 6 + 2 * np.pi * k for k in (-1, 0, 1)]
    for e in ends:
        assert min(abs(e - w) for w in exact) < 1e-10


def test_adaptive_precision_refuses_and_says_what_to_do_instead():
    H, target = sine_homotopy()
    with pytest.raises(ValueError) as caught:
        bertini.HomotopySolver(H, start_points(), target, mptype='adaptive')

    message = str(caught.value)
    assert 'not polynomial' in message
    assert "mptype='double'" in message                 # the fixed-precision way out
    assert 'amp_config=' in message                     # the hand-set way out


def test_deriving_an_amp_config_from_an_analytic_system_refuses():
    """The C++ derivation refuses too, for anyone who reaches past the seam."""
    _, target = sine_homotopy()
    with pytest.raises(RuntimeError):
        bertini.tracking.amp_config_from(target)


def test_adaptive_precision_works_when_the_bounds_are_set_by_hand():
    """The escape hatch the refusal names: nothing about adaptive precision is broken for an
    analytic system, only the recipe that derives its bounds from a degree.  So choose the two
    error bounds yourself and hand them over."""
    H, target = sine_homotopy()

    config = bertini.tracking.AMPConfig()
    config.jacobian_eval_error_bound = 8
    config.function_eval_error_bound = 4
    config.linear_solve_error_bound = 1

    solver = bertini.HomotopySolver(H, start_points(), target,
                                    mptype='adaptive', amp_config=config)
    solver.solve()

    ends = [complex(s[0]) for s in solver.all_solutions() if len(s)]
    assert len(ends) == len(ZEROS_OF_SINE)
    exact = [np.pi / 6 + 2 * np.pi * k for k in (-1, 0, 1)] + \
            [5 * np.pi / 6 + 2 * np.pi * k for k in (-1, 0, 1)]
    for e in ends:
        assert min(abs(e - w) for w in exact) < 1e-10


def test_the_endgame_finds_the_cycle_number_at_a_transcendental_branch_point():
    """sin(x) = 1 has a double root at pi/2.  The power series endgame's assumption survives
    for an analytic family (Weierstrass preparation), and it reports cycle number 2."""
    H, target = sine_homotopy('1')
    points = [np.array([bertini.multiprec.complex_mp(v)])
              for v in ('0', '3.14159265358979323846264338327950288')]
    solver = bertini.HomotopySolver(H, points, target, mptype='multiple')
    solver.solve()

    assert all(m.cycle_num == 2 for m in solver.solution_metadata())
    for s in solver.all_solutions():
        assert abs(complex(s[0]) - np.pi / 2) < 1e-9
