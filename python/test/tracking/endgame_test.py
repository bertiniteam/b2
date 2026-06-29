# This file is part of Bertini 2.
#
# python/test/endgame_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/endgame_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/endgame_test.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
#
#   silviana amethyst
#   University of Notre Dame
#
#  silviana amethyst
#  UWEC
#  Spring 2018
#


__author__ = 'ofloveandhate'


import pytest
import numpy as np

from bertini import *
from bertini.function_tree.symbol import *
from bertini.function_tree.root import *
from bertini.function_tree import *
from bertini.tracking import *
from bertini.endgame import *

import bertini.system.start_system as ss
import bertini.multiprec as mp
from bertini.multiprec import Float as mpfr_float
from bertini.multiprec import Complex as mpfr_complex


AMBIENT_PRECISION = 50


def test_using_total_degree_ss():
    default_precision(AMBIENT_PRECISION)

    x = Variable("x")
    y = Variable("y")
    t = Variable("t")

    sys = System()

    var_grp = VariableGroup()
    var_grp.append(x)
    var_grp.append(y)

    sys.add_variable_group(var_grp)

    sys.add_function((x-1)**3)
    sys.add_function((y-1)**2)

    sys.homogenize()
    sys.auto_patch()

    assert sys.is_patched() == 1
    assert sys.is_homogeneous() == 1

    td = ss.TotalDegree(sys)

    assert td.is_patched() == 1
    assert td.is_homogeneous() == 1

    gamma = Rational.rand()

    final_system = (1-t)*sys + gamma*t*td
    final_system.add_path_variable(t)

    prec_config = AMPConfig(final_system)

    stepping_pref = SteppingConfig()
    newton_pref = NewtonConfig()

    tracker = AMPTracker(final_system)

    tracker.setup(Predictor.RK4, 1e-5, 1e5, stepping_pref, newton_pref)
    tracker.precision_setup(prec_config)

    num_paths_to_track = td.num_start_points()
    n = int(str(num_paths_to_track))  # this line sucks, wtf.

    t_start = mpfr_complex(1)
    t_endgame_boundary = mpfr_complex("0.1")
    t_final = mpfr_complex(0)

    bdry_points = []

    for i in range(n):
        default_precision(AMBIENT_PRECISION)
        final_system.precision(AMBIENT_PRECISION)
        start_point = td.start_point_mp(i)

        bdry_pt = np.array(np.zeros((3)).astype(np.int64), dtype=mpfr_complex)

        track_success_code = tracker.track_path(bdry_pt, t_start, t_endgame_boundary, start_point)
        bdry_points.append(bdry_pt)

        assert track_success_code == SuccessCode.Success

    tracker.setup(Predictor.HeunEuler, 1e-6, 1e5, stepping_pref, newton_pref)
    my_endgame = AMPCauchyEG(tracker, t_endgame_boundary)

    final_homogenized_solutions = [np.empty(dtype=mpfr_complex, shape=(3,)) for i in range(n)]

    for i in range(n):
        final_system.precision(bdry_points[i][0].precision)

        track_success_code = my_endgame.run(bdry_points[i])

        final_homogenized_solutions[i] = my_endgame.final_approximation()

        assert track_success_code == SuccessCode.Success

    dehomogenized_solns = [sys.dehomogenize_point(soln) for soln in final_homogenized_solutions]

    exact_soln = np.array([mpfr_complex(1), mpfr_complex(1)])

    for soln in dehomogenized_solns:
        diff = exact_soln - soln
        # NB: np.sum / np.prod / np.mean over multiprecision (mpfr/mpc) arrays can raise
        # SystemError on some numpy + eigenpy builds -- numpy cannot construct the reduction
        # identity element for these custom dtypes.  See the "Known gotchas" page in the docs.
        # np.dot(diff, diff) == sum(diff_i**2) (numpy's dot does not conjugate) and goes through
        # the dtype's dot slot, which works everywhere; it preserves this assertion exactly.
        assert mp.abs(np.sqrt(np.dot(diff, diff))) < 1e-10
