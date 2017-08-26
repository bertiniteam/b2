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
#  Copyright(C) 2016 by Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
#
#   Dani Brake
#   University of Notre Dame
#


__author__ = 'ofloveandhate'

from pybertini import *
from pybertini.function_tree.symbol import *
from pybertini.function_tree.root import *
from pybertini.function_tree import *
from pybertini.tracking import *
from pybertini.tracking.config import *
from pybertini.endgame import *
from pybertini.endgame.config import *

import unittest
import numpy as np
import pdb


class EndgameTest(unittest.TestCase):
    def setUp(self):
        self.ambient_precision = 50;

    def test_using_total_degree_ss(self):
        default_precision(self.ambient_precision);

        x = Variable("x");
        y = Variable("y");
        t = Variable("t");
        #
        sys = System();
        #
        var_grp = VariableGroup();
        var_grp.append(x);
        var_grp.append(y);

        sys.add_variable_group(var_grp);

        sys.add_function((x-1)**3)
        sys.add_function((y-1)**2)

        sys.homogenize();
        sys.auto_patch();

        self.assertEqual(sys.is_patched(), 1)
        self.assertEqual(sys.is_homogeneous(), 1)

        td = TotalDegree(sys);

        self.assertEqual(td.is_patched(), 1)
        self.assertEqual(td.is_homogeneous(), 1)

        gamma = Rational.rand();


        final_system = (1-t)*sys + gamma*t*td;
        final_system.add_path_variable(t);

        print(final_system)
        prec_config = AMPConfig(final_system);

        stepping_pref = SteppingConfig();
        newton_pref = NewtonConfig();



        tracker = AMPTracker(final_system);

        tracker.setup(Predictor.RK4, 1e-5, 1e5, stepping_pref, newton_pref);
        tracker.precision_setup(prec_config);

        num_paths_to_track = td.num_start_points();
        n = int(str(num_paths_to_track));

        t_start = mpfr_complex(1);
        t_endgame_boundary = mpfr_complex("0.1");
        t_final = mpfr_complex(0);

        bdry_points = [VectorXmp() for i in range(n)]
        for i in range(n):
            default_precision(self.ambient_precision);
            final_system.precision(self.ambient_precision);
            start_point = td.start_pointmp(i);

            bdry_pt = VectorXmp();
            track_success_code = tracker.track_path(bdry_pt,t_start, t_endgame_boundary, start_point);
            bdry_points[i] = bdry_pt;

            self.assertEqual(track_success_code, SuccessCode.Success)



        tracker.setup(Predictor.HeunEuler, 1e-6, 1e5, stepping_pref, newton_pref);
        my_endgame = AMPCauchyEG(tracker);


        final_homogenized_solutions = [VectorXmp() for i in range(n)]
        for i in range(n):
            default_precision(bdry_points[i][0].precision());
            final_system.precision(bdry_points[i][0].precision());
            track_success_code = my_endgame.run(mpfr_complex(t_endgame_boundary),bdry_points[i]);
            final_homogenized_solutions[i] = my_endgame.final_approximation();
            print(final_system.dehomogenize_point(final_homogenized_solutions[i]));
            self.assertEqual(track_success_code, SuccessCode.Success)


if __name__ == '__main__':
    unittest.main();
