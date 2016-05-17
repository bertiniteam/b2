# This file is part of Bertini 2.
# 
# python/test/system_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/test/system_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/test/system_test.py.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Copyright(C) 2016 by Bertini2 Development Team
# 
#  See <http://www.gnu.org/licenses/> for a copy of the license, 
#  as well as COPYING.  Bertini2 is provided with permitted 
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
# 
#   James Collins
#   West Texas A&M University
#   Spring 2016
# 


__author__ = 'James Collins'

from pybertini import *
import unittest
import numpy as np
import pdb


class AMPTrackingTest(unittest.TestCase):
    def setUp(self):
        self.toldbl = 1e-15;
        self.x = Variable("x");
        self.y = Variable("y");
        self.z = Variable("z");
        self.t = Variable("t");
        self.a = Float("4.897", "1.23")
        #
        self.f = Function(self.x*self.y);
        self.g = Function(pow(self.x,2)*self.y - self.a*self.z*self.x);



    def test_tracker_linear(self):
        x = self.x;  y = self.y; t = self.t;
        s = System();

        vars = VariableGroup();
        vars.append(y);
        s.add_function(y-t);
        s.add_path_variable(t);
        s.add_variable_group(vars);

        ampconfig = amp_config_from(s);

        tracker = AMPTracker(s);

        stepping_pref = Stepping_mp();
        newton_pref = Newton();

        tracker.setup(Predictor.Euler, mpfr_float("1e-5"), mpfr_float("1e5"), stepping_pref, newton_pref);
        tracker.precision_setup(ampconfig);

        t_start = mpfr_complex(1)
        t_end = mpfr_complex(0)

        y_start = VectorXmp([mpfr_complex(1)]);

        y_end = VectorXmp();

        tracker.track_path(y_end, t_start, t_end, y_start);

        self.assertEqual(y_end.rows(), 1)
        self.assertLessEqual(norm(y_end[0]-mpfr_complex(0)), mpfr_float("1e-5"))




    def test_tracker_quad(self):
        x = self.x;  y = self.y; t = self.t;
        s = System();

        vars = VariableGroup();
        vars.append(y);
        s.add_function(y-t**2);
        s.add_path_variable(t);
        s.add_variable_group(vars);

        ampconfig = amp_config_from(s);

        tracker = AMPTracker(s);

        stepping_pref = Stepping_mp();
        newton_pref = Newton();

        tracker.setup(Predictor.Euler, mpfr_float("1e-5"), mpfr_float("1e5"), stepping_pref, newton_pref);
        tracker.precision_setup(ampconfig);

        t_start = mpfr_complex(1)
        t_end = mpfr_complex(-1)

        y_start = VectorXmp([mpfr_complex(1)]);

        y_end = VectorXmp();

        tracker.track_path(y_end, t_start, t_end, y_start);

        self.assertEqual(y_end.rows(), 1)
        self.assertLessEqual(norm(y_end[0]-mpfr_complex(1)), mpfr_float("1e-5"))



    def test_tracker_sqrt(self):
        x = self.x;  y = self.y; t = self.t;
        s = System();

        vars = VariableGroup();
        vars.append(y); vars.append(x);
        s.add_function(x-t);
        s.add_function(y**2 - x)
        s.add_path_variable(t);
        s.add_variable_group(vars);

        ampconfig = amp_config_from(s);

        tracker = AMPTracker(s);

        stepping_pref = Stepping_mp();
        newton_pref = Newton();

        tracker.setup(Predictor.Euler, mpfr_float("1e-5"), mpfr_float("1e5"), stepping_pref, newton_pref);
        tracker.precision_setup(ampconfig);

        t_start = mpfr_complex(1)
        t_end = mpfr_complex(0)

        y_start = VectorXmp([mpfr_complex(1), mpfr_complex(1)]);

        y_end = VectorXmp();

        tracker.track_path(y_end, t_start, t_end, y_start);

        self.assertEqual(y_end.rows(), 2)
        self.assertLessEqual(norm(y_end[0]-mpfr_complex(0)), mpfr_float("1e-5"))
        self.assertLessEqual(norm(y_end[1]-mpfr_complex(0)), mpfr_float("1e-5"))

        y_start = VectorXmp([mpfr_complex(1), mpfr_complex(-1)]);

        tracker.track_path(y_end, t_start, t_end, y_start);

        self.assertEqual(y_end.rows(), 2)
        self.assertLessEqual(norm(y_end[0]-mpfr_complex(0)), mpfr_float("1e-5"))
        self.assertLessEqual(norm(y_end[1]-mpfr_complex(0)), mpfr_float("1e-5"))


        y_start = VectorXmp([mpfr_complex(-1), mpfr_complex(-1)]);

        tracker.track_path(y_end, t_start, t_end, y_start);

        self.assertEqual(y_end.rows(), 2)
        self.assertLessEqual(norm(y_end[0]-mpfr_complex(0)), mpfr_float("1e-5"))
        self.assertLessEqual(norm(y_end[1]-mpfr_complex(0)), mpfr_float("1e-5"))


        y_start = VectorXmp([mpfr_complex(-1), mpfr_complex(0,1)]);

        tracker.track_path(y_end, t_start, t_end, y_start);

        self.assertEqual(y_end.rows(), 2)
        self.assertLessEqual(norm(y_end[0]-mpfr_complex(0)), mpfr_float("1e-5"))
        self.assertLessEqual(norm(y_end[1]-mpfr_complex(0)), mpfr_float("1e-5"))


if __name__ == '__main__':
    unittest.main();


