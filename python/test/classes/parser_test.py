# This file is part of Bertini 2.
# 
# python/test/parser_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/test/parser_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/test/parser_test.py.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Copyright(C) 2016-2018 by Bertini2 Development Team
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
#  Danielle Brake
#  UWEC
#  Spring 2018
#



__author__ = 'jcollins'

from pybertini import *
import unittest
import numpy as np
import pdb

import pybertini.parse as pp
import pybertini.minieigen as mi


class ParserTest(unittest.TestCase):
    def setUp(self):
        self.tol_d = 1e-15;



    def test_create_system(self):
        tol_d = self.tol_d;
        input = 'function f, g; variable_group x,y,z; f = 3*x*y*z; g = x^2 + y^2 + z^2 - 1;';
        sys = pp.system(input);
        #
        vals = mi.VectorXd((complex(-2.43,.21 ),complex(4.84, -1.94),complex(-6.48, -.731)))
        sysEval = sys.eval(vals);
        #
        self.assertLessEqual(np.abs(sysEval[0].real / (233.2850778)-1), tol_d)
        self.assertLessEqual(np.abs(sysEval[0].imag / (-86.5039806)-1), tol_d)
        self.assertLessEqual(np.abs(sysEval[1].real / (65.978839)-1), tol_d)
        self.assertLessEqual(np.abs(sysEval[1].imag / (-10.32604)-1), tol_d)
        #
        sys.differentiate()
        sysJac = sys.jacobian(vals)


if __name__ == '__main__':
    unittest.main();

