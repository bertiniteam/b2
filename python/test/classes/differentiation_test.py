# This file is part of Bertini 2.
# 
# python/test/differentiation_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/test/differentiation_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/test/differentiation_test.py.  If not, see <http://www.gnu.org/licenses/>.
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


from pybertini import *
from pybertini.function_tree.symbol import *
from pybertini.function_tree.root import *
from pybertini.function_tree import *
import numpy as np;
import unittest
import pdb

import pybertini.multiprec as mp
from pybertini.multiprec import Float as mpfr_float
from pybertini.multiprec import Complex as mpfr_complex


class DiffTest(unittest.TestCase):
    def setUp(self):
        default_precision(30);
        self.x = Variable("x")
        self.x.set_current_value(complex(-2.43,.21 ))
        self.y = Variable("y")
        self.y.set_current_value(complex(4.84, -1.94))
        self.z = Variable('z')
        self.z.set_current_value(complex(-6.48, -.731))
        self.p = Variable('p')
        self.p.set_current_value(complex(-.321, -.72))
        self.a = Float("3.12", ".612")
        self.b = Float("-.823", "2.62")
        self.tol_d = float(1e-14);
        #
        self.x.set_current_value(mpfr_complex("-2.43",".21" ))
        self.y.set_current_value(mpfr_complex("4.84", "-1.94"))
        self.z.set_current_value(mpfr_complex("-6.48", "-.731"))
        self.p.set_current_value(mpfr_complex("-.321", "-.72"))
        self.a = Float(mpfr_complex("3.12", ".612"))
        self.b = Float(mpfr_complex("-.823", "2.62"))
        self.tol_mp = mpfr_float("1e-27");


    def test_sum_rule(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = x+y+a;
        df = f.differentiate();
        #
        self.assertLessEqual(np.abs(df.eval_d(x).real/(1.0)-1), tol_d)
        self.assertLessEqual(np.abs(df.eval_d(x).imag-0), tol_d)
        #
        self.assertLessEqual(mp.abs(df.eval_mp(x).real/mpfr_float("1.0")-1), tol_mp)
        self.assertLessEqual(mp.abs(df.eval_mp(x).imag-mpfr_float("0")), tol_mp)
        #
        df.reset()
        self.assertLessEqual(np.abs(df.eval_d(y).real/(1.0)-1), tol_d)
        self.assertLessEqual(np.abs(df.eval_d(y).imag-0), tol_d)
        #
        self.assertLessEqual(mp.abs(df.eval_mp(y).real/mpfr_float("1.0")-1), tol_mp)
        self.assertLessEqual(mp.abs(df.eval_mp(y).imag-mpfr_float("0")), tol_mp)
        #
        df.reset()
        self.assertLessEqual(np.abs(df.eval_d(z).real-0), tol_d)
        self.assertLessEqual(np.abs(df.eval_d(z).imag-0), tol_d)
        #
        self.assertLessEqual(mp.abs(df.eval_mp(z).real-mpfr_float("0.0")), tol_mp)
        self.assertLessEqual(mp.abs(df.eval_mp(z).imag-mpfr_float("0")), tol_mp)

    def test_power_rule(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = x**2 + y**3;
        df = f.differentiate();
        #
        self.assertLessEqual(np.abs(df.eval_d(x).real / (-4.86)-1), tol_d)
        self.assertLessEqual(np.abs(df.eval_d(x).imag / (0.42)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(df.eval_mp(x).real / mpfr_float("-4.86")-1), tol_mp)
        self.assertLessEqual(mp.abs(df.eval_mp(x).imag / mpfr_float("0.42")-1), tol_mp)
        #
        df.reset()
        self.assertLessEqual(np.abs(df.eval_d(y).real / (58.9860)-1), tol_d)
        self.assertLessEqual(np.abs(df.eval_d(y).imag / (-56.3376)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(df.eval_mp(y).real / mpfr_float("58.9860")-1), tol_mp)
        self.assertLessEqual(mp.abs(df.eval_mp(y).imag / mpfr_float("-56.3376")-1), tol_mp)


    def test_prod_rule(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = x**2*y**4 - a*x*y*z**2;
        df = f.differentiate();
        #
        self.assertLessEqual(np.abs(df.eval_d(x).real / (-559.28968169592)-1), tol_d)
        self.assertLessEqual(np.abs(df.eval_d(x).imag / (3577.05276993648)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(df.eval_mp(x).real / mpfr_float("-559.28968169592")-1), tol_mp)
        self.assertLessEqual(mp.abs(df.eval_mp(x).imag / mpfr_float("3577.05276993648")-1), tol_mp)
        #
        df.reset()
        self.assertLessEqual(np.abs(df.eval_d(y).real / (1161.85042980828)-1), tol_d)
        self.assertLessEqual(np.abs(df.eval_d(y).imag / (-3157.24325320476)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(df.eval_mp(y).real / mpfr_float("1161.85042980828")-1), tol_mp)
        self.assertLessEqual(mp.abs(df.eval_mp(y).imag / mpfr_float("-3157.24325320476")-1), tol_mp)
        #
        df.reset()
        self.assertLessEqual(np.abs(df.eval_d(z).real / (-520.5265859088)-1), tol_d)
        self.assertLessEqual(np.abs(df.eval_d(z).imag / (84.7479679056)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(df.eval_mp(z).real / mpfr_float("-520.5265859088")-1), tol_mp)
        self.assertLessEqual(mp.abs(df.eval_mp(z).imag / mpfr_float("84.7479679056")-1), tol_mp)


    def test_trancendental(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = sin(x*y) + exp(z*y) - log(x*x);
        df = f.differentiate();
        #
        self.assertLessEqual(np.abs(df.eval_d(x).real / (-17.648420086229721902138620795382021306411662490)-1), 9e-13)
        self.assertLessEqual(np.abs(df.eval_d(x).imag / (-803.11883403426275105632833868183320319093878729)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(df.eval_mp(x).real / mpfr_float("-17.648420086229721902138620795382021306411662490")-1), tol_mp)
        self.assertLessEqual(mp.abs(df.eval_mp(x).imag / mpfr_float("-803.11883403426275105632833868183320319093878729")-1), tol_mp)
        #
        df.reset()
        self.assertLessEqual(np.abs(df.eval_d(y).real / (-100.97157179433748763552280062599971478593963953)-1), tol_d)
        self.assertLessEqual(np.abs(df.eval_d(y).imag / (361.98093991820979266721712882115615553425318528)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(df.eval_mp(y).real / mpfr_float("-100.97157179433748763552280062599971478593963953")-1), tol_mp)
        self.assertLessEqual(mp.abs(df.eval_mp(y).imag / mpfr_float("361.98093991820979266721712882115615553425318528")-1), tol_mp)
        #
        df.reset()
        self.assertLessEqual(np.abs(df.eval_d(z).real / (-2.1642907643013779167501866500194314960002972412e-14)-1), tol_d)
        self.assertLessEqual(np.abs(df.eval_d(z).imag / (2.1105887207247540399884720817624768568595288922e-14)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(df.eval_mp(z).real / mpfr_float("-2.1642907643013779167501866500194314960002972412e-14")-1), tol_mp)
        self.assertLessEqual(mp.abs(df.eval_mp(z).imag / mpfr_float("2.1105887207247540399884720817624768568595288922e-14")-1), tol_mp)


if __name__ == '__main__':
    unittest.main();


    
