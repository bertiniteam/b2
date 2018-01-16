# This file is part of Bertini 2.
# 
# python/test/function_tree_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/test/function_tree_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/test/function_tree_test.py.  If not, see <http://www.gnu.org/licenses/>.
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

import pybertini.multiprec as mp
from pybertini.multiprec import Float as mpfr_float
from pybertini.multiprec import Complex as mpfr_complex

class SymbolTest(unittest.TestCase):
    def setUp(self):
        default_precision(30);
        self.x_d = complex(-2.43,.21 );
        self.y_d = complex(4.84, -1.94);
        self.z_d = complex(-6.48, -.731);
        self.p_d = complex(-.321, -.72);
        self.tol_d = float(1e-15);
        #
        self.x_mp = mpfr_complex("-2.43",".21" );
        self.y_mp = mpfr_complex("4.84", "-1.94");
        self.z_mp = mpfr_complex("-6.48", "-.731");
        self.p_mp = mpfr_complex("-.321", "-.72");
        self.tol_mp = mpfr_float("1e-27");


    def test_Float_construct(self):
        x_d = self.x_d; y_d = self.y_d; z_d = self.z_d; p_d = self.p_d; tol_d = self.tol_d;
        x_mp = self.x_mp; y_mp = self.y_mp; z_mp = self.z_mp; p_mp = self.p_mp; tol_mp = self.tol_mp;
        x = Float("4.3", "-9e-3");
        x = Float(y_mp);
        x = Float(mpfr_float("9.3"), mpfr_float("-3"));
        x = Float("9.2", "-43.2e2");


    def test_Float_eval(self):
        x_d = self.x_d; y_d = self.y_d; z_d = self.z_d; p_d = self.p_d; tol_d = self.tol_d;
        x_mp = self.x_mp; y_mp = self.y_mp; z_mp = self.z_mp; p_mp = self.p_mp; tol_mp = self.tol_mp;
        x = Float(x_mp); y = Float(y_mp); z = Float(z_mp);
        #
        self.assertLessEqual(np.abs(x.eval_d().real/(-2.43)-1), tol_d);
        self.assertLessEqual(np.abs(x.eval_d().imag/(.21)-1), tol_d);
        #
        self.assertLessEqual(mp.abs(y.eval_mp().real/mpfr_float("4.84")-1), tol_mp);
        self.assertLessEqual(mp.abs(y.eval_mp().imag/mpfr_float("-1.94")-1), tol_mp);

    def test_Float_funcs(self):
        x_d = self.x_d; y_d = self.y_d; z_d = self.z_d; p_d = self.p_d; tol_d = self.tol_d;
        x_mp = self.x_mp; y_mp = self.y_mp; z_mp = self.z_mp; p_mp = self.p_mp; tol_mp = self.tol_mp;
        x = Float(x_mp); y = Float(y_mp); z = Float(z_mp);
        #
        self.assertEqual(x.degree(), 0)
        d = y.differentiate();
        self.assertLessEqual(mp.abs(d.eval_mp().real-mpfr_float("0")), tol_mp)
        self.assertLessEqual(mp.abs(d.eval_mp().imag-mpfr_float("0")), tol_mp)
        self.assertTrue(y.is_homogeneous());
        self.assertTrue(y.is_polynomial());





    def test_Variable_construct(self):
        x_d = self.x_d; y_d = self.y_d; z_d = self.z_d; p_d = self.p_d; tol_d = self.tol_d;
        x_mp = self.x_mp; y_mp = self.y_mp; z_mp = self.z_mp; p_mp = self.p_mp; tol_mp = self.tol_mp;
        x = Variable("x");


    def test_Variable_eval(self):
        x_d = self.x_d; y_d = self.y_d; z_d = self.z_d; p_d = self.p_d; tol_d = self.tol_d;
        x_mp = self.x_mp; y_mp = self.y_mp; z_mp = self.z_mp; p_mp = self.p_mp; tol_mp = self.tol_mp;
        x = Variable("x"); y = Variable("y");
        x.set_current_value(x_d); x.set_current_value(x_mp)
        #
        self.assertLessEqual(np.abs(x.eval_d().real/(-2.43)-1), tol_d)
        self.assertLessEqual(np.abs(x.eval_d().imag/(.21)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(x.eval_mp().real/mpfr_float("-2.43")-1), tol_mp)
        self.assertLessEqual(mp.abs(x.eval_mp().imag/mpfr_float(".21")-1), tol_mp)

    def test_Variable_funcs(self):
        x_d = self.x_d; y_d = self.y_d; z_d = self.z_d; p_d = self.p_d; tol_d = self.tol_d;
        x_mp = self.x_mp; y_mp = self.y_mp; z_mp = self.z_mp; p_mp = self.p_mp; tol_mp = self.tol_mp;
        x = Variable("x"); y = Variable("y")
        #
        self.assertEqual(x.degree(), 1)
        self.assertEqual(x.degree(x), 1)
        self.assertEqual(x.degree(y), 0)
        d = y.differentiate()
        self.assertTrue(y.is_homogeneous())
        self.assertTrue(y.is_polynomial())



    def test_Pi_construct(self):
        x_d = self.x_d; y_d = self.y_d; z_d = self.z_d; p_d = self.p_d; tol_d = self.tol_d;
        x_mp = self.x_mp; y_mp = self.y_mp; z_mp = self.z_mp; p_mp = self.p_mp; tol_mp = self.tol_mp;
        x = Pi()
        y = make_pi()
        #
        self.assertLessEqual(np.abs(x.eval_d().real/(3.1415926535897932384626433832795028841971693994)-1), tol_d)
        self.assertLessEqual(np.abs(x.eval_d().imag - (0)), tol_d)
        #
        self.assertLessEqual(mp.abs(x.eval_mp().real/mpfr_float("3.1415926535897932384626433832795028841971693994")-1), tol_mp)
        self.assertLessEqual(mp.abs(x.eval_mp().imag - mpfr_float("0")), tol_mp)
        #
        self.assertLessEqual(np.abs(y.eval_d().real/(3.1415926535897932384626433832795028841971693994)-1), tol_d)
        self.assertLessEqual(np.abs(y.eval_d().imag - (0)), tol_d)
        #
        self.assertLessEqual(mp.abs(y.eval_mp().real/mpfr_float("3.1415926535897932384626433832795028841971693994")-1), tol_mp)
        self.assertLessEqual(mp.abs(y.eval_mp().imag - mpfr_float("0")), tol_mp)


    def test_E_construct(self):
        x_d = self.x_d; y_d = self.y_d; z_d = self.z_d; p_d = self.p_d; tol_d = self.tol_d;
        x_mp = self.x_mp; y_mp = self.y_mp; z_mp = self.z_mp; p_mp = self.p_mp; tol_mp = self.tol_mp;
        x = E()
        y = make_e()
        #
        self.assertLessEqual(np.abs(x.eval_d().real/(2.7182818284590452353602874713526624977572470937)-1), tol_d)
        self.assertLessEqual(np.abs(x.eval_d().imag - (0)), tol_d)
        #
        self.assertLessEqual(mp.abs(x.eval_mp().real - mpfr_float("2.7182818284590452353602874713526624977572470937")), tol_mp)
        self.assertLessEqual(mp.abs(x.eval_mp().imag - mpfr_float("0")), tol_mp)
        #
        self.assertLessEqual(np.abs(y.eval_d().real/(2.7182818284590452353602874713526624977572470937)-1), tol_d)
        self.assertLessEqual(np.abs(y.eval_d().imag - (0)), tol_d)
        #
        self.assertLessEqual(mp.abs(y.eval_mp().real/mpfr_float("2.7182818284590452353602874713526624977572470937")-1), tol_mp)
        self.assertLessEqual(mp.abs(y.eval_mp().imag - mpfr_float("0")), tol_mp)


    def test_I_construct(self):
        x_d = self.x_d; y_d = self.y_d; z_d = self.z_d; p_d = self.p_d; tol_d = self.tol_d;
        x_mp = self.x_mp; y_mp = self.y_mp; z_mp = self.z_mp; p_mp = self.p_mp; tol_mp = self.tol_mp;
        y = make_i()
        #
        self.assertLessEqual(np.abs(y.eval_d().real - (0)), tol_d)
        self.assertLessEqual(np.abs(y.eval_d().imag/(1.0)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(y.eval_mp().real - mpfr_float("0")), tol_mp)
        self.assertLessEqual(mp.abs(y.eval_mp().imag/mpfr_float("1.0")-1), tol_mp)





class OperatorTest(unittest.TestCase):
    def setUp(self):
        default_precision(30);
        self.x = Variable("x")
        self.x.set_current_value(complex(-2.43,.21 ))
        self.y = Variable("y")
        self.y.set_current_value(complex(4.84, -1.94))
        self.z = Variable("z")
        self.z.set_current_value(complex(-6.48, -.731))
        self.p = Variable("p")
        self.p.set_current_value(complex(-.321, -.72))
        self.a = Float(mpfr_complex("3.12", ".612"))
        self.b = Float(mpfr_complex("-.823", "2.62"))
        self.tol_d = float(9e-14);
        #
        self.x.set_current_value(mpfr_complex("-2.43",".21" ))
        self.y.set_current_value(mpfr_complex("4.84", "-1.94"))
        self.z.set_current_value(mpfr_complex("-6.48", "-.731"))
        self.p.set_current_value(mpfr_complex("-.321", "-.72"))
        self.a = Float(mpfr_complex("3.12", ".612"))
        self.b = Float(mpfr_complex("-.823", "2.62"))
        self.tol_mp = mpfr_float("1e-27");

    def test_plus(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = Function(x+y+a)
        self.assertLessEqual(np.abs(f.eval_d().real/(5.53)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag/(-1.118)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real/mpfr_float("5.53")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag/mpfr_float("-1.118")-1), tol_mp)
        #
        f = Function(x+Float("3.87"))
        self.assertLessEqual(np.abs(f.eval_d().real/(1.44)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag/(0.21)-1), tol_d)
        #
        f = Function(x+mpfr_complex("3.87", "-2.1"))
        self.assertLessEqual(np.abs(f.eval_d().real/(1.44)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag/(-1.89)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real/mpfr_float("1.44")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag/mpfr_float("-1.89")-1), tol_mp)
        #
        f = Function(x+(-5))
        self.assertLessEqual(np.abs(f.eval_d().real/(-7.43)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag/(0.21)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real/mpfr_float("-7.43")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag/mpfr_float("0.21")-1), tol_mp)
        #
        f = Function(x); f += y; f += a;
        self.assertLessEqual(np.abs(f.eval_d().real/(5.53)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag/(-1.118)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real/mpfr_float("5.53")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag/mpfr_float("-1.118")-1), tol_mp)
        #
        f = Function(x); f += Float("3.87");
        self.assertLessEqual(np.abs(f.eval_d().real/(1.44)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag/(0.21)-1), tol_d)



    def test_sub(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = Function(x-y-a)
        self.assertLessEqual(np.abs(f.eval_d().real/(-10.39)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag/(1.538)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real/mpfr_float("-10.39")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag/mpfr_float("1.538")-1), tol_mp)
        #
        f = Function(y-Float("3.87"))
        self.assertLessEqual(np.abs(f.eval_d().real/(0.97)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag/(-1.94)-1), tol_d)
        #
        f = Function(y-Float("3.87","0"))
        self.assertLessEqual(np.abs(f.eval_d().real / (0.97)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-1.94)-1), tol_d)
        #
        #
        f = Function(y-mpfr_complex("3.87", "-2.1"))
        self.assertLessEqual(np.abs(f.eval_d().real / (0.97)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (0.16)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("0.97")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("0.16")-1), tol_mp)
        #
        f = Function(y-(-5))
        self.assertLessEqual(np.abs(f.eval_d().real / (9.84)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-1.94)-1), tol_d)
        #
        f = Function(x); f -= y; f -= a;
        self.assertLessEqual(np.abs(f.eval_d().real / (-10.39)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (1.538)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("-10.39")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("1.538")-1), tol_mp)
        #
        f = Function(y); f -= Float("3.87");
        self.assertLessEqual(np.abs(f.eval_d().real / (0.97)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-1.94)-1), tol_d)


    def test_num_times_var(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = Function(a*x*b*y)
        #
        self.assertLessEqual(np.abs(f.eval_d().real / (3.4011196056)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-110.9953448712)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("3.4011196056")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("-110.9953448712")-1), tol_mp)

    def test_var_times_var(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = Function(x*y*z)
        self.assertLessEqual(np.abs(f.eval_d().real / (77.7616926)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-28.8346602)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("77.7616926")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("-28.8346602")-1), tol_mp)

    def test_var_div_var(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = Function(x/y)
        self.assertLessEqual(np.abs(f.eval_d().real / (-.44755270475041560619657805305047592426404601827)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-.13600253041648890000441351713180233327939034617)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("-.44755270475041560619657805305047592426404601827")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("-.13600253041648890000441351713180233327939034617")-1), tol_mp)

    def test_trans_funcs(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = Function(sin(x))
        self.assertLessEqual(np.abs(f.eval_d().real/(-.66749329633668695550441899166308616328986315948)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag/(-.16020928942503633132090203927960650380076680938)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("-.66749329633668695550441899166308616328986315948")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("-.16020928942503633132090203927960650380076680938")-1), tol_mp)
        #
        f = Function(cos(y))
        self.assertLessEqual(np.abs(f.eval_d().real / (0.45194679593300564730917329070452033759984813611)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-3.3798161097977088705360399142708324234265626016)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("0.45194679593300564730917329070452033759984813611")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("-3.3798161097977088705360399142708324234265626016")-1), tol_mp)
        #
        f = Function(tan(z))
        self.assertLessEqual(np.abs(f.eval_d().real / (-.11998086808607765336591715593714295443402911227)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-.63859741450762243500349270264429166927927928889)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("-.11998086808607765336591715593714295443402911227")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("-.63859741450762243500349270264429166927927928889")-1), tol_mp)
        #
        f = Function(asin(x))
        self.assertLessEqual(np.abs(f.eval_d().real / (-1.4763431474004472804452143435221887167393328861)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (1.5406263884278099750127157814537559611048741005)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("-1.4763431474004472804452143435221887167393328861")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("1.5406263884278099750127157814537559611048741005")-1), tol_mp)
        #
        f = Function(acos(y))
        self.assertLessEqual(np.abs(f.eval_d().real / (0.38769800860408664087229892623614567735197135529)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (2.3379037587834289977359318611458042347281923566)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("0.38769800860408664087229892623614567735197135529")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("2.3379037587834289977359318611458042347281923566")-1), tol_mp)
        #
        f = Function(atan(z))
        self.assertLessEqual(np.abs(f.eval_d().real / (-1.4195347801361539102032503530226060969949192059)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-0.16801358511827150554928904776095870747673962940e-1)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("-1.4195347801361539102032503530226060969949192059")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("-0.16801358511827150554928904776095870747673962940e-1")-1), tol_mp)
        #
        f = Function(exp(y))
        self.assertLessEqual(np.abs(f.eval_d().real / (-45.639359208255772966298371983389382308765171859)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-117.94721623715960520658000231550940351946595854)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("-45.639359208255772966298371983389382308765171859")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("-117.94721623715960520658000231550940351946595854")-1), tol_mp)
        #
        f = Function(log(z))
        self.assertLessEqual(np.abs(f.eval_d().real / (1.8750432590213669716781046977781508038070552297)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-3.0292589170775161726973168096174940982043177322)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("1.8750432590213669716781046977781508038070552297")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("-3.0292589170775161726973168096174940982043177322")-1), tol_mp)

    def test_power(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = Function(y**3);
        self.assertLessEqual(np.abs(f.eval_d().real / (58.732432)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-129.035608)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("58.732432")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("-129.035608")-1), tol_mp)
        #
        f = Function(pow(y,3));
        self.assertLessEqual(np.abs(f.eval_d().real / (58.732432)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-129.035608)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("58.732432")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("-129.035608")-1), tol_mp)
        #
        f = Function(x**p);
        self.assertLessEqual(np.abs(f.eval_d().real / (-.35190932545709434788093164550270669097948909024)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-6.7687858345625791466707575744042177964518271087)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("-.35190932545709434788093164550270669097948909024")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("-6.7687858345625791466707575744042177964518271087")-1), tol_mp)
        #
        f = Function(pow(x,p));
        self.assertLessEqual(np.abs(f.eval_d().real / (-.35190932545709434788093164550270669097948909024)-1), tol_d)
        self.assertLessEqual(np.abs(f.eval_d().imag / (-6.7687858345625791466707575744042177964518271087)-1), tol_d)
        #
        self.assertLessEqual(mp.abs(f.eval_mp().real / mpfr_float("-.35190932545709434788093164550270669097948909024")-1), tol_mp)
        self.assertLessEqual(mp.abs(f.eval_mp().imag / mpfr_float("-6.7687858345625791466707575744042177964518271087")-1), tol_mp)


    def test_Operator_degree(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = Function(x**3-9*y*y*z*pow(x,4));
        w = VariableGroup(); w.append(x); w.append(y);
        self.assertEqual(f.degree(),7)
        self.assertEqual(f.degree(x),4)
        self.assertEqual(f.degree(y),2)
        self.assertEqual(f.degree(z),1)
        self.assertEqual(f.degree(w),6)
        w = VariableGroup(); w.append(y); w.append(z);
        self.assertEqual(f.degree(w),3)

    def test_Operator_ishom(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = Function(y**2 - 9*y*z + 3*x*x - y*82*x - pow(z,2));
        self.assertTrue(f.is_homogeneous())
        f = Function(y**2 - 9*y*z + 3*x*x - y*82*x - pow(z,4));
        self.assertFalse(f.is_homogeneous())
        f = Function(y**2 - 9*y*z + 3*x*x - 5);
        self.assertFalse(f.is_homogeneous())

    def test_Operator_ispoly(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        f = Function(y**2 - 9*y*z + 3*x*x - y*82*x - pow(z,2));
        self.assertTrue(f.is_polynomial());
        #
        f = Function(x**2*y - 9 + sin(z));
        self.assertFalse(f.is_polynomial())

    def test_Homogenize(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;
        #
        h = Variable("h");
        f = Function(x**2 + y**2 + z**2 - 1);
        #
        vars = VariableGroup();
        vars.append(y); vars.append(x); vars.append(z);
        #
        f.homogenize(vars,h);
        self.assertEqual(f.degree(h),2)
        self.assertTrue(f.is_homogeneous())
        #
        self.assertFalse(f.is_homogeneous(x))
        self.assertFalse(f.is_homogeneous(y))
        self.assertFalse(f.is_homogeneous(z))
        self.assertFalse(f.is_homogeneous(h))
        #
        vars.append(h);
        self.assertTrue(f.is_homogeneous(vars))





if __name__ == '__main__':
    unittest.main();


