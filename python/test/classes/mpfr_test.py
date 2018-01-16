# This file is part of Bertini 2.
# 
# python/test/mpfr_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/test/mpfr_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/test/mpfr_test.py.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Copyright(C) 2016-2018 by Bertini2 Development Team
# 
#  See <http://www.gnu.org/licenses/> for a copy of the license, 
#  as well as COPYING.  Bertini2 is provided with permitted 
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
# 
#   Danielle Brake
#   University of Wisconsin - Eau Claire
#   Fall 2017, Spring 2018
#  
#   James Collins
#   West Texas A&M University
#   Spring 2016
# 

import pybertini as pb

from pybertini import multiprec as mp

import unittest
import pdb


dbltol = 1e-15;
class MPFRFloat(unittest.TestCase):
    def setUp(self):
        mp.default_precision(30);
        self.x = mp.Float("4.23")
        self.y = mp.Float("-3.86")
        self.z = mp.Float("1.1495")
        self.p = mp.Float(".34")
        self.tol = mp.Float("1e-27");

    def test_arith_int(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        self.assertLessEqual(mp.abs((x+8) - mp.Float("12.23")), tol)
        self.assertLessEqual(mp.abs((y-2) - mp.Float("-5.86")), tol)
        self.assertLessEqual(mp.abs((8+x) - mp.Float("12.23")), tol)
        self.assertLessEqual(mp.abs((2-y) - mp.Float("5.86")), tol)
        self.assertLessEqual(mp.abs((z*6) - mp.Float("6.897")), tol)
        self.assertLessEqual(mp.abs((6*z) - mp.Float("6.897")), tol)
        self.assertLessEqual(mp.abs((y/3) - mp.Float("-1.2866666666666666666666666666666667")), tol)
        self.assertLessEqual(mp.abs((3/y) - mp.Float("-.77720207253886010362694300518134714")), tol)
        self.assertLessEqual(mp.abs((x**3) - mp.Float("75.686967")), tol)
        #
        result = mp.Float(x);
        result += 8;
        self.assertLessEqual(mp.abs(result - mp.Float("12.23")), tol)
        result = mp.Float(y);
        result -= 2;
        self.assertLessEqual(mp.abs(result - mp.Float("-5.86")), tol)
        result = mp.Float(z);
        result *= 6;
        self.assertLessEqual(mp.abs(result - mp.Float("6.897")), tol)
        result = mp.Float(y);
        result /= 3;
        self.assertLessEqual(mp.abs(result - mp.Float("-1.2866666666666666666666666666666667")), tol)
        #
        self.assertLessEqual(mp.abs((-z) - mp.Float("-1.1495")), tol)

    def test_arith_mpfr(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        self.assertLessEqual(mp.abs((x+y) - mp.Float("0.37")), tol)
        self.assertLessEqual(mp.abs((z-y) - mp.Float("5.0095")), tol)
        self.assertLessEqual(mp.abs((z*y) - mp.Float("-4.437070")), tol)
        self.assertLessEqual(mp.abs((y/x) - mp.Float("-.91252955082742316784869976359338061")), tol)
        self.assertLessEqual(mp.abs((x**y) - mp.Float("0.0038223124228935822000384505727705508")), tol)
        self.assertLessEqual(mp.abs((-z) - mp.Float("-1.1495")), tol)
        #
        result = mp.Float(x);
        result += y;
        self.assertLessEqual(mp.abs(result - mp.Float("0.37")), tol)
        result = mp.Float(z);
        result -= y;
        self.assertLessEqual(mp.abs(result - mp.Float("5.0095")), tol)
        result = mp.Float(z);
        result *= y;
        self.assertLessEqual(mp.abs(result - mp.Float("-4.437070")), tol)
        result = mp.Float(y);
        result /= x;
        self.assertLessEqual(mp.abs(result - mp.Float("-.91252955082742316784869976359338061")), tol)


    def test_trancendentals(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        self.assertLessEqual(mp.abs((mp.exp(x)) - mp.Float("68.717232173846461408252914213396109")), tol)
        self.assertLessEqual(mp.abs((mp.log(z)) - mp.Float("0.13932706522109918666170810230684295")), tol)
        self.assertLessEqual(mp.abs((mp.sqrt(z)) - mp.Float("1.0721473779289860297522254519889560")), tol)
        #
        self.assertLessEqual(mp.abs((mp.sin(x)) - mp.Float("-.88588921129660245121088859729926237")), tol)
        self.assertLessEqual(mp.abs((mp.cos(y)) - mp.Float("-.75285494656729525719980460936483635")), tol)
        self.assertLessEqual(mp.abs((mp.tan(z)) - mp.Float("2.2315038042849919118711153687209483")), tol)
        #
        self.assertLessEqual(mp.abs((mp.asin(p)) - mp.Float("0.34691689752716170922069696210451452")), tol)
        self.assertLessEqual(mp.abs((mp.acos(p)) - mp.Float("1.2238794292677349100106247295352369")), tol)
        self.assertLessEqual(mp.abs((mp.atan(z)) - mp.Float("0.85483739856328448882289109284144652")), tol)
        #
        self.assertLessEqual(mp.abs((mp.sinh(x)) - mp.Float("34.351339891649022639414777866662100")), tol)
        self.assertLessEqual(mp.abs((mp.cosh(y)) - mp.Float("23.743209684188284295743755381842167")), tol)
        self.assertLessEqual(mp.abs((mp.tanh(z)) - mp.Float("0.81758837109637920976170104688035086")), tol)

    def test_change_prec(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        mp.default_precision(40);
        tol = mp.Float("1e-37");
        t = mp.Float("4.23")
        self.assertLessEqual(mp.abs(t**(-2) - mp.Float("0.055888089689206333238323580861682566828183245868")), tol)
        mp.default_precision(30);
        tol = mp.Float("1e-27");










class MPFRComplex(unittest.TestCase):
    def setUp(self):
        mp.default_precision(30);
        self.x = mp.Complex("-2.43",".21" )
        self.y = mp.Complex("4.84", "-1.94")
        self.z = mp.Complex("-6.48", "-.731")
        self.p = mp.Complex("-.321", "-.72")
        self.tol = mp.Float("1e-27");

    def test_construct(self):
        t = mp.Complex(3.452)
        t = mp.Complex(mp.Float("-5.6"))
        t = mp.Complex("3.89")
        t = mp.Complex(mp.Float("2.98"), mp.Float("-1e-4"))
        t = mp.Complex(3.4, 3.5)
        t = mp.Complex("6e2", mp.Float("4.32"))
        t = mp.Complex(mp.Float("4.32"), "6e2")

    def test_arith_mp_float(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        a = mp.Float("3.12"); b = mp.Float("-5.92")
        res = mp.Complex(x+a)
        self.assertLessEqual(mp.abs(res.real - mp.Float("0.69")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("0.21")), tol)
        res = mp.Complex(y-b)
        self.assertLessEqual(mp.abs(res.real - mp.Float("10.76")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-1.94")), tol)
        res = mp.Complex(a+x)
        self.assertLessEqual(mp.abs(res.real - mp.Float("0.69")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("0.21")), tol)
        res = mp.Complex(b-y)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-10.76")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("1.94")), tol)
        res = mp.Complex(z*a)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-20.2176")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-2.28072")), tol)
        res = mp.Complex(a*z)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-20.2176")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-2.28072")), tol)
        res = mp.Complex(y/b)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-.81756756756756756756756756756756756756756756757")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float(".3277027027027027027027027027027027027027027027")), tol)
        res = mp.Complex(b/y)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-1.0538301972842157915642975887484736586585850264")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-.42240301296102864372618539714298324334662292381")), tol)
        res = mp.Complex(x**a)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-16.054376621961088182387920766649714821973863952")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-1.7411284591111236754359685247799914985638458821")), tol)
        #
        #
        #
        res = mp.Complex(x);
        res += a;
        self.assertLessEqual(mp.abs(res.real - mp.Float("0.69")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("0.21")), tol)
        res = mp.Complex(y);
        res -= b;
        self.assertLessEqual(mp.abs(res.real - mp.Float("10.76")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-1.94")), tol)
        res = mp.Complex(z);
        res *= a;
        self.assertLessEqual(mp.abs(res.real - mp.Float("-20.2176")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-2.28072")), tol)
        res = mp.Complex(y);
        res /= b;
        self.assertLessEqual(mp.abs(res.real - mp.Float("-.81756756756756756756756756756756756756756756757")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float(".3277027027027027027027027027027027027027027027")), tol)
        #
        res = mp.Complex(x**4);
        self.assertLessEqual(mp.abs(res.real - mp.Float("33.30735228")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-11.96306496")), tol)
        #
        res = mp.Complex(-z);
        self.assertLessEqual(mp.abs(res.real - mp.Float("6.48")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float(".731")), tol)



    def test_arith_mp_complex(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        #
        res = mp.Complex(x+y)
        self.assertLessEqual(mp.abs(res.real - mp.Float("2.41")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-1.73")), tol)
        res = mp.Complex(y-z)
        self.assertLessEqual(mp.abs(res.real - mp.Float("11.32")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-1.209")), tol)
        res = mp.Complex(y+x)
        self.assertLessEqual(mp.abs(res.real - mp.Float("2.41")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-1.73")), tol)
        res = mp.Complex(z-y)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-11.32")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("1.209")), tol)
        res = mp.Complex(z*x)
        self.assertLessEqual(mp.abs(res.real - mp.Float("15.89991")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("0.41553")), tol)
        res = mp.Complex(x*z)
        self.assertLessEqual(mp.abs(res.real - mp.Float("15.89991")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float(".41553")), tol)
        res = mp.Complex(y/x)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-2.0454866364094805849722642460917801311144730207")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("0.62158345940494200706001008572869389813414019163")), tol)
        res = mp.Complex(x/y)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-.44755270475041560619657805305047592426404601827")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-.13600253041648890000441351713180233327939034617")), tol)
        res = mp.Complex(y**z)
        self.assertLessEqual(mp.abs(res.real - mp.Float("0.0000051612634484879218649489640888954160904291899461")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("0.16242051733741136410199105656393042124100116889e-4")), tol)
        #
        #
        #
        res = mp.Complex(x);
        res += y;
        self.assertLessEqual(mp.abs(res.real - mp.Float("2.41")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-1.73")), tol)
        res = mp.Complex(y);
        res -= z;
        self.assertLessEqual(mp.abs(res.real - mp.Float("11.32")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-1.209")), tol)
        res = mp.Complex(z);
        res *= x;
        self.assertLessEqual(mp.abs(res.real - mp.Float("15.89991")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float(".41553")), tol)
        res = mp.Complex(y);
        res /= x;
        self.assertLessEqual(mp.abs(res.real - mp.Float("-2.0454866364094805849722642460917801311144730207")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("0.62158345940494200706001008572869389813414019163")), tol)


    def test_trancendentals(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        #
        res = mp.exp(x)
        self.assertLessEqual(mp.abs(res.real - mp.Float("0.086102743899954532232498058731947255067424332219")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("0.018352149302889219131202317785160051395400089327")), tol)
        res = mp.log(y)
        self.assertLessEqual(mp.abs(res.real - mp.Float("1.6514099178148475691128039277241118340531698491")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-0.38121862770417378405072154507774424569831993182")), tol)
        res = mp.sqrt(z)
        self.assertLessEqual(mp.abs(res.real - mp.Float("0.14335482322754515813189359093523204816185445814")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-2.5496177371015053245565185485769652617478797292")), tol)
        res = mp.sin(x)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-0.66749329633668695550441899166308616328986315948")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-0.16020928942503633132090203927960650380076680938")), tol)
        res = mp.cos(y)
        self.assertLessEqual(mp.abs(res.real - mp.Float("0.45194679593300564730917329070452033759984813611")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-3.3798161097977088705360399142708324234265626016")), tol)
        res = mp.tan(z)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-0.11998086808607765336591715593714295443402911227")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-0.63859741450762243500349270264429166927927928889")), tol)
        res = mp.asin(x)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-1.4763431474004472804452143435221887167393328861")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("1.5406263884278099750127157814537559611048741005")), tol)
        res = mp.acos(y)
        self.assertLessEqual(mp.abs(res.real - mp.Float("0.38769800860408664087229892623614567735197135529")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("2.3379037587834289977359318611458042347281923566")), tol)
        res = mp.atan(z)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-1.4195347801361539102032503530226060969949192059")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-0.016801358511827150554928904776095870747673962940")), tol)
        res = mp.sinh(x)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-5.5116175435238027338707341903303682792175349461")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("1.1931117850318239903301857156967336540973461581")), tol)
        res = mp.cosh(y)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-22.821106324812396153305984541517740047129741299")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-58.969920999917202046449299163531439999586349377")), tol)
        res = mp.tanh(z)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-.99999948909538256503828034023523935671055767287")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-0.0000046773288796255165542679839050497228616710086412")), tol)
        res = mp.square(z)
        self.assertLessEqual(mp.abs(res.real - mp.Float("41.456039")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("9.47376")), tol)
        res = mp.cube(z)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-261.70981416")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("-91.694329309")), tol)





    def test_misc_funcs(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        #
        res = mp.norm(x)
        self.assertLessEqual(mp.abs(res - mp.Float("5.949")), tol)
        res = mp.abs2(x)
        self.assertLessEqual(mp.abs(res - mp.Float("5.949")), tol)
        res = mp.conj(x)
        self.assertLessEqual(mp.abs(res.real - x.real), tol)
        self.assertLessEqual(mp.abs(res.imag - (-x.imag)), tol)
        res = mp.polar(mp.Float("3.21"), mp.Float("-5.62"))
        self.assertLessEqual(mp.abs(res.real - mp.Float("2.5295931897050156212406076422629449344206513531")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("1.9761726378527774897831544771425943545375239972")), tol)
        res = mp.arg(y)
        self.assertLessEqual(mp.abs(res - mp.Float("-.38121862770417378405072154507774424569831993182")), tol)
        res = mp.inverse(z)
        self.assertLessEqual(mp.abs(res.real - mp.Float("-0.15238180880075963272315628064317633672297417498")), tol)
        self.assertLessEqual(mp.abs(res.imag - mp.Float("0.017189984912554828938368401412062021935878722517")), tol)

    def test_change_prec(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        mp.default_precision(45);
        tol = mp.Float("1e-37");
        t = mp.Complex("-2.43",".21")
        t = t**(-2);
        self.assertLessEqual(mp.abs(t.real - mp.Float("0.16560329111110602501494676510297183141930819429")), tol)
        self.assertLessEqual(mp.abs(t.imag - mp.Float("0.028838165251841866149715852522538399390278791818")), tol)
        mp.default_precision(30);
        tol = mp.Float("1e-27");



if __name__ == '__main__':
    unittest.main();
