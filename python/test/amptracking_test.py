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


__author__ = 'jcollins'

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
        self.a = Float("4.897", "1.23")
        #
        self.f = Function(self.x*self.y);
        self.g = Function(pow(self.x,2)*self.y - self.a*self.z*self.x);

    def test_system_create(self):
        self.x = Variable("x");
        self.y = Variable("y");
        self.f = Function(self.x*self.y);
        #
        s = System();
        s.add_ungrouped_variable(self.x);
        s.add_ungrouped_variable(self.y);
        s.add_function(self.f)


    def test_system_eval(self):
        exact_real = (-32.841085, -150.5480559)
        exact_imag = (-26.66705, -258.97936865)
        #
        s = System();
        s.add_ungrouped_variable(self.x);
        s.add_ungrouped_variable(self.y);
        s.add_ungrouped_variable(self.z);
        s.add_function(self.f)
        s.add_function(self.g)
        #
        v = VectorXd.Zero(3);
        v[0] = complex(3.5,2.89); v[1] = complex(-9.32,.0765); v[2] = complex(5.4,-2.13);
        #
        e = s.eval(v)
        #
        self.assertTrue(np.abs(e[0].real - exact_real[0]) < self.toldbl*np.abs(exact_real[0]));
        self.assertTrue(np.abs(e[0].imag - exact_imag[0]) < self.toldbl*np.abs(exact_imag[0]));
        self.assertTrue(np.abs(e[1].real - exact_real[1]) < self.toldbl*np.abs(exact_real[1]));
        self.assertTrue(np.abs(e[1].imag - exact_imag[1]) < self.toldbl*np.abs(exact_imag[1]));
        #
        #
        s = parse_system('function f1, f2; variable_group x,y,z; f1 = x*y; f2 = x^2*y - z*x;')
        self.toldbl = mpfr_float('1e-27');
        exact_real = (mpfr_float('-32.841085'), mpfr_float('-62.9317230'))
        exact_imag = (mpfr_float('-26.66705'), mpfr_float('-196.39641065'))
        self.a = mpfr_complex('4.897', '1.23')
        v = VectorXmp((mpfr_complex('3.5', '2.89'), mpfr_complex('-9.32', '.0765'), mpfr_complex('5.4', '-2.13')));
        #
        e = s.eval(v)
        #
        self.assertLessEqual(abs(e[0].real / exact_real[0]-1) , self.toldbl);
        self.assertLessEqual(abs(e[0].imag / exact_imag[0]-1) , self.toldbl);
        self.assertLessEqual(abs(e[1].real / exact_real[1]-1) , self.toldbl);
        self.assertLessEqual(abs(e[1].imag / exact_imag[1]-1) , self.toldbl);





    def test_system_Jac(self):
        exact_real = ((-9.32, 3.5, 0), \
                      (-94.745870,3.8979,-13.5848))
        exact_imag = ((.0765, 2.89, 0),\
                      (-49.54549,20.230,-18.45733))
        #
        s = System();
        s.add_ungrouped_variable(self.x);
        s.add_ungrouped_variable(self.y);
        s.add_ungrouped_variable(self.z);
        s.add_function(self.f)
        s.add_function(self.g)
        #
        v = VectorXd.Zero(3);
        v[0] = complex(3.5,2.89); v[1] = complex(-9.32,.0765); v[2] = complex(5.4,-2.13);
        #
        s.differentiate();
        e = s.jacobian(v)
        #
        self.assertTrue(np.abs(e[0][0].real - exact_real[0][0]) <= self.toldbl*np.abs(exact_real[0][0]));
        self.assertTrue(np.abs(e[0][0].imag - exact_imag[0][0]) <= self.toldbl*np.abs(exact_imag[0][0]));
        self.assertTrue(np.abs(e[0][1].real - exact_real[0][1]) <= self.toldbl*np.abs(exact_real[0][1]));
        self.assertTrue(np.abs(e[0][1].imag - exact_imag[0][1]) <= self.toldbl*np.abs(exact_imag[0][1]));
        self.assertTrue(np.abs(e[0][2].real - exact_real[0][2]) <= self.toldbl*np.abs(exact_real[0][2]));
        self.assertTrue(np.abs(e[0][2].imag - exact_imag[0][2]) <= self.toldbl*np.abs(exact_imag[0][2]));
        self.assertTrue(np.abs(e[1][0].real - exact_real[1][0]) <= self.toldbl*np.abs(exact_real[1][0]));
        self.assertTrue(np.abs(e[1][0].imag - exact_imag[1][0]) <= self.toldbl*np.abs(exact_imag[1][0]));
        self.assertTrue(np.abs(e[1][1].real - exact_real[1][1]) <= self.toldbl*np.abs(exact_real[1][1]));
        self.assertTrue(np.abs(e[1][1].imag - exact_imag[1][1]) <= self.toldbl*np.abs(exact_imag[1][1]));
        self.assertTrue(np.abs(e[1][2].real - exact_real[1][2]) <= self.toldbl*np.abs(exact_real[1][2]));
        self.assertTrue(np.abs(e[1][2].imag - exact_imag[1][2]) <= self.toldbl*np.abs(exact_imag[1][2]));
        #
        #
        #
        #
        s = parse_system('function f1, f2; variable_group x,y,z; f1 = x*y; f2 = x^2*y - z*x;')
        self.toldbl = mpfr_float('1e-27');
        exact_real = ((mpfr_float('-9.32'), mpfr_float('3.5'), mpfr_float('0')), \
                      (mpfr_float('-71.082170'),mpfr_float('3.8979'),mpfr_float('-3.5')))
        exact_imag = ((mpfr_float('.0765'), mpfr_float('2.89'), mpfr_float('0')),\
                      (mpfr_float('-51.20410'),mpfr_float('20.230'),mpfr_float('-2.89')))
        v = VectorXmp((mpfr_complex('3.5', '2.89'), mpfr_complex('-9.32', '.0765'), mpfr_complex('5.4', '-2.13')));
        #
        s.differentiate();
        e = s.jacobian(v);
        #
        self.assertLessEqual(abs(e[0][0].real / exact_real[0][0]-1) , self.toldbl);
        self.assertLessEqual(abs(e[0][0].imag / exact_imag[0][0]-1) , self.toldbl);
        self.assertLessEqual(abs(e[0][1].real / exact_real[0][1]-1) , self.toldbl);
        self.assertLessEqual(abs(e[0][1].imag / exact_imag[0][1]-1) , self.toldbl);
        self.assertLessEqual(abs(e[0][2].real ) , self.toldbl);
        self.assertLessEqual(abs(e[0][2].imag ) , self.toldbl);
        self.assertLessEqual(abs(e[1][0].real / exact_real[1][0]-1) , self.toldbl);
        self.assertLessEqual(abs(e[1][0].imag / exact_imag[1][0]-1) , self.toldbl);
        self.assertLessEqual(abs(e[1][1].real / exact_real[1][1]-1) , self.toldbl);
        self.assertLessEqual(abs(e[1][1].imag / exact_imag[1][1]-1) , self.toldbl);
        self.assertLessEqual(abs(e[1][2].real / exact_real[1][2]-1) , self.toldbl);
        self.assertLessEqual(abs(e[1][2].imag / exact_imag[1][2]-1) , self.toldbl);



    def test_add_systems(self):
        x = self.x;  y = self.y;
        s1 = System(); s2 = System();
        #
        vars = VariableGroup();
        vars.append(x); vars.append(y);
        #
        s1.add_variable_group(vars)
        s1.add_function(y+1)
        s1.add_function(x*y)
        #
        s2.add_variable_group(vars)
        s2.add_function(-y-1)
        s2.add_function(-x*y)
        #
        s1 += s2;
        values = VectorXd((2,3))
        v = s1.eval(values)
        #
        self.assertEqual(v[0], 0.0)
        self.assertEqual(v[1], 0.0)
        #
        deg = s1.degrees()
        self.assertEqual(len(deg),2)
        #
        self.assertEqual(deg[0], 1)
        self.assertEqual(deg[1], 2)


    def test_mult_system_node(self):
        tol_d = self.toldbl;
        sys = parse_system('function f1, f2; variable_group x,y,z; f1 = x+2; f2 = y*y;')
        #
        z = Variable("z");
        sys *= Integer(2);
        #
        vals = VectorXd((complex(-2.43,.21 ),complex(4.84, -1.94),complex(-6.48, -.731)))
        sysEval = sys.eval(vals);
        #
        self.assertLessEqual(np.abs(sysEval[0].real / (-.86)-1), tol_d)
        self.assertLessEqual(np.abs(sysEval[0].imag / (0.42)-1), tol_d)
        self.assertLessEqual(np.abs(sysEval[1].real / (39.3240)-1), tol_d)
        self.assertLessEqual(np.abs(sysEval[1].imag / (-37.5584)-1), tol_d)



if __name__ == '__main__':
    unittest.main();


