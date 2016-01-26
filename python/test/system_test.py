__author__ = 'jcollins'

from libpybertini import *
from minieigen import *
import unittest
import numpy as np
import cmath


class SystemTest(unittest.TestCase):
    def setUp(self):
        self.toldbl = 1e-15;
        self.x = Variable("x");
        self.y = Variable("y");
        self.z = Variable("z");
        self.a = Float(4.897, 1.23)

        self.f = Function(self.x*self.y);
        self.g = Function(pow(self.x,2)*self.y - self.a*self.z*self.x);

    def test_system_create(self):
        self.x = Variable("x");
        self.y = Variable("y");
        self.f = Function(self.x*self.y);

        s = System();
        s.add_ungrouped_variable(self.x);
        s.add_ungrouped_variable(self.y);
        s.add_function(self.f)


    def test_system_eval(self):
        exact_real = (-32.841085, -150.5480559)
        exact_imag = (-26.66705, -258.97936865)

        s = System();
        s.add_ungrouped_variable(self.x);
        s.add_ungrouped_variable(self.y);
        s.add_ungrouped_variable(self.z);
        s.add_function(self.f)
        s.add_function(self.g)

        v = VectorXc.Zero(3);
        v[0] = complex(3.5,2.89); v[1] = complex(-9.32,.0765); v[2] = complex(5.4,-2.13);

        e = s.eval(v)

        self.assertTrue(np.abs(e[0].real - exact_real[0]) < self.toldbl*np.abs(exact_real[0]));
        self.assertTrue(np.abs(e[0].imag - exact_imag[0]) < self.toldbl*np.abs(exact_imag[0]));
        self.assertTrue(np.abs(e[1].real - exact_real[1]) < self.toldbl*np.abs(exact_real[1]));
        self.assertTrue(np.abs(e[1].imag - exact_imag[1]) < self.toldbl*np.abs(exact_imag[1]));

    def test_system_Jac(self):
        exact_real = ((-9.32, 3.5, 0), \
                      (-94.745870,3.8979,-13.5848))
        exact_imag = ((.0765, 2.89, 0),\
                      (-49.54549,20.230,-18.45733))

        s = System();
        s.add_ungrouped_variable(self.x);
        s.add_ungrouped_variable(self.y);
        s.add_ungrouped_variable(self.z);
        s.add_function(self.f)
        s.add_function(self.g)

        v = VectorXc.Zero(3);
        v[0] = complex(3.5,2.89); v[1] = complex(-9.32,.0765); v[2] = complex(5.4,-2.13);

        e = s.jacobian(v)

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



# if __name__ == '__main__':
#     unittest.main()