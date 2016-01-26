__author__ = 'jcollins'

from libpybertini import *
from minieigen import *
import unittest
import numpy as np
import cmath


class SystemTest(unittest.TestCase):
    def test_system_create(self):
        x = Variable("x");
        y = Variable("y");
        f = Function(x*y);

        s = System();
        s.add_ungrouped_variable(x);
        s.add_ungrouped_variable(y);
        s.add_function(f)


    def test_system_eval(self):
        toldbl = 1e-15;
        x = Variable("x");
        y = Variable("y");
        z = Variable("z");
        a = Float(4.897, 1.23)

        f = Function(x*y);
        g = Function(pow(x,2)*y - a*z*x);
        exact_real = (-32.841085, -150.5480559)
        exact_imag = (-26.66705, -258.97936865)

        s = System();
        s.add_ungrouped_variable(x);
        s.add_ungrouped_variable(y);
        s.add_ungrouped_variable(z);
        s.add_function(f)
        s.add_function(g)

        v = VectorXc.Zero(3);
        v[0] = complex(3.5,2.89); v[1] = complex(-9.32,.0765); v[2] = complex(5.4,-2.13);

        e = s.eval(v)

        self.assertTrue(np.abs(e[0].real - exact_real[0]) < toldbl*np.abs(exact_real[0]));
        self.assertTrue(np.abs(e[0].imag - exact_imag[0]) < toldbl*np.abs(exact_imag[0]));
        self.assertTrue(np.abs(e[1].real - exact_real[1]) < toldbl*np.abs(exact_real[1]));
        self.assertTrue(np.abs(e[1].imag - exact_imag[1]) < toldbl*np.abs(exact_imag[1]));





unittest.main()