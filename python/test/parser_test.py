__author__ = 'jcollins'

from libpybertini import *
import unittest
import numpy as np
import pdb


class ParserTest(unittest.TestCase):
    def setUp(self):
        self.tol_d = 1e-15;



    def test_create_system(self):
        tol_d = self.tol_d;
        input = 'function f, g; variable_group x,y,z; f = 3*x*y*z; g = x^2 + y^2 + z^2 - 1;';
        sys = parse_system(input);

        vals = VectorXd((complex(-2.43,.21 ),complex(4.84, -1.94),complex(-6.48, -.731)))
        sysEval = sys.eval(vals);

        self.assertLessEqual(np.abs(sysEval[0].real / (233.2850778)-1), tol_d)
        self.assertLessEqual(np.abs(sysEval[0].imag / (-86.5039806)-1), tol_d)
        self.assertLessEqual(np.abs(sysEval[1].real / (65.978839)-1), tol_d)
        self.assertLessEqual(np.abs(sysEval[1].imag / (-10.32604)-1), tol_d)

        sys.differentiate()
        sysJac = sys.jacobian(vals)