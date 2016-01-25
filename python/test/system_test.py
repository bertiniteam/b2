__author__ = 'jcollins'

from libpybertini import *
import unittest
import cmath


class SystemTest(unittest.TestCase):
    def test_system_create(self):
        x = Variable("x");
        y = Variable("y");
        f = Function(x*y);

        s = System;
        s.add_ungrouped_variables(x);
        s.add_ungrouped_variables(y);
        s.add_function(f)

        