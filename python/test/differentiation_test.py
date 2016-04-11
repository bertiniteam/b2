from pybertini import *
import numpy as np;
import unittest
import pdb


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

        f = x+y+a;
        df = f.differentiate();

        self.assertLessEqual(np.abs(df.evald(x).real/(1.0)-1), tol_d)
        self.assertLessEqual(np.abs(df.evald(x).imag-0), tol_d)

        self.assertLessEqual(abs(df.evalmp(x).real/mpfr_float("1.0")-1), tol_mp)
        self.assertLessEqual(abs(df.evalmp(x).imag-mpfr_float("0")), tol_mp)

        df.reset()
        self.assertLessEqual(np.abs(df.evald(y).real/(1.0)-1), tol_d)
        self.assertLessEqual(np.abs(df.evald(y).imag-0), tol_d)

        self.assertLessEqual(abs(df.evalmp(y).real/mpfr_float("1.0")-1), tol_mp)
        self.assertLessEqual(abs(df.evalmp(y).imag-mpfr_float("0")), tol_mp)

        df.reset()
        self.assertLessEqual(np.abs(df.evald(z).real-0), tol_d)
        self.assertLessEqual(np.abs(df.evald(z).imag-0), tol_d)

        self.assertLessEqual(abs(df.evalmp(z).real-mpfr_float("0.0")), tol_mp)
        self.assertLessEqual(abs(df.evalmp(z).imag-mpfr_float("0")), tol_mp)

    def test_power_rule(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;

        f = x**2 + y**3;
        df = f.differentiate();

        self.assertLessEqual(np.abs(df.evald(x).real / (-4.86)-1), tol_d)
        self.assertLessEqual(np.abs(df.evald(x).imag / (0.42)-1), tol_d)

        self.assertLessEqual(abs(df.evalmp(x).real / mpfr_float("-4.86")-1), tol_mp)
        self.assertLessEqual(abs(df.evalmp(x).imag / mpfr_float("0.42")-1), tol_mp)

        df.reset()
        self.assertLessEqual(np.abs(df.evald(y).real / (58.9860)-1), tol_d)
        self.assertLessEqual(np.abs(df.evald(y).imag / (-56.3376)-1), tol_d)

        self.assertLessEqual(abs(df.evalmp(y).real / mpfr_float("58.9860")-1), tol_mp)
        self.assertLessEqual(abs(df.evalmp(y).imag / mpfr_float("-56.3376")-1), tol_mp)


    def test_prod_rule(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;

        f = x**2*y**4 - a*x*y*z**2;
        df = f.differentiate();

        self.assertLessEqual(np.abs(df.evald(x).real / (-559.28968169592)-1), tol_d)
        self.assertLessEqual(np.abs(df.evald(x).imag / (3577.05276993648)-1), tol_d)

        self.assertLessEqual(abs(df.evalmp(x).real / mpfr_float("-559.28968169592")-1), tol_mp)
        self.assertLessEqual(abs(df.evalmp(x).imag / mpfr_float("3577.05276993648")-1), tol_mp)

        df.reset()
        self.assertLessEqual(np.abs(df.evald(y).real / (1161.85042980828)-1), tol_d)
        self.assertLessEqual(np.abs(df.evald(y).imag / (-3157.24325320476)-1), tol_d)

        self.assertLessEqual(abs(df.evalmp(y).real / mpfr_float("1161.85042980828")-1), tol_mp)
        self.assertLessEqual(abs(df.evalmp(y).imag / mpfr_float("-3157.24325320476")-1), tol_mp)

        df.reset()
        self.assertLessEqual(np.abs(df.evald(z).real / (-520.5265859088)-1), tol_d)
        self.assertLessEqual(np.abs(df.evald(z).imag / (84.7479679056)-1), tol_d)

        self.assertLessEqual(abs(df.evalmp(z).real / mpfr_float("-520.5265859088")-1), tol_mp)
        self.assertLessEqual(abs(df.evalmp(z).imag / mpfr_float("84.7479679056")-1), tol_mp)


    def test_trancendental(self):
        x = self.x; y = self.y; z = self.z; p = self.p; a = self.a; b = self.b;
        tol_d = self.tol_d; tol_mp = self.tol_mp;

        f = sin(x*y) + exp(z*y) - log(x*x);
        df = f.differentiate();

        self.assertLessEqual(np.abs(df.evald(x).real / (-17.648420086229721902138620795382021306411662490)-1), 9e-13)
        self.assertLessEqual(np.abs(df.evald(x).imag / (-803.11883403426275105632833868183320319093878729)-1), tol_d)

        self.assertLessEqual(abs(df.evalmp(x).real / mpfr_float("-17.648420086229721902138620795382021306411662490")-1), tol_mp)
        self.assertLessEqual(abs(df.evalmp(x).imag / mpfr_float("-803.11883403426275105632833868183320319093878729")-1), tol_mp)

        df.reset()
        self.assertLessEqual(np.abs(df.evald(y).real / (-100.97157179433748763552280062599971478593963953)-1), tol_d)
        self.assertLessEqual(np.abs(df.evald(y).imag / (361.98093991820979266721712882115615553425318528)-1), tol_d)

        self.assertLessEqual(abs(df.evalmp(y).real / mpfr_float("-100.97157179433748763552280062599971478593963953")-1), tol_mp)
        self.assertLessEqual(abs(df.evalmp(y).imag / mpfr_float("361.98093991820979266721712882115615553425318528")-1), tol_mp)

        df.reset()
        self.assertLessEqual(np.abs(df.evald(z).real / (-2.1642907643013779167501866500194314960002972412e-14)-1), tol_d)
        self.assertLessEqual(np.abs(df.evald(z).imag / (2.1105887207247540399884720817624768568595288922e-14)-1), tol_d)

        self.assertLessEqual(abs(df.evalmp(z).real / mpfr_float("-2.1642907643013779167501866500194314960002972412e-14")-1), tol_mp)
        self.assertLessEqual(abs(df.evalmp(z).imag / mpfr_float("2.1105887207247540399884720817624768568595288922e-14")-1), tol_mp)
