from libpybertini import *
import unittest


dbltol = 1e-15;
class MPFRGeneral(unittest.TestCase):
    def setUp(self):
        default_precision(30);
        self.x = mpfr_float("4.23")
        self.y = mpfr_float("-3.86")
        self.z = mpfr_float("1.1495")
        self.p = mpfr_float(".34")
        self.tol = mpfr_float("1e-26");

class MPFRFloatTest(MPFRGeneral):

    def test_arith_int(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        self.assertLessEqual(abs((x+8) - mpfr_float("12.23")), tol)
        self.assertLessEqual(abs((y-2) - mpfr_float("-5.86")), tol)
        self.assertLessEqual(abs((8+x) - mpfr_float("12.23")), tol)
        self.assertLessEqual(abs((2-y) - mpfr_float("5.86")), tol)
        self.assertLessEqual(abs((z*6) - mpfr_float("6.897")), tol)
        self.assertLessEqual(abs((6*z) - mpfr_float("6.897")), tol)
        self.assertLessEqual(abs((y/3) - mpfr_float("-1.2866666666666666666666666666666667")), tol)
        self.assertLessEqual(abs((3/y) - mpfr_float("-.77720207253886010362694300518134714")), tol)
        self.assertLessEqual(abs((x**3) - mpfr_float("75.686967")), tol)

        result = mpfr_float(x);
        result += 8;
        self.assertLessEqual(abs(result - mpfr_float("12.23")), tol)
        result = mpfr_float(y);
        result -= 2;
        self.assertLessEqual(abs(result - mpfr_float("-5.86")), tol)
        result = mpfr_float(z);
        result *= 6;
        self.assertLessEqual(abs(result - mpfr_float("6.897")), tol)
        result = mpfr_float(y);
        result /= 3;
        self.assertLessEqual(abs(result - mpfr_float("-1.2866666666666666666666666666666667")), tol)

        self.assertLessEqual(abs((z**3) - mpfr_float("1.518892112375")), tol)

    def test_arith_mpfr(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        self.assertLessEqual(abs((x+y) - mpfr_float("0.37")), tol)
        self.assertLessEqual(abs((z-y) - mpfr_float("5.0095")), tol)
        self.assertLessEqual(abs((z*y) - mpfr_float("-4.437070")), tol)
        self.assertLessEqual(abs((y/x) - mpfr_float("-.91252955082742316784869976359338061")), tol)
        self.assertLessEqual(abs((x**y) - mpfr_float("0.0038223124228935822000384505727705508")), tol)
        self.assertLessEqual(abs((-z) - mpfr_float("-1.1495")), tol)

        result = mpfr_float(x);
        result += y;
        self.assertLessEqual(abs(result - mpfr_float("0.37")), tol)
        result = mpfr_float(z);
        result -= y;
        self.assertLessEqual(abs(result - mpfr_float("5.0095")), tol)
        result = mpfr_float(z);
        result *= y;
        self.assertLessEqual(abs(result - mpfr_float("-4.437070")), tol)
        result = mpfr_float(y);
        result /= x;
        self.assertLessEqual(abs(result - mpfr_float("-.91252955082742316784869976359338061")), tol)


    def test_transendentals(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        self.assertLessEqual(abs((exp(x)) - mpfr_float("68.717232173846461408252914213396109")), tol)
        self.assertLessEqual(abs((log(z)) - mpfr_float("0.13932706522109918666170810230684295")), tol)
        self.assertLessEqual(abs((sqrt(z)) - mpfr_float("1.0721473779289860297522254519889560")), tol)

        self.assertLessEqual(abs((sin(x)) - mpfr_float("-.88588921129660245121088859729926237")), tol)
        self.assertLessEqual(abs((cos(y)) - mpfr_float("-.75285494656729525719980460936483635")), tol)
        self.assertLessEqual(abs((tan(z)) - mpfr_float("2.2315038042849919118711153687209483")), tol)

        self.assertLessEqual(abs((asin(p)) - mpfr_float("0.34691689752716170922069696210451452")), tol)
        self.assertLessEqual(abs((acos(p)) - mpfr_float("1.2238794292677349100106247295352369")), tol)
        self.assertLessEqual(abs((atan(z)) - mpfr_float("0.85483739856328448882289109284144652")), tol)

        self.assertLessEqual(abs((sinh(x)) - mpfr_float("34.351339891649022639414777866662100")), tol)
        self.assertLessEqual(abs((cosh(y)) - mpfr_float("23.743209684188284295743755381842167")), tol)
        self.assertLessEqual(abs((tanh(z)) - mpfr_float("0.81758837109637920976170104688035086")), tol)
