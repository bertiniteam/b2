from libpybertini import *
import unittest
import pdb


dbltol = 1e-15;
class MPFRFloat(unittest.TestCase):
    def setUp(self):
        default_precision(30);
        self.x = mpfr_float("4.23")
        self.y = mpfr_float("-3.86")
        self.z = mpfr_float("1.1495")
        self.p = mpfr_float(".34")
        self.tol = mpfr_float("1e-27");

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

        self.assertLessEqual(abs((-z) - mpfr_float("-1.1495")), tol)

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


    def test_trancendentals(self):
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

    def test_change_prec(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        default_precision(40);
        tol = mpfr_float("1e-37");
        t = mpfr_float("4.23")
        self.assertLessEqual(abs(t**(-2) - mpfr_float("0.055888089689206333238323580861682566828183245868")), tol)
        default_precision(30);
        tol = mpfr_float("1e-27");










class MPFRComplex(unittest.TestCase):
    def setUp(self):
        default_precision(30);
        self.x = mpfr_complex("-2.43",".21" )
        self.y = mpfr_complex("4.84", "-1.94")
        self.z = mpfr_complex("-6.48", "-.731")
        self.p = mpfr_complex("-.321", "-.72")
        self.tol = mpfr_float("1e-27");

    def test_construct(self):
        t = mpfr_complex(3.452)
        t = mpfr_complex(mpfr_float("-5.6"))
        t = mpfr_complex("3.89")
        t = mpfr_complex(mpfr_float("2.98"), mpfr_float("-1e-4"))
        t = mpfr_complex(3.4, 3.5)
        t = mpfr_complex("6e2", mpfr_float("4.32"))
        t = mpfr_complex(mpfr_float("4.32"), "6e2")

    def test_arith_mp_float(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        a = mpfr_float("3.12"); b = mpfr_float("-5.92")
        res = mpfr_complex(x+a)
        self.assertLessEqual(abs(res.real - mpfr_float("0.69")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("0.21")), tol)
        res = mpfr_complex(y-b)
        self.assertLessEqual(abs(res.real - mpfr_float("10.76")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-1.94")), tol)
        res = mpfr_complex(a+x)
        self.assertLessEqual(abs(res.real - mpfr_float("0.69")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("0.21")), tol)
        res = mpfr_complex(b-y)
        self.assertLessEqual(abs(res.real - mpfr_float("-10.76")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("1.94")), tol)
        res = mpfr_complex(z*a)
        self.assertLessEqual(abs(res.real - mpfr_float("-20.2176")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-2.28072")), tol)
        res = mpfr_complex(a*z)
        self.assertLessEqual(abs(res.real - mpfr_float("-20.2176")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-2.28072")), tol)
        res = mpfr_complex(y/b)
        self.assertLessEqual(abs(res.real - mpfr_float("-.81756756756756756756756756756756756756756756757")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float(".3277027027027027027027027027027027027027027027")), tol)
        res = mpfr_complex(b/y)
        self.assertLessEqual(abs(res.real - mpfr_float("-1.0538301972842157915642975887484736586585850264")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-.42240301296102864372618539714298324334662292381")), tol)
        res = mpfr_complex(x**a)
        self.assertLessEqual(abs(res.real - mpfr_float("-16.054376621961088182387920766649714821973863952")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-1.7411284591111236754359685247799914985638458821")), tol)



        res = mpfr_complex(x);
        res += a;
        self.assertLessEqual(abs(res.real - mpfr_float("0.69")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("0.21")), tol)
        res = mpfr_complex(y);
        res -= b;
        self.assertLessEqual(abs(res.real - mpfr_float("10.76")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-1.94")), tol)
        res = mpfr_complex(z);
        res *= a;
        self.assertLessEqual(abs(res.real - mpfr_float("-20.2176")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-2.28072")), tol)
        res = mpfr_complex(y);
        res /= b;
        self.assertLessEqual(abs(res.real - mpfr_float("-.81756756756756756756756756756756756756756756757")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float(".3277027027027027027027027027027027027027027027")), tol)

        res = mpfr_complex(x**4);
        self.assertLessEqual(abs(res.real - mpfr_float("33.30735228")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-11.96306496")), tol)

        res = mpfr_complex(-z);
        self.assertLessEqual(abs(res.real - mpfr_float("6.48")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float(".731")), tol)



    def test_arith_mp_complex(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;

        res = mpfr_complex(x+y)
        self.assertLessEqual(abs(res.real - mpfr_float("2.41")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-1.73")), tol)
        res = mpfr_complex(y-z)
        self.assertLessEqual(abs(res.real - mpfr_float("11.32")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-1.209")), tol)
        res = mpfr_complex(y+x)
        self.assertLessEqual(abs(res.real - mpfr_float("2.41")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-1.73")), tol)
        res = mpfr_complex(z-y)
        self.assertLessEqual(abs(res.real - mpfr_float("-11.32")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("1.209")), tol)
        res = mpfr_complex(z*x)
        self.assertLessEqual(abs(res.real - mpfr_float("15.89991")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("0.41553")), tol)
        res = mpfr_complex(x*z)
        self.assertLessEqual(abs(res.real - mpfr_float("15.89991")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float(".41553")), tol)
        res = mpfr_complex(y/x)
        self.assertLessEqual(abs(res.real - mpfr_float("-2.0454866364094805849722642460917801311144730207")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("0.62158345940494200706001008572869389813414019163")), tol)
        res = mpfr_complex(x/y)
        self.assertLessEqual(abs(res.real - mpfr_float("-.44755270475041560619657805305047592426404601827")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-.13600253041648890000441351713180233327939034617")), tol)
        res = mpfr_complex(y**z)
        self.assertLessEqual(abs(res.real - mpfr_float("0.0000051612634484879218649489640888954160904291899461")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("0.16242051733741136410199105656393042124100116889e-4")), tol)



        res = mpfr_complex(x);
        res += y;
        self.assertLessEqual(abs(res.real - mpfr_float("2.41")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-1.73")), tol)
        res = mpfr_complex(y);
        res -= z;
        self.assertLessEqual(abs(res.real - mpfr_float("11.32")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-1.209")), tol)
        res = mpfr_complex(z);
        res *= x;
        self.assertLessEqual(abs(res.real - mpfr_float("15.89991")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float(".41553")), tol)
        res = mpfr_complex(y);
        res /= x;
        self.assertLessEqual(abs(res.real - mpfr_float("-2.0454866364094805849722642460917801311144730207")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("0.62158345940494200706001008572869389813414019163")), tol)


    def test_trancendentals(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;

        res = exp(x)
        self.assertLessEqual(abs(res.real - mpfr_float("0.086102743899954532232498058731947255067424332219")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("0.018352149302889219131202317785160051395400089327")), tol)
        res = log(y)
        self.assertLessEqual(abs(res.real - mpfr_float("1.6514099178148475691128039277241118340531698491")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-0.38121862770417378405072154507774424569831993182")), tol)
        res = sqrt(z)
        self.assertLessEqual(abs(res.real - mpfr_float("0.14335482322754515813189359093523204816185445814")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-2.5496177371015053245565185485769652617478797292")), tol)
        res = sin(x)
        self.assertLessEqual(abs(res.real - mpfr_float("-0.66749329633668695550441899166308616328986315948")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-0.16020928942503633132090203927960650380076680938")), tol)
        res = cos(y)
        self.assertLessEqual(abs(res.real - mpfr_float("0.45194679593300564730917329070452033759984813611")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-3.3798161097977088705360399142708324234265626016")), tol)
        res = tan(z)
        self.assertLessEqual(abs(res.real - mpfr_float("-0.11998086808607765336591715593714295443402911227")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-0.63859741450762243500349270264429166927927928889")), tol)
        res = asin(x)
        self.assertLessEqual(abs(res.real - mpfr_float("-1.4763431474004472804452143435221887167393328861")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("1.5406263884278099750127157814537559611048741005")), tol)
        res = acos(y)
        self.assertLessEqual(abs(res.real - mpfr_float("0.38769800860408664087229892623614567735197135529")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("2.3379037587834289977359318611458042347281923566")), tol)
        res = atan(z)
        self.assertLessEqual(abs(res.real - mpfr_float("-1.4195347801361539102032503530226060969949192059")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-0.016801358511827150554928904776095870747673962940")), tol)
        res = sinh(x)
        self.assertLessEqual(abs(res.real - mpfr_float("-5.5116175435238027338707341903303682792175349461")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("1.1931117850318239903301857156967336540973461581")), tol)
        res = cosh(y)
        self.assertLessEqual(abs(res.real - mpfr_float("-22.821106324812396153305984541517740047129741299")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-58.969920999917202046449299163531439999586349377")), tol)
        res = tanh(z)
        self.assertLessEqual(abs(res.real - mpfr_float("-.99999948909538256503828034023523935671055767287")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-0.0000046773288796255165542679839050497228616710086412")), tol)
        res = square(z)
        self.assertLessEqual(abs(res.real - mpfr_float("41.456039")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("9.47376")), tol)
        res = cube(z)
        self.assertLessEqual(abs(res.real - mpfr_float("-261.70981416")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("-91.694329309")), tol)





    def test_misc_funcs(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;

        res = norm(x)
        self.assertLessEqual(abs(res - mpfr_float("5.949")), tol)
        res = abs2(x)
        self.assertLessEqual(abs(res - mpfr_float("5.949")), tol)
        res = conj(x)
        self.assertLessEqual(abs(res.real - x.real), tol)
        self.assertLessEqual(abs(res.imag - (-x.imag)), tol)
        res = polar(mpfr_float("3.21"), mpfr_float("-5.62"))
        self.assertLessEqual(abs(res.real - mpfr_float("2.5295931897050156212406076422629449344206513531")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("1.9761726378527774897831544771425943545375239972")), tol)
        res = arg(y)
        self.assertLessEqual(abs(res - mpfr_float("-.38121862770417378405072154507774424569831993182")), tol)
        res = inverse(z)
        self.assertLessEqual(abs(res.real - mpfr_float("-0.15238180880075963272315628064317633672297417498")), tol)
        self.assertLessEqual(abs(res.imag - mpfr_float("0.017189984912554828938368401412062021935878722517")), tol)

    def test_change_prec(self):
        x = self.x; y = self.y; z = self.z; p = self.p; tol = self.tol;
        default_precision(45);
        tol = mpfr_float("1e-37");
        t = mpfr_complex("-2.43",".21")
        t = t**(-2);
        self.assertLessEqual(abs(t.real - mpfr_float("0.16560329111110602501494676510297183141930819429")), tol)
        self.assertLessEqual(abs(t.imag - mpfr_float("0.028838165251841866149715852522538399390278791818")), tol)
        default_precision(30);
        tol = mpfr_float("1e-27");
