//This file is part of Bertini 2.
//
//complex_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//complex_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with complex_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire



#include "bertini2/num_traits.hpp"
#include <boost/test/unit_test.hpp>
#include <fstream>


#include "eigen_extensions.hpp"
#include <Eigen/Dense>

#include "externs.hpp"


BOOST_AUTO_TEST_SUITE(complex_multiprecision_class)

using mpfr_float = bertini::mpfr_float;
using bertini::DefaultPrecision;

BOOST_AUTO_TEST_CASE(complex_create_default_constructor)
{
	bertini::complex z;
}


BOOST_AUTO_TEST_CASE(complex_create_two_input_constructor)
{
	bertini::complex z("0.1","1.2");
}


BOOST_AUTO_TEST_CASE(complex_addition)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	bertini::complex z("0.1","1.2");
	bertini::complex v("0.2","1.3");
	
	bertini::complex r = z+v;
	BOOST_CHECK(abs(r.real()-bertini::mpfr_float("0.3")) < threshold_clearance_mp);
	BOOST_CHECK(abs(r.imag()-bertini::mpfr_float("2.5")) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(complex_subtraction)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	bertini::complex z("0.1","1.2");
	bertini::complex v("0.2","1.3");
	
	bertini::complex r = z-v;
	BOOST_CHECK(abs(r.real()-bertini::mpfr_float("-0.1")) < threshold_clearance_mp);
	BOOST_CHECK(abs(r.imag()-bertini::mpfr_float("-0.1")) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(complex_negation)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	bertini::complex z("0.1","1.2");

	bertini::complex r = -z;

	BOOST_CHECK_EQUAL(r.real(),-bertini::mpfr_float("0.1"));
	BOOST_CHECK_EQUAL(r.imag(),-bertini::mpfr_float("1.2"));
}

BOOST_AUTO_TEST_CASE(complex_multiplication)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex z("0.1000","1.2000");
	bertini::complex v("0.2000","1.3000");
	bertini::complex r = z*v;
	
	BOOST_CHECK(abs(r.real()-bertini::mpfr_float("-1.54")) < threshold_clearance_mp);
	BOOST_CHECK(abs(r.imag()-bertini::mpfr_float("0.37"))< threshold_clearance_mp);
	
}


BOOST_AUTO_TEST_CASE(complex_division)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex z("1.5","2.25");
	bertini::complex v("-3.1","5.1");
	bertini::complex r = z/v;
	
	BOOST_CHECK(abs(r.real()-bertini::mpfr_float("0.191605839416058394160583941605839416058394160583941605839416")) < threshold_clearance_mp);
	BOOST_CHECK(abs(r.imag()-bertini::mpfr_float("-0.410583941605839416058394160583941605839416058394160583941606"))< threshold_clearance_mp);
	
}





BOOST_AUTO_TEST_CASE(complex_inverse)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = inverse(z);
	
	BOOST_CHECK( abs(real(w) - mpfr_float("0.205128205128205128205128205128205128205128205128205128205128")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK( abs(imag(w) - mpfr_float("-0.307692307692307692307692307692307692307692307692307692307692")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}





BOOST_AUTO_TEST_CASE(complex_sqrt)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = sqrt(z);
	
	BOOST_CHECK(abs(real(w)- mpfr_float("1.44985576120488481677238036203436657121811982450524518214799")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("0.775939255547105301187773883610045740361564192945364951422177")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}




BOOST_AUTO_TEST_CASE(complex_log)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = log(z);
	
	BOOST_CHECK(abs(real(w)- mpfr_float("0.994792606278987440587524714788831870899124261482342501702857")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("0.982793723247329067985710611014666014496877453631628556761425")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}





BOOST_AUTO_TEST_CASE(complex_exponential)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25"), v("-3.1", "5.1");
	bertini::complex w = pow(z,v);
	
	BOOST_CHECK(abs(real(w)- mpfr_float("-0.000134184252069350633500916903057767602709221538143686248955454")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("0.000273589335652288672638569400353054443664478211013806804752627")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}



/////////////////trig


BOOST_AUTO_TEST_CASE(complex_sin)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = sin(z);
	
	BOOST_CHECK(abs(real(w)- mpfr_float("4.78455206454183468376293313858001326644589827567133864164966")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("0.331840118511466433278036564370915664859189055829310017561679")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}

BOOST_AUTO_TEST_CASE(complex_cos)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = cos(z);
	
	
	BOOST_CHECK(abs(real(w)- mpfr_float("0.339295764714918536764798162132262895147653901636178423640742")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("-4.67941686644942009239504224985289863706184473303907201576996")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}


BOOST_AUTO_TEST_CASE(complex_tan)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = tan(z);
	
	BOOST_CHECK(abs(real(w)- mpfr_float("0.0032055151478664658165188420647201618843131312777936688172756")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("1.02223515949711949776179466571527189987176458555122839197245")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}











//////////hyperbolic trig

BOOST_AUTO_TEST_CASE(complex_sinh)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = sinh(z);
	
	BOOST_CHECK(abs(real(w)- mpfr_float("-1.33755718909601135383265663229426431626814744207128061926332")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("1.83034686972219823005238120890074929109631502976918245359956")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}

BOOST_AUTO_TEST_CASE(complex_cosh)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = cosh(z);
	
	
	BOOST_CHECK(abs(real(w)- mpfr_float("-1.47772167013515546574129058780587714414840064875386904640065")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("1.65673527269339558987730878420158930664270667966911313867276")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}


BOOST_AUTO_TEST_CASE(complex_tanh)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = tanh(z);
	
	BOOST_CHECK(abs(real(w)- mpfr_float("1.01633467755934300564918114729269683730547835877531900704421")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("-0.0991725055603753242720729815306256466293467541150851175179255")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}








////////////////////// arc trig functions


BOOST_AUTO_TEST_CASE(complex_asin)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = asin(z);
	
	BOOST_CHECK(abs(real(w)- mpfr_float("0.557728098701882481783572689617979376112207607697432876677129")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("1.70220126917599043453132270948063371483983368727808240752544")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}

BOOST_AUTO_TEST_CASE(complex_acos)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = acos(z);
	
	BOOST_CHECK(abs(real(w)- mpfr_float("1.01306822809301413744774900202177206598637709199012003381034")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("-1.70220126917599043453132270948063371483983368727808240752544")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}


BOOST_AUTO_TEST_CASE(complex_atan)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = atan(z);

	BOOST_CHECK(abs(real(w)- mpfr_float("1.34897118928106882765457017435994248505893277188774702053495")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("0.303034028741274232429033400709495120396582520146639573148632")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}




////////////////////// arc hyperbolic trig functions


BOOST_AUTO_TEST_CASE(complex_asinh)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = asinh(z);
	
	BOOST_CHECK(abs(real(w)- mpfr_float("1.67614738083746762411099043767860402156031800041110901006856")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("0.950053155715849776576329694197937868331910565261475162022272")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}

BOOST_AUTO_TEST_CASE(complex_acosh)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = acosh(z);
	
	BOOST_CHECK(abs(real(w)- mpfr_float("1.70220126917599043453132270948063371483983368727808240752544")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("1.01306822809301413744774900202177206598637709199012003381034")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}


BOOST_AUTO_TEST_CASE(complex_atanh)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.5", "2.25");
	bertini::complex w = atanh(z);
	
	BOOST_CHECK(abs(real(w)- mpfr_float("0.188961443693877322997198271883423688031632018662054076093211")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- mpfr_float("1.26114018722767258645714375731198514259936666747537259674077")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}








BOOST_AUTO_TEST_CASE(complex_absolute_value)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	BOOST_CHECK_EQUAL(abs(bertini::complex("1.0","1.0")), sqrt(mpfr_float("2.0")));
	
}

BOOST_AUTO_TEST_CASE(complex_abs_of_i_is_one)
{
	bertini::complex i("0.0","1.0");
	bertini::mpfr_float one = abs(i);
	BOOST_CHECK_EQUAL(one,bertini::mpfr_float("1.0"));
}




BOOST_AUTO_TEST_CASE(complex_eleventh_power)
{
	using mpfr_float = bertini::mpfr_float;
	bertini::complex z("2.0","0.5"), w;
	
	w = pow(z,11);
	BOOST_CHECK_EQUAL(real(w), mpfr_float("-1319867") / mpfr_float("512") );
	BOOST_CHECK_EQUAL(imag(w),  mpfr_float("2529647") / mpfr_float("2048"));
	
	
}




BOOST_AUTO_TEST_CASE(complex_conjugation)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("2.0", sqrt(mpfr_float("2.0")));
	bertini::complex w = conj(z);
	
	BOOST_CHECK_EQUAL(real(w),mpfr_float("2.0"));
	BOOST_CHECK_EQUAL(imag(w),-sqrt(mpfr_float("2.0")));
}






BOOST_AUTO_TEST_CASE(complex_argument)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex z("1.0","1.0");
	
	BOOST_CHECK_EQUAL(arg(z), acos(mpfr_float("-1.0"))/mpfr_float("4.0"));
}





BOOST_AUTO_TEST_CASE(complex_make_random_50)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);  
	
	bertini::complex z = bertini::complex::rand();

}



BOOST_AUTO_TEST_CASE(complex_make_random_100)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(100);  
	bertini::complex z = bertini::complex::rand();


}



BOOST_AUTO_TEST_CASE(complex_serialization)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	
	bertini::complex z("1.23456", acos(mpfr_float("-1")));
	
	{
		std::ofstream fout("serialization_test_complex");
		
		boost::archive::text_oarchive oa(fout);
		
		// write class instance to archive
		oa << z;
	}
	
	//	bertini::complex w;
	bertini::complex w;
	{
		std::ifstream fin("serialization_test_complex");
		
		boost::archive::text_iarchive ia(fin);
		// read class state from archive
		ia >> w;
	}
	
	BOOST_CHECK_EQUAL(real(z),real(w));
	BOOST_CHECK_EQUAL(imag(z),imag(w));
	
}





BOOST_AUTO_TEST_CASE(complex_precision_predictable)
{
	DefaultPrecision(30);
	bertini::complex a(1,2);
	BOOST_CHECK_EQUAL(Precision(a),30);

	DefaultPrecision(50);
	bertini::complex b(3,4);
	BOOST_CHECK_EQUAL(Precision(b),50);

	DefaultPrecision(70);
	bertini::complex c(5,6);
	BOOST_CHECK_EQUAL(Precision(c),70);

	a = b;
	BOOST_CHECK_EQUAL(Precision(a),50);

	a = b+c;
	BOOST_CHECK_EQUAL(Precision(a),70);


	DefaultPrecision(90);
	a = b+c;
	BOOST_CHECK_EQUAL(Precision(a),90);

	a = b;
	BOOST_CHECK_EQUAL(Precision(a),50);
}





BOOST_AUTO_TEST_CASE(complex_get_from_stream_parens_with_comma)
{
	std::stringstream ss;
	ss << "(0.1234,-4.12397)";
	
	bertini::complex z;
	ss >> z;
	
	BOOST_CHECK_EQUAL(z.real(),bertini::mpfr_float("0.1234"));
	BOOST_CHECK_EQUAL(z.imag(),bertini::mpfr_float("-4.12397"));
}


BOOST_AUTO_TEST_CASE(complex_get_from_stream_parens_no_comma)
{
	std::stringstream ss;
	ss << "(-3.651263418976498712e-2)";
	
	bertini::complex z;
	ss >> z;
	
	BOOST_CHECK_EQUAL(z.real(),bertini::mpfr_float("-3.651263418976498712e-2"));
	BOOST_CHECK_EQUAL(z.imag(),bertini::mpfr_float("0"));
}


BOOST_AUTO_TEST_CASE(complex_get_from_stream_no_parens)
{
	std::stringstream ss;
	ss << "-3.651263418976498712e-2";
	
	bertini::complex z;
	ss >> z;
	
	BOOST_CHECK_EQUAL(z.real(),bertini::mpfr_float("-3.651263418976498712e-2"));
	BOOST_CHECK_EQUAL(z.imag(),bertini::mpfr_float("0"));
}




BOOST_AUTO_TEST_CASE(precision_of_mpfr_is_16)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex a("1.23124");
	BOOST_CHECK_EQUAL(bertini::Precision(a), CLASS_TEST_MPFR_DEFAULT_DIGITS);
}

BOOST_AUTO_TEST_CASE(precision_of_mpfr_complex_is_16)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex a("1.23124","-0.6789124678912394");
	BOOST_CHECK_EQUAL(bertini::Precision(a), CLASS_TEST_MPFR_DEFAULT_DIGITS);
}


BOOST_AUTO_TEST_CASE(precision_random_real_default)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	auto a = bertini::complex::RandomReal();
	BOOST_CHECK_EQUAL(bertini::Precision(a), CLASS_TEST_MPFR_DEFAULT_DIGITS);
	BOOST_CHECK_EQUAL(CLASS_TEST_MPFR_DEFAULT_DIGITS, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(precision_random_real_highest)
{	using namespace bertini;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex a;
	RandomReal(a,1000);
	BOOST_CHECK_EQUAL(bertini::Precision(a), 1000);
	BOOST_CHECK_EQUAL(CLASS_TEST_MPFR_DEFAULT_DIGITS, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(precision_random_unit_default)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	auto a = bertini::complex::RandomUnit();
	BOOST_CHECK_EQUAL(bertini::Precision(a), CLASS_TEST_MPFR_DEFAULT_DIGITS);
	BOOST_CHECK_EQUAL(CLASS_TEST_MPFR_DEFAULT_DIGITS, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(precision_random_unit_highest)
{	using namespace bertini;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex a;
	bertini::RandomUnit(a,1000);
	BOOST_CHECK_EQUAL(bertini::Precision(a), 1000);
	BOOST_CHECK_EQUAL(CLASS_TEST_MPFR_DEFAULT_DIGITS, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(precision_random_default)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	auto a = bertini::complex::rand();
	BOOST_CHECK_EQUAL(bertini::Precision(a), CLASS_TEST_MPFR_DEFAULT_DIGITS);
	BOOST_CHECK_EQUAL(CLASS_TEST_MPFR_DEFAULT_DIGITS, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(precision_random_highest)
{	using namespace bertini;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex a;
	bertini::rand(a,1000);
	BOOST_CHECK_EQUAL(Precision(a), 1000);
	BOOST_CHECK_EQUAL(CLASS_TEST_MPFR_DEFAULT_DIGITS, DefaultPrecision());
}


BOOST_AUTO_TEST_SUITE_END()









BOOST_AUTO_TEST_SUITE(complex_eigen_compatibility)


using mpfr_float = bertini::mpfr_float;
using bertini::DefaultPrecision;

BOOST_AUTO_TEST_CASE(mpfr_complex_eigen_norm_of_vector)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	Eigen::Matrix<bertini::complex, Eigen::Dynamic, Eigen::Dynamic> A(1,3);
	A << bertini::complex("1.0","1.0"), bertini::complex("1.0","1.0"), bertini::complex("1.0","1.0");
	bertini::mpfr_float n = A.norm();

	BOOST_CHECK(abs(n/sqrt(bertini::mpfr_float("6"))-mpfr_float("1"))<threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(mpfr_complex_eigen_negative_of_vector)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	Eigen::Matrix<bertini::complex, Eigen::Dynamic, Eigen::Dynamic> A(1,3);
	A << bertini::complex("1.0","1.0"), bertini::complex("1.0","1.0"), bertini::complex("1.0","1.0");
	auto B = -A;

	BOOST_CHECK_EQUAL(B(0), -A(0));
	BOOST_CHECK_EQUAL(B(1), -A(1));
	BOOST_CHECK_EQUAL(B(2), -A(2));
}


BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(miscellaneous_complex_tests_out_of_place)

using mpfr_float = bertini::mpfr_float;
using bertini::DefaultPrecision;

BOOST_AUTO_TEST_CASE(mpfr_float_serialization)
{
	
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	
	mpfr_float q = acos( mpfr_float("-1.0") );
	
	
	{
		std::ofstream fout("serialization_test_mpfr");
		
		boost::archive::text_oarchive oa(fout);
		
		// write class instance to archive
		oa << q;
	}
	
	mpfr_float w;
	{
		std::ifstream fin("serialization_test_mpfr");
		
		boost::archive::text_iarchive ia(fin);
		// read class state from archive
		ia >> w;
	}
	
	BOOST_CHECK_EQUAL(q,w);
	
}


BOOST_AUTO_TEST_CASE(mpfr_float_serialization2)
{
	using mpfr_float = bertini::mpfr_float;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	
	mpfr_float q = acos( mpfr_float("-1.0") );
	
	
	{
		std::ofstream fout("serialization_test_mpfr");
		
		boost::archive::text_oarchive oa(fout);
		
		// write class instance to archive
		oa << q;
	}
	
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS*2);

	mpfr_float w;
	{
		std::ifstream fin("serialization_test_mpfr");
		
		boost::archive::text_iarchive ia(fin);
		// read class state from archive
		ia >> w;
	}
	
	BOOST_CHECK_EQUAL(q,w);
	
}

BOOST_AUTO_TEST_SUITE_END()

