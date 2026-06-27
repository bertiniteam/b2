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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire



#include "bertini2/num_traits.hpp"
#include <boost/test/unit_test.hpp>
#include <fstream>


#include "bertini2/eigen_extensions.hpp"
#include <Eigen/Dense>

#include "externs.hpp"


BOOST_AUTO_TEST_SUITE(complex_multiprecision_class)

using bertini::Precision;
using real_mp = bertini::real_mp;
using bertini::DefaultPrecision;

BOOST_AUTO_TEST_CASE(complex_create_default_constructor)
{
	bertini::complex_mp z;
}


BOOST_AUTO_TEST_CASE(complex_create_two_input_constructor)
{
	bertini::complex_mp z("0.1","1.2");
}


BOOST_AUTO_TEST_CASE(complex_addition)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	bertini::complex_mp z("0.1","1.2");
	bertini::complex_mp v("0.2","1.3");
	
	bertini::complex_mp r = z+v;
	BOOST_CHECK(abs(r.real()-bertini::real_mp("0.3")) < threshold_clearance_mp);
	BOOST_CHECK(abs(r.imag()-bertini::real_mp("2.5")) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(complex_subtraction)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	bertini::complex_mp z("0.1","1.2");
	bertini::complex_mp v("0.2","1.3");
	
	bertini::complex_mp r = z-v;
	BOOST_CHECK(abs(r.real()-bertini::real_mp("-0.1")) < threshold_clearance_mp);
	BOOST_CHECK(abs(r.imag()-bertini::real_mp("-0.1")) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(complex_negation)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	bertini::complex_mp z("0.1","1.2");

	bertini::complex_mp r = -z;

	BOOST_CHECK_EQUAL(r.real(),-bertini::real_mp("0.1"));
	BOOST_CHECK_EQUAL(r.imag(),-bertini::real_mp("1.2"));
}

BOOST_AUTO_TEST_CASE(complex_multiplication)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex_mp z("0.1000","1.2000");
	bertini::complex_mp v("0.2000","1.3000");
	bertini::complex_mp r = z*v;
	
	BOOST_CHECK(abs(r.real()-bertini::real_mp("-1.54")) < threshold_clearance_mp);
	BOOST_CHECK(abs(r.imag()-bertini::real_mp("0.37"))< threshold_clearance_mp);
	
}


BOOST_AUTO_TEST_CASE(complex_division)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex_mp z("1.5","2.25");
	bertini::complex_mp v("-3.1","5.1");
	bertini::complex_mp r = z/v;
	
	BOOST_CHECK(abs(r.real()-bertini::real_mp("0.191605839416058394160583941605839416058394160583941605839416")) < threshold_clearance_mp);
	BOOST_CHECK(abs(r.imag()-bertini::real_mp("-0.410583941605839416058394160583941605839416058394160583941606"))< threshold_clearance_mp);
	
}



// commented out as boost::multiprecision::mpc doesn't offer `inverse`

// BOOST_AUTO_TEST_CASE(complex_inverse)
// {
// 	using real_mp = bertini::real_mp;
// 	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
// 	bertini::complex_mp z("1.5", "2.25");
// 	bertini::complex_mp w = inverse(z);
	
// 	BOOST_CHECK( abs(real(w) - real_mp("0.205128205128205128205128205128205128205128205128205128205128")) < threshold_clearance_mp);
// 	// this value computed with matlab's vpa.
// 	BOOST_CHECK( abs(imag(w) - real_mp("-0.307692307692307692307692307692307692307692307692307692307692")) < threshold_clearance_mp);
// 	// this value computed with matlab's vpa.
	
// }





BOOST_AUTO_TEST_CASE(complex_sqrt)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = sqrt(z);
	
	BOOST_CHECK(abs(real(w)- real_mp("1.44985576120488481677238036203436657121811982450524518214799")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("0.775939255547105301187773883610045740361564192945364951422177")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}




BOOST_AUTO_TEST_CASE(complex_log)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = log(z);
	
	BOOST_CHECK(abs(real(w)- real_mp("0.994792606278987440587524714788831870899124261482342501702857")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("0.982793723247329067985710611014666014496877453631628556761425")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}


BOOST_AUTO_TEST_CASE(complex_pow_expressionofreal)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	real_mp v("-3.1");
	bertini::complex_mp w = pow(z,2*v);
		
	//	0.0020583559161721337087583598028377697763165142184117220013738 + 0.000395572868694696378274143631981503385520868787960481823572936
	BOOST_CHECK(abs(real(w)- real_mp("0.0020583559161721337087583598028377697763165142184117220013738")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("0.000395572868694696378274143631981503385520868787960481823572936")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
}



BOOST_AUTO_TEST_CASE(complex_exponential)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25"), v("-3.1", "5.1");
	bertini::complex_mp w = pow(z,v);
	
	BOOST_CHECK(abs(real(w)- real_mp("-0.000134184252069350633500916903057767602709221538143686248955454")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("0.000273589335652288672638569400353054443664478211013806804752627")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}



/////////////////trig


BOOST_AUTO_TEST_CASE(complex_sin)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = sin(z);
	
	BOOST_CHECK(abs(real(w)- real_mp("4.78455206454183468376293313858001326644589827567133864164966")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("0.331840118511466433278036564370915664859189055829310017561679")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}

BOOST_AUTO_TEST_CASE(complex_cos)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = cos(z);
	
	
	BOOST_CHECK(abs(real(w)- real_mp("0.339295764714918536764798162132262895147653901636178423640742")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("-4.67941686644942009239504224985289863706184473303907201576996")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}


BOOST_AUTO_TEST_CASE(complex_tan)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = tan(z);
	
	BOOST_CHECK(abs(real(w)- real_mp("0.0032055151478664658165188420647201618843131312777936688172756")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("1.02223515949711949776179466571527189987176458555122839197245")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}











//////////hyperbolic trig

BOOST_AUTO_TEST_CASE(complex_sinh)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = sinh(z);
	
	BOOST_CHECK(abs(real(w)- real_mp("-1.33755718909601135383265663229426431626814744207128061926332")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("1.83034686972219823005238120890074929109631502976918245359956")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}

BOOST_AUTO_TEST_CASE(complex_cosh)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = cosh(z);
	
	
	BOOST_CHECK(abs(real(w)- real_mp("-1.47772167013515546574129058780587714414840064875386904640065")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("1.65673527269339558987730878420158930664270667966911313867276")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}


BOOST_AUTO_TEST_CASE(complex_tanh)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = tanh(z);
	
	BOOST_CHECK(abs(real(w)- real_mp("1.01633467755934300564918114729269683730547835877531900704421")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("-0.0991725055603753242720729815306256466293467541150851175179255")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}








////////////////////// arc trig functions


BOOST_AUTO_TEST_CASE(complex_asin)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = asin(z);
	
	BOOST_CHECK(abs(real(w)- real_mp("0.557728098701882481783572689617979376112207607697432876677129")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("1.70220126917599043453132270948063371483983368727808240752544")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}

BOOST_AUTO_TEST_CASE(complex_acos)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = acos(z);
	
	BOOST_CHECK(abs(real(w)- real_mp("1.01306822809301413744774900202177206598637709199012003381034")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("-1.70220126917599043453132270948063371483983368727808240752544")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}


BOOST_AUTO_TEST_CASE(complex_atan)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = atan(z);

	BOOST_CHECK(abs(real(w)- real_mp("1.34897118928106882765457017435994248505893277188774702053495")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("0.303034028741274232429033400709495120396582520146639573148632")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}




////////////////////// arc hyperbolic trig functions


BOOST_AUTO_TEST_CASE(complex_asinh)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = asinh(z);
	
	BOOST_CHECK(abs(real(w)- real_mp("1.67614738083746762411099043767860402156031800041110901006856")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("0.950053155715849776576329694197937868331910565261475162022272")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}

BOOST_AUTO_TEST_CASE(complex_acosh)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = acosh(z);
	
	BOOST_CHECK(abs(real(w)- real_mp("1.70220126917599043453132270948063371483983368727808240752544")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("1.01306822809301413744774900202177206598637709199012003381034")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}


BOOST_AUTO_TEST_CASE(complex_atanh)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.5", "2.25");
	bertini::complex_mp w = atanh(z);
	
	BOOST_CHECK(abs(real(w)- real_mp("0.188961443693877322997198271883423688031632018662054076093211")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	BOOST_CHECK(abs(imag(w)- real_mp("1.26114018722767258645714375731198514259936666747537259674077")) < threshold_clearance_mp);
	// this value computed with matlab's vpa.
	
}








BOOST_AUTO_TEST_CASE(complex_absolute_value)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	BOOST_CHECK_EQUAL(abs(bertini::complex_mp("1.0","1.0")), sqrt(real_mp("2.0")));
	
}

BOOST_AUTO_TEST_CASE(complex_abs_of_i_is_one)
{
	bertini::complex_mp i("0.0","1.0");
	bertini::real_mp one = abs(i);
	BOOST_CHECK_EQUAL(one,bertini::real_mp("1.0"));
}




BOOST_AUTO_TEST_CASE(complex_eleventh_power)
{
	using real_mp = bertini::real_mp;
	bertini::complex_mp z("2.0","0.5"), w;
	
	w = pow(z,11);
	BOOST_CHECK_EQUAL(real(w), real_mp("-1319867") / real_mp("512") );
	BOOST_CHECK_EQUAL(imag(w),  real_mp("2529647") / real_mp("2048"));
	
	
}




BOOST_AUTO_TEST_CASE(complex_conjugation)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("2.0", sqrt(real_mp("2.0")));
	bertini::complex_mp w = conj(z);
	
	BOOST_CHECK_EQUAL(real(w),real_mp("2.0"));
	BOOST_CHECK_EQUAL(imag(w),-sqrt(real_mp("2.0")));
}






BOOST_AUTO_TEST_CASE(complex_argument)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	bertini::complex_mp z("1.0","1.0");
	
	BOOST_CHECK_EQUAL(arg(z), acos(real_mp("-1.0"))/real_mp("4.0"));
}





BOOST_AUTO_TEST_CASE(complex_make_random_50)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);  
	
	bertini::complex_mp z = bertini::multiprecision::rand();

}



BOOST_AUTO_TEST_CASE(complex_make_random_100)
{
	DefaultPrecision(100);  
	bertini::complex_mp z = bertini::multiprecision::rand();
}


#ifndef B2_FORBID_MIXED_ARITHMETIC
BOOST_AUTO_TEST_CASE(interoperability_with_rational)
{
	complex_mp z;
	bertini::mpq_rational r;

	complex_mp w = z*r;
}
#endif


BOOST_AUTO_TEST_CASE(complex_serialization)
{
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	
	bertini::complex_mp z("1.23456", acos(real_mp("-1")));
	
	{
		std::ofstream fout("serialization_test_complex");
		
		boost::archive::text_oarchive oa(fout);
		
		// write class instance to archive
		oa << z;
	}
	
	//	bertini::complex_mp w;
	bertini::complex_mp w;
	{
		std::ifstream fin("serialization_test_complex");
		
		boost::archive::text_iarchive ia(fin);
		// read class state from archive
		ia >> w;
	}
	
	BOOST_CHECK_EQUAL(real(z),real(w));
	BOOST_CHECK_EQUAL(imag(z),imag(w));
	
}





BOOST_AUTO_TEST_CASE(complex_precision_predictable_add)
{

	bertini::scoped_mpfr_precision_options_this_thread precision_options(boost::multiprecision::variable_precision_options::preserve_source_precision);


	DefaultPrecision(30);
	bertini::complex_mp a(1,2);

	DefaultPrecision(50);
	bertini::complex_mp b(3,4);

	DefaultPrecision(70);
	bertini::complex_mp c(5,6);

	a = b+c;
	BOOST_CHECK_EQUAL(Precision(a),std::max(Precision(b),Precision(c)));

	DefaultPrecision(90);
	a = b+c;
	BOOST_CHECK_EQUAL(Precision(a),90);
}


BOOST_AUTO_TEST_CASE(complex_precision_predictable_sub)
{
	bertini::scoped_mpfr_precision_options_this_thread precision_options(boost::multiprecision::variable_precision_options::preserve_source_precision);


	DefaultPrecision(30);
	bertini::complex_mp a(1,2);
	BOOST_CHECK_EQUAL(Precision(a),30);

	DefaultPrecision(50);
	bertini::complex_mp b(3,4);
	BOOST_CHECK_EQUAL(Precision(b),50);

	DefaultPrecision(70);
	bertini::complex_mp c(5,6);
	BOOST_CHECK_EQUAL(Precision(c),70);

	a = b-c;
	BOOST_CHECK_EQUAL(Precision(a),std::max(Precision(b),Precision(c)));

	DefaultPrecision(90);
	a = b-c;
	BOOST_CHECK_EQUAL(Precision(a),90);
}


BOOST_AUTO_TEST_CASE(complex_precision_predictable_mul)
{
	bertini::scoped_mpfr_precision_options_this_thread precision_options(boost::multiprecision::variable_precision_options::preserve_source_precision);

	DefaultPrecision(30);
	bertini::complex_mp a(1,2);
	BOOST_CHECK_EQUAL(Precision(a),30);

	DefaultPrecision(50);
	bertini::complex_mp b(3,4);
	BOOST_CHECK_EQUAL(Precision(b),50);

	DefaultPrecision(70);
	bertini::complex_mp c(5,6);
	BOOST_CHECK_EQUAL(Precision(c),70);

	a = b*c;
	BOOST_CHECK_EQUAL(Precision(a),std::max(Precision(b),Precision(c)));

	DefaultPrecision(90);
	a = b*c;
	BOOST_CHECK_EQUAL(Precision(a),90);
}



BOOST_AUTO_TEST_CASE(complex_precision_predictable_div)
{
	bertini::scoped_mpfr_precision_options_this_thread precision_options(boost::multiprecision::variable_precision_options::preserve_source_precision);

	DefaultPrecision(30);
	bertini::complex_mp a(1,2);
	BOOST_CHECK_EQUAL(Precision(a),30);

	DefaultPrecision(50);
	bertini::complex_mp b(3,4);
	BOOST_CHECK_EQUAL(Precision(b),50);

	DefaultPrecision(70);
	bertini::complex_mp c(5,6);
	BOOST_CHECK_EQUAL(Precision(c),70);

	a = b/c;
	BOOST_CHECK_EQUAL(Precision(a),std::max(Precision(b),Precision(c)));

	DefaultPrecision(90);
	a = b/c;
	BOOST_CHECK_EQUAL(Precision(a),90);
}


BOOST_AUTO_TEST_CASE(complex_precision_predictable_pow)
{
	bertini::scoped_mpfr_precision_options_this_thread precision_options(boost::multiprecision::variable_precision_options::preserve_source_precision);

	DefaultPrecision(30);
	bertini::complex_mp a(1,2);
	BOOST_CHECK_EQUAL(Precision(a),30);

	DefaultPrecision(50);
	bertini::complex_mp b(3,4);
	BOOST_CHECK_EQUAL(Precision(b),50);

	DefaultPrecision(70);
	bertini::complex_mp c(5,6);
	BOOST_CHECK_EQUAL(Precision(c),70);

	a = pow(b,c);
	BOOST_CHECK_EQUAL(Precision(a),std::max(Precision(b),Precision(c)));

	DefaultPrecision(90);
	a = pow(b,c);
	BOOST_CHECK_EQUAL(Precision(a),90);
}

BOOST_AUTO_TEST_CASE(complex_precision_predictable_trig)
{
	bertini::scoped_mpfr_precision_options_this_thread precision_options(boost::multiprecision::variable_precision_options::preserve_source_precision);

	DefaultPrecision(30);
	bertini::complex_mp a(1,2);
	BOOST_CHECK_EQUAL(Precision(a),30);

	DefaultPrecision(50);
	bertini::complex_mp b(3,4);
	BOOST_CHECK_EQUAL(Precision(b),50);

	DefaultPrecision(70);
	bertini::complex_mp c(5,6);
	BOOST_CHECK_EQUAL(Precision(c),DefaultPrecision());

	a = tan(b);
	BOOST_CHECK_EQUAL(Precision(a),DefaultPrecision());

	DefaultPrecision(90);
	a = tan(b);
	BOOST_CHECK_EQUAL(Precision(a),DefaultPrecision());
}


BOOST_AUTO_TEST_CASE(complex_precision_predictable_arg)
{
	bertini::scoped_mpfr_precision_options_this_thread precision_options(boost::multiprecision::variable_precision_options::preserve_source_precision);

	DefaultPrecision(30);
	bertini::complex_mp a(1,2);
	BOOST_CHECK_EQUAL(Precision(a),30);

	DefaultPrecision(50);
	bertini::complex_mp b(3,4);
	BOOST_CHECK_EQUAL(Precision(b),50);

	DefaultPrecision(70);
	bertini::complex_mp c(5,6);
	BOOST_CHECK_EQUAL(Precision(c),70);

	a = arg(b);
	BOOST_CHECK_EQUAL(Precision(a),DefaultPrecision());

	DefaultPrecision(90);
	a = arg(b);
	BOOST_CHECK_EQUAL(Precision(a),DefaultPrecision());
}




BOOST_AUTO_TEST_CASE(complex_precision_predictable2)
{
	bertini::scoped_mpfr_precision_options_this_thread precision_options(boost::multiprecision::variable_precision_options::preserve_source_precision);

	complex_mp::default_precision(50);
	real_mp::default_precision(50);
	complex_mp a;
	BOOST_CHECK_EQUAL(a.real().precision(),50);
	complex_mp::default_precision(40);
	real_mp::default_precision(40);
	BOOST_CHECK_EQUAL(a.real().precision(),50);

	complex_mp::default_precision(50);
	real_mp::default_precision(50);
	complex_mp b;
	BOOST_CHECK_EQUAL(b.real().precision(),50);
	real_mp::default_precision(40);
	BOOST_CHECK_EQUAL(b.real().precision(),50);



	complex_mp::default_precision(50);
	real_mp::default_precision(50);
	complex_mp c;
	BOOST_CHECK_EQUAL(c.real().precision(),50);
	complex_mp::default_precision(40);
	BOOST_CHECK_EQUAL(c.real().precision(),50);
}


BOOST_AUTO_TEST_CASE(complex_precision_predictable3)
{
	bertini::scoped_mpfr_precision_options_this_thread precision_options(bertini::DefaultPrecisionPolicy());

	DefaultPrecision(50);
	complex_mp a;

	DefaultPrecision(40);
	complex_mp b(a);

	complex_mp c = a;

	BOOST_CHECK_EQUAL(c.precision(), 50);
	BOOST_CHECK_EQUAL(b.precision(), 50);
	BOOST_CHECK_EQUAL(a.precision(), 50);
}

BOOST_AUTO_TEST_CASE(complex_precision_move)
{
	bertini::scoped_mpfr_precision_options_this_thread precision_options(bertini::DefaultPrecisionPolicy());

	DefaultPrecision(50);
	complex_mp a;

	DefaultPrecision(100);
	complex_mp b;

	DefaultPrecision(50);
	a = std::move(b);

	BOOST_CHECK_EQUAL(a.precision(),100);
}

BOOST_AUTO_TEST_CASE(complex_precision_move2)
{
	bertini::scoped_mpfr_precision_options_this_thread precision_options(bertini::DefaultPrecisionPolicy());

	DefaultPrecision(50);
	complex_mp a;

	DefaultPrecision(100);
	complex_mp b;

	a = std::move(b);
	BOOST_CHECK_EQUAL(a.precision(),100);
}

BOOST_AUTO_TEST_CASE(real_precision_increase)
{
	bertini::scoped_mpfr_precision_options_this_thread precision_options(bertini::DefaultPrecisionPolicy());

	DefaultPrecision(50);
	real_mp a(2);
	a.precision(100);
	BOOST_CHECK(abs(real_mp(2)-a) < 1e-50);
}

BOOST_AUTO_TEST_CASE(complex_precision_increase)
{
	bertini::scoped_mpfr_precision_options_this_thread precision_options(bertini::DefaultPrecisionPolicy());

	DefaultPrecision(50);
	complex_mp a(2,3);
	a.precision(100);
	BOOST_CHECK(abs(complex_mp(2,3)-a) < 1e-50);
}

BOOST_AUTO_TEST_CASE(construct_from_nondefault_reals)
{
	bertini::scoped_mpfr_precision_options_this_thread precision_options(bertini::DefaultPrecisionPolicy());

	DefaultPrecision(100);
	real_mp a(1);
	real_mp b(2);


	DefaultPrecision(50);

	complex_mp z(a,b);

	BOOST_CHECK_EQUAL(a.precision(), 100);
	BOOST_CHECK_EQUAL(b.precision(), 100);
	BOOST_CHECK_EQUAL(z.precision(), 100);
}


BOOST_AUTO_TEST_CASE(complex_get_from_stream_parens_with_comma)
{
	std::stringstream ss;
	ss << "(0.1234,-4.12397)";
	
	bertini::complex_mp z;
	ss >> z;
	
	BOOST_CHECK_EQUAL(z.real(),bertini::real_mp("0.1234"));
	BOOST_CHECK_EQUAL(z.imag(),bertini::real_mp("-4.12397"));
}


BOOST_AUTO_TEST_CASE(complex_get_from_stream_parens_no_comma)
{
	std::stringstream ss;
	ss << "(-3.651263418976498712e-2)";
	
	bertini::complex_mp z;
	ss >> z;
	
	BOOST_CHECK_EQUAL(z.real(),bertini::real_mp("-3.651263418976498712e-2"));
	BOOST_CHECK_EQUAL(z.imag(),bertini::real_mp("0"));
}


BOOST_AUTO_TEST_CASE(complex_get_from_stream_no_parens)
{
	std::stringstream ss;
	ss << "-3.651263418976498712e-2";
	
	bertini::complex_mp z;
	ss >> z;
	
	BOOST_CHECK_EQUAL(z.real(),bertini::real_mp("-3.651263418976498712e-2"));
	BOOST_CHECK_EQUAL(z.imag(),bertini::real_mp("0"));
}



BOOST_AUTO_TEST_CASE(float_get_from_stream_parens_with_comma)
{
	std::stringstream ss;
	ss << "0.1234";
	
	bertini::real_mp z;
	ss >> z;
	
	BOOST_CHECK_EQUAL(z,bertini::real_mp("0.1234"));
}






BOOST_AUTO_TEST_CASE(precision_of_mpfr_is_16)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex_mp a("1.23124");
	BOOST_CHECK_EQUAL(bertini::Precision(a), CLASS_TEST_MPFR_DEFAULT_DIGITS);
}

BOOST_AUTO_TEST_CASE(precision_of_mpfr_complex_is_16)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex_mp a("1.23124","-0.6789124678912394");
	BOOST_CHECK_EQUAL(bertini::Precision(a), CLASS_TEST_MPFR_DEFAULT_DIGITS);
}


BOOST_AUTO_TEST_CASE(precision_random_real_default)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	auto a = bertini::multiprecision::RandomReal();
	BOOST_CHECK_EQUAL(bertini::Precision(a), CLASS_TEST_MPFR_DEFAULT_DIGITS);
	BOOST_CHECK_EQUAL(CLASS_TEST_MPFR_DEFAULT_DIGITS, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(precision_random_real_highest)
{	using namespace bertini;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex_mp a;
	bertini::multiprecision::RandomRealAssign (a,1000);
	BOOST_CHECK_EQUAL(bertini::Precision(a), 1000);
	BOOST_CHECK_EQUAL(CLASS_TEST_MPFR_DEFAULT_DIGITS, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(precision_random_unit_default)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	auto a = bertini::multiprecision::RandomUnit();
	BOOST_CHECK_EQUAL(bertini::Precision(a), CLASS_TEST_MPFR_DEFAULT_DIGITS);
	BOOST_CHECK_EQUAL(CLASS_TEST_MPFR_DEFAULT_DIGITS, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(precision_random_unit_highest)
{	using namespace bertini;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex_mp a;
	bertini::multiprecision::RandomUnitAssign(a,1000);
	BOOST_CHECK_EQUAL(bertini::Precision(a), 1000);
	BOOST_CHECK_EQUAL(CLASS_TEST_MPFR_DEFAULT_DIGITS, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(precision_random_default)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	auto a = bertini::multiprecision::rand();
	BOOST_CHECK_EQUAL(bertini::Precision(a), CLASS_TEST_MPFR_DEFAULT_DIGITS);
	BOOST_CHECK_EQUAL(CLASS_TEST_MPFR_DEFAULT_DIGITS, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(precision_random_highest)
{	using namespace bertini;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	bertini::complex_mp a;
	bertini::multiprecision::rand_assign(a,1000);
	BOOST_CHECK_EQUAL(Precision(a), 1000);
	BOOST_CHECK_EQUAL(CLASS_TEST_MPFR_DEFAULT_DIGITS, DefaultPrecision());
	BOOST_CHECK_EQUAL(a,a);
}


BOOST_AUTO_TEST_CASE(precision_through_swap)
{	using namespace bertini;
	DefaultPrecision(50);
	complex_mp a;
	DefaultPrecision(40);
	complex_mp b;

	swap(a,b);
	BOOST_CHECK_EQUAL(a.precision(), 40);
	BOOST_CHECK_EQUAL(b.precision(), 50);
}

BOOST_AUTO_TEST_CASE(precision_through_swap2)
{	using namespace bertini;
	DefaultPrecision(50);
	complex_mp a;
	DefaultPrecision(40);
	complex_mp b;

	a.swap(b);
	BOOST_CHECK_EQUAL(a.precision(), 40);
	BOOST_CHECK_EQUAL(b.precision(), 50);
}


BOOST_AUTO_TEST_CASE(precision_equality)
{	using namespace bertini;
	DefaultPrecision(1000);
	complex_mp a("1.324719827394086120398419082734980126734089612309871092830981236748901273498071240986123094861246981263481263489016238947147129807419028748901273409127349087124612576129076541203975704195690418570914657910465091256016501650916509165097164509164509761409561097561097650791650971465097165097162059761209561029756019265019726509126509172650971625097162450971309756104975610274650917825018740981274098127409182375701465172340923847120836540491320467127043127893281461230951097260126309812374091265091824981231236409851274",
		"-0.80743891267394610982659071452346156102764312401571972642394120395608291471029347812645125986123123904123471209381289471230512983491286102875870192091283712396550981723409812740981263471230498715096104897123094710923879065981740928740981271801391209238470129560941870129387409812883437894183883841283700483832883218128438938184289148239164079329657861209381892037483468937489237419236509823723705612893489712412306531274812364980127304981648712483248732");

	BOOST_CHECK_EQUAL(a,a);
}

BOOST_AUTO_TEST_CASE(precision_equality_default_differs)
{	using namespace bertini;
	DefaultPrecision(1000);
	complex_mp a("1.324719827394086120398419082734980126734089612309871092830981236748901273498071240986123094861246981263481263489016238947147129807419028748901273409127349087124612576129076541203975704195690418570914657910465091256016501650916509165097164509164509761409561097561097650791650971465097165097162059761209561029756019265019726509126509172650971625097162450971309756104975610274650917825018740981274098127409182375701465172340923847120836540491320467127043127893281461230951097260126309812374091265091824981231236409851274",
		"-0.80743891267394610982659071452346156102764312401571972642394120395608291471029347812645125986123123904123471209381289471230512983491286102875870192091283712396550981723409812740981263471230498715096104897123094710923879065981740928740981271801391209238470129560941870129387409812883437894183883841283700483832883218128438938184289148239164079329657861209381892037483468937489237419236509823723705612893489712412306531274812364980127304981648712483248732");
	DefaultPrecision(40);

	BOOST_CHECK_EQUAL(a,a);
}


// Regression coverage for the macos-15-intel SIGABRT inside System.eval.
// Boost.Multiprecision >= 1.87 default-constructs mpfr temporaries at
// thread_default_precision; if that thread-local is 0, mpfr_init2 aborts.
// These tests pin down our contract: bertini::DefaultPrecision keeps
// the static and thread-local defaults aligned, and default construction
// yields useful precision (not 0).

BOOST_AUTO_TEST_CASE(default_precision_aligns_thread_local_with_static)
{
	using namespace bertini;
	DefaultPrecision(57);
	BOOST_CHECK_EQUAL(real_mp::default_precision(),         57u);
	BOOST_CHECK_EQUAL(complex_mp::default_precision(),       57u);
	BOOST_CHECK_EQUAL(real_mp::thread_default_precision(),  57u);
	BOOST_CHECK_EQUAL(complex_mp::thread_default_precision(),57u);

	// restore so subsequent tests aren't surprised
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
}


BOOST_AUTO_TEST_CASE(default_constructed_mpfr_float_inherits_default_precision)
{
	using namespace bertini;
	DefaultPrecision(63);
	real_mp f;                            // no explicit precision
	BOOST_CHECK_EQUAL(f.precision(), 63u);   // must not be 0

	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
}


BOOST_AUTO_TEST_CASE(default_constructed_mpfr_complex_inherits_default_precision)
{
	using namespace bertini;
	DefaultPrecision(71);
	complex_mp z;                          // no explicit precision
	BOOST_CHECK_EQUAL(z.precision(), 71u);   // must not be 0

	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
}


BOOST_AUTO_TEST_CASE(default_precision_recovers_thread_local_from_zero)
{
	using namespace bertini;
	// Simulate a freshly-spawned thread on Boost 1.87: thread_default is 0.
	real_mp::thread_default_precision(0);
	complex_mp::thread_default_precision(0);

	// The fix: DefaultPrecision must repair both, not just the static side.
	// Without this, the next default-constructed mpfr temporary aborts inside
	// mpfr_init2 (this is what surfaced on macos-15-intel via boost::python).
	DefaultPrecision(50);
	BOOST_CHECK_EQUAL(real_mp::thread_default_precision(),  50u);
	BOOST_CHECK_EQUAL(complex_mp::thread_default_precision(),50u);

	real_mp  f;
	complex_mp z;
	BOOST_CHECK_EQUAL(f.precision(), 50u);
	BOOST_CHECK_EQUAL(z.precision(), 50u);

	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
}


BOOST_AUTO_TEST_SUITE_END()









BOOST_AUTO_TEST_SUITE(complex_eigen_compatibility)


using real_mp = bertini::real_mp;
using bertini::DefaultPrecision;

BOOST_AUTO_TEST_CASE(mpfr_complex_eigen_norm_of_vector)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	Eigen::Matrix<bertini::complex_mp, Eigen::Dynamic, Eigen::Dynamic> A(1,3);
	A << bertini::complex_mp("1.0","1.0"), bertini::complex_mp("1.0","1.0"), bertini::complex_mp("1.0","1.0");
	bertini::real_mp n = A.norm();

	BOOST_CHECK(abs(n/sqrt(bertini::real_mp("6"))-real_mp("1"))<threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(mpfr_complex_eigen_negative_of_vector)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	Eigen::Matrix<bertini::complex_mp, Eigen::Dynamic, Eigen::Dynamic> A(1,3);
	A << bertini::complex_mp("1.0","1.0"), bertini::complex_mp("1.0","1.0"), bertini::complex_mp("1.0","1.0");
	auto B = -A;

	BOOST_CHECK_EQUAL(B(0), -A(0));
	BOOST_CHECK_EQUAL(B(1), -A(1));
	BOOST_CHECK_EQUAL(B(2), -A(2));
}


BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(miscellaneous_complex_tests_out_of_place)

using real_mp = bertini::real_mp;
using bertini::DefaultPrecision;

BOOST_AUTO_TEST_CASE(mpfr_float_serialization)
{
	
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	
	real_mp q = acos( real_mp("-1.0") );
	
	
	{
		std::ofstream fout("serialization_test_mpfr");
		
		boost::archive::text_oarchive oa(fout);
		
		// write class instance to archive
		oa << q;
	}
	
	real_mp w;
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
	using real_mp = bertini::real_mp;
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	
	real_mp q = acos( real_mp("-1.0") );
	
	
	{
		std::ofstream fout("serialization_test_mpfr");
		
		boost::archive::text_oarchive oa(fout);
		
		// write class instance to archive
		oa << q;
	}
	
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS*2);

	real_mp w;
	{
		std::ifstream fin("serialization_test_mpfr");
		
		boost::archive::text_iarchive ia(fin);
		// read class state from archive
		ia >> w;
	}
	
	BOOST_CHECK_EQUAL(q,w);
	
}

BOOST_AUTO_TEST_SUITE_END()

