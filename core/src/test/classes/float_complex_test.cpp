

#include <boost/test/unit_test.hpp>
#include <complex>



#include "float.hpp"


BOOST_AUTO_TEST_SUITE(use_float_in_std_complex)

BOOST_AUTO_TEST_CASE(Float_make_std_complex_number_i_squared_is_minus_one_high_precision)
{
	
	boost::multiprecision::mpfr_float::default_precision(1000);
	
	std::complex<bertini::Float> c(0,0),d(0,1),a;
	
	a = d*d;
	BOOST_CHECK(a.real()==bertini::Float(-1));
	BOOST_CHECK(a.imag()==bertini::Float(0));
}



BOOST_AUTO_TEST_CASE(Float_std_complex_number_abs_and_arg)
{

	boost::multiprecision::mpfr_float::default_precision(100);
	
	
	std::complex<bertini::Float> c(sqrt(bertini::Float(2.0)),bertini::Float(1.0)),d(0.0,1.0);
	
	std::cout.precision(100);
	std::cout << c << std::endl;
	std::cout << std::abs(c) << "\n" << sqrt(bertini::Float(2.0)) << std::endl;
	std::cout << abs(c)-sqrt(bertini::Float(2.0))  << std::endl;
	
	
	BOOST_CHECK( abs(abs(c)-sqrt(bertini::Float(2.0))) < bertini::Float("1e-999") );
	
	BOOST_CHECK_EQUAL(std::arg(c),bertini::Float::pi()/bertini::Float("4.0"));
	
	
	BOOST_CHECK_EQUAL(std::abs(d),bertini::Float(1));
	BOOST_CHECK_EQUAL(std::arg(d),bertini::Float::pi()/bertini::Float("2.0"));

	
}



BOOST_AUTO_TEST_SUITE_END()





