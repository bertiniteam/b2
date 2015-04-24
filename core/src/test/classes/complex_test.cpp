
#include "complex.hpp"
#include <boost/test/unit_test.hpp>

#define private public

BOOST_AUTO_TEST_SUITE(complex_multiprecision_class)

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
	bertini::complex z("0.1","1.2");
	bertini::complex v("0.2","1.3");
	
	bertini::complex r = z+v;
	BOOST_CHECK_EQUAL(r.real(),boost::multiprecision::mpfr_float("0.3"));
	BOOST_CHECK_EQUAL(r.imag(),boost::multiprecision::mpfr_float("2.5"));
}



BOOST_AUTO_TEST_CASE(complex_multiplication)
{
	boost::multiprecision::mpfr_float::default_precision(50);
	bertini::complex z("0.1000","1.2000");
	bertini::complex v("0.2000","1.3000");
	bertini::complex r = z*v;
	
	BOOST_CHECK(abs(r.real()-boost::multiprecision::mpfr_float("-1.54")) < 1e-50);
	BOOST_CHECK(abs(r.imag()-boost::multiprecision::mpfr_float("0.37"))< 1e-50);
	
}






BOOST_AUTO_TEST_CASE(complex_abs_of_i_is_one)
{
	bertini::complex i("0.0","1.0");
	boost::multiprecision::mpfr_float one = abs(i);
	std::cout << one << std::endl;
	BOOST_CHECK_EQUAL(one,boost::multiprecision::mpfr_float("1.0"));
}



BOOST_AUTO_TEST_CASE(complex_write_to_screen)
{
	bertini::complex i("0.0","1.0");
	
	std::cout << i.real() << " " << i.imag() << std::endl;
	std::cout << i << std::endl;
}








BOOST_AUTO_TEST_CASE(complex_get_from_stream_parens_with_comma)
{
	std::stringstream ss;
	ss << "(0.1234,-4.12397)";
	
	bertini::complex z;
	ss >> z;
	
	BOOST_CHECK_EQUAL(z.real(),boost::multiprecision::mpfr_float("0.1234"));
	BOOST_CHECK_EQUAL(z.imag(),boost::multiprecision::mpfr_float("-4.12397"));
}


BOOST_AUTO_TEST_CASE(complex_get_from_stream_parens_no_comma)
{
	std::stringstream ss;
	ss << "(-3.651263418976498712e-2)";
	
	bertini::complex z;
	ss >> z;
	
	BOOST_CHECK_EQUAL(z.real(),boost::multiprecision::mpfr_float("-3.651263418976498712e-2"));
	BOOST_CHECK_EQUAL(z.imag(),boost::multiprecision::mpfr_float("0"));
}


BOOST_AUTO_TEST_CASE(complex_get_from_stream_no_parens)
{
	std::stringstream ss;
	ss << "-3.651263418976498712e-2";
	
	bertini::complex z;
	ss >> z;
	
	BOOST_CHECK_EQUAL(z.real(),boost::multiprecision::mpfr_float("-3.651263418976498712e-2"));
	BOOST_CHECK_EQUAL(z.imag(),boost::multiprecision::mpfr_float("0"));
}



BOOST_AUTO_TEST_CASE(complex_eleventh_power)
{

	bertini::complex z("2.0","0.5"), w;
	
	w = pow(z,11);
	std::cout << "z^11 = " << w << std::endl;
	
}




BOOST_AUTO_TEST_SUITE_END()


