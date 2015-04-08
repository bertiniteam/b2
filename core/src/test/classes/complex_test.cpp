
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
	bertini::complex z("0.1","1.2");
	bertini::complex v("0.2","1.3");
	
	bertini::complex r = z*v;
	BOOST_CHECK_EQUAL(r.real(),boost::multiprecision::mpfr_float("-1.540"));
	BOOST_CHECK_EQUAL(r.imag(),boost::multiprecision::mpfr_float("0.370"));
}






BOOST_AUTO_TEST_CASE(complex_abs_of_i_is_one)
{
	bertini::complex i("0.0","1.0");
	boost::multiprecision::mpfr_float one = abs(i);
	std::cout << one << std::endl;
	BOOST_CHECK_EQUAL(one,boost::multiprecision::mpfr_float("1.0"));
}



BOOST_AUTO_TEST_SUITE_END()


