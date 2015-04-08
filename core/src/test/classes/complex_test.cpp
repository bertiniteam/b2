
#include "complex.hpp"
#include <boost/test/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(complex_multiprecision_class)

BOOST_AUTO_TEST_CASE(complex_create_default_constructor)
{
	bertini::complex z;
}


BOOST_AUTO_TEST_CASE(complex_create_two_input_constructor)
{
	bertini::complex z("0.1","1.2");
}



BOOST_AUTO_TEST_SUITE_END()


