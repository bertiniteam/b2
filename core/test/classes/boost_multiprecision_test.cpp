#include <boost/test/unit_test.hpp>
#include <boost/multiprecision/mpc.hpp>

using mpfr_float = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>, boost::multiprecision::et_on>;
using mpz_int = boost::multiprecision::number<boost::multiprecision::backends::gmp_int, boost::multiprecision::et_on>;
using mpq_rational = boost::multiprecision::number<boost::multiprecision::backends::gmp_rational, boost::multiprecision::et_on>;
using mpc_complex = boost::multiprecision::number<boost::multiprecision::backends::mpc_complex_backend<0>, boost::multiprecision::et_on>;

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_SUITE(boost_multiprecision)

BOOST_AUTO_TEST_CASE(precision_complex_rational_add)
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = a+b;

	BOOST_CHECK_EQUAL(c.precision(),30);
}


BOOST_AUTO_TEST_CASE(precision_complex_rational_add_other_order, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = b+a;

	BOOST_CHECK_EQUAL(c.precision(),30);
}



BOOST_AUTO_TEST_CASE(precision_complex_rational_sub, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = a-b;

	BOOST_CHECK_EQUAL(c.precision(),30);
}


BOOST_AUTO_TEST_CASE(precision_complex_rational_sub_other_order, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = b-a;

	BOOST_CHECK_EQUAL(c.precision(),30);
}



BOOST_AUTO_TEST_CASE(precision_complex_rational_mul, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = a*b;

	BOOST_CHECK_EQUAL(c.precision(),30);
}

BOOST_AUTO_TEST_CASE(precision_complex_rational_mul_other_order, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = b*a;

	BOOST_CHECK_EQUAL(c.precision(),30);
}


BOOST_AUTO_TEST_CASE(precision_complex_rational_div, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = a/b;

	BOOST_CHECK_EQUAL(c.precision(),30);
}


BOOST_AUTO_TEST_CASE(precision_complex_rational_div_other_order, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = b/a;

	BOOST_CHECK_EQUAL(c.precision(),30);
}



BOOST_AUTO_TEST_SUITE_END() // numtraits tests