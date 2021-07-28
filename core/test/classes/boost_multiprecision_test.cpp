#include <boost/test/unit_test.hpp>
#include <boost/multiprecision/mpc.hpp>

using mpfr_float = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>, boost::multiprecision::et_on>;
using mpz_int = boost::multiprecision::number<boost::multiprecision::backends::gmp_int, boost::multiprecision::et_on>;
using mpq_rational = boost::multiprecision::number<boost::multiprecision::backends::gmp_rational, boost::multiprecision::et_on>;
using mpc_complex = boost::multiprecision::number<boost::multiprecision::backends::mpc_complex_backend<0>, boost::multiprecision::et_on>;

namespace utf = boost::unit_test;

// BOOST_AUTO_TEST_SUITE(boost_multiprecision_tutorial)


// using mp_t = boost::multiprecision::debug_adaptor_t<boost::multiprecision::mpfr_float>;

// struct scoped_mpfr_precision
// {
// 	unsigned saved_digits10;
// 	scoped_mpfr_precision(unsigned digits10) : saved_digits10(mp_t::thread_default_precision())
// 	{
// 		mp_t::thread_default_precision(digits10);
// 	}

// 	~scoped_mpfr_precision()
// 	{
// 		mp_t::thread_default_precision(saved_digits10);
// 	}

// 	void reset(unsigned digits10)
// 	{
// 		mp_t::thread_default_precision(digits10);
// 	}

// 	void reset()
// 	{
// 		mp_t::thread_default_precision(saved_digits10);
// 	}
// };

// struct scoped_mpfr_precision_options
// {

// 	boost::multiprecision::variable_precision_options saved_options;

// 	scoped_mpfr_precision_options(boost::multiprecision::variable_precision_options opts) : saved_options(mp_t::thread_default_variable_precision_options())
// 	{
// 		mp_t::thread_default_variable_precision_options(opts);
// 	}

// 	~scoped_mpfr_precision_options()
// 	{
// 		mp_t::thread_default_variable_precision_options(saved_options);
// 	}

// 	void reset(boost::multiprecision::variable_precision_options opts)
// 	{
// 		mp_t::thread_default_variable_precision_options(opts);
// 	}
// };



// BOOST_AUTO_TEST_CASE(precision_2_to_the_1000_minus_1_preserve_source)
// {
// 	scoped_mpfr_precision scope_1(mp_t::thread_default_precision());
// 	scoped_mpfr_precision_options scope_2(boost::multiprecision::variable_precision_options::preserve_source_precision);


// 	// calculate 2^1000 - 1:
// 	mpz_int i(1);
// 	i = i << 1000;
// 	i -= 1;
// 	std::cout << i << std::endl;
// 	mp_t f = i;

// 	BOOST_CHECK(f.precision()>mp_t::thread_default_precision());
// }


// BOOST_AUTO_TEST_CASE(precision_2_to_the_1000_minus_1_preserve_target)
// {
// 	scoped_mpfr_precision scope_1(mp_t::thread_default_precision());
// 	scoped_mpfr_precision_options scope_2(boost::multiprecision::variable_precision_options::preserve_target_precision);


// 	// calculate 2^1000 - 1:
// 	mpz_int i(1);
// 	i = i << 1000;
// 	i -= 1;

// 	mp_t f(0);
// 	f = i;

// 	BOOST_CHECK_EQUAL(f.precision(),mp_t::default_precision());
// }


// BOOST_AUTO_TEST_CASE(precision_2_to_the_1000_minus_1_preserve_all)
// {
// 	scoped_mpfr_precision scope_1(mp_t::thread_default_precision());
// 	scoped_mpfr_precision_options scope_2(boost::multiprecision::variable_precision_options::preserve_all_precision);

// 	// calculate 2^1000 - 1:
// 	mpz_int i(1);
// 	i = i << 1000;
// 	i -= 1;

// 	mp_t f = i;

// 	BOOST_CHECK(f.precision()>mp_t::thread_default_precision());
// }

// BOOST_AUTO_TEST_SUITE_END() // boost_multiprecision_tutorial












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



BOOST_AUTO_TEST_CASE(precision_complex_longlong_div_set_float_16digits)
{
	mpfr_float::default_precision(16);

	mpc_complex b(1,0);

	unsigned long long d{13};

	mpc_complex result{b/d};
	BOOST_CHECK_EQUAL(result.precision(),16);
}

BOOST_AUTO_TEST_CASE(precision_complex_longlong_div_set_complex_16digits, *utf::depends_on("boost_multiprecision/precision_complex_longlong_div_set_float_16digits"))
{
	mpc_complex::default_precision(16);

	mpc_complex b(1,0);

	unsigned long long d{13};

	mpc_complex result{b/d};
	BOOST_CHECK_EQUAL(result.precision(),16);
}



BOOST_AUTO_TEST_CASE(precision_complex_longlong_div_set_both_16digits, *utf::depends_on("boost_multiprecision/precision_complex_longlong_div_set_complex_16digits"))
{
	mpfr_float::default_precision(16);
	mpc_complex::default_precision(16);

	mpc_complex b(1,0);

	unsigned long long d{13};

	mpc_complex result{b/d};
	BOOST_CHECK_EQUAL(result.precision(),16);
}

BOOST_AUTO_TEST_CASE(precision_complex_longlong_div_set_float20_complex16, *utf::depends_on("boost_multiprecision/precision_complex_longlong_div_set_both_16digits"))
{
	mpfr_float::default_precision(20);
	mpc_complex::default_precision(16);

	mpc_complex b(1,0);

	unsigned long long d{13};

	mpc_complex result{b/d};
	BOOST_CHECK_EQUAL(result.precision(),16);
}




BOOST_AUTO_TEST_CASE(precision_complex_longlong_div_set_float16_complex20, *utf::depends_on("boost_multiprecision/precision_complex_longlong_div_set_float20_complex16"))
{
	mpfr_float::default_precision(16);
	mpc_complex::default_precision(20);

	mpc_complex b(1,0);

	unsigned long long d{13};

	mpc_complex result{b/d};
	BOOST_CHECK_EQUAL(result.precision(),20);
}



BOOST_AUTO_TEST_CASE(precision_complex_longlong_div_set_float20_complex20)
{
	mpfr_float::default_precision(20);
	mpc_complex::default_precision(20);

	mpc_complex b(1,0);

	unsigned long long d{13};

	mpc_complex result{b/d};
	BOOST_CHECK_EQUAL(result.precision(),20);
}



BOOST_AUTO_TEST_SUITE_END() // boost_multiprecision tests