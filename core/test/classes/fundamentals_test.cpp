//This file is part of Bertini 2.
//
//fundamentals_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//fundamentals_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with fundamentals_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire


#include <boost/test/unit_test.hpp>
#include <boost/multiprecision/mpfr.hpp>

#include <iostream>
#include <vector>

#include "bertini2/double_extensions.hpp"
#include "bertini2/num_traits.hpp"









BOOST_AUTO_TEST_SUITE(super_fundamentals)

using real_mp = bertini::real_mp;
using mpq_rational = bertini::mpq_rational;
using complex_dbl = bertini::complex_dbl;
using bertini::DefaultPrecision;
#include <limits>











BOOST_AUTO_TEST_CASE(complex_pow)
{
	auto x = bertini::rand_complex();

	using bertini::pow;
	BOOST_CHECK_EQUAL(pow(x, 2), x*x);
}


BOOST_AUTO_TEST_CASE(complex_pow_on_real_stays_real)
{
	using bertini::pow;
	
	for (int ii=0; ii<100; ++ii)
	{
		auto x = complex_dbl(bertini::RandReal());
		auto result = pow(x, 2);
		BOOST_CHECK_EQUAL(result, x*x);
		BOOST_CHECK_EQUAL(imag(result), 0);
	}
}


BOOST_AUTO_TEST_CASE(complex_double_nans)
{
	complex_dbl x(0., std::numeric_limits<double>::quiet_NaN());
	using bertini::isnan;
	BOOST_CHECK(isnan(x));
}


BOOST_AUTO_TEST_CASE(complex_double_nans2)
{
	complex_dbl x(std::numeric_limits<double>::quiet_NaN(), 0.);
	using bertini::isnan;
	BOOST_CHECK(isnan(x));
}


BOOST_AUTO_TEST_CASE(complex_double_nans3)
{
	complex_dbl x(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
	using bertini::isnan;
	BOOST_CHECK(isnan(x));
}


BOOST_AUTO_TEST_CASE(complex_double_nans4)
{
	complex_dbl x(0., 0.);
	using bertini::isnan;
	BOOST_CHECK(!isnan(x));
}

BOOST_AUTO_TEST_CASE(complex_double_random_subsequent_not_equal)
{
	using T = complex_dbl;

	auto x = bertini::RandomUnit<T>();
	auto y = bertini::RandomUnit<T>();

	BOOST_CHECK(x!=y);
}


BOOST_AUTO_TEST_CASE(mpfr_float_can_be_nan)
{
	DefaultPrecision(50);
	real_mp a(std::numeric_limits<double>::quiet_NaN());
	BOOST_CHECK(isnan(a));
}


BOOST_AUTO_TEST_CASE(constructing_mpfr_from_double)
{
	DefaultPrecision(50);
	real_mp from_double(0.1);
	real_mp from_string("0.1");

	BOOST_CHECK(abs(from_string - from_double) < std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(construct_rational_from_integers)
{
	mpq_rational p(1,2);
	mpq_rational q(-1,2);

	BOOST_CHECK_EQUAL(p+q, mpq_rational(0));
}

// commented out code below checks construction of mpq_rationals from other types.  they're broken.  :(
// BOOST_AUTO_TEST_CASE(construct_rational_from_mpfr)
// {
// 	DefaultPrecision(50);
// 	real_mp p("1.1");
// 	mpq_rational q(p);
// 	BOOST_CHECK_EQUAL(mpq_rational(11,10), q);
// }


// BOOST_AUTO_TEST_CASE(construct_rational_from_string)
// {
// 	mpq_rational q("1.1");
// 	BOOST_CHECK_EQUAL(mpq_rational(11,10), q);
// }


// this commented out test will fail... because of mixed precision arithmetic
// BOOST_AUTO_TEST_CASE(multiple_mpfr_by_double)
// {
// 	DefaultPrecision(50);
// 	real_mp a("0.1");
// 	double factor = 0.1;

// 	real_mp result = a*factor;
// 	real_mp expected("0.01");

// 	BOOST_CHECK_CLOSE(expected, result, 1e-50);
// }



BOOST_AUTO_TEST_CASE(making_mpfr_from_pow_int_base)
{
	DefaultPrecision(50);

	real_mp result = pow(real_mp(10), -5);
	real_mp expected("1e-5");

	BOOST_CHECK_CLOSE(expected, result, 1e-50);
}

BOOST_AUTO_TEST_CASE(making_mpfr_from_pow_str_base)
{
	DefaultPrecision(50);

	real_mp result = pow(real_mp("10"), -5);
	real_mp expected("1e-5");

	BOOST_CHECK_CLOSE(expected, result, 1e-50);
}


BOOST_AUTO_TEST_CASE(making_mpfr_from_pow_doub_exp)
{
	DefaultPrecision(50);

using boost::multiprecision::pow;

	real_mp result = pow(real_mp(10), -5);
	real_mp expected("1e-5");

	BOOST_CHECK_CLOSE(expected, result, 1e-50);
}


BOOST_AUTO_TEST_CASE(making_mpfr_from_pow_int_base_mpfr_exp)
{
	DefaultPrecision(50);

	real_mp result = pow(10, real_mp(-5));
	real_mp expected("1e-5");

	BOOST_CHECK_CLOSE(expected, result, 1e-50);
}


BOOST_AUTO_TEST_CASE(make_rational_from_double)
{
	mpq_rational result(0.1);
	mpq_rational expected(1,10);
	BOOST_CHECK_CLOSE(result, expected, 1e-16);
}

BOOST_AUTO_TEST_CASE(make_random_mpfr_float_50)
{	
	using mpfr_50 = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<50>, boost::multiprecision::et_off>;

	using namespace boost::multiprecision;
	using namespace boost::random;
	DefaultPrecision(50);
	uniform_real_distribution<mpfr_50> ur(0,1);
	independent_bits_engine<mt19937, 50L*1000L/301L, mpz_int> gen;

	auto c = ur(gen);

	mpfr_50 r(c);
	mpfr_50 s(ur(gen));
}


BOOST_AUTO_TEST_CASE(make_random_mpfr_float_100)
{	
	using mpfr_100 = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<100>, boost::multiprecision::et_off>;

	using namespace boost::multiprecision;
	using namespace boost::random;
	DefaultPrecision(100);
	uniform_real_distribution<mpfr_100> ur(0,1);
	independent_bits_engine<mt19937, 100L*1000L/301L, mpz_int> gen;

	auto c = ur(gen);

	mpfr_100 r(c);
	mpfr_100 s(ur(gen));
}


BOOST_AUTO_TEST_CASE(RandomMP_default_precision_50)
{
	using namespace bertini;
	DefaultPrecision(50);
	auto a = RandomMp(real_mp(-1),real_mp(1));
	BOOST_CHECK_EQUAL(a.precision(), DefaultPrecision());
	BOOST_CHECK_EQUAL(50, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(RandomMP_default_precision_100)
{
	using namespace bertini;
	DefaultPrecision(100);
	auto a = RandomMp(real_mp(-1),real_mp(1));
	BOOST_CHECK_EQUAL(a.precision(), DefaultPrecision());
	BOOST_CHECK_EQUAL(100, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(RandomMP_nondefault_precision_100)
{
	using namespace bertini;
	DefaultPrecision(100);
	real_mp a;

	RandomMpAssign(a,500);

	BOOST_CHECK_EQUAL(a.precision(), 500);
	BOOST_CHECK_EQUAL(100, DefaultPrecision());
}

// The multiprecision RandomMp generator now draws from the single per-thread
// engine (ThreadEngine) shared by every random type, so SetGlobalSeed controls
// it.  Before this it used a private, unseedable independent_bits_engine.  These
// tests pin that contract: same seed -> identical mp draws; different seed ->
// different draws.
BOOST_AUTO_TEST_CASE(RandomMP_honors_global_seed)
{
	using namespace bertini;
	DefaultPrecision(50);

	auto draw_five = []{
		std::vector<real_mp> v;
		for (int ii = 0; ii < 5; ++ii)
			v.push_back(RandomMp(50));
		return v;
	};

	SetGlobalSeed(1234u);
	auto first = draw_five();

	SetGlobalSeed(1234u);
	auto second = draw_five();

	for (int ii = 0; ii < 5; ++ii)
		BOOST_CHECK_EQUAL(first[static_cast<size_t>(ii)], second[static_cast<size_t>(ii)]);

	SetGlobalSeed(5678u);
	auto third = draw_five();

	bool any_different = false;
	for (int ii = 0; ii < 5; ++ii)
		if (third[static_cast<size_t>(ii)] != first[static_cast<size_t>(ii)])
			any_different = true;
	BOOST_CHECK(any_different);
}

BOOST_AUTO_TEST_CASE(max_et_on)
{
	real_mp a(1), b(2), c(4);
	// auto d = max(a,b*b+c);

}




























struct scoped_mpfr_precision_options_all_threads
{
   boost::multiprecision::variable_precision_options saved_options;
   scoped_mpfr_precision_options_all_threads(boost::multiprecision::variable_precision_options opts) : saved_options(real_mp::default_variable_precision_options())
   {
      real_mp::default_variable_precision_options(opts);
   }
   ~scoped_mpfr_precision_options_all_threads()
   {
      real_mp::default_variable_precision_options(saved_options);
   }
   void reset(boost::multiprecision::variable_precision_options opts)
   {
      real_mp::default_variable_precision_options(opts);
   }
};







struct scoped_mpfr_precision_options_this_thread
{
   boost::multiprecision::variable_precision_options saved_options;
   scoped_mpfr_precision_options_this_thread(boost::multiprecision::variable_precision_options opts) : saved_options(real_mp::thread_default_variable_precision_options())
   {
      real_mp::thread_default_variable_precision_options(opts);
   }
   ~scoped_mpfr_precision_options_this_thread()
   {
      real_mp::thread_default_variable_precision_options(saved_options);
   }
   void reset(boost::multiprecision::variable_precision_options opts)
   {
      real_mp::thread_default_variable_precision_options(opts);
   }
};










BOOST_AUTO_TEST_CASE(precision_through_arithemetic)
{
	DefaultPrecision(50);


	scoped_mpfr_precision_options_this_thread scoped_opts1(boost::multiprecision::variable_precision_options::preserve_related_precision);

	real_mp x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);

// https://github.com/boostorg/multiprecision/issues/60
// 
// "Copying or move-assignment copies the precision of the source.
// Assignment keeps the precision of the target."

	DefaultPrecision(30);
	real_mp y = pow(x,2);
	BOOST_CHECK_EQUAL(y.precision(), 50);
	

	real_mp z = x;
	real_mp q(x);
	BOOST_CHECK_EQUAL(z.precision(), 50);
	BOOST_CHECK_EQUAL(q.precision(), 50);

	BOOST_CHECK(fabs(z - x) < 1e-50);


	DefaultPrecision(70);

	z = x;

	BOOST_CHECK_EQUAL(z.precision(),50);
	BOOST_CHECK(fabs(z - x) < 1e-50);


	y.precision(70);
	z.precision(30);

	BOOST_CHECK_EQUAL(y.precision(),70);
	BOOST_CHECK_EQUAL(z.precision(),30);
	BOOST_CHECK_EQUAL(x.precision(),50);



	scoped_mpfr_precision_options_this_thread scoped_opts2(boost::multiprecision::variable_precision_options::preserve_target_precision);

	y = z*x;


	BOOST_CHECK_EQUAL(z.precision(),30);
	BOOST_CHECK_EQUAL(x.precision(),50);
	BOOST_CHECK_EQUAL(y.precision(), 70);

}




BOOST_AUTO_TEST_CASE(precision_in_construction)
{
	DefaultPrecision(50);

	real_mp x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);

	DefaultPrecision(30);
	
	real_mp a = x;
	real_mp b(x);
	BOOST_CHECK_EQUAL(a.precision(), 50);
	BOOST_CHECK_EQUAL(b.precision(), 50);
	BOOST_CHECK_EQUAL(a,x);
	BOOST_CHECK_EQUAL(b,x);

	DefaultPrecision(70);

	real_mp c = x;
	real_mp d(x);
	BOOST_CHECK_EQUAL(c.precision(), 50);
	BOOST_CHECK_EQUAL(d.precision(), 50);
	BOOST_CHECK_EQUAL(c,x);
	BOOST_CHECK_EQUAL(d,x);
}


BOOST_AUTO_TEST_CASE(precision_through_arithemetic2)
{
	DefaultPrecision(50);
	real_mp a(1);

	DefaultPrecision(400);
	real_mp b(2);

	DefaultPrecision(600);
	real_mp c(3);

	a = b;
	BOOST_CHECK_EQUAL(a.precision(),400); // the precision of source

	a = b+c;
	BOOST_CHECK_EQUAL(a.precision(),600); // the bigger of 400,600 is 600
}



BOOST_AUTO_TEST_CASE(precision_mpfr_constructed_from_string)
{
	DefaultPrecision(30);
	real_mp x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(),30);
}


BOOST_AUTO_TEST_CASE(precision_of_double_is_16)
{
	double a(1.23124);
	BOOST_CHECK_EQUAL(bertini::Precision(a), 16);
}

BOOST_AUTO_TEST_CASE(precision_of_complex_double_is_16)
{
	std::complex<double> a(1.23124, -0.12345679);
	BOOST_CHECK_EQUAL(bertini::Precision(a), 16);
}


BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(numtraits)

using bertini::DefaultPrecision;

template<typename T>
using NumTraits = bertini::NumTraits<T>;

BOOST_AUTO_TEST_CASE(num_digits_double)
{
	using T = double;
	
	DefaultPrecision(16);
	BOOST_CHECK_EQUAL(NumTraits<T>::NumDigits(), 16);

	DefaultPrecision(30);
	BOOST_CHECK_EQUAL(NumTraits<T>::NumDigits(), 16);

	DefaultPrecision(100);
	BOOST_CHECK_EQUAL(NumTraits<T>::NumDigits(), 16);
}

BOOST_AUTO_TEST_CASE(num_digits_complex_double)
{
	using T = std::complex<double>;
	
	DefaultPrecision(16);
	BOOST_CHECK_EQUAL(NumTraits<T>::NumDigits(), 16);

	DefaultPrecision(30);
	BOOST_CHECK_EQUAL(NumTraits<T>::NumDigits(), 16);

	DefaultPrecision(100);
	BOOST_CHECK_EQUAL(NumTraits<T>::NumDigits(), 16);
}


BOOST_AUTO_TEST_CASE(num_digits_mpfr_float)
{
	using T = bertini::real_mp;
	
	DefaultPrecision(16);
	BOOST_CHECK_EQUAL(NumTraits<T>::NumDigits(), 16);

	DefaultPrecision(30);
	BOOST_CHECK_EQUAL(NumTraits<T>::NumDigits(), 30);

	DefaultPrecision(100);
	BOOST_CHECK_EQUAL(NumTraits<T>::NumDigits(), 100);
}


BOOST_AUTO_TEST_CASE(num_digits_mpfr_complex)
{
	using T = bertini::complex_mp;
	
	DefaultPrecision(16);
	BOOST_CHECK_EQUAL(NumTraits<T>::NumDigits(), 16);

	DefaultPrecision(30);
	BOOST_CHECK_EQUAL(NumTraits<T>::NumDigits(), 30);

	DefaultPrecision(100);
	BOOST_CHECK_EQUAL(NumTraits<T>::NumDigits(), 100);
}


BOOST_AUTO_TEST_SUITE_END() // numtraits tests
