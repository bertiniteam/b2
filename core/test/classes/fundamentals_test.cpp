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
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame


#include <boost/test/unit_test.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>
#include <iostream>

#include "bertini2/limbo.hpp"
#include "bertini2/num_traits.hpp"

#include <boost/random.hpp>

using mpfr_float = bertini::mpfr_float;

using bertini::DefaultPrecision;

BOOST_AUTO_TEST_SUITE(super_fundamentals)



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
	auto a = RandomMp(mpfr_float(-1),mpfr_float(1));
	BOOST_CHECK_EQUAL(a.precision(), DefaultPrecision());
	BOOST_CHECK_EQUAL(50, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(RandomMP_default_precision_100)
{
	using namespace bertini;
	DefaultPrecision(100);
	auto a = RandomMp(mpfr_float(-1),mpfr_float(1));
	BOOST_CHECK_EQUAL(a.precision(), DefaultPrecision());
	BOOST_CHECK_EQUAL(100, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(RandomMP_nondefault_precision_100)
{
	using namespace bertini;
	DefaultPrecision(100);
	mpfr_float a;

	RandomMp(a,500);

	BOOST_CHECK_EQUAL(a.precision(), 500);
	BOOST_CHECK_EQUAL(100, DefaultPrecision());
}

BOOST_AUTO_TEST_CASE(max_et_on)
{
	mpfr_float a(1), b(2), c(4);
	// auto d = max(a,b*b+c);

}

BOOST_AUTO_TEST_CASE(precision_through_arithemetic)
{
	DefaultPrecision(50);

	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);

	DefaultPrecision(30);
	mpfr_float y = pow(x,2);
	BOOST_CHECK_EQUAL(y.precision(), 30);
	

	mpfr_float z = x;
	BOOST_CHECK_EQUAL(z.precision(), 30);

	BOOST_CHECK(fabs(z - mpfr_float("0.012345678901234567890123456789")) < 1e-30);



	DefaultPrecision(70);

	z = x;

	BOOST_CHECK_EQUAL(z.precision(),50);
	BOOST_CHECK(fabs(z - mpfr_float("0.01234567890123456789012345678901234567890123456789")) < 1e-50);


	y.precision(70);
	z.precision(30);

	BOOST_CHECK_EQUAL(y.precision(),70);
	BOOST_CHECK_EQUAL(z.precision(),30);
	BOOST_CHECK_EQUAL(x.precision(),50);

	y = z*x;

	// y is of precision 70 because it's a pre-existing variable.
	BOOST_CHECK_EQUAL(y.precision(), 70);
}

BOOST_AUTO_TEST_CASE(precision_mpfr_constructed_from_string)
{
	DefaultPrecision(30);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(),30);
}

	
BOOST_AUTO_TEST_CASE(index_and_subscript_generation1)
{

	std::vector<size_t> dimensions{2,2};
	std::vector<size_t> v;

	std::vector<size_t> solution{0,0};
	v = bertini::IndexToSubscript(0ul,dimensions);
	BOOST_CHECK(v==solution);

	solution[0] = 1; solution[1] = 0;
	v = bertini::IndexToSubscript(1ul,dimensions);
	BOOST_CHECK(v==solution);

	solution[0] = 0; solution[1] = 1;
	v = bertini::IndexToSubscript(2ul,dimensions);
	BOOST_CHECK(v==solution);

	solution[0] = 1; solution[1] = 1;
	v = bertini::IndexToSubscript(3ul,dimensions);
	BOOST_CHECK(v==solution);

}



BOOST_AUTO_TEST_CASE(index_and_subscript_generation2)
{

	size_t index = 20;

	std::vector<size_t> dimensions{2,3,4,5};

	std::vector<size_t> v = bertini::IndexToSubscript(index,dimensions);

	std::vector<size_t> solution{0,1,3,0};
	BOOST_CHECK(v==solution);
}

BOOST_AUTO_TEST_CASE(index_and_subscript_generation3)
{

	size_t index = 119;

	std::vector<size_t> dimensions{2,3,4,5};

	std::vector<size_t> v = bertini::IndexToSubscript(index,dimensions);

	std::vector<size_t> solution{1,2,3,4};
	BOOST_CHECK(v==solution);
}



BOOST_AUTO_TEST_CASE(index_and_subscript_generation_out_of_range)
{


	std::vector<size_t> dimensions{2,3,4,5};
	BOOST_CHECK_THROW(bertini::IndexToSubscript(120ul,dimensions),std::out_of_range);
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




