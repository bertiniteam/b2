//This file is part of Bertini 2.
//
//test/generating/mpfr_complex.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/generating/mpfr_complex.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/generating/mpfr_complex.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

/**
\file test/generating/mpfr_complex.cpp  Tests the generating of text from mpfr_complex's.
*/


#include "bertini2/io/generators.hpp"
#include "bertini2/io/parsers.hpp"

#include <boost/test/unit_test.hpp>



using mpfr = bertini::mpfr;

BOOST_AUTO_TEST_SUITE(mpfr_complex_generating)

BOOST_AUTO_TEST_SUITE(round_trip)



BOOST_AUTO_TEST_CASE(zero)
{	
	bertini::DefaultPrecision(30);
	mpfr z(0);


	std::string result;
    std::back_insert_iterator<std::string> sink(result);

    BOOST_CHECK(bertini::generators::Classic::generate(sink, z));

    mpfr rt;
    BOOST_CHECK(bertini::parsers::Classic::parse(result.begin(), result.end(), rt));
    BOOST_CHECK_EQUAL(rt,z);
}

BOOST_AUTO_TEST_CASE(one)
{	
	bertini::DefaultPrecision(30);
	mpfr z(1);


	std::string result;
    std::back_insert_iterator<std::string> sink(result);

    BOOST_CHECK(bertini::generators::Classic::generate(sink, z));

    mpfr rt;
    BOOST_CHECK(bertini::parsers::Classic::parse(result.begin(), result.end(), rt));
    BOOST_CHECK_EQUAL(rt,z);
}

BOOST_AUTO_TEST_CASE(sqrt_2)
{   
    bertini::DefaultPrecision(30);
    auto z = sqrt(mpfr(2));


    std::string result;
    std::back_insert_iterator<std::string> sink(result);

    BOOST_CHECK(bertini::generators::Classic::generate(sink, z));

    mpfr rt;
    BOOST_CHECK(bertini::parsers::Classic::parse(result.begin(), result.end(), rt));
    BOOST_CHECK_EQUAL(rt,z);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
