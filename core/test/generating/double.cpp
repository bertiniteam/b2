//This file is part of Bertini 2.
//
//test/generating/double.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/generating/double.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/generating/double.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

/**
\file test/generating/double.cpp  Tests the generating of text from double's.
*/


#include "bertini2/io/generators.hpp"
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(double_generation)

BOOST_AUTO_TEST_CASE(zero)
{	
	bertini::DefaultPrecision(30);
	double z(0);


	std::string result;
    std::back_insert_iterator<std::string> sink(result);

    BOOST_CHECK(bertini::generators::Classic::generate(sink, z));

    BOOST_CHECK_EQUAL(std::stod(result),z);
}


BOOST_AUTO_TEST_CASE(one)
{   
    bertini::DefaultPrecision(30);
    double z(1);


    std::string result;
    std::back_insert_iterator<std::string> sink(result);

    BOOST_CHECK(bertini::generators::Classic::generate(sink, z));

    BOOST_CHECK_EQUAL(std::stod(result),z);
}


BOOST_AUTO_TEST_CASE(sqrt_2)
{   
    bertini::DefaultPrecision(30);
    double z = sqrt(2);


    std::string result;
    std::back_insert_iterator<std::string> sink(result);

    BOOST_CHECK(bertini::generators::Classic::generate(sink, z));

    BOOST_CHECK_EQUAL(std::stod(result),z);
}

BOOST_AUTO_TEST_CASE(check_135_477005)
{   
    bertini::DefaultPrecision(30);
    double z = 135.477005;


    std::string result;
    std::back_insert_iterator<std::string> sink(result);

    BOOST_CHECK(bertini::generators::Classic::generate(sink, z));

    BOOST_CHECK_EQUAL(std::stod(result),z);
}


BOOST_AUTO_TEST_SUITE_END()
