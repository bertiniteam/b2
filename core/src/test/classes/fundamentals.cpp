//This file is part of Bertini 2.0.
//
//fundamentals.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//fundamentals.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with fundamentals.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  fundamentals.cpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015


#include <boost/test/unit_test.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <iostream>


using mpfr = boost::multiprecision::mpfr_float;

BOOST_AUTO_TEST_SUITE(super_fundamentals)


BOOST_AUTO_TEST_CASE(precision_through_arithemetic)
{
	mpfr::default_precision(50);

	mpfr x("0.01234567890123456789012345678901234567890123456789");

	mpfr::default_precision(30);
	mpfr y = pow(x,2);

	BOOST_CHECK_EQUAL(y.precision(), 30);

	mpfr z = x;

	BOOST_CHECK_EQUAL(z.precision(), 30);

	BOOST_CHECK(fabs(z - mpfr("0.012345678901234567890123456789")) < 1e-30);



	mpfr::default_precision(70);

	z = x;

	BOOST_CHECK_EQUAL(z.precision(),50);
	BOOST_CHECK(fabs(z - mpfr("0.01234567890123456789012345678901234567890123456789")) < 1e-50);


	y.precision(70);
	z.precision(30);

	BOOST_CHECK_EQUAL(y.precision(),70);
	BOOST_CHECK_EQUAL(z.precision(),30);
	BOOST_CHECK_EQUAL(x.precision(),50);

	y = z*x;

	// y is of precision 70 because it's a pre-existing variable.
	BOOST_CHECK_EQUAL(y.precision(), 70);
}





BOOST_AUTO_TEST_SUITE_END()




