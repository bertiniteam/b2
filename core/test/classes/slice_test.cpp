//This file is part of Bertini 2.
//
//slice_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//slice_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with slice_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

//  slice_test.cpp
//

//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring 2016

/**
\file slice_test.cpp Unit testing for slicing
*/

#include <boost/test/unit_test.hpp>

#include "bertini2/system/slice.hpp"

BOOST_AUTO_TEST_SUITE(linear_slicing)

using namespace bertini;
using Var = std::shared_ptr<bertini::node::Variable>;

using bertini::DefaultPrecision;

BOOST_AUTO_TEST_CASE(slice_basic_complex)
{
	Var x = MakeVariable("x"), y = MakeVariable("y"), z = MakeVariable("z");

	VariableGroup vars{x,y,z};

	auto s = LinearSlice::RandomComplex(vars,2);

	BOOST_CHECK_EQUAL(s.Dimension(),2);
	BOOST_CHECK_EQUAL(s.NumVariables(),3);
}


BOOST_AUTO_TEST_CASE(slice_basic_crazy_overslice)
{
	Var x = MakeVariable("x"), y = MakeVariable("y"), z = MakeVariable("z");

	VariableGroup vars{x,y,z};

	auto s = LinearSlice::RandomComplex(vars,6);

	BOOST_CHECK_EQUAL(s.Dimension(),6);
	BOOST_CHECK_EQUAL(s.NumVariables(),3);
}


BOOST_AUTO_TEST_CASE(slice_basic_real)
{
	Var x = MakeVariable("x"), y = MakeVariable("y"), z = MakeVariable("z");

	VariableGroup vars{x,y,z};

	auto s = LinearSlice::RandomReal(vars,2);

	BOOST_CHECK_EQUAL(s.Dimension(),2);
	BOOST_CHECK_EQUAL(s.NumVariables(),3);
}


BOOST_AUTO_TEST_SUITE_END()

